! Coarse-grained tight-binding projection after Liu et al. (2022).
! Implemented by various AIs, coordinated by Anh-Luan Phan (Sep+Oct 2026)
!
! ======================================================================
! OVERVIEW (read this first)
! ======================================================================
! Mid-layer between orbital Hamiltonian assembly (upt%ham, CSR) and the
! eigensolvers (Lanczos, JD, LAPACK, SLEPc). Builds a smaller reduced
! operator H_red = Q^H H Q, solves on it, then lifts eigenvectors back.
!
! Single implementation: ICGN (controlled by upt%icgn_* parameters).
!
!   selfenergy_order < 0  — plain projection (Liu CG / ICG without Sigma):
!       partition atoms into blocks (METIS, or graph-BFS fallback),
!       diagonalize each block, keep states in the core+buffer window
!       [core_emin - bottom_buffer, core_emax + top_buffer], project
!       inter-block couplings; only round-off-level entries (|V| < 1e-10 eV) are dropped.
!       epsilon (outer_coupling_threshold) only decides which discarded states are
!       promoted to outers (|V_ij|^2 > epsilon*|E_i-E_j| with a core state i).
!
!   selfenergy_order >= 0 — same projection plus Neumann self-energy:
!       Sigma = H_PQ (E0 - D)^-1 [W (E0 - D)^-1]^order H_QP
!       (P = retained, Q = discarded), merged into H_red.
!       order = 0, 1, 2, ... is the series truncation.
!
! Optional core-state outers (add_core_outers > 0): also keep discarded
! states that couple strongly to core states (level-1 only).
!
! Life cycle (uptight.F90):
!   1. After H is built: icgn_set_original, then icgn_prepare fills
!      upt%icgn_ham / icgn_U and upt%icgn_blocks (rows, eval, q, ...).
!   2. icgn_swap_in points upt%ham at the reduced operator; solvers see
!      only that matrix (no CG awareness required).
!   3. After solve: icgn_finalize_eigenvectors lifts (physical = Q * reduced)
!      and restores the original operator via icgn_swap_out.
!
! Parallelism (OpenMP, prepare only): block diagonalization and pair
! projections over blocks/pairs. Thread budget = omp_get_max_threads()
! (capped by UPT_CG_THREADS); nesting flattened so BLAS/LAPACK never
! oversubscribes.
!
! CSR convention: Mi row pointers, Mj columns, M values; sparse_fmt
! 'F'/'U'/'L'. Reduced operator always 'F', diagonal present, columns
! strictly ascending per row. Input may be F/U/L (see stored_entry).
! Orbitals of one atom are contiguous (n_spin * n_st(atom)); blocks are
! whole atoms.
! ======================================================================
module coarse_grain
   use, intrinsic :: iso_c_binding, only : c_int, c_double, c_char
  use precision, only : dp
  use upt_param, only : OUPT, ICGNBlock
  use sparse_matrix, only : CSR, create_matrix, destroy_matrix
      use mpi_globals, only : num_procs, id0, id, shift_init, shift_end, &
         shift_init_Mi, shift_end_Mi
  implicit none
  private

  ! Set to .true. to print detailed coarse-graining statistics (per-block
  ! progress, timings, METIS/BFS partition notes...). They are routed through
  ! icgn_log_progress. Messages that matter to the user (dimension summary,
  ! fall-backs, Neumann diagnostics) use icgn_log_important and are always shown.
  logical, parameter :: ICGN_LOG_STATS = .false.

  ! A ICGNPair holds the data of one pair (a,b), a<b, of blocks.
  type ICGNPair
     integer :: a = 0, b = 0
     ! v = transformed coupling; h = physical-space coupling collected sparsely
     complex(dp), dimension(:,:), pointer :: v => null()
     complex(dp), dimension(:,:), pointer :: h => null()
     ! t = H_AB * Q_B accumulated directly from the global CSR Hamiltonian.
     ! Used by the ICGN projection path to avoid materializing dense H_AB.
     complex(dp), dimension(:,:), pointer :: t => null()
     ! g = full (all block states) projected coupling S_a^H H_ab S_b, kept only on
     ! the outer path so the Neumann correction need not project again.
     complex(dp), dimension(:,:), pointer :: g => null()
  end type ICGNPair

  ! ---- Block-sparse Neumann correction (see icgn_sigma_blocksparse) --------
  ! Q-space data of one block: local indices of the discarded (Q) states and the
  ! resolvent R = 1/(E0 - E) of each of them.
  type ICGNQBlock
     integer :: np = 0, nq = 0
     integer, pointer :: qidx(:) => null()
     real(dp), pointer :: res(:) => null()
  end type ICGNQBlock

  ! Dense (filtered) coupling blocks of one block pair (A,B), A<B:
  !   gab = g(P_A, Q_B), gba = g(P_B, Q_A)  (already the conjugate-transposed
  !   orientation used by the edge lists),  wab = W(Q_A, Q_B),  W(B,A) = wab^H.
  type ICGNQPair
     complex(dp), pointer :: gab(:,:) => null()
     complex(dp), pointer :: gba(:,:) => null()
     complex(dp), pointer :: wab(:,:) => null()
  end type ICGNQPair

  ! Pointer to a matrix (array-of-pointers helper)
  type ICGNMat
     complex(dp), pointer :: m(:,:) => null()
  end type ICGNMat

  ! ---- Geometry of the nested (cube -> row -> plane) partition ---------------
  ! Cube id ic = ix + nx*(iy-1) + nx*ny*(iz-1) (all 1-based), where ix/iy/iz run
  ! along nested_axes(1:3).  The cubes of one row (fixed iy,iz) and the rows of
  ! one plane (fixed iz) therefore have CONSECUTIVE ids.  perm(1:n) lists the
  ! orbitals sorted by (cube id, orbital index): every cube, row and plane owns
  ! one contiguous slice of perm (cstart marks the cube slices), and pos is the
  ! inverse of perm.  The position of an orbital inside any unit (cube/row/plane)
  ! is pos(orbital) - (start of that unit's slice) + 1; Q matrices of every level
  ! use exactly this row order.
  type ICGNHier
     integer :: nx = 0, ny = 0, nz = 0
     integer, allocatable :: cube_of_orb(:), perm(:), pos(:), cstart(:)
  end type ICGNHier

  ! Retained-fraction statistics of the nested tiers (1 = cubes, 2 = rows,
  ! 3 = planes).  Updated from OpenMP loops under a critical section.
  integer(8), save :: nst_cnt(3) = 0_8, nst_tot(3) = 0_8, nst_ret(3) = 0_8
  real(dp), save :: nst_fsum(3) = 0.0_dp

  ! Round-off filter: couplings with |V| < icgn_zero_tol (eV) are stored as zero.
  ! Chosen well above accumulated round-off of the projections (~1e-14 eV) and far
  ! below any coupling that matters physically.
  real(dp), parameter :: icgn_zero_tol = 1.0e-10_dp

  ! Block-sparse storage of the Neumann self-energy Sigma (P-P space).  Row a of
  ! blocks holds the blocks Sigma(P_a, P_blk(t)) for blk(t) >= a only; the part
  ! below the block diagonal follows from Hermiticity (see icgn_csr_add_sigma).
  type ICGNSigRow
     integer :: n = 0
     integer, pointer :: blk(:) => null()
     type(ICGNMat), pointer :: mat(:) => null()
  end type ICGNSigRow

  public :: icgn_configure, icgn_prepare, icgn_clear, icgn_active, icgn_lift
  public :: icgn_log_progress, icgn_log_important
  public :: icgn_get_active, icgn_lift_active
  public :: icgn_get_info
  public :: icgn_swap_in, icgn_swap_out, icgn_set_original, icgn_forget, icgn_destroy_operator, &
           icgn_finalize_eigenvectors

  interface
       subroutine upt_log_message(message, length) bind(C, name='upt_log_message')
          import :: c_char, c_int
          character(kind=c_char), intent(in) :: message(*)
          integer(c_int), value :: length
       end subroutine upt_log_message

       integer(c_int) function icgn_metis_partition(nvtxs, xadj, adjncy, vwgt, &
          adjwgt, nparts, ufactor, seed, part) bind(C, name='upt_icgn_metis_partition')
       import :: c_int
       integer(c_int), value :: nvtxs, nparts, ufactor, seed
       integer(c_int), intent(in) :: xadj(*), adjncy(*), vwgt(*), adjwgt(*)
       integer(c_int), intent(out) :: part(*)
     end function icgn_metis_partition
  end interface

  ! Saved state of the reduced operator while it is swapped into upt%ham, see
  ! icgn_swap_in()/icgn_finalize_eigenvectors() below. The original operator is
  ! remembered separately by icgn_set_original() at the moment it is built, so the
  ! restore can never pick up the reduced matrix by mistake.
  type(CSR), target  :: icgn_original_ham          ! value snapshot, NOT an alias
  type(CSR), target  :: icgn_original_u
  logical            :: icgn_original_valid = .false.
  type(CSR), pointer :: icgn_saved_ham => null()
  type(CSR), pointer :: icgn_saved_u   => null()
  integer :: icgn_saved_shift_init = 0
  integer :: icgn_saved_shift_end = 0
  integer :: icgn_saved_shift_init_mi = 0
  integer :: icgn_saved_shift_end_mi = 0
  integer :: icgn_saved_n_spin = 0
  ! .true. when icgn_U is Q^T U Q (Kramers shortcut valid on reduced space).
  ! .false. when identity fallback is used → swap_in forces n_spin=1.
  logical :: icgn_u_projected = .false.

contains

! ===========================================================================
! SECTION A — Lifecycle & public API
! Swap/restore, configure, query, clear, prepare, lift.
! Solvers only see H via swap_in; they never know about coarse-graining.
! ===========================================================================

  ! Called by the orchestrator right after the orbital operator has been built
  ! and before the reduced one is prepared.
  !
  ! Ownership: this is a SHALLOW snapshot of the CSR header (pointer addresses
  ! only, not a deep copy of the matrix data). The caller MUST keep the original
  ! Hamiltonian heap alive until icgn_swap_out / icgn_forget / icgn_destroy_operator.
  ! Destroying or reallocating upt%ham between set_original and restore is a
  ! contract violation and leaves dangling pointers in the snapshot.
  subroutine icgn_set_original(upt)
    type(OUPT), intent(inout), target :: upt

    ! A previous solve may have exited before its normal finalize hook. If the
    ! original snapshot is still valid, restore it before taking a fresh snapshot;
    ! otherwise fail closed and keep the current swap bookkeeping intact.
    if (associated(icgn_saved_ham)) then
       if (.not. icgn_original_valid) then
          call icgn_log_important(upt, &
               '(icgn) set_original: cannot refresh snapshot while original operator is invalid')
          return
       end if
       if (.not. associated(icgn_original_ham%M)) then
          call icgn_log_important(upt, &
               '(icgn) set_original: cannot refresh snapshot; original Hamiltonian storage is unavailable')
          return
       end if
       if (.not. associated(icgn_original_ham%Mi)) then
          call icgn_log_important(upt, &
               '(icgn) set_original: cannot refresh snapshot; original row pointers are unavailable')
          return
       end if
       call icgn_swap_out(upt)
    end if

    icgn_original_ham = upt%ham    ! value copy of struct header (pointer addresses, not heap data)
    icgn_original_u   = upt%U
    icgn_original_valid = associated(icgn_original_ham%M) .and. associated(icgn_original_ham%Mi)

  end subroutine icgn_set_original

  !----------------------------------------------------------------------------!
  ! Single entry point for using a reduced operator from outside uptight.
  !
  ! Solvers do not know whether they are working on the original or on a
  ! reduced operator: the swap is done here, once, so that upt%ham (and the
  ! dimension dependent shift bookkeeping every driver relies on) already
  ! describes whatever operator has to be solved. icgn_finalize_eigenvectors()
  ! brings the eigenvectors back to the orbital basis and restores the original
  ! operator; icgn_swap_out() only restores.
  !----------------------------------------------------------------------------!

  subroutine icgn_swap_in(upt, active)
    type(OUPT), intent(inout), target :: upt
    logical, intent(out) :: active
    type(CSR), pointer :: reduced_ham, reduced_u
    logical :: got
    integer :: nred

    active = .false.
    if (associated(icgn_saved_ham)) then
      active = .true.
      return
    end if

    ! Never guess which matrix the original one was. icgn_set_original() is the only
    ! thing allowed to record it, and it runs right after the orbital operator has
    ! been built. Without that record, coarse-graining simply does not activate.
    if (.not. icgn_original_valid) return

    call icgn_get_active(upt, reduced_ham, reduced_u, got)
    if (.not. got) return

    nred = reduced_ham%nrow

    icgn_saved_ham => icgn_original_ham
    icgn_saved_u   => icgn_original_u
    icgn_saved_shift_init = shift_init
    icgn_saved_shift_end = shift_end
    icgn_saved_shift_init_mi = shift_init_Mi(id)
    icgn_saved_shift_end_mi = shift_end_Mi(id)
    icgn_saved_n_spin = upt%n_spin

    upt%ham = reduced_ham
    upt%U   = reduced_u
    ! When U_red = Q^T U Q was built, keep n_spin so Lanczos/JD can still use the
    ! Kramers shortcut on the reduced space (same contract as direct mode).
    ! Identity fallback means the shortcut is not valid → force n_spin=1.
    if (.not. icgn_u_projected) upt%n_spin = 1

    ! The result arrays of the previous solve belong to the previous operator and
    ! are kept across k-points so that band tracking can pick them up. Their row
    ! count belongs to a different basis now, so they must not be handed to the
    ! solver: upt_set_state and the drivers both refuse a mismatch. States found
    ! in the orbital basis cannot seed a search in the reduced one anyway.
    if (associated(upt%eigen_vectors)) deallocate(upt%eigen_vectors)
    if (associated(upt%eigen_values))  deallocate(upt%eigen_values)
    if (associated(upt%particles))     deallocate(upt%particles)

    shift_init = 1
    shift_end = nred
    shift_init_Mi(id) = 1
    shift_end_Mi(id) = nred

    active = .true.
  end subroutine icgn_swap_in

  subroutine icgn_swap_out(upt)
    type(OUPT), intent(inout), target :: upt

    ! No-op unless a reduced operator is currently in place. Without this guard
    ! a call made after the Hamiltonian has been destroyed would put the dangling
    ! icgn_original_ham pointer back into upt%ham.
    if (.not. associated(icgn_saved_ham)) return

    ! Restore only when the shallow snapshot still looks alive. If the original
    ! heap was freed against the set_original contract, skip the assignment
    ! rather than planting a dangling CSR into upt%ham.
    if (icgn_original_valid .and. associated(icgn_original_ham%M) .and. &
        associated(icgn_original_ham%Mi)) then
      upt%ham = icgn_original_ham
      upt%U   = icgn_original_u
    else if (icgn_original_valid) then
      write(*,'(a)') '(icgn) swap_out: original Hamiltonian snapshot is no longer associated; skip restore'
      icgn_original_valid = .false.
    end if

    shift_init = icgn_saved_shift_init
    shift_end = icgn_saved_shift_end
    shift_init_Mi(id) = icgn_saved_shift_init_mi
    shift_end_Mi(id) = icgn_saved_shift_end_mi
    upt%n_spin = icgn_saved_n_spin
    nullify(icgn_saved_ham, icgn_saved_u)
  end subroutine icgn_swap_out

  ! Brings the eigenvectors back to the orbital basis and restores the original
  ! operator. This is the only thing a caller has to do after a solve that ran on
  ! a reduced operator: the solvers themselves are unaware of it.
  subroutine icgn_finalize_eigenvectors(upt)
    type(OUPT), intent(inout), target :: upt
    complex(dp), allocatable :: reduced(:,:), lifted(:,:)
    integer :: nfull, nred, num_ev
    character(len=256) :: log_message

    if (.not. associated(icgn_saved_ham)) return
    if (.not. associated(upt%eigen_vectors)) then
      call icgn_swap_out(upt)
      return
    end if

    nred  = size(upt%eigen_vectors, 1)
    num_ev = size(upt%eigen_vectors, 2)

    ! Take the orbital dimension from what the coarse-graining layer recorded
    ! when it built the reduced operator. icgn_saved_ham only says which matrix was
    ! in place when the swap happened, and it is not a reliable source for the
    ! original size once the Hamiltonian has been rebuilt.
    nfull = icgn_orbital_dim(upt)
    ! nfull == nred is legitimate (every state retained): the vectors still have
    ! to be rotated back from the block eigenbasis, so only nfull < nred, or a
    ! row count that differs from the recorded reduced dimension, is an error.
    if (nfull < nred .or. nred /= upt%icgn_reduced_dim) then
      write(log_message,'(a,i0,a,i0,a,i0,a)') &
           '(icgn) finalize: inconsistent dimensions (nred=', nred, &
           ', reduced_dim=', upt%icgn_reduced_dim, ', nfull=', nfull, &
           '), keeping reduced vectors'
      call icgn_log_important(upt, trim(log_message))
      call icgn_swap_out(upt)
      return
    end if

    write(log_message,'(a,i0,a,i0,a,i0)') '(icgn) finalize: lifting nred=', nred, &
         ' -> nfull=', nfull, ', num_ev=', num_ev
    call icgn_log_progress(upt, trim(log_message))

    allocate(reduced(nred, num_ev))
    reduced = upt%eigen_vectors(1:nred, 1:num_ev)

    allocate(lifted(nfull, num_ev))
    ! Map the eigenvectors of the reduced operator back to the physical basis.
    call icgn_lift_active(upt, reduced, lifted)

    deallocate(upt%eigen_vectors)
    allocate(upt%eigen_vectors(nfull, num_ev))
    upt%eigen_vectors = lifted

    deallocate(reduced, lifted)

    call icgn_swap_out(upt)
  end subroutine icgn_finalize_eigenvectors

  ! Orbital dimension of the original Hamiltonian, as recorded when the reduced
  ! operator was built. Zero when coarse-graining is not active.
  function icgn_orbital_dim(upt) result(n)
    type(OUPT), intent(in) :: upt
    integer :: n
    logical :: ready, conv
    integer :: od, rd, nb
    real(dp) :: cf, s2

    n = 0
    if (icgn_active(upt)) then
      call icgn_get_info(upt, ready, od, rd, nb, cf, s2, conv)
      if (ready) n = od
    end if
  end function icgn_orbital_dim

  ! Drops every pointer into the reduced operator. Must be called before the
  ! operator is destroyed, otherwise icgn_original_ham/icgn_original_u and the
  ! swap bookkeeping keep pointing at freed memory.
  subroutine icgn_forget(upt)
    type(OUPT), intent(inout) :: upt
    nullify(icgn_saved_ham, icgn_saved_u)
    icgn_original_valid = .false.
    ! Null out pointer components to avoid stale pointer values (defensive)
    nullify(icgn_original_ham%M, icgn_original_ham%Mi, icgn_original_ham%Mj)
    icgn_original_ham%nrow = 0
    icgn_original_ham%nnz  = 0
    nullify(icgn_original_u%M, icgn_original_u%Mi, icgn_original_u%Mj)
    icgn_original_u%nrow = 0
    icgn_original_u%nnz  = 0
    icgn_saved_shift_init = 0
    icgn_saved_shift_end = 0
    icgn_saved_shift_init_mi = 0
    icgn_saved_shift_end_mi = 0
    icgn_saved_n_spin = 0
    icgn_u_projected = .false.
  end subroutine icgn_forget

  ! Drops every pointer into the reduced operator. Only for teardown. Must not touch the record of
  ! the original operator.
  subroutine icgn_destroy_operator(upt)
    type(OUPT), intent(inout), target :: upt
    integer :: i
    ! Do not destroy icgn_ham while upt%ham still aliases it. If the original
    ! snapshot has gone stale, fail closed and preserve the live operator.
    if (associated(icgn_saved_ham)) then
       if (.not. icgn_original_valid .or. .not. associated(icgn_original_ham%M) .or. &
           .not. associated(icgn_original_ham%Mi)) then
          write(*,'(a)') '(icgn) destroy_operator: original Hamiltonian unavailable; keeping swapped operator intact'
          return
       end if
       call icgn_swap_out(upt)
    end if
    call icgn_forget(upt)

    if (associated(upt%icgn_ham%M)) then

       call destroy_matrix(upt%icgn_ham)

    end if
    if (associated(upt%icgn_U%M)) call destroy_matrix(upt%icgn_U)
    if (associated(upt%icgn_blocks)) then
       do i = 1, size(upt%icgn_blocks)
          if (associated(upt%icgn_blocks(i)%rows)) deallocate(upt%icgn_blocks(i)%rows)
          if (associated(upt%icgn_blocks(i)%eval)) deallocate(upt%icgn_blocks(i)%eval)
          if (associated(upt%icgn_blocks(i)%q)) deallocate(upt%icgn_blocks(i)%q)
          if (associated(upt%icgn_blocks(i)%evals_full)) deallocate(upt%icgn_blocks(i)%evals_full)
          if (associated(upt%icgn_blocks(i)%S_full)) deallocate(upt%icgn_blocks(i)%S_full)
          if (associated(upt%icgn_blocks(i)%retained_idx)) deallocate(upt%icgn_blocks(i)%retained_idx)
       end do
       deallocate(upt%icgn_blocks)

    end if
    upt%icgn_ready = .false.
    upt%icgn_original_dim = 0
    upt%icgn_reduced_dim = 0
    upt%icgn_cut_fraction = 0.0_dp
  end subroutine icgn_destroy_operator

  !
  ! icgn_prepare: partition atoms → diagonalize blocks → select retained
  ! states (core + buffer + optional level-1 outers) → project couplings →
  ! build icgn_ham CSR. When selfenergy_order >= 0, add Neumann Sigma in
  ! block-sparse form and merge into the CSR. When selfenergy_order < 0 the
  ! correction is skipped (plain projection).
  ! ==========================================================================

  subroutine icgn_configure(upt, enabled, nblocks, core_emin, core_emax, &
     top_buffer, bottom_buffer, epsilon, add_core_outers, selfenergy_order, E0, metis_imbalance, &
       check_convergence, pi_maxiter, pi_tol)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled, check_convergence
    integer, intent(in) :: nblocks, add_core_outers, selfenergy_order, pi_maxiter
   real(dp), intent(in) :: core_emin, core_emax, top_buffer, bottom_buffer, epsilon, E0, metis_imbalance, pi_tol
    upt%icgn_enabled = enabled
    upt%icgn_num_blocks = nblocks
    upt%icgn_core_emin = core_emin
    upt%icgn_core_emax = core_emax
   upt%icgn_top_buffer = top_buffer
   upt%icgn_bottom_buffer = bottom_buffer
    upt%icgn_epsilon = epsilon
    upt%icgn_add_core_outers = max(0, add_core_outers)
    upt%icgn_selfenergy_order = selfenergy_order
    upt%icgn_E0 = E0
    ! Nested tiers default to the scalar values; UPT_set_coarse_graining_nested_tiers
    ! may override them afterwards (it must be called after this routine).
    upt%icgn_tier_emin = core_emin
    upt%icgn_tier_emax = core_emax
    upt%icgn_tier_top  = top_buffer
    upt%icgn_tier_bot  = bottom_buffer
    upt%icgn_tier_eps  = epsilon
    upt%icgn_tier_E0   = E0
    upt%icgn_nested_outers = max(0, add_core_outers)
    upt%icgn_nested_orders = selfenergy_order
    upt%icgn_metis_imbalance = metis_imbalance
    ! Temporarily disabled Neumann convergence checker; preserve the old
    ! assignments below as comments for future reconsideration.
    ! upt%icgn_check_convergence = check_convergence
    ! upt%icgn_pi_maxiter = pi_maxiter
    ! upt%icgn_pi_tol     = pi_tol
    upt%icgn_check_convergence = .false.
    upt%icgn_pi_maxiter = 1000
    upt%icgn_pi_tol = 1.0e-3_dp
    upt%icgn_sigma_T2   = -1.0_dp
    upt%icgn_pi_converged = .true.
    upt%icgn_ready = .false.
  end subroutine icgn_configure

  logical function icgn_active(upt)
    type(OUPT), intent(in) :: upt
    icgn_active = upt%icgn_enabled .and. upt%icgn_ready
  end function icgn_active

  subroutine icgn_get_info(upt, ready, orig_dim, red_dim, nblocks, cut_frac, &
       sigma_T2, pi_converged)
    type(OUPT), intent(in) :: upt
    logical, intent(out) :: ready, pi_converged
    integer, intent(out) :: orig_dim, red_dim, nblocks
    real(dp), intent(out) :: cut_frac, sigma_T2
    ready        = upt%icgn_ready
    orig_dim     = upt%icgn_original_dim
    red_dim      = upt%icgn_reduced_dim
    nblocks      = upt%icgn_num_blocks
    cut_frac     = upt%icgn_cut_fraction
    sigma_T2     = upt%icgn_sigma_T2
    pi_converged = upt%icgn_pi_converged
  end subroutine icgn_get_info

  subroutine icgn_get_active(upt, active_ham, active_u, active)
      type(OUPT), intent(in), target :: upt
    type(CSR), pointer, intent(out) :: active_ham, active_u
    logical, intent(out) :: active
    active = .false.
    nullify(active_ham, active_u)
    if (icgn_active(upt)) then
       active_ham => upt%icgn_ham
       active_u => upt%icgn_U
       active = .true.
    end if
  end subroutine icgn_get_active

  subroutine icgn_lift_active(upt, reduced, physical)
    type(OUPT), intent(in) :: upt
    complex(dp), intent(in) :: reduced(:,:)
    complex(dp), intent(out) :: physical(:,:)
    if (icgn_active(upt)) then
       call icgn_lift(upt, reduced, physical)
    else
       physical = reduced
    end if
  end subroutine icgn_lift_active

  subroutine icgn_clear(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    ! Never free the reduced matrix while upt%ham aliases it. Restore the
    ! original first; if its snapshot is no longer valid, leave the operator
    ! untouched rather than create a dangling pointer.
    if (associated(icgn_saved_ham)) then
       if (.not. icgn_original_valid .or. .not. associated(icgn_original_ham%M) .or. &
           .not. associated(icgn_original_ham%Mi)) then
          write(*,'(a)') '(icgn) clear: original Hamiltonian unavailable; keeping swapped operator intact'
          return
       end if
       call icgn_swap_out(upt)
    end if

    if (associated(upt%icgn_ham%M)) then

       call destroy_matrix(upt%icgn_ham)

    end if

    if (associated(upt%icgn_U%M)) then

       call destroy_matrix(upt%icgn_U)

    end if
    if (associated(upt%icgn_blocks)) then

       do i = 1, size(upt%icgn_blocks)
          if (associated(upt%icgn_blocks(i)%rows))         deallocate(upt%icgn_blocks(i)%rows)
          if (associated(upt%icgn_blocks(i)%eval))         deallocate(upt%icgn_blocks(i)%eval)
          if (associated(upt%icgn_blocks(i)%q))            deallocate(upt%icgn_blocks(i)%q)
          if (associated(upt%icgn_blocks(i)%evals_full))   deallocate(upt%icgn_blocks(i)%evals_full)
          if (associated(upt%icgn_blocks(i)%S_full))       deallocate(upt%icgn_blocks(i)%S_full)
          if (associated(upt%icgn_blocks(i)%retained_idx)) deallocate(upt%icgn_blocks(i)%retained_idx)
       end do
       deallocate(upt%icgn_blocks)

    end if
    upt%icgn_ready        = .false.
    upt%icgn_original_dim = 0
    upt%icgn_reduced_dim  = 0
    upt%icgn_cut_fraction = 0.0_dp

  end subroutine icgn_clear

  ! ---------------------------------------------------------------------------
  ! icgn_prepare: build the reduced operator and optional Neumann correction.
  ! pairs hold inter-block couplings (physical h and/or projected v/g) reused
  ! for H_PP assembly and, when selfenergy_order >= 0, for H_PQ / H_QQ.
  ! ---------------------------------------------------------------------------
  subroutine icgn_prepare(upt, ierr)
    character(len=*), parameter :: core_outer_warning = &
         'Including core-state outers higher than level 1 is not yet implemented --> ' // &
         'Fall back to level-1 outers only.'
    type(OUPT), intent(inout) :: upt
    integer, intent(out) :: ierr

    integer :: n, na, i, j, k, p, status, nedge, maxedge
    integer :: r, c, br, bc, npair, total_ret, pos
    integer :: hash_size, slot, min_atom, max_atom
    integer(8), allocatable :: edge_key(:)
    integer, allocatable :: atom_of(:), label(:), row_of(:)
    integer, allocatable :: counts(:), cursor(:), offsets(:), bsize(:)
    integer, allocatable :: blk_of(:), ierr_blk(:)
    integer(c_int), allocatable :: xadj(:), adjncy(:), vwgt(:), adjwgt(:), part(:)
    real(dp), allocatable :: edge_weight(:)
    real(dp) :: max_weight, all_weight, cut_weight
    ! Item 8 experiment temporarily disabled for performance comparison:
    ! real(dp), parameter :: energy_margin = 1.0e-3_dp
    ! real(dp) :: emin_keep, emax_keep, best_dist, boundary_eval
    ! integer :: idx_lo, idx_hi
    type(ICGNPair), allocatable :: pairs(:)

    logical, allocatable :: is_core(:,:), keep_mask(:,:)

    integer :: ia, ib, nred, nb
    integer(8) :: npq, nqq
    integer :: nth, saved_levels
    integer, allocatable :: roff_blk(:), local_ret_idx(:,:)
    type(ICGNSigRow), allocatable :: sig(:)
    real(dp) :: max_sigma, E0_used, t_all, t_stage, t_mark
    logical :: need_qq
    ! logical :: fast_path  ! item 8 experiment temporarily disabled
    ! Block-sparse Neumann correction (P-Q / Q-Q blocks + Sigma via ICGNSigRow)
    type(ICGNQBlock), allocatable :: bq(:)
    type(ICGNQPair), allocatable :: qp(:)
    integer, allocatable :: nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, allocatable :: nbr_isA(:)
    integer(8) :: npq_nz, nqq_nz

    ierr = 0
    npq = 0; nqq = 0

    call icgn_clear(upt)

    if (.not. upt%icgn_enabled) return
   call icgn_log_progress(upt, 'mode=icgn preparation started')
    if (num_procs /= 1) then
       ierr = 1; write(*,*) '(icgn) MPI not supported'; return
    end if

    na = upt%basis%n_basis
    n  = upt%ham%nrow
    if (na < 1) then
       ierr = 2; write(*,*) '(icgn) invalid n_atoms'; return
    end if
    if (upt%icgn_num_blocks < 0) then
       ierr = 2; write(*,*) '(icgn) num_blocks must be non-negative (0 = auto)'; return
    end if
    if (upt%icgn_partition_type /= 2) call icgn_resolve_num_blocks(upt, na, n)
   if (upt%icgn_core_emin >= upt%icgn_core_emax .or. upt%icgn_top_buffer < 0.0_dp .or. &
      upt%icgn_bottom_buffer < 0.0_dp) then
       ierr = 3; write(*,*) '(icgn) invalid core window or buffer'; return
    end if

    if (.not. associated(upt%ham%M)) then
       ierr = 4; write(*,*) '(icgn) Hamiltonian not initialized'; return
    end if

    ! ---- Build atom→orbital mapping (identical to icgn_prepare) ------------

    allocate(atom_of(n), offsets(na+1), bsize(na))

    pos = 1

    do i = 1, na
       offsets(i) = pos
       bsize(i)   = upt%n_spin * upt%basis%n_st(i)
       do j = 1, bsize(i)
          atom_of(pos) = i; pos = pos + 1
       end do
    end do

    offsets(na+1) = pos
    if (pos-1 /= n) then
       ierr = 5; write(*,*) '(icgn) atom/orbital mapping inconsistent'
       call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
       return
    end if

    ! ---- Nested hierarchical CG (cube→row→plane→bulk) -----------------------
    if (upt%icgn_partition_type == 2) then
       ! Per-tier add_core_outers / neumann_order live in icgn_nested_outers/orders.
       call icgn_nested_prepare(upt, na, n, offsets, bsize, atom_of, ierr)
       if (allocated(atom_of)) deallocate(atom_of)
       if (allocated(offsets)) deallocate(offsets)
       if (allocated(bsize)) deallocate(bsize)
       return
    end if

    ! ---- Block partition: geometry (1D PL) or coupling (METIS/graph) --------
    allocate(label(na), counts(na))
    if (upt%icgn_partition_type == 1) then
       allocate(part(na))
       call icgn_geometry_partition(upt, na, bsize, part, status)
       if (status /= 0) then
          ierr = status
          call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
          return
       end if
       do i = 1, na; label(i) = int(part(i)) + 1; end do
       deallocate(part)
       call icgn_log_progress(upt, 'mode=icgn geometry partition complete')
    else
       maxedge = max(1, upt%ham%nnz)
       ! Aggregate all orbital couplings into one weighted undirected edge per atom
       ! pair. The open-addressed hash avoids passing duplicate edges to METIS.
       hash_size = max(17, 2*maxedge + 1)
       allocate(edge_key(hash_size), edge_weight(hash_size))
       edge_key = 0_8; edge_weight = 0.0_dp
       nedge = 0; max_weight = 0.0_dp; counts = 0
       do r = 1, n
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k)
             if (atom_of(r) == atom_of(c)) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             br = atom_of(r); bc = atom_of(c)
             min_atom = min(br, bc); max_atom = max(br, bc)
             slot = int(mod(int(min_atom-1,8)*int(na,8) + int(max_atom,8), int(hash_size,8))) + 1
             do while (edge_key(slot) /= 0_8 .and. edge_key(slot) /= int(min_atom-1,8)*int(na,8)+int(max_atom,8))
                slot = slot + 1
                if (slot > hash_size) slot = 1
             end do
             if (edge_key(slot) == 0_8) then
                edge_key(slot) = int(min_atom-1,8)*int(na,8) + int(max_atom,8)
                nedge = nedge + 1
                counts(min_atom) = counts(min_atom) + 1
                counts(max_atom) = counts(max_atom) + 1
             end if
             edge_weight(slot) = edge_weight(slot) + abs(upt%ham%M(k))**2
          end do
       end do
       do slot = 1, hash_size
          if (edge_key(slot) /= 0_8) max_weight = max(max_weight, edge_weight(slot))
       end do
       if (nedge == 0 .or. max_weight == 0.0_dp) then
          ierr = 7; write(*,*) '(icgn) atom graph has no couplings'
          call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
          return
       end if
       allocate(xadj(na+1), cursor(na), adjncy(2*nedge), adjwgt(2*nedge), vwgt(na), part(na))
       xadj(1) = 0_c_int
       do i = 1, na
          xadj(i+1) = xadj(i) + int(counts(i), c_int)
          cursor(i)  = int(xadj(i)) + 1
          vwgt(i)    = int(bsize(i), c_int)
       end do
       do slot = 1, hash_size
          if (edge_key(slot) == 0_8) cycle
          br = int((edge_key(slot)-1_8) / int(na,8)) + 1
          bc = int(mod(edge_key(slot)-1_8, int(na,8))) + 1
          p = max(1, nint(edge_weight(slot) / max_weight * 1000000.0_dp))
          adjncy(cursor(br)) = int(bc-1,c_int); adjwgt(cursor(br)) = int(p,c_int); cursor(br)=cursor(br)+1
          adjncy(cursor(bc)) = int(br-1,c_int); adjwgt(cursor(bc)) = int(p,c_int); cursor(bc)=cursor(bc)+1
       end do
       deallocate(edge_key)
       ! METIS UFACTOR is (imbalance - 1) * 1000 (e.g. 1.03 -> 30),
       ! not imbalance * 1000.
       status = icgn_metis_partition(int(na,c_int), xadj, adjncy, vwgt, adjwgt, &
            int(upt%icgn_num_blocks,c_int), &
            int(max(1, nint(1000.0_dp*(upt%icgn_metis_imbalance-1.0_dp))),c_int), 42_c_int, part)
       if (status /= 0) then
          call icgn_log_important(upt, 'METIS unavailable or failed; using built-in graph-BFS fallback partitioning')
          call icgn_graph_partition(na, upt%icgn_num_blocks, vwgt, xadj, adjncy, adjwgt, part)
       else
          call icgn_log_progress(upt, 'METIS partitioning done')
       end if
       do i = 1, na; label(i) = int(part(i)) + 1; end do
       call icgn_log_progress(upt, 'mode=icgn graph partition complete')
    end if

    all_weight = 0.0_dp; cut_weight = 0.0_dp
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_weight = all_weight + abs(upt%ham%M(k))**2
          if (label(atom_of(r)) /= label(atom_of(c))) &
               cut_weight = cut_weight + abs(upt%ham%M(k))**2
          ! Geometry PL must be block-tridiagonal (or circulant under PBC):
          ! no H coupling between non-adjacent blocks. With periodic BC along
          ! the partition axis, blocks 1 and nb are neighbours (ring).
          if (upt%icgn_partition_type == 1 .and. abs(upt%ham%M(k)) > 0.0_dp) then
             block
               integer :: lb, lc, dblk, nb_chk
               lb = label(atom_of(r)); lc = label(atom_of(c))
               nb_chk = upt%icgn_num_blocks
               dblk = abs(lb - lc)
               if (upt%basis%periodic .and. nb_chk > 1) &
                    dblk = min(dblk, nb_chk - dblk)
               if (dblk > 1) then
                  ierr = 24
                  write(*,'(a,i0,a,i0,a,i0,a)') &
                       '(icgn) ERROR: geometry partition has non-adjacent block coupling |i-j|=', &
                       dblk, ' (blocks ', lb, ' and ', lc, &
                       '). PL too thin vs ETB interaction range; increase layer content or reduce num_blocks.'
                  call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
                  return
               end if
             end block
          end if
       end do
    end do
    if (all_weight > 0.0_dp) then
       upt%icgn_cut_fraction = cut_weight / all_weight
    else
       upt%icgn_cut_fraction = 0.0_dp
    end if

    ! ---- Allocate icgn_blocks and fill rows arrays --------------------------
    deallocate(counts)
    nb = upt%icgn_num_blocks
    allocate(upt%icgn_blocks(nb), counts(nb))
    counts = 0
    do i = 1, na; counts(label(i)) = counts(label(i)) + bsize(i); end do
    do i = 1, nb
       upt%icgn_blocks(i)%nrow = counts(i)
       allocate(upt%icgn_blocks(i)%rows(counts(i)))
    end do
    ! cursor is only allocated on the METIS/graph path; geometry needs it too.
    if (allocated(cursor)) deallocate(cursor)
    allocate(cursor(nb))
    cursor = 0
    do i = 1, na
       br = label(i)
       do j = offsets(i), offsets(i+1)-1
          cursor(br) = cursor(br) + 1
          upt%icgn_blocks(br)%rows(cursor(br)) = j
       end do
    end do

    ! ---- Diagonalize each block fully (store S_full, evals_full) ------------
    ! Blocks are independent.  When there are at least as many blocks as threads,
    ! each thread owns whole blocks and LAPACK runs single-threaded inside the
    ! region (nesting is flattened), so the CPU is never oversubscribed.  With
    ! fewer blocks than threads, blocks run one after another and the threaded
    ! LAPACK/BLAS uses the full thread budget instead.
    allocate(row_of(n), blk_of(n), ierr_blk(nb)); row_of = 0; blk_of = 0; ierr_blk = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
          blk_of(upt%icgn_blocks(i)%rows(j)) = i
       end do
    end do
    nth = icgn_thread_budget()
    if (nb < nth) nth = 1
    nth = icgn_outer_threads(nb, nth, maxval(upt%icgn_blocks(1:nb)%nrow))
    t_all = icgn_wall(); t_stage = t_all
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel do default(shared) private(i) num_threads(nth) schedule(dynamic,1) if(nth > 1)
    do i = 1, nb
       !$omp critical (icgn_log)
       call icgn_log_block(upt, 'icgn', i, upt%icgn_blocks(i)%nrow, 'processing')
       !$omp end critical (icgn_log)
       call icgn_diagonalize_block(upt, i, row_of, blk_of, ierr_blk(i))
    end do
    !$omp end parallel do
    call icgn_restore_nesting(saved_levels)
    call icgn_report_time('diagonalize blocks', t_stage)
    block
       character(len=8) :: tm
       integer :: ltm, stm
       call get_environment_variable('UPT_CG_TIMING', tm, length=ltm, status=stm)
       if (ltm > 0) write(*,'(a,i0,a,i0,a,i0,a,i0)') '[ICGN-TIME] diag info: nblocks=', nb, &
            ' max_block_dim=', maxval(upt%icgn_blocks(1:nb)%nrow), ' outer_threads=', nth, ' budget=', icgn_thread_budget()
    end block
    do i = 1, nb
       if (ierr_blk(i) /= 0) then
          ierr = ierr_blk(i)
          call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
          return
       end if
    end do
    deallocate(ierr_blk, blk_of)

    ! ---- Core + buffer selection, then boundary-degeneracy expansion --------
    ! Nominal window: [core_emin - bottom_buffer, core_emax + top_buffer].
    ! Greedy expand at both edges: if consecutive (energy-sorted) states have
    ! |E_{k+1}-E_k| <= delta, keep the outside one too (all-or-nothing multiplet
    ! at the cut).  Fast-path ZHEEVR uses a 1e-3 guard so those near-edge states
    ! exist; dummy uncomputed columns are skipped.
    allocate(is_core(nb, maxval(counts)), keep_mask(nb, maxval(counts)))
    is_core = .false.; keep_mask = .false.

    !$omp parallel do default(shared) private(i, j) schedule(dynamic,1)
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          if (.not. associated(upt%icgn_blocks(i)%evals_full)) cycle
          if (.not. associated(upt%icgn_blocks(i)%S_full)) cycle
          ! Skip ZHEEVR placeholder columns (zero vectors, sentinel evals).
          if (maxval(abs(upt%icgn_blocks(i)%S_full(:, j))) <= 0.0_dp) cycle
          associate(e => upt%icgn_blocks(i)%evals_full(j))
            if (e >= upt%icgn_core_emin .and. e <= upt%icgn_core_emax) then
               is_core(i,j)   = .true.
               keep_mask(i,j) = .true.
            else if (e >= upt%icgn_core_emin - upt%icgn_bottom_buffer .and. &
                     e <= upt%icgn_core_emax + upt%icgn_top_buffer) then
               keep_mask(i,j) = .true.
            end if
          end associate
       end do
       call icgn_expand_keep_boundary(upt%icgn_blocks(i), keep_mask(i, :), &
            upt%icgn_core_emin - upt%icgn_bottom_buffer, &
            upt%icgn_core_emax + upt%icgn_top_buffer, 1.0e-6_dp)
    end do
    !$omp end parallel do

    row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
    end do

    ! Level-1 outer: collect physical block couplings, then project
    ! each pair once at matrix level.  The zero-outer/no-Neumann case
    ! defers pair construction until q is available and uses the ICGN-style
    ! sparse projection path instead of materializing dense physical H_ab.
    if (upt%icgn_add_core_outers > 0 .or. upt%icgn_selfenergy_order >= 0) then
       ! collect_pairs allocates pairs to the exact number of edged pairs.
       t_stage = icgn_wall()
       call icgn_collect_pairs(upt, atom_of, label, row_of, pairs, npair, ierr)
       if (ierr /= 0) then
          call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
          return
       end if
       call icgn_report_time('collect pair couplings H_ab', t_stage)
       if (icgn_timing_on()) call icgn_report_boundary(pairs, npair)
    else
       ! Plain path: pairs filled later by sparse fast-path in build.
       allocate(pairs(0)); npair = 0
    end if
    if (upt%icgn_add_core_outers > 0) then
       ! Pair projections are independent (each touches only pairs(i)).  The full
       ! projected coupling pairs(i)%v and the physical block pairs(i)%h are kept:
       ! build_icgn_reduced_hamiltonian only slices v, and the Neumann correction
       ! reuses h, so nothing is recollected or reprojected later.
       nth = min(icgn_thread_budget(), max(1, npair))
       nth = icgn_outer_threads(npair, nth, maxval(upt%icgn_blocks(1:nb)%nrow))
       t_stage = icgn_wall()
       call icgn_flatten_nesting(saved_levels)
       !$omp parallel do default(shared) private(i, ia, ib) num_threads(nth) schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          call project_pair(pairs(i), upt%icgn_blocks(ia)%S_full, upt%icgn_blocks(ib)%S_full, .true.)
       end do
       !$omp end parallel do
       call icgn_restore_nesting(saved_levels)
       call icgn_report_time('project outer pairs', t_stage)
    end if

    ! Outer selection is optional. Level 1 is the currently implemented
    ! core-state outer level; higher requested levels fall back to it.
    t_stage = icgn_wall()
    if (upt%icgn_add_core_outers > 1) then
       call icgn_log_important(upt, core_outer_warning)
    end if
    if (upt%icgn_add_core_outers > 0) then
       do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       do j = 1, upt%icgn_blocks(ia)%nrow
          if (is_core(ia, j)) then
             do k = 1, upt%icgn_blocks(ib)%nrow
                if (.not. keep_mask(ib, k)) then
                   if (icgn_outer_passes(pairs(i)%v(j,k), upt%icgn_blocks(ia)%evals_full(j), &
                       upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) &
                       keep_mask(ib, k) = .true.
                end if
             end do
          end if
       end do
       do k = 1, upt%icgn_blocks(ib)%nrow
          if (is_core(ib, k)) then
             do j = 1, upt%icgn_blocks(ia)%nrow
                if (.not. keep_mask(ia, j)) then
                   if (icgn_outer_passes(pairs(i)%v(j,k), upt%icgn_blocks(ib)%evals_full(k), &
                       upt%icgn_blocks(ia)%evals_full(j), upt%icgn_epsilon)) &
                       keep_mask(ia, j) = .true.
                end if
             end do
          end if
       end do
    end do
       end if
    call icgn_report_time('outer selection (keep_mask, serial)', t_stage)
    ! The projected pairs (v, h) are kept: build_icgn_reduced_hamiltonian slices v
    ! (reuse_full) and the Neumann correction reuses h.  Freed later.

    t_mark = icgn_wall()
    ! ---- Apply keep_mask to nret, q, eval, retained_idx --------------------
    total_ret = 0
    do i = 1, nb
       upt%icgn_blocks(i)%nret = count(keep_mask(i, 1:upt%icgn_blocks(i)%nrow))
       total_ret = total_ret + upt%icgn_blocks(i)%nret
    end do
    if (total_ret == 0) then
       ierr = 9; write(*,*) '(icgn) no states retained'
       call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
       return
    end if
    upt%icgn_original_dim = n
    upt%icgn_reduced_dim  = total_ret
    ! nred is used below by icgn_U also when the Neumann correction is off
    ! (selfenergy_order < 0); it must not depend on that branch.
    nred = total_ret

    do i = 1, nb
       if (upt%icgn_blocks(i)%nret == 0) cycle
       allocate(upt%icgn_blocks(i)%eval(upt%icgn_blocks(i)%nret))
       allocate(upt%icgn_blocks(i)%q(upt%icgn_blocks(i)%nrow, upt%icgn_blocks(i)%nret))
       allocate(upt%icgn_blocks(i)%retained_idx(upt%icgn_blocks(i)%nret))
       j = 0
       do k = 1, upt%icgn_blocks(i)%nrow
          if (.not. keep_mask(i, k)) cycle
          j = j + 1
          upt%icgn_blocks(i)%eval(j)         = upt%icgn_blocks(i)%evals_full(k)
          upt%icgn_blocks(i)%q(:, j)         = upt%icgn_blocks(i)%S_full(:, k)
          upt%icgn_blocks(i)%retained_idx(j) = k
       end do
    end do
    deallocate(is_core, keep_mask)

    ! ---- Rebuild row_of for reduced-ham build -------------------------------
    row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
    end do

    call icgn_report_time('apply keep_mask (q/eval/retained_idx)', t_mark)
    ! ---- Build reduced Hamiltonian (keeping pairs for Sigma) ----------------
    t_stage = icgn_wall()
    call build_icgn_reduced_hamiltonian(upt, atom_of, label, row_of, pairs, npair, &
         upt%icgn_add_core_outers > 0, ierr)
    if (ierr /= 0) then
       call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
       return
    end if
    call icgn_report_time('build reduced H (pair projection)', t_stage)

    ! ====================================================================
    ! Self-energy correction via Neumann series (skipped when neumann_order < 0)
    ! ====================================================================
    if (upt%icgn_selfenergy_order >= 0) then
       ! Neumann correction is applied as requested, but its convergence is
       ! deliberately NOT checked. The user is responsible for validating it.
       call icgn_log_important(upt, &
            'WARNING: Uptight does NOT check for Neumann-series convergence. Use it with care!')
       ! H_PP already in upt%icgn_ham (full sorted CSR). Build Sigma in block-sparse
       ! form and merge; no dense nred×nred intermediate.
       t_mark = icgn_wall()
       nred = total_ret

       allocate(roff_blk(nb+1)); roff_blk(1) = 1
       do i = 1, nb; roff_blk(i+1) = roff_blk(i) + upt%icgn_blocks(i)%nret; end do

       allocate(local_ret_idx(nb, maxval([(upt%icgn_blocks(i)%nrow, i=1,nb)])))
       local_ret_idx = 0
       do i = 1, nb
          do j = 1, upt%icgn_blocks(i)%nret
             local_ret_idx(i, upt%icgn_blocks(i)%retained_idx(j)) = j
          end do
       end do
       call icgn_report_time('Neumann setup (index tables)', t_mark)

       E0_used = upt%icgn_E0
       ! Only fall back when E0 lies outside the open core window.
       ! E0 = 0 is a valid choice whenever 0 is inside (core_min, core_max).
       if (E0_used <= upt%icgn_core_emin .or. E0_used >= upt%icgn_core_emax) then
          E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp
          call icgn_log_important(upt, &
               'E0 outside (core_min, core_max); using window midpoint instead')
       end if
       need_qq = (upt%icgn_selfenergy_order >= 1)
       t_stage = icgn_wall()
       call icgn_build_neumann_blocks(upt, upt%icgn_blocks, nb, pairs, npair, local_ret_idx, E0_used, need_qq, &
            bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, npq_nz, nqq_nz)
       ! Keep nonzero counts in int64; converting to default INTEGER here could
       ! wrap negative and silently skip Sigma or the convergence diagnostic.
       npq = npq_nz; nqq = nqq_nz
       call icgn_report_time('Neumann: block P-Q/Q-Q blocks', t_stage)
       t_mark = icgn_wall()

       call icgn_free_block_eigensystems(upt)
       call destroy_pairs(pairs)

       upt%icgn_sigma_T2     = -1.0_dp
       upt%icgn_pi_converged = .true.
       ! Temporarily disabled: do not run the Neumann convergence checker.
       ! if (upt%icgn_check_convergence .and. nqq > 0) then
       !    t_stage = icgn_wall()
       !    call icgn_power_iteration_blk(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, &
       !         upt%icgn_pi_maxiter, upt%icgn_pi_tol, upt%icgn_sigma_T2, upt%icgn_pi_converged)
       !    call icgn_report_convergence(upt)
       !    call icgn_report_time('power iteration (block matvec)', t_stage)
       ! end if
       call icgn_report_time('free S_full/pairs', t_mark)

       allocate(sig(nb))
       t_stage = icgn_wall()
       if (npq > 0) call icgn_sigma_blocksparse(nb, upt%icgn_selfenergy_order, &
            maxval(upt%icgn_blocks(1:nb)%nrow), bq, qp, &
            nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, sig)
       call icgn_report_time('Sigma (block-sparse zgemm)', t_stage)
       call icgn_free_neumann_blocks(bq, qp, npair)
       deallocate(nbr_ptr, nbr_blk, nbr_pair, nbr_rev, nbr_isA)

       max_sigma = 0.0_dp
       do i = 1, nb
          do j = 1, sig(i)%n
             if (associated(sig(i)%mat(j)%m)) &
                  max_sigma = max(max_sigma, maxval(abs(sig(i)%mat(j)%m)))
          end do
       end do
       if (icgn_timing_on()) then
          write(*,'(a,i0,a,i0,a,i0,a,es10.3)') '[ICGN-TIME] stats: nred=', nred, &
               ' npq=', npq, ' nqq=', nqq, ' max|Sigma|=', max_sigma
       end if
       t_mark = icgn_wall()
       deallocate(local_ret_idx)

       call icgn_csr_add_sigma(upt%icgn_ham, roff_blk, sig, nred, ierr)
       if (ierr /= 0) then
          call icgn_log_important(upt, &
               'ERROR: Sigma merge aborted because CSR entry count exceeds default-integer capacity')
          call icgn_free_sigrows(sig)
          deallocate(roff_blk)
          call icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
          ierr = 12
          return
       end if
       call icgn_free_sigrows(sig)
       deallocate(roff_blk)
       call icgn_report_time('merge Sigma into CSR (block-sparse)', t_mark)
    end if  ! icgn_selfenergy_order >= 0 (Neumann correction)

    ! When the Neumann branch is skipped (selfenergy_order < 0) nothing above has
    ! released the pair matrices (v/h/t) nor the full block eigensystems.
    ! Free them here; both calls are no-ops if the Neumann branch already did it.
    if (allocated(pairs)) call destroy_pairs(pairs)
    call icgn_free_block_eigensystems(upt)

    ! Build icgn_U: Q^T U Q when original U exists (Kramers-capable), else identity.
    ! pairs already freed above; U build only needs blocks%q and upt%U.
    call icgn_build_reduced_U(upt, nred, ierr)
    if (ierr /= 0) return

    ! Last timing line: wall from start of block diagonalization through U_red.
    call icgn_report_time('TOTAL ICGN PREPARATION TIME', t_all)

    upt%icgn_ready = .true.
       call icgn_log_info(upt, 'icgn', n, total_ret, &
          upt%icgn_num_blocks, upt%icgn_cut_fraction, upt%icgn_sigma_T2, &
          upt%icgn_pi_converged, &
          upt%icgn_check_convergence .and. upt%icgn_selfenergy_order >= 0 .and. nqq > 0)

    ! Cleanup (geometry path never allocates the METIS/graph temporaries)
    if (allocated(atom_of)) deallocate(atom_of)
    if (allocated(offsets)) deallocate(offsets)
    if (allocated(bsize)) deallocate(bsize)
    if (allocated(label)) deallocate(label)
    if (allocated(row_of)) deallocate(row_of)
    if (allocated(counts)) deallocate(counts)
    if (allocated(cursor)) deallocate(cursor)
    if (allocated(edge_weight)) deallocate(edge_weight)
    if (allocated(xadj)) deallocate(xadj)
    if (allocated(adjncy)) deallocate(adjncy)
    if (allocated(adjwgt)) deallocate(adjwgt)
    if (allocated(vwgt)) deallocate(vwgt)
    if (allocated(part)) deallocate(part)

  end subroutine icgn_prepare

  ! Lift eigenvectors from reduced basis back to orbital (physical) space:
  ! physical(rows_of_block,:) = q_block * reduced(offset : offset+nret-1, :).
  ! One zgemm per block (same helper used elsewhere in this module).
  subroutine icgn_lift(upt, reduced, physical)
    type(OUPT), intent(in) :: upt
    complex(dp), intent(in) :: reduced(:,:)
    complex(dp), intent(out) :: physical(:,:)
    integer :: i, j, off, nn, nr, nev
    complex(dp), allocatable :: tmp(:,:)
    physical = (0.0_dp, 0.0_dp)
    off = 1
    nev = size(reduced, 2)
    do i = 1, size(upt%icgn_blocks)
       nr = upt%icgn_blocks(i)%nret
       if (nr > 0) then
          nn = upt%icgn_blocks(i)%nrow
          allocate(tmp(nn, nev))
          call icgn_zgemm('N', upt%icgn_blocks(i)%q, reduced(off:off+nr-1, :), tmp)
          do j = 1, nn
             physical(upt%icgn_blocks(i)%rows(j), :) = tmp(j, :)
          end do
          deallocate(tmp)
       end if
       off = off + nr
    end do
  end subroutine icgn_lift

  ! Build upt%icgn_U = Q^T U Q (block-diagonal) so Lanczos/JD keep the Kramers
  ! shortcut on the reduced space.  The solver applies partner = conjg(U*v)
  ! (antilinear).  With real U (TB onsite ±I), matching the physical partner
  ! Q^H conjg(U Q v) requires transpose, not conjugate-transpose:
  !   U_red = Q^T U Q
  ! Falls back to identity when original U is missing or n_spin /= 2; then
  ! icgn_u_projected stays .false. and swap_in forces n_spin=1.
  !
  ! Parallelism: dense U_red_ii per CG block (same OpenMP pattern as diag).
  ! CSR assembly is serial.
  subroutine icgn_build_reduced_U(upt, nred, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: nred
    integer, intent(out) :: ierr
    integer :: nb, i, jr, k, r, c, lc, nn, nr, oa, nnz, pos, jcol, nfull
    integer :: nth, saved_levels, max_nn
    integer, allocatable :: roff(:), locmap(:), rowcount(:), next(:)
    complex(dp), allocatable :: u_loc(:,:), work(:,:)
    type(ICGNMat), allocatable :: ured(:)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    real(dp), parameter :: drop = 1.0e-14_dp
    real(dp) :: t_mark

    ierr = 0
    icgn_u_projected = .false.
    if (associated(upt%icgn_U%M)) call destroy_matrix(upt%icgn_U)

    ! Fallback: identity (no Kramers on reduced space).
    if (.not. associated(upt%U%M) .or. .not. associated(upt%U%Mi) .or. &
        upt%n_spin /= 2 .or. nred <= 0) then
       call create_matrix(upt%icgn_U, nred, nred, max(nred, 1))
       upt%icgn_U%sparse_fmt = 'F'
       upt%icgn_U%Mi(1) = 1
       do i = 1, nred
          upt%icgn_U%Mj(i) = i
          upt%icgn_U%M(i)  = (1.0_dp, 0.0_dp)
          upt%icgn_U%Mi(i+1) = i + 1
       end do
       upt%icgn_U%nnz = nred
       return
    end if

    nb = upt%icgn_num_blocks
    nfull = upt%U%nrow
    allocate(roff(nb+1))
    roff(1) = 1
    do i = 1, nb
       roff(i+1) = roff(i) + upt%icgn_blocks(i)%nret
    end do

    ! ---- Parallel: project U onto each retained block basis ----------------
    allocate(ured(nb))
    do i = 1, nb
       nullify(ured(i)%m)
    end do

    max_nn = 1
    do i = 1, nb
       if (upt%icgn_blocks(i)%nret > 0) max_nn = max(max_nn, upt%icgn_blocks(i)%nrow)
    end do
    nth = icgn_thread_budget()
    if (nb < nth) nth = 1
    nth = icgn_outer_threads(nb, nth, max_nn)
    t_mark = icgn_wall()
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel do default(shared) private(i, jr, k, r, c, lc, nn, nr, locmap, u_loc, work) &
    !$omp& num_threads(nth) schedule(dynamic,1) if(nth > 1)
    do i = 1, nb
       nr = upt%icgn_blocks(i)%nret
       if (nr == 0) cycle
       if (.not. associated(upt%icgn_blocks(i)%q)) cycle
       nn = upt%icgn_blocks(i)%nrow
       allocate(locmap(nfull))
       locmap = 0
       do jr = 1, nn
          locmap(upt%icgn_blocks(i)%rows(jr)) = jr
       end do
       allocate(u_loc(nn, nn))
       u_loc = czero
       do jr = 1, nn
          r = upt%icgn_blocks(i)%rows(jr)
          do k = upt%U%Mi(r), upt%U%Mi(r+1) - 1
             c = upt%U%Mj(k)
             if (c < 1 .or. c > nfull) cycle
             lc = locmap(c)
             if (lc == 0) cycle
             u_loc(jr, lc) = upt%U%M(k)
          end do
       end do
       allocate(work(nn, nr), ured(i)%m(nr, nr))
       call icgn_zgemm('N', u_loc, upt%icgn_blocks(i)%q, work)  ! U_loc * q
       call icgn_zgemm('T', upt%icgn_blocks(i)%q, work, ured(i)%m)  ! q^T *(U_loc q)
       deallocate(u_loc, work, locmap)
    end do
    !$omp end parallel do
    call icgn_restore_nesting(saved_levels)
    call icgn_report_time('build U_red (Q^T U Q, parallel)', t_mark)

    ! ---- Serial: count nnz and assemble CSR --------------------------------
    allocate(rowcount(nred))
    rowcount = 0
    do i = 1, nb
       if (.not. associated(ured(i)%m)) cycle
       nr = size(ured(i)%m, 1)
       oa = roff(i)
       do jr = 1, nr
          do jcol = 1, nr
             if (abs(ured(i)%m(jr, jcol)) > drop) &
                  rowcount(oa + jr - 1) = rowcount(oa + jr - 1) + 1
          end do
       end do
    end do

    nnz = sum(rowcount)
    if (nnz <= 0) then
       do i = 1, nb
          if (associated(ured(i)%m)) deallocate(ured(i)%m)
       end do
       deallocate(ured, roff, rowcount)
       call create_matrix(upt%icgn_U, nred, nred, max(nred, 1))
       upt%icgn_U%sparse_fmt = 'F'
       upt%icgn_U%Mi(1) = 1
       do i = 1, nred
          upt%icgn_U%Mj(i) = i
          upt%icgn_U%M(i)  = (1.0_dp, 0.0_dp)
          upt%icgn_U%Mi(i+1) = i + 1
       end do
       upt%icgn_U%nnz = nred
       return
    end if

    call create_matrix(upt%icgn_U, nred, nred, nnz)
    upt%icgn_U%sparse_fmt = 'F'
    upt%icgn_U%Mi(1) = 1
    do i = 1, nred
       upt%icgn_U%Mi(i+1) = upt%icgn_U%Mi(i) + rowcount(i)
    end do
    upt%icgn_U%nnz = nnz
    allocate(next(nred))
    do i = 1, nred
       next(i) = upt%icgn_U%Mi(i)
    end do

    do i = 1, nb
       if (.not. associated(ured(i)%m)) cycle
       nr = size(ured(i)%m, 1)
       oa = roff(i)
       do jr = 1, nr
          do jcol = 1, nr
             if (abs(ured(i)%m(jr, jcol)) <= drop) cycle
             pos = next(oa + jr - 1)
             upt%icgn_U%Mj(pos) = oa + jcol - 1
             upt%icgn_U%M(pos)  = ured(i)%m(jr, jcol)
             next(oa + jr - 1)  = pos + 1
          end do
       end do
       deallocate(ured(i)%m)
    end do

    call icgn_sort_csr_rows(upt%icgn_U)

    deallocate(ured, roff, rowcount, next)
    icgn_u_projected = .true.
    call icgn_log_progress(upt, '(icgn) U_red = Q^T U Q built (Kramers shortcut enabled on reduced space)')
  end subroutine icgn_build_reduced_U

! ===========================================================================
! SECTION B — Prepare helpers (partition, diagonalize, pairs, H_red)
! ===========================================================================

  ! Resolve upt%icgn_num_blocks: 0 => auto heuristic; always clamp to [1, na].
  ! Auto aims for ~target_orb orbitals per block (balance diag cost ~ N^3/B^2
  ! against parallelism) and prefers B >= thread budget so outer parallel-do
  ! over blocks can fill the CPUs.
  subroutine icgn_resolve_num_blocks(upt, na, n_orb)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: na, n_orb
    integer :: nb, T, requested
    character(len=256) :: message
    integer, parameter :: target_orb = 200
    requested = upt%icgn_num_blocks
    if (requested == 0) then
       T = icgn_thread_budget()
       nb = max(1, (max(1, n_orb) + target_orb - 1) / target_orb)  ! ceil(n_orb/target)
       nb = max(nb, T)
       write(message,'(a,i0,a,i0,a,i0,a,i0,a)') &
            '(icgn) auto num_blocks=', nb, ' (n_orb=', n_orb, &
            ', n_atoms=', na, ', threads=', T, ')'
       call icgn_log_important(upt, trim(message))
    else
       nb = requested
    end if
    if (nb > na) then
       write(message,'(a,i0,a,i0,a)') &
            '(icgn) WARNING: num_blocks=', nb, ' exceeds n_atoms=', na, &
            ' — falling back to num_blocks = n_atoms'
       call icgn_log_important(upt, trim(message))
       nb = na
    end if
    if (nb < 1) nb = 1
    upt%icgn_num_blocks = nb
  end subroutine icgn_resolve_num_blocks

  ! ============================================================================
  ! icgn_graph_partition: connectivity-aware fallback partition when METIS is
  ! unavailable.
  !
  ! Algorithm: weighted greedy graph growing (BFS-seeded, priority-queue-free),
  ! with dynamic target rebalancing and best-fit neighbour selection.
  ! Guarantees:
  !   (1) Every atom in a block is reachable from the block seed via bonds —
  !       so no block is a disjoint set of atoms (a block may still consist of
  !       several BFS components if the seed's own component runs out before
  !       reaching the target weight; each such component is itself internally
  !       connected).
  !   (2) All orbitals of an atom go to the same block (works at the atom level,
  !       so the orbital → block mapping is done after by the caller).
  !   (3) Balance: after each block is closed, the target weight for remaining
  !       blocks is recomputed from the remaining unassigned weight, so any
  !       overshoot/undershoot in one block is spread over the rest instead of
  !       accumulating into the last block.
  !
  ! Inputs:
  !   na         — number of atoms
  !   nblocks    — desired number of blocks
  !   vwgt(na)   — vertex weight of each atom (= number of orbitals, bsize(i))
  !   xadj(na+1) — CSR row pointers of atom adjacency graph (0-based)
  !   adjncy(*)  — CSR column indices (0-based atom indices)
  !   adjwgt(*)  — CSR edge weights (integer, Hamiltonian coupling strength)
  !
  ! Output:
  !   part(na)   — block index (0-based, in [0, nblocks-1]) for each atom
  ! ============================================================================
  subroutine icgn_graph_partition(na, nblocks, vwgt, xadj, adjncy, adjwgt, part)
    use, intrinsic :: iso_c_binding, only : c_int
    integer,            intent(in)  :: na, nblocks
    integer(c_int),     intent(in)  :: vwgt(na)
    integer(c_int),     intent(in)  :: xadj(na+1), adjncy(*), adjwgt(*)
    integer(c_int),     intent(out) :: part(na)

    ! --- local ---
    integer :: i, j, atom, nb_atom, blk, seed
    integer :: target_wt, cur_wt, total_wt, remaining_wt, remaining_blocks
    integer :: qhead, qtail, qsize
    integer, allocatable :: queue(:)     ! BFS queue (atom indices, 1-based)
    logical, allocatable :: visited(:)  ! atom already assigned?
    integer :: best_deg, deg, adj_start, adj_end
    integer :: n_unvisited
    integer :: best_fit_atom, best_fit_wt, best_fit_gap
    logical :: found_fit

    ! ----- compute total weight -----
    total_wt = 0
    do i = 1, na
       total_wt = total_wt + int(vwgt(i))
    end do

    allocate(visited(na), queue(na))
    visited = .false.
    part    = int(nblocks - 1, c_int)   ! default: last block (catches unvisited atoms)

    remaining_wt     = total_wt
    remaining_blocks = nblocks
    n_unvisited      = na

    do blk = 0, nblocks - 2    ! assign blocks 0 .. nblocks-2; last gets remainder

       if (n_unvisited == 0) exit

       ! ---- recompute target dynamically from what's left, so overshoot/
       ! undershoot in earlier blocks doesn't accumulate into later ones ----
       target_wt = (remaining_wt + remaining_blocks - 1) / remaining_blocks   ! ceiling

       cur_wt = 0
       qhead  = 1
       qtail  = 0
       qsize  = 0

       ! ---- grow the block, possibly across several connected components ----
       do while (cur_wt < target_wt .and. n_unvisited > 0)

          ! ---- need a new BFS seed (either first seed of this block, or a
          ! new component because the previous one ran dry) ----
          if (qhead > qtail) then
             seed     = -1
             best_deg = -1
             do i = 1, na
                if (visited(i)) cycle
                deg = int(xadj(i+1) - xadj(i))   ! number of inter-atom bonds
                if (deg > best_deg) then
                   best_deg = deg
                   seed     = i
                end if
             end do
             if (seed == -1) exit    ! no unvisited atoms left (shouldn't happen, n_unvisited>0)

             visited(seed) = .true.
             part(seed)    = int(blk, c_int)
             cur_wt        = cur_wt + int(vwgt(seed))
             n_unvisited   = n_unvisited - 1
             qtail         = qtail + 1
             queue(qtail)  = seed
             if (cur_wt >= target_wt) exit
          end if

          atom = queue(qhead); qhead = qhead + 1

          adj_start = int(xadj(atom)) + 1    ! convert 0-based xadj to 1-based
          adj_end   = int(xadj(atom+1))

          ! Repeatedly pick the best-fit unvisited neighbour of this atom:
          ! prefer the one whose weight fills up to target_wt without
          ! exceeding it (minimises overshoot); if every remaining neighbour
          ! would overshoot, take the smallest one (minimises the overshoot
          ! amount rather than picking arbitrarily / by edge weight alone).
          ! Cost: O(degree^2) per atom — degree is small for tight-binding.
          do while (cur_wt < target_wt)
             found_fit     = .false.
             best_fit_atom = -1
             best_fit_wt   = -1
             best_fit_gap  = huge(0)   ! smallest (target - cur - w) >= 0 seen so far
             do j = adj_start, adj_end
                i = int(adjncy(j)) + 1    ! 0-based → 1-based
                if (.not. visited(i)) then
                   block
                     integer :: wv, gap
                     wv = int(vwgt(i))
                     gap = target_wt - cur_wt - wv
                     if (gap >= 0) then
                        ! fits without overshoot: keep the tightest fit
                        if (.not. found_fit .or. gap < best_fit_gap) then
                           found_fit     = .true.
                           best_fit_gap  = gap
                           best_fit_atom = i
                           best_fit_wt   = wv
                        end if
                     else if (.not. found_fit) then
                        ! would overshoot: among these, keep the smallest atom
                        ! (minimises how far over target_wt we go)
                        if (best_fit_atom == -1 .or. wv < best_fit_wt) then
                           best_fit_atom = i
                           best_fit_wt   = wv
                        end if
                     end if
                   end block
                end if
             end do
             if (best_fit_atom == -1) exit     ! no more unvisited neighbours of this atom

             nb_atom          = best_fit_atom
             visited(nb_atom) = .true.
             part(nb_atom)    = int(blk, c_int)
             cur_wt           = cur_wt + int(vwgt(nb_atom))
             n_unvisited      = n_unvisited - 1
             qtail            = qtail + 1
             queue(qtail)     = nb_atom
          end do

       end do  ! grow block (possibly multi-component)

       remaining_wt     = remaining_wt - cur_wt
       remaining_blocks = remaining_blocks - 1

    end do  ! blk loop

    ! ---- any unvisited atom (leftover / disconnected remainder) → last block ----
    ! (part already initialised to nblocks-1)

    deallocate(visited, queue)
  end subroutine icgn_graph_partition

  ! ============================================================================
  ! Geometry partition: 1D principal layers (PL) along x/y/z.
  ! 1) Detect atomic layers from projected coordinates (gap clustering).
  ! 2) Pack consecutive layers into K balanced blocks (DP on layer sizes).
  ! part(1:na) is 0-based block id.  On failure ierr /= 0.
  ! ============================================================================
  subroutine icgn_geometry_partition(upt, na, bsize, part, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: na
    integer, intent(in) :: bsize(na)
    integer, intent(out) :: part(na)
    integer, intent(out) :: ierr
    integer :: dir, i, j, L, K, a, best_dir, n_gap
    real(dp) :: extent(3), u, tol, gap, max_small, min_large, jump, best_jump
    real(dp), allocatable :: u_atom(:), u_sorted(:), gaps(:)
    integer, allocatable :: order(:), layer_of_atom(:), layer_size(:), layer_natom(:)
    integer, allocatable :: cum(:), cut(:), dp_prev(:)
    real(dp), allocatable :: cost(:,:)
    real(dp) :: target, mean, dev, best_cost, c
    integer :: s, t, prev, blk, first, last
    character(len=256) :: msg

    ierr = 0
    part = 0
    if (.not. associated(upt%basis%coord)) then
       ierr = 20; write(*,*) '(icgn) geometry partition: basis coordinates missing'; return
    end if
    if (size(upt%basis%coord, 2) < na) then
       ierr = 20; write(*,*) '(icgn) geometry partition: coord size mismatch'; return
    end if

    ! Direction: 1/2/3 = x/y/z; 0 = auto (longest extent)
    dir = upt%icgn_partition_dir
    if (dir < 0 .or. dir > 3) dir = 0
    do j = 1, 3
       extent(j) = maxval(upt%basis%coord(j, 1:na)) - minval(upt%basis%coord(j, 1:na))
    end do
    if (dir == 0) then
       best_dir = 1
       do j = 2, 3
          if (extent(j) > extent(best_dir)) best_dir = j
       end do
       dir = best_dir
       write(msg,'(a,i0,a,3(f10.4,1x))') &
            '(icgn) geometry partition: auto direction=', dir, ' extents=', extent
       call icgn_log_important(upt, trim(msg))
    end if

    allocate(u_atom(na), order(na))
    do i = 1, na
       u_atom(i) = upt%basis%coord(dir, i)
       order(i) = i
    end do
    ! Insertion sort atoms by projected coordinate (na is modest for layer detect).
    do i = 2, na
       a = order(i); u = u_atom(a); j = i - 1
       do while (j >= 1)
          if (u_atom(order(j)) <= u) exit
          order(j+1) = order(j); j = j - 1
       end do
       order(j+1) = a
    end do

    allocate(u_sorted(na), gaps(max(1, na-1)))
    do i = 1, na
       u_sorted(i) = u_atom(order(i))
    end do
    n_gap = na - 1
    do i = 1, n_gap
       gaps(i) = u_sorted(i+1) - u_sorted(i)
    end do

    tol = upt%icgn_layer_tol
    if (tol <= 0.0_dp) then
       ! Auto tol: largest jump between small (intra-layer jitter) and large (inter-layer) gaps.
       if (n_gap < 1) then
          L = 1
          allocate(layer_of_atom(na), layer_size(1), layer_natom(1))
          layer_of_atom = 1; layer_size(1) = sum(bsize); layer_natom(1) = na
       else
          block
            integer, allocatable :: gord(:)
            integer :: gi, gj, gt
            allocate(gord(n_gap))
            do gi = 1, n_gap; gord(gi) = gi; end do
            do gi = 2, n_gap
               gt = gord(gi); u = gaps(gt); gj = gi - 1
               do while (gj >= 1)
                  if (gaps(gord(gj)) <= u) exit
                  gord(gj+1) = gord(gj); gj = gj - 1
               end do
               gord(gj+1) = gt
            end do
            best_jump = -1.0_dp; tol = -1.0_dp
            do gi = 1, n_gap - 1
               jump = gaps(gord(gi+1)) - gaps(gord(gi))
               if (jump > best_jump) then
                  best_jump = jump
                  tol = 0.5_dp * (gaps(gord(gi)) + gaps(gord(gi+1)))
               end if
            end do
            deallocate(gord)
          end block
          if (tol <= 0.0_dp) then
             ierr = 21
             write(*,*) '(icgn) geometry partition: cannot auto-detect layer tolerance (no clear gap clusters)'
             write(*,*) '(icgn) set layer_tolerance explicitly (Angstrom)'
             deallocate(u_atom, order, u_sorted, gaps)
             return
          end if
          write(msg,'(a,f10.6,a)') '(icgn) geometry partition: auto layer_tolerance=', tol, ' Ang'
          call icgn_log_important(upt, trim(msg))
       end if
    end if

    if (.not. allocated(layer_of_atom)) then
       allocate(layer_of_atom(na))
       L = 1
       layer_of_atom(order(1)) = 1
       do i = 2, na
          if (gaps(i-1) > tol) L = L + 1
          layer_of_atom(order(i)) = L
       end do
       allocate(layer_size(L), layer_natom(L))
       layer_size = 0; layer_natom = 0
       do i = 1, na
          layer_size(layer_of_atom(i)) = layer_size(layer_of_atom(i)) + bsize(i)
          layer_natom(layer_of_atom(i)) = layer_natom(layer_of_atom(i)) + 1
       end do
    end if

    write(msg,'(a,i0,a,i0)') '(icgn) geometry partition: atomic layers=', L, ', atoms=', na
    call icgn_log_important(upt, trim(msg))

    K = upt%icgn_num_blocks
    ! Geometry PL: each block needs >= 1 atomic layer. Clamp when auto heuristic
    ! (or an explicit request) exceeds L; never hard-fail on this alone.
    if (K > L) then
       write(msg,'(a,i0,a,i0,a)') &
            '(icgn) geometry: num_blocks=', K, ' > atomic layers=', L, &
            ' — clamping to n_layers (min(L, heuristic/request))'
       call icgn_log_important(upt, trim(msg))
       K = L
       upt%icgn_num_blocks = K
    end if
    if (K < 1) K = 1
    if (K == 1) then
       part = 0
       deallocate(u_atom, order, u_sorted, gaps, layer_of_atom, layer_size, layer_natom)
       return
    end if

    ! DP: pack L consecutive layers into K groups, minimise sum of squared
    ! deviations of group orbital counts from mean N/K.
    mean = real(sum(layer_size), dp) / real(K, dp)
    allocate(cum(0:L))
    cum(0) = 0
    do i = 1, L
       cum(i) = cum(i-1) + layer_size(i)
    end do
    allocate(cost(0:L, 0:K))
    cost = huge(1.0_dp)
    cost(0, 0) = 0.0_dp
    allocate(dp_prev(L*K)); dp_prev = 0
    do t = 1, K
       do j = t, L   ! at least one layer per remaining group: j-s >= K-t, s >= t-1
          best_cost = huge(1.0_dp); prev = t - 1
          do s = t - 1, j - 1
             if (cost(s, t-1) >= huge(1.0_dp) * 0.5_dp) cycle
             dev = real(cum(j) - cum(s), dp) - mean
             c = cost(s, t-1) + dev * dev
             if (c < best_cost) then
                best_cost = c; prev = s
             end if
          end do
          cost(j, t) = best_cost
          dp_prev((t-1)*L + j) = prev
       end do
    end do
    if (cost(L, K) >= huge(1.0_dp) * 0.5_dp) then
       ierr = 23; write(*,*) '(icgn) geometry partition: DP failed'; return
    end if

    ! Backtrack cuts: cut(t) = last layer index of group t
    allocate(cut(0:K))
    cut(K) = L
    do t = K, 1, -1
       cut(t-1) = dp_prev((t-1)*L + cut(t))
    end do
    cut(0) = 0

    ! Map layers -> blocks, atoms -> part (0-based)
    do t = 1, K
       first = cut(t-1) + 1
       last = cut(t)
       do i = 1, na
          if (layer_of_atom(i) >= first .and. layer_of_atom(i) <= last) part(i) = t - 1
       end do
    end do

    write(msg,'(a,i0,a,f12.3)') '(icgn) geometry partition: blocks=', K, ', balance cost=', cost(L, K)
    call icgn_log_important(upt, trim(msg))

    deallocate(u_atom, order, u_sorted, gaps, layer_of_atom, layer_size, layer_natom, cum, cost, dp_prev, cut)
  end subroutine icgn_geometry_partition

  ! ============================================================================
  ! NESTED (MULTI-LEVEL) COARSE-GRAINING:  cube -> row -> plane -> bulk H_CG
  !
  ! The supercell is cut into nx*ny*nz cubes (grids along nested_axes).
  !   level 0: every cube is diagonalized (window states kept)            -> cube units
  !   level 1: cubes of one ROW  -> H_CG(row)   -> diagonalize -> row unit
  !   level 2: rows  of one PLANE-> H_CG(plane) -> diagonalize -> plane unit
  !   level 3: planes -> H_CG(bulk) (NOT diagonalized) handed to the solver.
  !
  ! Generic building block icgn_cg_assemble (no diagonalization): given a group of
  ! units, each with an orthonormal retained basis Q_u (rows = the unit's orbitals
  ! in the basis of H_in) and energies eval_u, it returns the CSR H_CG of the group:
  ! diag(eval) on the unit diagonals, Q_a^H H_ab Q_b on the pairs (a,b) (only the
  ! interface rows/columns of H_ab enter the products, as in the fast path).
  ! H_in is any CSR and the unit/row maps are supplied by ICGNHier, so the same
  ! routine serves every level.
  !
  ! Q bookkeeping: a unit's Q is always expressed in the ORBITAL basis, restricted
  ! to the unit's own orbitals.  Merging a group gives
  !     Q_group = blockdiag(Q_u) * Z        (Z = eigenvectors of H_CG(group))
  ! so the accumulated Q of a plane is supported on that plane's orbitals only:
  ! the final Q is block-diagonal over planes and is returned to uptight as
  ! icgn_blocks(plane) (rows, q, eval) - exactly the contract of the fast path.
  ! No other Q is kept: the Q of cubes/rows is freed as soon as its group is merged.
  !
  ! Memory/parallelism: planes are processed one after another (or T at a time when
  ! #planes >= threads), inside a plane the rows, inside a row the cubes; only the
  ! units of ONE group are alive at any time.  Exactly one loop level is parallel:
  ! the outermost level having at least as many items as threads (items of large
  ! dimension run one after another so threaded LAPACK works on them, see
  ! icgn_outer_threads).  Pair projections of a merge are parallel when the merge
  ! itself runs on the full thread budget (tb).
  ! ============================================================================

  ! Local unit index ul (0 = orbital c lies outside group g) and position lp of
  ! orbital c inside that unit.  lev=1: units = cubes of row g; lev=2: units = rows
  ! of plane g; lev=3: units = planes (g unused).
  pure subroutine icgn_hier_map(hh, lev, g, c, ul, lp)
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g, c
    integer, intent(out) :: ul, lp
    integer :: ic, ir, ip
    ic = hh%cube_of_orb(c)
    ir = (ic - 1) / hh%nx + 1
    ip = (ir - 1) / hh%ny + 1
    ul = 0; lp = 0
    select case (lev)
    case (1)
       if (ir /= g) return
       ul = ic - (g - 1) * hh%nx
       lp = hh%pos(c) - hh%cstart(ic) + 1
    case (2)
       if (ip /= g) return
       ul = ir - (g - 1) * hh%ny
       lp = hh%pos(c) - hh%cstart((ir - 1) * hh%nx + 1) + 1
    case default
       ul = ip
       lp = hh%pos(c) - hh%cstart((ip - 1) * hh%nx * hh%ny + 1) + 1
    end select
  end subroutine icgn_hier_map

  ! Slice [p1,p2] of perm owned by group g of level lev.
  pure subroutine icgn_hier_range(hh, lev, g, p1, p2)
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g
    integer, intent(out) :: p1, p2
    select case (lev)
    case (1)
       p1 = hh%cstart((g - 1) * hh%nx + 1)
       p2 = hh%cstart(g * hh%nx + 1) - 1
    case (2)
       p1 = hh%cstart((g - 1) * hh%nx * hh%ny + 1)
       p2 = hh%cstart(g * hh%nx * hh%ny + 1) - 1
    case default
       p1 = 1
       p2 = size(hh%perm)
    end select
  end subroutine icgn_hier_range

  ! First position in perm of unit u (local index) of group g, level lev.
  pure integer function icgn_hier_ustart(hh, lev, g, u) result(s)
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g, u
    select case (lev)
    case (1)
       s = hh%cstart((g - 1) * hh%nx + u)
    case (2)
       s = hh%cstart(((g - 1) * hh%ny + (u - 1)) * hh%nx + 1)
    case default
       s = hh%cstart((u - 1) * hh%nx * hh%ny + 1)
    end select
  end function icgn_hier_ustart

  ! Retained index range [i1,i2] among eigenvalues w(1:m) sorted ascending:
  ! nominal window [elo,ehi], then all-or-nothing growth of degenerate multiplets
  ! (|dE| <= delta) at both window edges (same rule as icgn_expand_keep_boundary).
  ! i1 > i2 means nothing is retained.
  pure subroutine icgn_window_range(w, m, elo, ehi, delta, i1, i2)
    real(dp), intent(in) :: w(:)
    integer, intent(in) :: m
    real(dp), intent(in) :: elo, ehi, delta
    integer, intent(out) :: i1, i2
    integer :: j
    i1 = m + 1; i2 = 0
    do j = 1, m
       if (w(j) >= elo .and. w(j) <= ehi) then
          if (i1 > m) i1 = j
          i2 = j
       end if
    end do
    if (i2 < i1) then
       i1 = 1; i2 = 0
       return
    end if
    do while (i2 < m)
       if (w(i2 + 1) - w(i2) > delta) exit
       i2 = i2 + 1
    end do
    do while (i1 > 1)
       if (w(i1) - w(i1 - 1) > delta) exit
       i1 = i1 - 1
    end do
  end subroutine icgn_window_range

  ! Nominal keep window [emin - bottom_buffer, emax + top_buffer] of nested tier t.
  pure subroutine icgn_tier_window(upt, t, elo, ehi)
    type(OUPT), intent(in) :: upt
    integer, intent(in) :: t
    real(dp), intent(out) :: elo, ehi
    elo = upt%icgn_tier_emin(t) - upt%icgn_tier_bot(t)
    ehi = upt%icgn_tier_emax(t) + upt%icgn_tier_top(t)
  end subroutine icgn_tier_window

  ! Window diagonalization of a dense Hermitian matrix h (both triangles filled)
  ! with the window of nested tier `tier`.
  ! Same policy as the fast path: ZHEEVR on [window -+ 1e-3 guard] (ZHEEVD fallback),
  ! then keep the nominal core+buffer window plus degenerate multiplets at its edges.
  ! w_out/z_out are allocated (size 0 when nothing is retained) unless ierr /= 0.
  subroutine icgn_window_diag(upt, tier, h, w_out, z_out, m, ierr)
    type(OUPT), intent(in) :: upt
    integer, intent(in) :: tier
    complex(dp), intent(in) :: h(:,:)
    real(dp), allocatable, intent(out) :: w_out(:)
    complex(dp), allocatable, intent(out) :: z_out(:,:)
    integer, intent(out) :: m, ierr
    real(dp), parameter :: margin = 1.0e-3_dp, delta = 1.0e-6_dp
    integer :: nn, mm, i1, i2
    real(dp) :: elo, ehi
    real(dp), allocatable :: w(:)
    complex(dp), allocatable :: a(:,:), z(:,:)

    nn = size(h, 1)
    m = 0; ierr = 0
    if (nn < 1) then
       allocate(w_out(0), z_out(0, 0))
       return
    end if
    call icgn_tier_window(upt, tier, elo, ehi)
    allocate(a(nn, nn), w(nn), z(nn, nn))
    a = h
    call dense_eigh_window(a, elo - margin, ehi + margin, w, z, mm, ierr)
    if (ierr /= 0) then
       a = h
       call coarse_eigh(a, w, ierr)
       if (ierr /= 0) then
          deallocate(a, w, z)
          return
       end if
       z = a
       mm = nn
    end if
    call icgn_window_range(w, mm, elo, ehi, delta, i1, i2)
    m = max(0, i2 - i1 + 1)
    allocate(w_out(m), z_out(nn, m))
    if (m > 0) then
       w_out = w(i1:i2)
       z_out = z(:, i1:i2)
    end if
    deallocate(a, w, z)
  end subroutine icgn_window_diag

  ! ---- retained-fraction statistics per nested tier ---------------------------
  subroutine icgn_stat_reset()
    nst_cnt = 0_8; nst_tot = 0_8; nst_ret = 0_8; nst_fsum = 0.0_dp
  end subroutine icgn_stat_reset

  ! Record one diagonalized unit of tier t: ntot states, nret retained.
  subroutine icgn_stat_add(t, ntot, nret)
    integer, intent(in) :: t, ntot, nret
    if (ntot <= 0 .or. t < 1 .or. t > 3) return
    !$omp critical (icgn_stats)
    nst_cnt(t) = nst_cnt(t) + 1_8
    nst_tot(t) = nst_tot(t) + int(ntot, 8)
    nst_ret(t) = nst_ret(t) + int(nret, 8)
    nst_fsum(t) = nst_fsum(t) + real(nret, dp) / real(ntot, dp)
    !$omp end critical (icgn_stats)
  end subroutine icgn_stat_add

  ! Serial: one log line per tier with the average retained fraction.
  subroutine icgn_stat_log(upt)
    type(OUPT), intent(in) :: upt
    character(len=12), parameter :: names(3) = (/ 'cubes       ', 'rows        ', 'planes      ' /)
    character(len=256) :: msg
    integer :: t
    do t = 1, 3
       if (nst_cnt(t) == 0) then
          write(msg,'(a,i0,a,a)') '(icgn) nested tier ', t, ' ', trim(names(t)) // &
               ': no unit diagonalized'
       else
          write(msg,'(a,i0,1x,a,a,i0,a,f7.3,a,f7.3,a,i0,a,i0,a)') '(icgn) nested tier ', t, &
               trim(names(t)), ': units=', nst_cnt(t), ' retained fraction: mean/unit=', &
               100.0_dp * nst_fsum(t) / real(nst_cnt(t), dp), '%  total=', &
               100.0_dp * real(nst_ret(t), dp) / real(max(1_8, nst_tot(t)), dp), '% (', &
               nst_ret(t), '/', nst_tot(t), ' states)'
       end if
       call icgn_log_important(upt, trim(msg))
    end do
  end subroutine icgn_stat_log

  subroutine icgn_free_unit(u)
    type(ICGNBlock), intent(inout) :: u
    if (associated(u%q)) deallocate(u%q)
    if (associated(u%eval)) deallocate(u%eval)
    if (associated(u%rows)) deallocate(u%rows)
    if (associated(u%S_full)) deallocate(u%S_full)
    if (associated(u%evals_full)) deallocate(u%evals_full)
    if (associated(u%retained_idx)) deallocate(u%retained_idx)
    nullify(u%q, u%eval, u%rows, u%S_full, u%evals_full, u%retained_idx)
    u%nret = 0
  end subroutine icgn_free_unit

  ! Number of states of the full eigensystem of a unit (0 if it has none).
  pure integer function icgn_unit_nfull(u) result(n)
    type(ICGNBlock), intent(in) :: u
    n = 0
    if (associated(u%evals_full)) n = size(u%evals_full)
  end function icgn_unit_nfull

  subroutine icgn_free_units(units)
    type(ICGNBlock), allocatable, intent(inout) :: units(:)
    integer :: i
    if (.not. allocated(units)) return
    do i = 1, size(units)
       call icgn_free_unit(units(i))
    end do
    deallocate(units)
  end subroutine icgn_free_units

  ! V = Qa^H H_ab Qb for the pair (a,b) whose H_ab entries are
  ! e_ord(m1:m2) -> (row e_la in a, column e_lb in b, value e_v).  Only the
  ! interface rows of a / columns of b take part: V = (Qa(R,:))^H (H(R,C) Qb(C,:)).
  subroutine icgn_project_entries(qa, qb, m1, m2, e_ord, e_la, e_lb, e_v, v)
    complex(dp), intent(in) :: qa(:,:), qb(:,:)
    integer, intent(in) :: m1, m2
    integer, intent(in) :: e_ord(:), e_la(:), e_lb(:)
    complex(dp), intent(in) :: e_v(:)
    complex(dp), pointer :: v(:,:)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp), czero = (0.0_dp, 0.0_dp)
    integer :: m, e, nra, ncb
    integer, allocatable :: rmap(:), cmap(:), ridx(:), cidx(:)
    complex(dp), allocatable :: hs(:,:), qas(:,:), qbs(:,:), t(:,:)

    allocate(v(size(qa, 2), size(qb, 2)))
    v = czero
    if (m2 < m1) return
    allocate(rmap(size(qa, 1)), cmap(size(qb, 1)), ridx(m2 - m1 + 1), cidx(m2 - m1 + 1))
    rmap = 0; cmap = 0; nra = 0; ncb = 0
    do m = m1, m2
       e = e_ord(m)
       if (rmap(e_la(e)) == 0) then
          nra = nra + 1; rmap(e_la(e)) = nra; ridx(nra) = e_la(e)
       end if
       if (cmap(e_lb(e)) == 0) then
          ncb = ncb + 1; cmap(e_lb(e)) = ncb; cidx(ncb) = e_lb(e)
       end if
    end do
    allocate(hs(nra, ncb))
    hs = czero
    do m = m1, m2
       e = e_ord(m)
       hs(rmap(e_la(e)), cmap(e_lb(e))) = hs(rmap(e_la(e)), cmap(e_lb(e))) + e_v(e)
    end do
    allocate(qas(nra, size(qa, 2)), qbs(ncb, size(qb, 2)), t(nra, size(qb, 2)))
    qas = qa(ridx(1:nra), :)
    qbs = qb(cidx(1:ncb), :)
    t = czero
    call icgn_zgemm_gen('N', 'N', cone, hs, qbs, czero, t)     ! (nra x ncb)(ncb x nb)
    call icgn_zgemm_gen('C', 'N', cone, qas, t, czero, v)      ! (na x nra)(nra x nb)
  end subroutine icgn_project_entries

  ! Make sure the columns of every CSR row are strictly ascending (cheap check
  ! first; only unsorted rows are insertion-sorted).
  subroutine icgn_csr_sort_raw(mi, mj, m, n)
    integer, intent(in) :: n
    integer, intent(in) :: mi(:)
    integer, intent(inout) :: mj(:)
    complex(dp), intent(inout) :: m(:)
    integer :: r, a, b, j, k, tmpj
    complex(dp) :: tmpv
    logical :: ok
    do r = 1, n
       a = mi(r); b = mi(r + 1) - 1
       ok = .true.
       do j = a + 1, b
          if (mj(j) <= mj(j - 1)) then
             ok = .false.; exit
          end if
       end do
       if (ok) cycle
       do j = a + 1, b
          tmpj = mj(j); tmpv = m(j); k = j - 1
          do while (k >= a)
             if (mj(k) <= tmpj) exit
             mj(k + 1) = mj(k); m(k + 1) = m(k); k = k - 1
          end do
          mj(k + 1) = tmpj; m(k + 1) = tmpv
       end do
    end do
  end subroutine icgn_csr_sort_raw


  ! ---------------------------------------------------------------------------
  ! Generic H_CG assembly of ONE group (no diagonalization).
  !   hin          : operator to project (any CSR; here always the orbital H)
  !   hh, lev, g   : unit/row maps of the group (see icgn_hier_map)
  !   units(1:nu)  : retained basis (q, eval, nret, nrow) of every unit
  !   tb           : thread budget of this call (1 inside an outer parallel loop)
  !   apply_filter : drop couplings with icgn_coupling_survives (final level only)
  ! Output (raw CSR, format 'F', diagonal present, columns ascending):
  !   rows/cols = retained states of the units concatenated in unit order.
  ! Diagonal blocks are diag(eval): valid because every unit basis diagonalizes
  ! H restricted to the unit (true at every level of the nested scheme).
  ! Parallel over unit pairs; only entries inside the group are scanned.
  ! ---------------------------------------------------------------------------
  subroutine icgn_cg_assemble(upt, hin, hh, lev, g, units, nu, tb, apply_filter, &
       mi, mj, mv, nnz, ierr, tier)
    type(OUPT), intent(in) :: upt
    type(CSR), intent(in) :: hin
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g, nu, tb
    type(ICGNBlock), intent(in) :: units(:)
    logical, intent(in) :: apply_filter
    integer, allocatable, intent(out) :: mi(:), mj(:)
    complex(dp), allocatable, intent(out) :: mv(:)
    integer, intent(out) :: nnz, ierr
    integer, intent(in), optional :: tier    ! nested tier whose window / epsilon apply
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    integer :: i, j, k, ia, ib, npair, nth, nn, u
    integer, allocatable :: pa(:), pb(:), slot_of(:,:)
    type(ICGNMat), allocatable :: vp(:)
    real(dp) :: elo, ehi, eps_t

    ierr = 0; nnz = 0
    nn = 0
    do u = 1, nu
       nn = nn + units(u)%nret
    end do
    if (nn < 1) then
       ierr = 34; write(*,*) '(icgn) assemble: no retained states in group'
       return
    end if
    if (present(tier)) then
       call icgn_tier_window(upt, tier, elo, ehi)
       eps_t = upt%icgn_tier_eps(tier)
    else
       elo = upt%icgn_core_emin - upt%icgn_bottom_buffer
       ehi = upt%icgn_core_emax + upt%icgn_top_buffer
       eps_t = upt%icgn_epsilon
    end if
    call icgn_group_pairs(hin, hh, lev, g, units, nu, tb, .false., pa, pb, npair, vp, slot_of, ierr)
    if (ierr /= 0) return
    if (apply_filter) then
       nth = max(1, min(tb, npair))
       !$omp parallel do default(shared) private(i, j, k, ia, ib) num_threads(nth) &
       !$omp& schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          ia = pa(i); ib = pb(i)
          do k = 1, units(ib)%nret
             do j = 1, units(ia)%nret
                if (.not. icgn_coupling_survives(vp(i)%m(j, k), units(ia)%eval(j), &
                     units(ib)%eval(k), eps_t, elo, ehi)) vp(i)%m(j, k) = czero
             end do
          end do
       end do
       !$omp end parallel do
    end if
    call icgn_csr_from_pairs(units, nu, slot_of, vp, npair, mi, mj, mv, nnz, ierr)
    do i = 1, npair
       if (associated(vp(i)%m)) deallocate(vp(i)%m)
    end do
    deallocate(vp, pa, pb, slot_of)
  end subroutine icgn_cg_assemble

  ! Pair couplings inside one group.  use_full = .false.: V = Q_a^H H_ab Q_b on
  ! the retained bases (q);  use_full = .true.: G = S_a^H H_ab S_b on the full
  ! eigenbases (S_full).  Only pairs of units that have a basis and at least one
  ! stored H entry are produced.  Parallel over pairs (budget tb).
  !   pa/pb(1:npair): unit indices (pa < pb);  vp(i)%m: dense block;
  !   slot_of(a,b), a<b: pair index or 0.
  subroutine icgn_group_pairs(hin, hh, lev, g, units, nu, tb, use_full, pa, pb, npair, vp, &
       slot_of, ierr)
    type(CSR), intent(in) :: hin
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g, nu, tb
    type(ICGNBlock), intent(in) :: units(:)
    logical, intent(in) :: use_full
    integer, allocatable, intent(out) :: pa(:), pb(:), slot_of(:,:)
    integer, intent(out) :: npair, ierr
    type(ICGNMat), allocatable, intent(out) :: vp(:)

    integer :: u, i, k, p, p1, p2, r, c, ne, e, nth, ur, uc, lr, lc, ia, ib
    logical, allocatable :: active(:)
    integer, allocatable :: e_la(:), e_lb(:), e_slot(:), p_start(:), p_fill(:), e_ord(:)
    complex(dp), allocatable :: e_v(:)

    ierr = 0
    allocate(active(nu))
    do u = 1, nu
       if (use_full) then
          active(u) = associated(units(u)%S_full)
          if (active(u)) active(u) = size(units(u)%S_full, 2) > 0
       else
          active(u) = units(u)%nret > 0
       end if
    end do
    call icgn_hier_range(hh, lev, g, p1, p2)

    ! ---- pass 1: count inter-unit entries (upper-triangle convention) ------
    ne = 0
    do p = p1, p2
       r = hh%perm(p)
       call icgn_hier_map(hh, lev, g, r, ur, lr)
       if (ur < 1) cycle
       if (.not. active(ur)) cycle
       do k = hin%Mi(r), hin%Mi(r + 1) - 1
          c = hin%Mj(k)
          call icgn_hier_map(hh, lev, g, c, uc, lc)
          if (uc < 1 .or. uc == ur) cycle
          if (.not. active(uc)) cycle
          if (.not. stored_entry(hin%sparse_fmt, r, c)) cycle
          ne = ne + 1
       end do
    end do

    ! ---- pass 2: store entries normalised to H_ab(la,lb), a < b ------------
    allocate(e_la(max(1, ne)), e_lb(max(1, ne)), e_slot(max(1, ne)), e_v(max(1, ne)))
    allocate(slot_of(nu, nu), pa(max(1, nu * (nu - 1) / 2)), pb(max(1, nu * (nu - 1) / 2)))
    slot_of = 0; npair = 0; e = 0
    do p = p1, p2
       r = hh%perm(p)
       call icgn_hier_map(hh, lev, g, r, ur, lr)
       if (ur < 1) cycle
       if (.not. active(ur)) cycle
       do k = hin%Mi(r), hin%Mi(r + 1) - 1
          c = hin%Mj(k)
          call icgn_hier_map(hh, lev, g, c, uc, lc)
          if (uc < 1 .or. uc == ur) cycle
          if (.not. active(uc)) cycle
          if (.not. stored_entry(hin%sparse_fmt, r, c)) cycle
          e = e + 1
          if (ur < uc) then
             ia = ur; ib = uc
             e_la(e) = lr; e_lb(e) = lc; e_v(e) = hin%M(k)
          else
             ia = uc; ib = ur
             e_la(e) = lc; e_lb(e) = lr; e_v(e) = conjg(hin%M(k))
          end if
          if (slot_of(ia, ib) == 0) then
             npair = npair + 1
             slot_of(ia, ib) = npair; pa(npair) = ia; pb(npair) = ib
          end if
          e_slot(e) = slot_of(ia, ib)
       end do
    end do

    ! ---- counting sort of the entries by pair ------------------------------
    allocate(p_start(npair + 1), p_fill(max(1, npair)), e_ord(max(1, ne)))
    p_start = 0
    do e = 1, ne
       p_start(e_slot(e) + 1) = p_start(e_slot(e) + 1) + 1
    end do
    do i = 1, npair
       p_start(i + 1) = p_start(i + 1) + p_start(i)
    end do
    p_fill(1:npair) = p_start(1:npair)
    do e = 1, ne
       i = e_slot(e)
       p_fill(i) = p_fill(i) + 1
       e_ord(p_fill(i)) = e
    end do
    deallocate(e_slot, p_fill)

    ! ---- pair projections (each pair owns its block: no locking) -----------
    allocate(vp(max(1, npair)))
    nth = max(1, min(tb, npair))
    !$omp parallel do default(shared) private(i, ia, ib) num_threads(nth) &
    !$omp& schedule(dynamic,1) if(nth > 1)
    do i = 1, npair
       ia = pa(i); ib = pb(i)
       if (use_full) then
          call icgn_project_entries(units(ia)%S_full, units(ib)%S_full, p_start(i) + 1, &
               p_start(i + 1), e_ord, e_la, e_lb, e_v, vp(i)%m)
       else
          call icgn_project_entries(units(ia)%q, units(ib)%q, p_start(i) + 1, &
               p_start(i + 1), e_ord, e_la, e_lb, e_v, vp(i)%m)
       end if
    end do
    !$omp end parallel do
    deallocate(e_la, e_lb, e_v, e_ord, p_start, active)
  end subroutine icgn_group_pairs

  ! CSR (format 'F', diagonal = eval of the retained states, both triangles of
  ! every non-zero pair block, columns ascending) from the pair blocks vp, whose
  ! sizes must be nret_a x nret_b.  Emission order (unit u: partners w<u, diagonal,
  ! partners v>u) already yields ascending columns.
  subroutine icgn_csr_from_pairs(units, nu, slot_of, vp, npair, mi, mj, mv, nnz, ierr)
    type(ICGNBlock), intent(in) :: units(:)
    integer, intent(in) :: nu, npair
    integer, intent(in) :: slot_of(:,:)
    type(ICGNMat), intent(in) :: vp(:)
    integer, allocatable, intent(out) :: mi(:), mj(:)
    complex(dp), allocatable, intent(out) :: mv(:)
    integer, intent(out) :: nnz, ierr
    integer :: u, v, i, j, k, r, nn, pos
    integer, allocatable :: coff(:), rowcount(:), next(:)
    integer(8) :: nnz8

    ierr = 0; nnz = 0
    allocate(coff(nu + 1))
    coff(1) = 1
    do u = 1, nu
       coff(u + 1) = coff(u) + units(u)%nret
    end do
    nn = coff(nu + 1) - 1
    allocate(rowcount(max(1, nn)))
    rowcount = 1
    do u = 1, nu
       do v = u + 1, nu
          i = slot_of(u, v)
          if (i == 0) cycle
          do k = 1, size(vp(i)%m, 2)
             do j = 1, size(vp(i)%m, 1)
                if (abs(vp(i)%m(j, k)) == 0.0_dp) cycle
                rowcount(coff(u) + j - 1) = rowcount(coff(u) + j - 1) + 1
                rowcount(coff(v) + k - 1) = rowcount(coff(v) + k - 1) + 1
             end do
          end do
       end do
    end do
    nnz8 = 0_8
    do r = 1, nn
       nnz8 = nnz8 + int(rowcount(r), kind=8)
    end do
    if (nnz8 > int(huge(nnz), kind=8)) then
       write(*,'(a,i0,a)') '(icgn) ERROR: nested H_CG nnz=', nnz8, &
            ' exceeds default INTEGER capacity; reduce the window or use more cubes'
       ierr = 12
       return
    end if
    nnz = int(nnz8)
    allocate(mi(nn + 1), mj(max(1, nnz)), mv(max(1, nnz)), next(max(1, nn)))
    mi(1) = 1
    do r = 1, nn
       mi(r + 1) = mi(r) + rowcount(r)
    end do
    if (nn > 0) next(1:nn) = mi(1:nn)
    do u = 1, nu
       do j = 1, units(u)%nret
          r = coff(u) + j - 1
          pos = next(r)
          mj(pos) = r
          mv(pos) = cmplx(units(u)%eval(j), 0.0_dp, kind=dp)
          next(r) = pos + 1
       end do
       do v = u + 1, nu
          i = slot_of(u, v)
          if (i == 0) cycle
          do j = 1, size(vp(i)%m, 1)
             do k = 1, size(vp(i)%m, 2)
                if (abs(vp(i)%m(j, k)) == 0.0_dp) cycle
                r = coff(u) + j - 1
                pos = next(r)
                mj(pos) = coff(v) + k - 1
                mv(pos) = vp(i)%m(j, k)
                next(r) = pos + 1
                r = coff(v) + k - 1
                pos = next(r)
                mj(pos) = coff(u) + j - 1
                mv(pos) = conjg(vp(i)%m(j, k))
                next(r) = pos + 1
             end do
          end do
       end do
    end do
    call icgn_csr_sort_raw(mi, mj, mv, nn)
    deallocate(rowcount, next, coff)
  end subroutine icgn_csr_from_pairs

  ! Level 0: diagonalize cube ic (dense H of its own orbitals) and keep the
  ! window states.  u%q rows follow the cube's slice of perm.
  subroutine icgn_nested_cube(upt, hh, ic, u, ierr)
    type(OUPT), intent(in) :: upt
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: ic
    type(ICGNBlock), intent(out) :: u
    integer, intent(out) :: ierr
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    integer :: cs, nn, jr, jc, r, c, k, m
    real(dp), allocatable :: w(:)
    complex(dp), allocatable :: z(:,:)
    complex(dp), pointer :: h(:,:)
    logical :: full

    ierr = 0
    cs = hh%cstart(ic)
    nn = hh%cstart(ic + 1) - cs
    u%nrow = nn; u%nret = 0
    if (nn == 0) return
    full = icgn_tier_full(upt, 1)
    if (full) then
       ! Outers / Neumann at the cube tier: the full eigensystem is kept for the
       ! row merge (which selects P / Q).  H is assembled directly in S_full and
       ! diagonalized in place (no second nn x nn copy).
       allocate(u%S_full(nn, nn))
       h => u%S_full
    else
       allocate(h(nn, nn))
    end if
    h = czero
    do jr = 1, nn
       r = hh%perm(cs + jr - 1)
       do k = upt%ham%Mi(r), upt%ham%Mi(r + 1) - 1
          c = upt%ham%Mj(k)
          if (hh%cube_of_orb(c) /= ic) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          jc = hh%pos(c) - cs + 1
          h(jr, jc) = upt%ham%M(k)
          if (r /= c) h(jc, jr) = conjg(upt%ham%M(k))
       end do
    end do
    if (full) then
       allocate(u%evals_full(nn))
       call coarse_eigh(h, u%evals_full, ierr)
       if (ierr /= 0) then
          deallocate(u%S_full, u%evals_full)
          nullify(u%S_full, u%evals_full)
       end if
       return
    end if
    call icgn_window_diag(upt, 1, h, w, z, m, ierr)
    deallocate(h)
    if (ierr /= 0) return
    u%nret = m
    call icgn_stat_add(1, nn, m)
    if (m > 0) then
       allocate(u%eval(m), u%q(nn, m))
       u%eval = w
       u%q = z(:, 1:m)
    end if
    deallocate(w, z)
  end subroutine icgn_nested_cube

  ! .true. when diagonalization tier t (1 = cubes, 2 = rows, 3 = planes) keeps
  ! its full spectrum for core-state outers and/or the Neumann self-energy.
  logical function icgn_tier_full(upt, t)
    type(OUPT), intent(in) :: upt
    integer, intent(in) :: t
    icgn_tier_full = (upt%icgn_nested_outers(t) > 0 .or. upt%icgn_nested_orders(t) >= 0)
  end function icgn_tier_full

  ! Merge the units of group g (level lev) into ONE unit.
  !   child tier (lev) window-only : H_CG = icgn_cg_assemble on the retained bases;
  !   child tier (lev) full        : icgn_full_group (outers, H_PP, Neumann Sigma)
  !                                  -> H_eff = H_PP + Sigma on the P bases.
  ! H_eff is diagonalized (window states only, or the FULL spectrum when the tier
  ! of the new unit, lev+1, needs outers / Neumann) and
  !     Q_group = blockdiag(Q_u^P) * Z
  ! is expressed on the group's own orbitals (perm slice).  The input units are
  ! freed here.
  subroutine icgn_nested_merge(upt, hh, lev, g, units, nu, tb, outu, ierr)
    type(OUPT), intent(in) :: upt
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g, nu, tb
    type(ICGNBlock), intent(inout) :: units(:)
    type(ICGNBlock), intent(out) :: outu
    integer, intent(out) :: ierr
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    integer :: u, i, k, p1, p2, nrow_g, nn, m, nact, ro, c1, c2, nnz, ncols
    integer, allocatable :: mi(:), mj(:), coff(:), roff(:)
    complex(dp), allocatable :: mv(:), h(:,:), z(:,:), tmp(:,:)
    real(dp), allocatable :: w(:)
    type(ICGNSigRow), allocatable :: sig(:)
    complex(dp), pointer :: qt(:,:)
    logical :: child_full, self_full

    ierr = 0
    child_full = icgn_tier_full(upt, lev)
    self_full = (lev < 3) .and. icgn_tier_full(upt, min(3, lev + 1))
    call icgn_hier_range(hh, lev, g, p1, p2)
    nrow_g = p2 - p1 + 1
    outu%nrow = nrow_g; outu%nret = 0
    nact = 0
    do u = 1, nu
       if (child_full) then
          if (icgn_unit_nfull(units(u)) > 0) nact = nact + 1
       else
          if (units(u)%nret > 0) nact = nact + 1
       end if
    end do
    if (nact == 0) go to 900
    if (nu == 1 .and. .not. child_full) then
       ! A single window-only unit: its basis already is the group basis and
       ! H_CG is diagonal in it, so no assembly / diagonalization is needed.
       if (self_full) then
          ! next tier wants the "full" spectrum: that is exactly the retained states
          outu%evals_full => units(1)%eval
          outu%S_full => units(1)%q
       else
          call icgn_stat_add(lev + 1, units(1)%nret, units(1)%nret)
          outu%nret = units(1)%nret
          outu%eval => units(1)%eval
          outu%q => units(1)%q
       end if
       nullify(units(1)%eval, units(1)%q)
       units(1)%nret = 0
       go to 900
    end if

    if (child_full) then
       call icgn_full_group(upt, hh, lev, g, units, nu, tb, mi, mj, mv, nnz, sig, roff, ierr)
    else
       call icgn_cg_assemble(upt, upt%ham, hh, lev, g, units, nu, tb, .false., &
            mi, mj, mv, nnz, ierr)
    end if
    if (ierr /= 0) go to 900
    nn = size(mi) - 1
    if (nn < 1) go to 900
    allocate(h(nn, nn))
    h = czero
    do i = 1, nn
       do k = mi(i), mi(i + 1) - 1
          h(i, mj(k)) = mv(k)
       end do
    end do
    deallocate(mi, mj, mv)
    if (allocated(sig)) then
       call icgn_dense_add_sigma(h, roff, sig)
       call icgn_free_sigrows(sig)
    end if

    if (self_full) then
       allocate(w(nn))
       call coarse_eigh(h, w, ierr)       ! h <- eigenvectors (all nn states)
       if (ierr /= 0) go to 900
       call move_alloc(h, z)
       ncols = nn
    else
       call icgn_window_diag(upt, lev + 1, h, w, z, m, ierr)
       deallocate(h)
       if (ierr /= 0) go to 900
       call icgn_stat_add(lev + 1, nn, m)
       ncols = m
    end if
    if (ncols == 0) go to 900

    if (self_full) then
       allocate(outu%S_full(nrow_g, ncols), outu%evals_full(ncols))
       outu%evals_full = w(1:ncols)
       qt => outu%S_full
    else
       allocate(outu%q(nrow_g, ncols), outu%eval(ncols))
       outu%eval = w(1:ncols)
       outu%nret = ncols
       qt => outu%q
    end if
    qt = czero
    allocate(coff(nu + 1))
    coff(1) = 1
    do u = 1, nu
       coff(u + 1) = coff(u) + units(u)%nret
    end do
    do u = 1, nu
       if (units(u)%nret == 0) cycle
       ro = icgn_hier_ustart(hh, lev, g, u) - p1
       c1 = coff(u); c2 = coff(u + 1) - 1
       allocate(tmp(units(u)%nrow, ncols))
       call icgn_zgemm('N', units(u)%q, z(c1:c2, :), tmp)
       qt(ro + 1:ro + units(u)%nrow, :) = tmp
       deallocate(tmp)
    end do
    deallocate(coff)

900 continue
    do u = 1, nu
       call icgn_free_unit(units(u))
    end do
    if (allocated(sig)) call icgn_free_sigrows(sig)
  end subroutine icgn_nested_merge

  ! Dense H_eff += Sigma (block-sparse upper blocks a2 >= a; the lower part follows
  ! from Hermiticity, diagonal blocks are Hermitianized).
  subroutine icgn_dense_add_sigma(h, roff, sig)
    complex(dp), intent(inout) :: h(:,:)
    integer, intent(in) :: roff(:)
    type(ICGNSigRow), intent(in) :: sig(:)
    integer :: a, t, a2, i, j
    do a = 1, size(sig)
       do t = 1, sig(a)%n
          a2 = sig(a)%blk(t)
          if (a2 == a) then
             do j = 1, size(sig(a)%mat(t)%m, 2)
                do i = 1, size(sig(a)%mat(t)%m, 1)
                   h(roff(a) + i - 1, roff(a) + j - 1) = h(roff(a) + i - 1, roff(a) + j - 1) &
                        + 0.5_dp * (sig(a)%mat(t)%m(i, j) + conjg(sig(a)%mat(t)%m(j, i)))
                end do
             end do
          else
             do j = 1, size(sig(a)%mat(t)%m, 2)
                do i = 1, size(sig(a)%mat(t)%m, 1)
                   h(roff(a) + i - 1, roff(a2) + j - 1) = h(roff(a) + i - 1, roff(a2) + j - 1) &
                        + sig(a)%mat(t)%m(i, j)
                   h(roff(a2) + j - 1, roff(a) + i - 1) = h(roff(a2) + j - 1, roff(a) + i - 1) &
                        + conjg(sig(a)%mat(t)%m(i, j))
                end do
             end do
          end if
       end do
    end do
  end subroutine icgn_dense_add_sigma

  ! Full-spectrum treatment of ONE group (the nested counterpart of the flat
  ! ICGN preparation), used when the tier of the units (tier = lev) has
  ! add_core_outers > 0 and/or neumann_order >= 0.  Every unit carries its full
  ! eigensystem (S_full on its own orbitals, evals_full).  Steps:
  !   1. nominal window -> keep / core flags (+ degenerate multiplets at the edges)
  !   2. full pair couplings G_ab = S_a^H H_ab S_b            (parallel over pairs)
  !   3. level-1 core-state outers: discarded states strongly coupled to core states
  !   4. retained bases (q, eval, retained_idx), S_full freed
  !   5. H_PP (sliced and filtered G) -> CSR
  !   6. Neumann Sigma = H_PQ R [W R]^n H_QP inside the group (block-sparse)
  ! Sigma stays local to the group: couplings through discarded states of OTHER
  ! groups are only seen at the next tier (through the retained states).
  ! On return: units(u)%q/eval/nret/retained_idx are the P basis, evals_full is still
  ! allocated (freed by the caller with the unit).  sig is allocated only if
  ! Sigma was computed; roff = offsets of the units in the P basis.
  subroutine icgn_full_group(upt, hh, lev, g, units, nu, tb, mi, mj, mv, nnz, sig, roff, ierr)
    type(OUPT), intent(in) :: upt
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: lev, g, nu, tb
    type(ICGNBlock), intent(inout) :: units(:)
    integer, allocatable, intent(out) :: mi(:), mj(:), roff(:)
    complex(dp), allocatable, intent(out) :: mv(:)
    integer, intent(out) :: nnz, ierr
    type(ICGNSigRow), allocatable, intent(out) :: sig(:)

    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    integer :: tier, outer_t, order_t, u, i, j, k, ia, ib, npair, nth, nn, maxnf, i1, i2
    integer, allocatable :: pa(:), pb(:), slot_of(:,:), local_ret_idx(:,:), nfu(:)
    integer, allocatable :: nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, allocatable :: keep(:,:), core(:,:), nbr_isA(:)
    type(ICGNMat), allocatable :: vp(:), vpp(:)
    type(ICGNPair), allocatable :: prs(:)
    type(ICGNQBlock), allocatable :: bq(:)
    type(ICGNQPair), allocatable :: qp(:)
    integer(8) :: npq, nqq
    real(dp) :: elo, ehi, E0, emin_t, emax_t, eps_t
    complex(dp) :: vv

    ierr = 0; nnz = 0
    tier = min(max(lev, 1), 3)
    outer_t = upt%icgn_nested_outers(tier)
    order_t = upt%icgn_nested_orders(tier)
    call icgn_tier_window(upt, tier, elo, ehi)
    emin_t = upt%icgn_tier_emin(tier)
    emax_t = upt%icgn_tier_emax(tier)
    eps_t = upt%icgn_tier_eps(tier)

    allocate(nfu(nu))
    do u = 1, nu
       nfu(u) = icgn_unit_nfull(units(u))
    end do
    maxnf = max(1, maxval(nfu))
    allocate(keep(maxnf, nu), core(maxnf, nu))
    keep = .false.; core = .false.
    do u = 1, nu
       if (nfu(u) == 0) cycle
       call icgn_window_range(units(u)%evals_full, nfu(u), elo, ehi, 1.0e-6_dp, i1, i2)
       if (i2 >= i1) keep(i1:i2, u) = .true.
       do j = 1, nfu(u)
          if (keep(j, u) .and. units(u)%evals_full(j) >= emin_t .and. &
               units(u)%evals_full(j) <= emax_t) core(j, u) = .true.
       end do
    end do

    ! ---- full pair couplings --------------------------------------------------
    call icgn_group_pairs(upt%ham, hh, lev, g, units, nu, tb, .true., pa, pb, npair, vp, &
         slot_of, ierr)
    if (ierr /= 0) return

    ! ---- level-1 core-state outers (serial: keep is updated while scanning) ---
    if (outer_t > 0) then
       do i = 1, npair
          ia = pa(i); ib = pb(i)
          do j = 1, nfu(ia)
             if (.not. core(j, ia)) cycle
             do k = 1, nfu(ib)
                if (keep(k, ib)) cycle
                if (icgn_outer_passes(vp(i)%m(j, k), units(ia)%evals_full(j), &
                     units(ib)%evals_full(k), eps_t)) keep(k, ib) = .true.
             end do
          end do
          do k = 1, nfu(ib)
             if (.not. core(k, ib)) cycle
             do j = 1, nfu(ia)
                if (keep(j, ia)) cycle
                if (icgn_outer_passes(vp(i)%m(j, k), units(ib)%evals_full(k), &
                     units(ia)%evals_full(j), eps_t)) keep(j, ia) = .true.
             end do
          end do
       end do
    end if

    ! ---- retained bases ---------------------------------------------------------
    nn = 0
    do u = 1, nu
       units(u)%nret = 0
       if (nfu(u) == 0) cycle
       units(u)%nret = count(keep(1:nfu(u), u))
       call icgn_stat_add(tier, nfu(u), units(u)%nret)
       nn = nn + units(u)%nret
       if (units(u)%nret == 0) cycle
       allocate(units(u)%retained_idx(units(u)%nret), units(u)%eval(units(u)%nret))
       allocate(units(u)%q(units(u)%nrow, units(u)%nret))
       k = 0
       do j = 1, nfu(u)
          if (.not. keep(j, u)) cycle
          k = k + 1
          units(u)%retained_idx(k) = j
          units(u)%eval(k) = units(u)%evals_full(j)
          units(u)%q(:, k) = units(u)%S_full(:, j)
       end do
    end do
    do u = 1, nu
       if (associated(units(u)%S_full)) deallocate(units(u)%S_full)
       nullify(units(u)%S_full)
    end do
    if (nn == 0) then
       do i = 1, npair
          if (associated(vp(i)%m)) deallocate(vp(i)%m)
       end do
       allocate(mi(1), mj(1), mv(1))
       mi(1) = 1
       return
    end if

    ! ---- H_PP: retained slice of G, filtered, -> CSR ----------------------------
    allocate(vpp(max(1, npair)))
    nth = max(1, min(tb, npair))
    !$omp parallel do default(shared) private(i, j, k, ia, ib, vv) num_threads(nth) &
    !$omp& schedule(dynamic,1) if(nth > 1)
    do i = 1, npair
       ia = pa(i); ib = pb(i)
       allocate(vpp(i)%m(units(ia)%nret, units(ib)%nret))
       do k = 1, units(ib)%nret
          do j = 1, units(ia)%nret
             vv = vp(i)%m(units(ia)%retained_idx(j), units(ib)%retained_idx(k))
             if (icgn_coupling_survives(vv, units(ia)%eval(j), units(ib)%eval(k), &
                  eps_t, elo, ehi)) then
                vpp(i)%m(j, k) = vv
             else
                vpp(i)%m(j, k) = czero
             end if
          end do
       end do
       if (order_t < 0) deallocate(vp(i)%m)
    end do
    !$omp end parallel do
    call icgn_csr_from_pairs(units, nu, slot_of, vpp, npair, mi, mj, mv, nnz, ierr)
    do i = 1, npair
       if (associated(vpp(i)%m)) deallocate(vpp(i)%m)
    end do
    deallocate(vpp)
    if (ierr /= 0) then
       do i = 1, npair
          if (associated(vp(i)%m)) deallocate(vp(i)%m)
       end do
       return
    end if

    ! ---- Neumann self-energy inside the group -----------------------------------
    if (order_t >= 0) then
       allocate(prs(max(1, npair)))
       do i = 1, npair
          prs(i)%a = pa(i); prs(i)%b = pb(i)
          prs(i)%g => vp(i)%m
       end do
       allocate(local_ret_idx(nu, maxnf))
       local_ret_idx = 0
       do u = 1, nu
          do j = 1, units(u)%nret
             local_ret_idx(u, units(u)%retained_idx(j)) = j
          end do
       end do
       allocate(roff(nu + 1))
       roff(1) = 1
       do u = 1, nu
          roff(u + 1) = roff(u) + units(u)%nret
       end do
       E0 = upt%icgn_tier_E0(tier)
       if (E0 <= emin_t .or. E0 >= emax_t) E0 = 0.5_dp * (emin_t + emax_t)
       call icgn_build_neumann_blocks(upt, units, nu, prs, npair, local_ret_idx, E0, &
            order_t >= 1, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, npq, nqq, tb, &
            elo, ehi, eps_t)
       do i = 1, npair
          if (associated(vp(i)%m)) deallocate(vp(i)%m)
          nullify(prs(i)%g)
       end do
       deallocate(prs, local_ret_idx)
       allocate(sig(nu))
       if (npq > 0) call icgn_sigma_blocksparse(nu, order_t, maxnf, bq, qp, nbr_ptr, nbr_blk, &
            nbr_pair, nbr_isA, nbr_rev, sig, tb)
       call icgn_free_neumann_blocks(bq, qp, npair)
       deallocate(nbr_ptr, nbr_blk, nbr_pair, nbr_rev, nbr_isA)
    end if
    deallocate(vp, pa, pb, slot_of, keep, core, nfu)
  end subroutine icgn_full_group

  ! Levels 0+1 for one row: diagonalize its cubes (parallel over cubes when the
  ! row runs on tb > 1 threads), then merge them into the row unit.
  subroutine icgn_nested_row(upt, hh, ir, tb, ru, ierr)
    type(OUPT), intent(in) :: upt
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: ir, tb
    type(ICGNBlock), intent(out) :: ru
    integer, intent(out) :: ierr
    integer :: nx, iu, ic1, ic, nth, cmax
    type(ICGNBlock), allocatable :: cub(:)
    integer, allocatable :: ierr_c(:)

    nx = hh%nx
    ic1 = (ir - 1) * nx + 1
    allocate(cub(nx), ierr_c(nx))
    ierr_c = 0
    cmax = 0
    do iu = 1, nx
       ic = ic1 + iu - 1
       cmax = max(cmax, hh%cstart(ic + 1) - hh%cstart(ic))
    end do
    nth = 1
    if (tb > 1) nth = icgn_outer_threads(nx, min(tb, nx), cmax)
    !$omp parallel do default(shared) private(iu) num_threads(nth) schedule(dynamic,1) if(nth > 1)
    do iu = 1, nx
       call icgn_nested_cube(upt, hh, ic1 + iu - 1, cub(iu), ierr_c(iu))
    end do
    !$omp end parallel do
    ierr = 0
    do iu = 1, nx
       if (ierr_c(iu) /= 0) ierr = ierr_c(iu)
    end do
    if (ierr == 0) then
       call icgn_nested_merge(upt, hh, 1, ir, cub, nx, tb, ru, ierr)
    else
       ru%nrow = hh%cstart(ic1 + nx) - hh%cstart(ic1)
    end if
    do iu = 1, nx
       call icgn_free_unit(cub(iu))
    end do
    deallocate(cub, ierr_c)
  end subroutine icgn_nested_row

  ! Levels 0-2 for one plane: build its rows (parallel over rows when
  ! ny >= tb, otherwise rows run one after another and parallelism goes to the
  ! cubes of each row), then merge the rows into the plane unit.
  subroutine icgn_nested_plane(upt, hh, ip, tb, pu, ierr)
    type(OUPT), intent(in) :: upt
    type(ICGNHier), intent(in) :: hh
    integer, intent(in) :: ip, tb
    type(ICGNBlock), intent(out) :: pu
    integer, intent(out) :: ierr
    integer :: nx, ny, iu, ir1, ic, nth, cmax
    type(ICGNBlock), allocatable :: rws(:)
    integer, allocatable :: ierr_r(:)

    nx = hh%nx; ny = hh%ny
    ir1 = (ip - 1) * ny + 1
    allocate(rws(ny), ierr_r(ny))
    ierr_r = 0
    cmax = 0
    do ic = (ip - 1) * nx * ny + 1, ip * nx * ny
       cmax = max(cmax, hh%cstart(ic + 1) - hh%cstart(ic))
    end do
    nth = 1
    if (tb > 1 .and. ny >= tb) nth = icgn_outer_threads(ny, tb, cmax)
    if (nth > 1) then
       !$omp parallel do default(shared) private(iu) num_threads(nth) schedule(dynamic,1)
       do iu = 1, ny
          call icgn_nested_row(upt, hh, ir1 + iu - 1, 1, rws(iu), ierr_r(iu))
       end do
       !$omp end parallel do
    else
       do iu = 1, ny
          call icgn_nested_row(upt, hh, ir1 + iu - 1, tb, rws(iu), ierr_r(iu))
          if (ierr_r(iu) /= 0) exit
       end do
    end if
    ierr = 0
    do iu = 1, ny
       if (ierr_r(iu) /= 0) ierr = ierr_r(iu)
    end do
    if (ierr == 0) then
       call icgn_nested_merge(upt, hh, 2, ip, rws, ny, tb, pu, ierr)
    else
       pu%nrow = hh%cstart(ip * nx * ny + 1) - hh%cstart((ip - 1) * nx * ny + 1)
    end if
    call icgn_free_units(rws)
    deallocate(ierr_r)
  end subroutine icgn_nested_plane

  ! Driver of the nested CG: cubes -> rows -> planes -> bulk H_CG.
  ! Fills upt%icgn_ham (H_CG of the plane units, never diagonalized here) and
  ! upt%icgn_blocks(1:nplanes) (rows, q, eval of every plane: the accumulated Q,
  ! block-diagonal over planes) so uptight's solve/lift/U_red need no changes.
  subroutine icgn_nested_run(upt, na, n_orb, offsets, bsize, atom_of, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: na, n_orb
    integer, intent(in) :: offsets(:), bsize(:), atom_of(:)
    integer, intent(out) :: ierr
    integer :: axes(3), grids(3), nx, ny, nz, nc, np
    integer :: i, g, a, o, ic, ip, nth_p, nred, nnz, r, k, c, p1, p2, T
    integer :: status, save_dir, save_nb, cmax, n_empty, pr
    integer, allocatable :: slab(:,:), part(:), cnt(:), cursor(:), ierr_p(:), mi(:), mj(:)
    complex(dp), allocatable :: mv(:)
    type(ICGNHier) :: hh
    type(ICGNBlock), allocatable :: pu(:)
    type(ICGNSigRow), allocatable :: sig(:)
    integer, allocatable :: roff(:)
    logical :: full3
    real(dp) :: t_all, t_stage, all_w, cut_w
    character(len=256) :: msg

    ierr = 0
    t_all = icgn_wall()
    axes = upt%icgn_nested_axes
    grids = upt%icgn_nested_grids
    do i = 1, 3
       if (axes(i) < 1 .or. axes(i) > 3) then
          ierr = 30; write(*,*) '(icgn) nested: nested_axes must be in {1,2,3}'; return
       end if
       if (grids(i) < 1) then
          ierr = 30; write(*,*) '(icgn) nested: nested_grids must be >= 1'; return
       end if
    end do
    if (axes(1) == axes(2) .or. axes(1) == axes(3) .or. axes(2) == axes(3)) then
       ierr = 30; write(*,*) '(icgn) nested: nested_axes must be a permutation of 1,2,3'; return
    end if
    if (.not. associated(upt%basis%coord)) then
       ierr = 30; write(*,*) '(icgn) nested: basis coordinates missing'; return
    end if
    if (size(upt%basis%coord, 2) < na) then
       ierr = 30; write(*,*) '(icgn) nested: coord size mismatch'; return
    end if
    nx = grids(1); ny = grids(2); nz = grids(3)
    nc = nx * ny * nz
    np = nz
    write(msg,'(a,3i4,a,3i4)') '(icgn) nested: axes=', axes, ' grids=', grids
    call icgn_log_important(upt, trim(msg))
    do i = 1, 3
       if (upt%icgn_tier_emin(i) >= upt%icgn_tier_emax(i) .or. upt%icgn_tier_top(i) < 0.0_dp .or. &
            upt%icgn_tier_bot(i) < 0.0_dp .or. upt%icgn_tier_eps(i) < 0.0_dp .or. &
            upt%icgn_nested_outers(i) < 0) then
          ierr = 3
          write(*,'(a,i0)') '(icgn) nested: invalid core window / buffer / epsilon / outers in tier ', i
          return
       end if
       write(msg,'(a,i0,a,2f10.4,a,2f8.4,a,es9.2,a,i0,a,i0,a,f10.4)') '(icgn) nested tier ', i, &
            ': core=', upt%icgn_tier_emin(i), upt%icgn_tier_emax(i), ' bot/top=', &
            upt%icgn_tier_bot(i), upt%icgn_tier_top(i), ' eps=', upt%icgn_tier_eps(i), &
            ' outers=', upt%icgn_nested_outers(i), ' neumann=', upt%icgn_nested_orders(i), &
            ' E0=', upt%icgn_tier_E0(i)
       call icgn_log_important(upt, trim(msg))
    end do
    if (all(upt%icgn_nested_outers == 0 .and. upt%icgn_nested_orders < 0)) &
         call icgn_log_important(upt, '(icgn) nested: fast path (window-only at every tier)')
    if (any(upt%icgn_nested_orders >= 0)) call icgn_log_important(upt, &
         'WARNING: Uptight does NOT check for Neumann-series convergence. Use it with care!')
    if (any(upt%icgn_nested_outers > 1)) call icgn_log_important(upt, &
         'Including core-state outers higher than level 1 is not yet implemented --> ' // &
         'Fall back to level-1 outers only.')
    full3 = icgn_tier_full(upt, 3)
    call icgn_stat_reset()

    ! ---- slabs along every axis (layer-aware, balanced; reuses geometry partition)
    allocate(slab(na, 3), part(na))
    slab = 0
    save_dir = upt%icgn_partition_dir; save_nb = upt%icgn_num_blocks
    do g = 1, 3
       if (grids(g) == 1) cycle
       upt%icgn_partition_dir = axes(g)
       upt%icgn_num_blocks = grids(g)
       call icgn_geometry_partition(upt, na, bsize, part, status)
       if (status == 0 .and. upt%icgn_num_blocks /= grids(g)) then
          status = 31
          write(*,'(a,i0,a)') '(icgn) nested: nested_grids(', g, ') exceeds the number of atomic layers'
       end if
       if (status /= 0) then
          upt%icgn_partition_dir = save_dir; upt%icgn_num_blocks = save_nb
          ierr = status
          return
       end if
       slab(:, g) = part
    end do
    upt%icgn_partition_dir = save_dir; upt%icgn_num_blocks = save_nb

    ! ---- orbital-level hierarchy maps ---------------------------------------
    hh%nx = nx; hh%ny = ny; hh%nz = nz
    allocate(hh%cube_of_orb(n_orb), hh%perm(n_orb), hh%pos(n_orb), hh%cstart(nc + 1))
    allocate(cnt(nc), cursor(nc))
    cnt = 0
    do a = 1, na
       ic = 1 + slab(a, 1) + nx * slab(a, 2) + nx * ny * slab(a, 3)
       do o = offsets(a), offsets(a + 1) - 1
          hh%cube_of_orb(o) = ic
       end do
       cnt(ic) = cnt(ic) + bsize(a)
    end do
    hh%cstart(1) = 1
    do ic = 1, nc
       hh%cstart(ic + 1) = hh%cstart(ic) + cnt(ic)
    end do
    cursor = hh%cstart(1:nc)
    do o = 1, n_orb
       ic = hh%cube_of_orb(o)
       hh%perm(cursor(ic)) = o
       hh%pos(o) = cursor(ic)
       cursor(ic) = cursor(ic) + 1
    end do
    cmax = maxval(cnt)
    n_empty = count(cnt == 0)
    deallocate(slab, part, cnt, cursor)
    if (n_empty > 0) then
       write(msg,'(a,i0,a,i0)') '(icgn) nested: WARNING empty cubes: ', n_empty, ' of ', nc
       call icgn_log_important(upt, trim(msg))
    end if

    ! ---- cubes -> rows -> planes: planes in parallel only if np >= threads ---
    T = icgn_thread_budget()
    nth_p = 1
    if (T > 1 .and. np >= T) nth_p = icgn_outer_threads(np, min(T, np), cmax)
    allocate(pu(np), ierr_p(np))
    ierr_p = 0
    t_stage = icgn_wall()
    if (nth_p > 1) then
       !$omp parallel do default(shared) private(ip) num_threads(nth_p) schedule(dynamic,1)
       do ip = 1, np
          call icgn_nested_plane(upt, hh, ip, 1, pu(ip), ierr_p(ip))
       end do
       !$omp end parallel do
    else
       do ip = 1, np
          call icgn_nested_plane(upt, hh, ip, T, pu(ip), ierr_p(ip))
          if (ierr_p(ip) /= 0) exit
       end do
    end if
    call icgn_report_time('nested: cubes -> rows -> planes', t_stage)
    do ip = 1, np
       if (ierr_p(ip) /= 0) ierr = ierr_p(ip)
    end do
    if (ierr /= 0) then
       call icgn_free_units(pu)
       return
    end if
    nred = 0
    do ip = 1, np
       nred = nred + pu(ip)%nret
    end do
    if (.not. full3 .and. nred < 1) then
       ierr = 9; write(*,*) '(icgn) no states retained'
       call icgn_free_units(pu)
       return
    end if

    ! ---- bulk H_CG: planes assembled with the generic routine ----------------
    t_stage = icgn_wall()
    if (full3) then
       ! plane tier with outers / Neumann: full-spectrum treatment of all planes
       call icgn_full_group(upt, hh, 3, 1, pu, np, T, mi, mj, mv, nnz, sig, roff, ierr)
       if (ierr == 0) then
          nred = 0
          do ip = 1, np
             nred = nred + pu(ip)%nret
          end do
          if (nred < 1) then
             ierr = 9; write(*,*) '(icgn) no states retained'
          end if
       end if
    else
       call icgn_cg_assemble(upt, upt%ham, hh, 3, 1, pu, np, T, .true., mi, mj, mv, nnz, ierr, tier=3)
    end if
    if (ierr /= 0) then
       call icgn_free_units(pu)
       return
    end if
    call icgn_report_time('nested: bulk H_CG (plane pairs)', t_stage)
    call icgn_stat_log(upt)
    if (associated(upt%icgn_ham%M)) call destroy_matrix(upt%icgn_ham)
    call create_matrix(upt%icgn_ham, nred, nred, nnz)
    upt%icgn_ham%sparse_fmt = 'F'
    upt%icgn_ham%Mi(1:nred + 1) = mi(1:nred + 1)
    upt%icgn_ham%Mj(1:nnz) = mj(1:nnz)
    upt%icgn_ham%M(1:nnz) = mv(1:nnz)
    upt%icgn_ham%nnz = nnz
    deallocate(mi, mj, mv)
    if (allocated(sig)) then
       t_stage = icgn_wall()
       call icgn_csr_add_sigma(upt%icgn_ham, roff, sig, nred, ierr)
       call icgn_free_sigrows(sig)
       if (ierr /= 0) then
          call icgn_log_important(upt, &
               'ERROR: Sigma merge aborted because CSR entry count exceeds default-integer capacity')
          ierr = 12
          call icgn_free_units(pu)
          return
       end if
       call icgn_report_time('nested: merge Sigma into bulk CSR', t_stage)
    end if

    ! ---- accumulated Q: one block per plane (moved, not copied) -------------
    allocate(upt%icgn_blocks(np))
    do ip = 1, np
       p1 = hh%cstart((ip - 1) * nx * ny + 1)
       p2 = hh%cstart(ip * nx * ny + 1) - 1
       upt%icgn_blocks(ip)%nrow = p2 - p1 + 1
       upt%icgn_blocks(ip)%nret = pu(ip)%nret
       allocate(upt%icgn_blocks(ip)%rows(p2 - p1 + 1))
       upt%icgn_blocks(ip)%rows = hh%perm(p1:p2)
       if (pu(ip)%nret > 0) then
          upt%icgn_blocks(ip)%eval => pu(ip)%eval
          upt%icgn_blocks(ip)%q => pu(ip)%q
          nullify(pu(ip)%eval, pu(ip)%q)
       end if
    end do
    call icgn_free_units(pu)

    ! cut fraction = weight of couplings between different planes
    all_w = 0.0_dp; cut_w = 0.0_dp
    do r = 1, n_orb
       pr = (hh%cube_of_orb(r) - 1) / (nx * ny)
       do k = upt%ham%Mi(r), upt%ham%Mi(r + 1) - 1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_w = all_w + abs(upt%ham%M(k))**2
          if ((hh%cube_of_orb(c) - 1) / (nx * ny) /= pr) cut_w = cut_w + abs(upt%ham%M(k))**2
       end do
    end do
    if (all_w > 0.0_dp) then
       upt%icgn_cut_fraction = cut_w / all_w
    else
       upt%icgn_cut_fraction = 0.0_dp
    end if
    upt%icgn_num_blocks = np
    upt%icgn_original_dim = n_orb
    upt%icgn_reduced_dim = nred

    call icgn_build_reduced_U(upt, nred, ierr)
    if (ierr /= 0) return
    call icgn_report_time('TOTAL ICGN PREPARATION TIME (nested)', t_all)
    upt%icgn_ready = .true.
    call icgn_log_info(upt, 'icgn', n_orb, nred, upt%icgn_num_blocks, upt%icgn_cut_fraction, &
         upt%icgn_sigma_T2, upt%icgn_pi_converged, .false.)
  end subroutine icgn_nested_run

  ! Entry point called by icgn_prepare when partition_type = 2.  Nesting is
  ! flattened for the whole run: exactly one parallel loop level is active, and
  ! BLAS/LAPACK is threaded only in the serial parts.
  subroutine icgn_nested_prepare(upt, na, n_orb, offsets, bsize, atom_of, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: na, n_orb
    integer, intent(in) :: offsets(:), bsize(:), atom_of(:)
    integer, intent(out) :: ierr
    integer :: saved_levels
    call icgn_flatten_nesting(saved_levels)
    call icgn_nested_run(upt, na, n_orb, offsets, bsize, atom_of, ierr)
    call icgn_restore_nesting(saved_levels)
  end subroutine icgn_nested_prepare

  ! Diagonalize block ib of icgn_blocks.  Thread-safe for distinct ib: it reads
  ! upt%ham and writes only upt%icgn_blocks(ib).  local(r) = position of row r

  ! inside its own block, blk_of(r) = block that owns row r (both cover all rows).
  subroutine icgn_diagonalize_block(upt, ib, local, blk_of, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, local(:), blk_of(:)
    integer, intent(out) :: ierr
    integer :: nn, k, r, c, jr, m, j
    complex(dp), allocatable :: h(:,:), z(:,:)
    real(dp), allocatable :: w(:)
    real(dp) :: vl, vu, margin
    logical :: fast_path
    ierr = 0; nn = upt%icgn_blocks(ib)%nrow
    if (nn == 0) then; upt%icgn_blocks(ib)%nret = 0; return; end if
    allocate(h(nn,nn), w(nn)); h = (0.0_dp, 0.0_dp)
    do jr = 1, nn
       r = upt%icgn_blocks(ib)%rows(jr)
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (blk_of(c) /= ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          h(local(r), local(c)) = upt%ham%M(k)
          if (r /= c) h(local(c), local(r)) = conjg(upt%ham%M(k))
       end do
    end do

    ! Fast path (no outer, no Neumann): only core+buffer states are kept later,
    ! so request that window from ZHEEVR. Outer / Neumann still need the full
    ! spectrum (discarded states enter P-Q and Q-Q couplings).
    fast_path = (upt%icgn_add_core_outers == 0 .and. upt%icgn_selfenergy_order < 0)
    if (fast_path) then
       vl = upt%icgn_core_emin - upt%icgn_bottom_buffer
       vu = upt%icgn_core_emax + upt%icgn_top_buffer
       ! Guard band outside the nominal keep window so ZHEEVR also returns
       ! near-edge states needed by boundary-degeneracy expansion (delta).
       margin = 1.0e-3_dp
       vl = vl - margin
       vu = vu + margin
       allocate(z(nn, nn))
       call dense_eigh_window(h, vl, vu, w, z, m, ierr)
       if (ierr /= 0) then
          ! Robust fallback: full ZHEEVD then same keep_mask filter as before.
          h = (0.0_dp, 0.0_dp)
          do jr = 1, nn
             r = upt%icgn_blocks(ib)%rows(jr)
             do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
                c = upt%ham%Mj(k)
                if (blk_of(c) /= ib) cycle
                if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
                h(local(r), local(c)) = upt%ham%M(k)
                if (r /= c) h(local(c), local(r)) = conjg(upt%ham%M(k))
             end do
          end do
          call coarse_eigh(h, w, ierr)
          if (ierr /= 0) then
             deallocate(h, w, z); return
          end if
          allocate(upt%icgn_blocks(ib)%evals_full(nn), upt%icgn_blocks(ib)%S_full(nn,nn))
          upt%icgn_blocks(ib)%evals_full = w
          upt%icgn_blocks(ib)%S_full     = h
          deallocate(h, w, z)
          return
       end if
       ! Store only the m window states as columns 1:m; remaining evals sit
       ! outside the keep window so the existing keep_mask logic is unchanged.
       allocate(upt%icgn_blocks(ib)%evals_full(nn), upt%icgn_blocks(ib)%S_full(nn,nn))
       upt%icgn_blocks(ib)%evals_full = vu + 1.0_dp
       upt%icgn_blocks(ib)%S_full     = (0.0_dp, 0.0_dp)
       do j = 1, m
          upt%icgn_blocks(ib)%evals_full(j) = w(j)
          upt%icgn_blocks(ib)%S_full(:, j)  = z(:, j)
       end do
       deallocate(h, w, z)
    else
       call coarse_eigh(h, w, ierr)
       if (ierr /= 0) then
          deallocate(h, w); return
       end if
       allocate(upt%icgn_blocks(ib)%evals_full(nn), upt%icgn_blocks(ib)%S_full(nn,nn))
       upt%icgn_blocks(ib)%evals_full = w
       upt%icgn_blocks(ib)%S_full     = h
       deallocate(h, w)
    end if
  end subroutine icgn_diagonalize_block

  ! Collect dense physical block couplings H_ab for every block pair that has
  ! at least one inter-block edge.  Allocates pairs to the exact pair count
  ! (no nb x nb map, no nb*(nb-1)/2).  Candidates collected from CSR, then
  ! sort+unique; dense h filled in parallel over pairs.
  subroutine icgn_collect_pairs(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(in) :: upt
    integer, intent(in) :: atom_of(:), label(:), local(:)
    type(ICGNPair), allocatable, intent(inout) :: pairs(:)
    integer, intent(out) :: npair, ierr
    integer :: r, k, c, ia, ib, slot, ne, e, m, i, nth, saved_levels, nc
    integer, allocatable :: cand_a(:), cand_b(:), e_slot(:), e_row(:), e_k(:), e_ord(:), p_start(:), p_fill(:)
    ierr = 0; npair = 0
    allocate(cand_a(max(1,upt%ham%nnz)), cand_b(max(1,upt%ham%nnz)))
    nc = 0
    ! Pass 1: collect every inter-block edge key (duplicates allowed).
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          if (upt%icgn_blocks(ia)%nrow == 0 .or. upt%icgn_blocks(ib)%nrow == 0) cycle
          nc = nc + 1
          cand_a(nc) = min(ia, ib); cand_b(nc) = max(ia, ib)
       end do
    end do
    call icgn_sort_pairs(cand_a, cand_b, nc)
    call icgn_unique_pairs(cand_a, cand_b, nc, npair)
    if (allocated(pairs)) call destroy_pairs(pairs)
    allocate(pairs(max(1, npair)))
    do i = 1, npair
       pairs(i)%a = cand_a(i); pairs(i)%b = cand_b(i)
    end do
    ! Pass 2: bucket CSR entries per pair (binary search for slot).
    allocate(e_slot(max(1,upt%ham%nnz)), e_row(max(1,upt%ham%nnz)), e_k(max(1,upt%ham%nnz)))
    ne = 0
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          if (upt%icgn_blocks(ia)%nrow == 0 .or. upt%icgn_blocks(ib)%nrow == 0) cycle
          slot = icgn_find_pair(cand_a, cand_b, npair, min(ia,ib), max(ia,ib))
          if (slot == 0) cycle
          ne = ne + 1
          if (ia < ib) then; e_slot(ne) = slot; else; e_slot(ne) = -slot; end if
          e_row(ne) = r; e_k(ne) = k
       end do
    end do
    deallocate(cand_a, cand_b)
    allocate(p_start(npair+1), p_fill(max(1,npair)), e_ord(max(1,ne)))
    p_start = 0
    do e = 1, ne
       p_start(abs(e_slot(e))+1) = p_start(abs(e_slot(e))+1) + 1
    end do
    do i = 1, npair
       p_start(i+1) = p_start(i+1) + p_start(i)
    end do
    p_fill(1:npair) = p_start(1:npair)
    do e = 1, ne
       i = abs(e_slot(e)); p_fill(i) = p_fill(i) + 1; e_ord(p_fill(i)) = e
    end do
    nth = min(icgn_thread_budget(), max(1, npair))
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel do default(shared) private(i, m, e, r, c, k) num_threads(nth) &
    !$omp& schedule(dynamic,1) if(nth > 1)
    do i = 1, npair
       allocate(pairs(i)%h(upt%icgn_blocks(pairs(i)%a)%nrow, upt%icgn_blocks(pairs(i)%b)%nrow))
       pairs(i)%h = (0.0_dp, 0.0_dp)
       do m = p_start(i)+1, p_start(i+1)
          e = e_ord(m); r = e_row(e); k = e_k(e); c = upt%ham%Mj(k)
          if (e_slot(e) > 0) then
             pairs(i)%h(local(r),local(c)) = pairs(i)%h(local(r),local(c)) + upt%ham%M(k)
          else
             pairs(i)%h(local(c),local(r)) = pairs(i)%h(local(c),local(r)) + conjg(upt%ham%M(k))
          end if
       end do
    end do
    !$omp end parallel do
    call icgn_restore_nesting(saved_levels)
    deallocate(e_slot, e_row, e_k, e_ord, p_start, p_fill)
  end subroutine icgn_collect_pairs

  ! Heapsort pairs (pa(i), pb(i)) by (a,b) ascending.  Used to unique inter-block
  ! pairs without a dense nb x nb map.
  pure subroutine icgn_sort_pairs(pa, pb, n)
    integer, intent(inout) :: pa(:), pb(:)
    integer, intent(in) :: n
    integer :: start, root, child, end_, ta, tb
    if (n <= 1) return
    start = n / 2
    do while (start >= 1)
       root = start
       do
          child = 2 * root
          if (child > n) exit
          if (child < n) then
             if (pa(child) < pa(child+1) .or. &
                  (pa(child) == pa(child+1) .and. pb(child) < pb(child+1))) child = child + 1
          end if
          if (pa(root) < pa(child) .or. &
               (pa(root) == pa(child) .and. pb(root) < pb(child))) then
             ta = pa(root); tb = pb(root)
             pa(root) = pa(child); pb(root) = pb(child)
             pa(child) = ta; pb(child) = tb
             root = child
          else
             exit
          end if
       end do
       start = start - 1
    end do
    end_ = n
    do while (end_ > 1)
       ta = pa(end_); tb = pb(end_)
       pa(end_) = pa(1); pb(end_) = pb(1)
       pa(1) = ta; pb(1) = tb
       end_ = end_ - 1
       root = 1
       do
          child = 2 * root
          if (child > end_) exit
          if (child < end_) then
             if (pa(child) < pa(child+1) .or. &
                  (pa(child) == pa(child+1) .and. pb(child) < pb(child+1))) child = child + 1
          end if
          if (pa(root) < pa(child) .or. &
               (pa(root) == pa(child) .and. pb(root) < pb(child))) then
             ta = pa(root); tb = pb(root)
             pa(root) = pa(child); pb(root) = pb(child)
             pa(child) = ta; pb(child) = tb
             root = child
          else
             exit
          end if
       end do
    end do
  end subroutine icgn_sort_pairs

  ! Compact sorted (pa,pb) in place; nu = number of unique pairs.
  pure subroutine icgn_unique_pairs(pa, pb, n, nu)
    integer, intent(inout) :: pa(:), pb(:)
    integer, intent(in) :: n
    integer, intent(out) :: nu
    integer :: i
    if (n <= 0) then; nu = 0; return; end if
    nu = 1
    do i = 2, n
       if (pa(i) /= pa(nu) .or. pb(i) /= pb(nu)) then
          nu = nu + 1
          pa(nu) = pa(i); pb(nu) = pb(i)
       end if
    end do
  end subroutine icgn_unique_pairs

  ! Binary search for (ia,ib) in sorted unique pair lists.  Returns slot or 0.
  pure integer function icgn_find_pair(pa, pb, n, ia, ib) result(slot)
    integer, intent(in) :: pa(:), pb(:), n, ia, ib
    integer :: lo, hi, mid
    lo = 1; hi = n; slot = 0
    do while (lo <= hi)
       mid = (lo + hi) / 2
       if (pa(mid) < ia .or. (pa(mid) == ia .and. pb(mid) < ib)) then
          lo = mid + 1
       else if (pa(mid) > ia .or. (pa(mid) == ia .and. pb(mid) > ib)) then
          hi = mid - 1
       else
          slot = mid; return
       end if
    end do
  end function icgn_find_pair

  ! ---------------------------------------------------------------------------
  ! OpenMP helpers for the ICGN preparation stages.
  !
  ! Thread budget = omp_get_max_threads() (already divided by the number of MPI
  ! processes in uptight.F90; under SLURM it follows OMP_NUM_THREADS /
  ! --cpus-per-task / the cpuset).  The optional environment variable
  ! UPT_CG_THREADS can only LOWER it (UPT_CG_THREADS=1 disables the outer
  ! parallelism, i.e. recovers the old serial-outer behaviour).  Inside an
  ! already-parallel region the budget is 1.  Without OpenMP it is always 1.
  ! ---------------------------------------------------------------------------
  subroutine project_pair(p, qa, qb, keep_h)
    type(ICGNPair), intent(inout) :: p
    complex(dp), intent(in) :: qa(:,:), qb(:,:)
    logical, intent(in), optional :: keep_h
    complex(dp), allocatable :: tmp(:,:)
    logical :: keep

    ! Fortran does not short-circuit .or./.and., so an absent optional argument
    ! must never appear in the same expression as its use.
    keep = .false.
    if (present(keep_h)) keep = keep_h

    ! First H_AB * Q_B, then Q_A^H * (H_AB * Q_B).
    ! The expensive arithmetic is now done at block level instead of once per
    ! physical Hamiltonian entry.
    allocate(p%v(size(qa,2), size(qb,2)))
    call icgn_project_sparse(p%h, qa, qb, p%v)
    if (.not. keep) deallocate(p%h)
  end subroutine project_pair

  subroutine project_accumulated_pair_icgn(p, qa)
    type(ICGNPair), intent(inout) :: p
    complex(dp), intent(in) :: qa(:,:)

    allocate(p%v(size(qa,2), size(p%t,2)))
    call icgn_zgemm('C', qa, p%t, p%v)
    deallocate(p%t)
  end subroutine project_accumulated_pair_icgn

  ! out = qa^H * h * qb, exploiting that the physical block coupling h has only a
  ! few nonzero rows and columns (the atoms next to the interface): only those
  ! rows of qa/h and columns of h/qb enter, so the dense n^3 products shrink to
  ! (#nonzero rows) x (#nonzero columns) sized ones.  Same result up to rounding.
  subroutine icgn_project_sparse(h, qa, qb, out)
    complex(dp), intent(in) :: h(:,:), qa(:,:), qb(:,:)
    complex(dp), intent(inout) :: out(:,:)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp), czero = (0.0_dp, 0.0_dp)
    integer :: nra, nrb, i, j, nr, nc
    integer, allocatable :: rz(:), cz(:)
    logical, allocatable :: rmask(:), cmask(:)
    complex(dp), allocatable :: hs(:,:), qas(:,:), qbs(:,:), tmp(:,:)

    nra = size(h,1); nrb = size(h,2)
    allocate(rmask(nra), cmask(nrb))
    rmask = .false.; cmask = .false.
    do j = 1, nrb
       do i = 1, nra
          if (h(i,j) /= czero) then
             rmask(i) = .true.; cmask(j) = .true.
          end if
       end do
    end do
    nr = count(rmask); nc = count(cmask)
    if (nr == 0 .or. nc == 0) then
       out = czero
       return
    end if
    allocate(rz(nr), cz(nc))
    nr = 0
    do i = 1, nra
       if (rmask(i)) then; nr = nr + 1; rz(nr) = i; end if
    end do
    nc = 0
    do j = 1, nrb
       if (cmask(j)) then; nc = nc + 1; cz(nc) = j; end if
    end do
    allocate(hs(nr, nc), qas(nr, size(qa,2)), qbs(nc, size(qb,2)), tmp(nr, size(qb,2)))
    hs  = h(rz, cz)
    qas = qa(rz, :)
    qbs = qb(cz, :)
    tmp = czero
    call icgn_zgemm_gen('N', 'N', cone, hs, qbs, czero, tmp)      ! (nr x nc)(nc x nB)
    call icgn_zgemm_gen('C', 'N', cone, qas, tmp, czero, out)     ! (nA x nr)(nr x nB)
  end subroutine icgn_project_sparse

  ! Write the coupling block p%v of the block pair (a,b) into the reduced CSR h.
  !   oa, ob     : offsets of the retained states of blocks a and b in the reduced basis
  !   fmt        : 'F' writes the block and its Hermitian conjugate, 'U' only
  !                the block (upper), 'L' only the conjugate (lower)
  !   next(:)    : running insertion position of every reduced row (updated)
  !   keep_zeros : if .true. exact zeros are written too. Default: zeros are
  !                skipped. The caller must size the CSR with the same rule
  !                (caller counts only non-zero entries).
  ! The number of slots written must match what the caller counted, otherwise
  ! the tail of M/Mj of a row is left uninitialized.
  subroutine emit_pair(h,p,oa,ob,fmt,next,keep_zeros)
    type(CSR),intent(inout)::h
    type(ICGNPair),intent(in)::p
    integer,intent(in)::oa,ob
    character(1),intent(in)::fmt
    integer,intent(inout)::next(:)
    logical,intent(in),optional::keep_zeros
    integer::i,j,k
    logical::skip_zero
    ! CSR sizing counts only non-zero entries, so
    ! zeros are skipped by default; with keep_zeros the caller must have
    ! sized the CSR from the dense blocks, otherwise the unfilled slots of
    ! M/Mj stay uninitialized.
    skip_zero = .true.
    if (present(keep_zeros)) skip_zero = .not. keep_zeros
    if(fmt/='L') then
       do i=1,size(p%v,1); do j=1,size(p%v,2)
          if (skip_zero .and. abs(p%v(i,j)) == 0.0_dp) cycle
          k=next(oa+i-1); h%Mj(k)=ob+j-1; h%M(k)=p%v(i,j); next(oa+i-1)=k+1
       end do; end do
    end if
    if(fmt=='F' .or. fmt=='L') then
       do j=1,size(p%v,2); do i=1,size(p%v,1)
          if (skip_zero .and. abs(p%v(i,j)) == 0.0_dp) cycle
          k=next(ob+j-1); h%Mj(k)=oa+i-1; h%M(k)=conjg(p%v(i,j)); next(ob+j-1)=k+1
       end do; end do
    end if
  end subroutine emit_pair

  ! Sort column indices of every CSR row into strictly ascending order
  ! (insertion sort; rows are short).  Required so the reduced operator
  ! always matches the documented output convention regardless of path.
  subroutine icgn_sort_csr_rows(h)
    type(CSR), intent(inout) :: h
    integer :: r, i, j, key_j
    complex(dp) :: key_m
    do r = 1, h%nrow
       do i = h%Mi(r)+1, h%Mi(r+1)-1
          key_j = h%Mj(i); key_m = h%M(i)
          j = i - 1
          do while (j >= h%Mi(r))
             if (h%Mj(j) <= key_j) exit
             h%Mj(j+1) = h%Mj(j); h%M(j+1) = h%M(j)
             j = j - 1
          end do
          h%Mj(j+1) = key_j; h%M(j+1) = key_m
       end do
    end do
  end subroutine icgn_sort_csr_rows

  ! Build reduced Hamiltonian H_PP (CSR) from retained states and pair couplings.
  subroutine build_icgn_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, reuse_full, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(ICGNPair), allocatable, intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    ! .true.: pairs already hold the full projected coupling v (and h) from
    ! icgn_prepare (outer path); only slicing to retained states is left.
    logical, intent(in) :: reuse_full
    integer, intent(out) :: ierr
    integer :: i, j, k, r, c, a, b, ia, ib, nnz, pos, slot, nred, nth, saved_levels
    integer :: ne, m, e
    integer, allocatable :: e_slot(:), e_row(:), e_k(:), e_ord(:), p_start(:), p_fill(:)
    real(dp) :: t_stage
    integer, allocatable :: roff(:), rowcount(:), next(:), cand_a(:), cand_b(:)
    complex(dp), allocatable :: g_full(:,:)
    integer :: nc
    ierr = 0; nred = upt%icgn_reduced_dim
    allocate(roff(upt%icgn_num_blocks+1)); roff(1) = 1
    do i = 1, upt%icgn_num_blocks; roff(i+1) = roff(i) + upt%icgn_blocks(i)%nret; end do

    if (upt%icgn_add_core_outers == 0 .and. upt%icgn_selfenergy_order < 0) then
       ! Fast path: accumulate T = H_AB * Q_B from global CSR (no dense H_AB).
       ! Discover unique pairs via sort+unique (no nb x nb map).
       npair = 0
       allocate(cand_a(max(1,upt%ham%nnz)), cand_b(max(1,upt%ham%nnz)))
       t_stage = icgn_wall()
       nc = 0
       do r = 1, upt%ham%nrow
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
             if (a == b) cycle
             if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             nc = nc + 1
             cand_a(nc) = min(a,b); cand_b(nc) = max(a,b)
          end do
       end do
       call icgn_sort_pairs(cand_a, cand_b, nc)
       call icgn_unique_pairs(cand_a, cand_b, nc, npair)
       if (allocated(pairs)) call destroy_pairs(pairs)
       allocate(pairs(max(1, npair)))
       do i = 1, npair
          pairs(i)%a = cand_a(i); pairs(i)%b = cand_b(i)
       end do
       allocate(e_slot(max(1,upt%ham%nnz)), e_row(max(1,upt%ham%nnz)), e_k(max(1,upt%ham%nnz)))
       ne = 0
       do r = 1, upt%ham%nrow
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
             if (a == b) cycle
             if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             ia = min(a,b); ib = max(a,b)
             slot = icgn_find_pair(cand_a, cand_b, npair, ia, ib)
             if (slot == 0) cycle
             ne = ne + 1
             if (a < b) then; e_slot(ne) = slot; else; e_slot(ne) = -slot; end if
             e_row(ne) = r; e_k(ne) = k
          end do
       end do
       deallocate(cand_a, cand_b)
       ! Counting sort of the entries by pair, keeping the scan order inside each
       ! pair so that the accumulation order (hence the result) is unchanged.
       allocate(p_start(npair+1), p_fill(max(1,npair)), e_ord(max(1,ne)))
       p_start = 0
       do e = 1, ne
          p_start(abs(e_slot(e))+1) = p_start(abs(e_slot(e))+1) + 1
       end do
       do i = 1, npair
          p_start(i+1) = p_start(i+1) + p_start(i)
       end do
       p_fill(1:npair) = p_start(1:npair)
       do e = 1, ne
          i = abs(e_slot(e)); p_fill(i) = p_fill(i) + 1; e_ord(p_fill(i)) = e
       end do
       call icgn_report_time('  fast path: discover pairs (serial)', t_stage)
       t_stage = icgn_wall()
       ! Phase B (parallel over pairs): each pair owns its T block, so the
       ! accumulation T = H_AB * Q_B and the left projection V_AB = Q_A^H T need no
       ! locking.
       nth = min(icgn_thread_budget(), max(1, npair))
       call icgn_flatten_nesting(saved_levels)
       !$omp parallel do default(shared) private(i, a, b, j, k, m, e, r, c) num_threads(nth) &
       !$omp& schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
          allocate(pairs(i)%t(upt%icgn_blocks(a)%nrow, upt%icgn_blocks(b)%nret))
          pairs(i)%t = (0.0_dp, 0.0_dp)
          do m = p_start(i)+1, p_start(i+1)
             e = e_ord(m); r = e_row(e); k = e_k(e); c = upt%ham%Mj(k)
             if (e_slot(e) > 0) then
                pairs(i)%t(local(r),:) = pairs(i)%t(local(r),:) + &
                     upt%ham%M(k) * upt%icgn_blocks(b)%q(local(c),:)
             else
                pairs(i)%t(local(c),:) = pairs(i)%t(local(c),:) + &
                     conjg(upt%ham%M(k)) * upt%icgn_blocks(b)%q(local(r),:)
             end if
          end do
          call project_accumulated_pair_icgn(pairs(i), upt%icgn_blocks(a)%q)
          do j = 1, upt%icgn_blocks(a)%nret
             do k = 1, upt%icgn_blocks(b)%nret
                if (.not. icgn_coupling_survives(pairs(i)%v(j,k), &
                    upt%icgn_blocks(a)%eval(j), upt%icgn_blocks(b)%eval(k), &
                    upt%icgn_epsilon, &
                    upt%icgn_core_emin - upt%icgn_bottom_buffer, &
                    upt%icgn_core_emax + upt%icgn_top_buffer)) pairs(i)%v(j,k) = (0.0_dp, 0.0_dp)
             end do
          end do
       end do
       !$omp end parallel do
       call icgn_restore_nesting(saved_levels)
       call icgn_report_time('  fast path: accumulate + project (parallel)', t_stage)
       deallocate(e_slot, e_row, e_k, e_ord, p_start, p_fill)
    else if (upt%icgn_add_core_outers == 0) then
       ! Zero-outer path with Neumann enabled: keep the full physical
       ! pair matrices because the correction needs P-Q and Q-Q couplings.
       ! Blocks with nret==0 have no retained basis (q unallocated); emit an
       ! empty v so the later CSR count/assembly is safe.
       nth = min(icgn_thread_budget(), max(1, npair))
       call icgn_flatten_nesting(saved_levels)
       !$omp parallel do default(shared) private(i, a, b, j, k) num_threads(nth) schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
          if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) then
             allocate(pairs(i)%v(upt%icgn_blocks(a)%nret, upt%icgn_blocks(b)%nret))
             cycle
          end if
          call project_pair(pairs(i), upt%icgn_blocks(a)%q, &
               upt%icgn_blocks(b)%q, .true.)
          do j = 1, upt%icgn_blocks(a)%nret
             do k = 1, upt%icgn_blocks(b)%nret
                if (.not. icgn_coupling_survives(pairs(i)%v(j,k), &
                    upt%icgn_blocks(a)%eval(j), upt%icgn_blocks(b)%eval(k), &
                    upt%icgn_epsilon, &
                    upt%icgn_core_emin - upt%icgn_bottom_buffer, &
                    upt%icgn_core_emax + upt%icgn_top_buffer)) pairs(i)%v(j,k) = (0.0_dp, 0.0_dp)
             end do
          end do
       end do
       !$omp end parallel do
       call icgn_restore_nesting(saved_levels)
    else
       if (.not. reuse_full) then
          ! collect_pairs allocates pairs to the exact edged-pair count.
          call icgn_collect_pairs(upt, atom_of, label, local, pairs, npair, ierr)
          if (ierr /= 0) return
       end if
       ! Per pair: project H_ab with the full eigenbasis (unless already done in
       ! icgn_prepare), then slice to the retained states.  Pairs are
       ! independent; g_full is a per-thread temporary.
       nth = min(icgn_thread_budget(), max(1, npair))
       call icgn_report_time('  outer path: start slice pairs', t_stage)
       t_stage = icgn_wall()
       call icgn_flatten_nesting(saved_levels)
       !$omp parallel do default(shared) private(i, a, b, j, k, g_full) num_threads(nth) &
       !$omp& schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
          ! After keep_mask some blocks may retain zero states. retained_idx/q/eval
          ! are not allocated for those blocks — must not touch them.
          if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) then
             if (associated(pairs(i)%v)) deallocate(pairs(i)%v)
             allocate(pairs(i)%v(upt%icgn_blocks(a)%nret, upt%icgn_blocks(b)%nret))
             if (associated(pairs(i)%h)) deallocate(pairs(i)%h)
             nullify(pairs(i)%g)
             cycle
          end if
          if (.not. associated(pairs(i)%v)) then
             ! Project path left this pair without a coupling matrix; emit empty v.
             allocate(pairs(i)%v(upt%icgn_blocks(a)%nret, upt%icgn_blocks(b)%nret))
             pairs(i)%v = (0.0_dp, 0.0_dp)
             if (associated(pairs(i)%h)) deallocate(pairs(i)%h)
             nullify(pairs(i)%g)
             cycle
          end if
          if (.not. reuse_full) call project_pair(pairs(i), upt%icgn_blocks(a)%S_full, &
               upt%icgn_blocks(b)%S_full, .true.)

          ! Slice to retained states using retained_idx.
          if (upt%icgn_selfenergy_order >= 0) then
             ! Neumann follows: keep the FULL projection in pairs(i)%g (it holds the
             ! P-Q and Q-Q couplings) so it is not computed a second time.
             pairs(i)%g => pairs(i)%v
             nullify(pairs(i)%v)
             allocate(pairs(i)%v(upt%icgn_blocks(a)%nret, upt%icgn_blocks(b)%nret))
             do j = 1, upt%icgn_blocks(b)%nret
                do k = 1, upt%icgn_blocks(a)%nret
                   pairs(i)%v(k,j) = pairs(i)%g(upt%icgn_blocks(a)%retained_idx(k), &
                                                 upt%icgn_blocks(b)%retained_idx(j))
                   if (.not. icgn_coupling_survives(pairs(i)%v(k,j), upt%icgn_blocks(a)%eval(k), &
                       upt%icgn_blocks(b)%eval(j), upt%icgn_epsilon, &
                    upt%icgn_core_emin - upt%icgn_bottom_buffer, &
                    upt%icgn_core_emax + upt%icgn_top_buffer)) &
                       pairs(i)%v(k,j) = (0.0_dp, 0.0_dp)
                end do
             end do
          else
             allocate(g_full(size(pairs(i)%v,1), size(pairs(i)%v,2)))
             g_full = pairs(i)%v
             deallocate(pairs(i)%v)
             allocate(pairs(i)%v(upt%icgn_blocks(a)%nret, upt%icgn_blocks(b)%nret))
             do j = 1, upt%icgn_blocks(b)%nret
                do k = 1, upt%icgn_blocks(a)%nret
                   pairs(i)%v(k,j) = g_full(upt%icgn_blocks(a)%retained_idx(k), &
                                             upt%icgn_blocks(b)%retained_idx(j))
                   if (.not. icgn_coupling_survives(pairs(i)%v(k,j), upt%icgn_blocks(a)%eval(k), &
                       upt%icgn_blocks(b)%eval(j), upt%icgn_epsilon, &
                    upt%icgn_core_emin - upt%icgn_bottom_buffer, &
                    upt%icgn_core_emax + upt%icgn_top_buffer)) &
                       pairs(i)%v(k,j) = (0.0_dp, 0.0_dp)
                end do
             end do
             deallocate(g_full)
          end if
          ! The physical block H_ab is no longer needed on this path.
          if (associated(pairs(i)%h)) deallocate(pairs(i)%h)
       end do
       !$omp end parallel do
       call icgn_restore_nesting(saved_levels)
       call icgn_report_time('  outer path: slice pairs done', t_stage)
    end if

    t_stage = icgn_wall()
    ! Build CSR in the canonical reduced format (always 'F', diagonal present,
    ! columns strictly ascending per row). Count both triangles for every
    ! non-zero coupling; zeros are omitted.
    if (nred <= 0) then
       ierr = 9
       write(*,*) '(icgn) build reduced H: nred <= 0'
       return
    end if
    allocate(rowcount(nred), next(nred)); rowcount = 1
    do i = 1, npair
       if (.not. associated(pairs(i)%v)) cycle
       a = pairs(i)%a; b = pairs(i)%b
       if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) cycle
       do j = 1, size(pairs(i)%v,1)
          do k = 1, size(pairs(i)%v,2)
             if (abs(pairs(i)%v(j,k)) == 0.0_dp) cycle
             rowcount(roff(a)+j-1) = rowcount(roff(a)+j-1) + 1
             rowcount(roff(b)+k-1) = rowcount(roff(b)+k-1) + 1
          end do
       end do
    end do
    ! Guard against 32-bit integer overflow of nnz (common when H_CG is dense).
    block
       integer(kind=8) :: nnz8
       nnz8 = 0_8
       do i = 1, nred
          nnz8 = nnz8 + int(rowcount(i), kind=8)
       end do
       if (nnz8 > int(huge(nnz), kind=8)) then
          write(*,'(a,i0,a)') '(icgn) ERROR: reduced CSR nnz=', nnz8, &
               ' exceeds default INTEGER capacity; tighten core window / reduce outers'
          ierr = 12
          deallocate(roff, rowcount, next)
          return
       end if
       nnz = int(nnz8)
    end block
    if (icgn_timing_on()) then
       write(*,'(a,i0,a,i0,a,i0)') '[ICGN-TIME] assemble prep: nred=', nred, &
            ' npair=', npair, ' nnz=', nnz
    end if
    call create_matrix(upt%icgn_ham, nred, nred, nnz)
    upt%icgn_ham%sparse_fmt = 'F'; upt%icgn_ham%Mi(1) = 1
    do i = 1, nred; upt%icgn_ham%Mi(i+1) = upt%icgn_ham%Mi(i) + rowcount(i); end do
    next = upt%icgn_ham%Mi(1:nred)
    do a = 1, upt%icgn_num_blocks
       do i = 1, upt%icgn_blocks(a)%nret
          pos = next(roff(a)+i-1)
          upt%icgn_ham%Mj(pos) = roff(a)+i-1
          upt%icgn_ham%M(pos)  = upt%icgn_blocks(a)%eval(i)
          next(roff(a)+i-1) = pos + 1
       end do
    end do
    do i = 1, npair
       if (.not. associated(pairs(i)%v)) cycle
       a = pairs(i)%a; b = pairs(i)%b
       if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) cycle
       call emit_pair(upt%icgn_ham, pairs(i), roff(a), roff(b), 'F', next)
    end do
    upt%icgn_ham%nnz = nnz
    ! Ensure column indices are strictly ascending inside every row.
    call icgn_sort_csr_rows(upt%icgn_ham)
    deallocate(roff, rowcount, next)
    call icgn_report_time('  assemble reduced CSR (serial)', t_stage)
  end subroutine build_icgn_reduced_hamiltonian

  subroutine destroy_pairs(pairs)
    type(ICGNPair),allocatable,intent(inout)::pairs(:)
    integer::i,nv,nh,nt
    if(.not.allocated(pairs)) then

       return
    end if
    nv=0; nh=0; nt=0
    do i=1,size(pairs)
       if(associated(pairs(i)%v)) nv=nv+1
       if(associated(pairs(i)%h)) nh=nh+1
       if(associated(pairs(i)%t)) nt=nt+1
    end do

         ! ' v=', nv, ' h=', nh, ' t=', nt
    do i=1,size(pairs)
       if(associated(pairs(i)%v)) deallocate(pairs(i)%v)
       if(associated(pairs(i)%h)) deallocate(pairs(i)%h)
       if(associated(pairs(i)%t)) deallocate(pairs(i)%t)
       if(associated(pairs(i)%g)) deallocate(pairs(i)%g)
    end do
    deallocate(pairs)

  end subroutine destroy_pairs

  ! Free full-block eigensystems left on upt%icgn_blocks after a failed prepare
  ! stage (diagonalize / collect / build).  Safe to call when blocks are absent.
  subroutine icgn_free_block_eigensystems(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    if (.not. associated(upt%icgn_blocks)) return
    do i = 1, size(upt%icgn_blocks)
       if (associated(upt%icgn_blocks(i)%S_full))     deallocate(upt%icgn_blocks(i)%S_full)
       if (associated(upt%icgn_blocks(i)%evals_full)) deallocate(upt%icgn_blocks(i)%evals_full)
       nullify(upt%icgn_blocks(i)%S_full, upt%icgn_blocks(i)%evals_full)
    end do
  end subroutine icgn_free_block_eigensystems

! ===========================================================================
! SECTION C — Neumann self-energy (Sigma)
! ===========================================================================

  ! Noise filter for couplings that are KEPT in a reduced / self-energy matrix.
  ! Only couplings whose magnitude is at the round-off level (|V| < icgn_zero_tol,
  ! in eV) are dropped, to save memory; every other coupling is stored unchanged.
  ! The physical threshold outer_coupling_threshold is NOT applied here: it only
  ! decides which discarded states are promoted to outers (icgn_outer_passes).
  ! (ea, eb, epsilon, win_emin, win_emax are kept only so that existing call
  ! sites need no change; they are ignored.)
  pure logical function icgn_coupling_survives(v, ea, eb, epsilon, win_emin, win_emax)
    complex(dp), intent(in) :: v
    real(dp), intent(in) :: ea, eb, epsilon, win_emin, win_emax
    real(dp) :: v2
    v2 = real(v,dp)**2 + aimag(v)**2
    icgn_coupling_survives = (v2 >= icgn_zero_tol * icgn_zero_tol)   ! false also for NaN
  end function icgn_coupling_survives

  ! Outer-state admission test (and nothing else): a discarded state j is promoted
  ! to the retained space when at least one coupling V_ij to a core state i
  ! satisfies   |V_ij|^2 > epsilon * |E_i - E_j|   (strict).  The product form
  ! avoids dividing by a vanishing gap; for degenerate states the right-hand
  ! side is 0, so any coupling above round-off passes.  Round-off couplings
  ! never pass, even for epsilon = 0.
  pure logical function icgn_outer_passes(v, ei, ej, epsilon)
    complex(dp), intent(in) :: v
    real(dp), intent(in) :: ei, ej, epsilon
    real(dp) :: v2
    v2 = real(v,dp)**2 + aimag(v)**2
    icgn_outer_passes = (v2 >= icgn_zero_tol * icgn_zero_tol) .and. &
         (v2 > max(epsilon, 0.0_dp) * abs(ei - ej))
  end function icgn_outer_passes


  ! Dense filtered P-Q / Q-Q coupling blocks of every block pair.  blk(1:nb) holds
  ! the full eigensystem (evals_full, S_full or the already projected pairs%g,
  ! retained_idx, nret) of every block; the number of states of a block is
  ! size(evals_full).  tb (optional): thread budget of this call; when present the
  ! caller owns the OpenMP nesting level (no flattening here, safe inside an
  ! outer parallel region).
  subroutine icgn_build_neumann_blocks(upt, blk, nb, pairs, npair, local_ret_idx, E0, need_w, &
       bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, npq_nz, nqq_nz, tb, &
       elo_o, ehi_o, eps_o)
    type(OUPT), intent(in) :: upt
    type(ICGNBlock), intent(in) :: blk(:)
    integer, intent(in) :: nb
    type(ICGNPair), intent(in) :: pairs(:)
    integer, intent(in) :: npair, local_ret_idx(:,:)
    real(dp), intent(in) :: E0
    logical, intent(in) :: need_w
    type(ICGNQBlock), allocatable, intent(out) :: bq(:)
    type(ICGNQPair), allocatable, intent(out) :: qp(:)
    integer, allocatable, intent(out) :: nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, allocatable, intent(out) :: nbr_isA(:)
    integer(8), intent(out) :: npq_nz, nqq_nz
    integer, intent(in), optional :: tb
    ! Optional filter window / epsilon (nested tiers); default: scalar ICGN values.
    real(dp), intent(in), optional :: elo_o, ehi_o, eps_o
    real(dp) :: elo_u, ehi_u, eps_u
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    real(dp), parameter :: thr = 1.0e-14_dp
    integer :: b, i, j, k, ia, ib, jq, kq, ra, rb, e, np_a, np_b, nq_a, nq_b, nfb
    integer :: nth, saved_levels, nra, nrb, n_near
    integer, allocatable :: cnt(:), fillp(:), eA(:), eB(:)
    complex(dp) :: v
    complex(dp), pointer :: g(:,:)
    complex(dp), allocatable :: tmp(:,:)
    logical :: own_g
    integer(8) :: nzp, nzq
    real(dp) :: gap, gap_min
    character(len=160) :: msg

    ! ---- per-block Q data ----
    elo_u = upt%icgn_core_emin - upt%icgn_bottom_buffer
    ehi_u = upt%icgn_core_emax + upt%icgn_top_buffer
    eps_u = upt%icgn_epsilon
    if (present(elo_o)) elo_u = elo_o
    if (present(ehi_o)) ehi_u = ehi_o
    if (present(eps_o)) eps_u = eps_o
    ! Discarded eigenvalues lie outside the core window, so E0 (inside the open
    ! core) cannot match them exactly; still warn if |E0-E| is tiny (numerical).
    allocate(bq(nb))
    n_near = 0
    gap_min = huge(1.0_dp)
    do b = 1, nb
       nfb = 0
       if (associated(blk(b)%evals_full)) nfb = size(blk(b)%evals_full)
       bq(b)%np = blk(b)%nret
       bq(b)%nq = nfb - blk(b)%nret
       allocate(bq(b)%qidx(max(0,bq(b)%nq)), bq(b)%res(max(0,bq(b)%nq)))
       kq = 0
       do k = 1, nfb
          if (local_ret_idx(b, k) == 0) then
             kq = kq + 1
             bq(b)%qidx(kq) = k
             gap = E0 - blk(b)%evals_full(k)
             if (abs(gap) < 1.0e-8_dp) then
                n_near = n_near + 1
                gap_min = min(gap_min, abs(gap))
             end if
             bq(b)%res(kq) = 1.0_dp / gap
          end if
       end do
    end do
    if (n_near > 0) then
       write(msg,'(a,i0,a,es12.5)') &
            'Neumann resolvent: |E0-Eq| < 1e-8 for ', n_near, &
            ' discarded state(s); min|gap|=', gap_min
       !$omp critical (icgn_log)
       call icgn_log_important(upt, trim(msg))
       !$omp end critical (icgn_log)
    end if

    ! ---- adjacency of the block graph (entries of block a: its neighbours) ----
    allocate(nbr_ptr(nb+1), cnt(nb))
    cnt = 0
    do i = 1, npair
       cnt(pairs(i)%a) = cnt(pairs(i)%a) + 1
       cnt(pairs(i)%b) = cnt(pairs(i)%b) + 1
    end do
    nbr_ptr(1) = 1
    do b = 1, nb
       nbr_ptr(b+1) = nbr_ptr(b) + cnt(b)
    end do
    allocate(nbr_blk(max(1,2*npair)), nbr_pair(max(1,2*npair)), nbr_isA(max(1,2*npair)), &
         nbr_rev(max(1,2*npair)), fillp(nb), eA(max(1,npair)), eB(max(1,npair)))
    fillp = nbr_ptr(1:nb)
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       e = fillp(ia); fillp(ia) = fillp(ia) + 1
       nbr_blk(e) = ib; nbr_pair(e) = i; nbr_isA(e) = .true.; eA(i) = e
       e = fillp(ib); fillp(ib) = fillp(ib) + 1
       nbr_blk(e) = ia; nbr_pair(e) = i; nbr_isA(e) = .false.; eB(i) = e
    end do
    do i = 1, npair
       nbr_rev(eA(i)) = eB(i); nbr_rev(eB(i)) = eA(i)
    end do
    deallocate(cnt, fillp, eA, eB)

    ! ---- dense filtered blocks of every pair (independent => parallel) ----
    allocate(qp(max(1,npair)))
    nzp = 0; nzq = 0
    if (present(tb)) then
       nth = max(1, min(tb, max(1, npair)))
    else
       nth = min(icgn_thread_budget(), max(1, npair))
       call icgn_flatten_nesting(saved_levels)
    end if
    !$omp parallel do default(shared) private(i, ia, ib, jq, kq, k, ra, rb, v, g, tmp, own_g, &
    !$omp&   np_a, np_b, nq_a, nq_b, nra, nrb) reduction(+:nzp, nzq) num_threads(nth) &
    !$omp&   schedule(dynamic,1) if(nth > 1)
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       np_a = bq(ia)%np; nq_a = bq(ia)%nq; np_b = bq(ib)%np; nq_b = bq(ib)%nq
       nra = blk(ia)%nrow; nrb = blk(ib)%nrow
       own_g = .false.
       if (associated(pairs(i)%g)) then
          g => pairs(i)%g                      ! full projection kept (outer path)
       else
          allocate(g(nra, nrb)); own_g = .true.
          call icgn_project_sparse(pairs(i)%h, blk(ia)%S_full, &
               blk(ib)%S_full, g)
       end if

       ! G(P_A, Q_B)
       allocate(qp(i)%gab(np_a, nq_b))
       do kq = 1, nq_b
          k = bq(ib)%qidx(kq)
          do jq = 1, np_a
             ra = blk(ia)%retained_idx(jq)
             v = g(ra, k)
             if (abs(v) < thr) then
                qp(i)%gab(jq, kq) = czero
             else if (.not. icgn_coupling_survives(v, blk(ia)%evals_full(ra), &
                  blk(ib)%evals_full(k), eps_u, elo_u, ehi_u)) then
                qp(i)%gab(jq, kq) = czero
             else
                qp(i)%gab(jq, kq) = v; nzp = nzp + 1
             end if
          end do
       end do
       ! G(P_B, Q_A) = conj of g(Q_A, P_B)
       allocate(qp(i)%gba(np_b, nq_a))
       do kq = 1, nq_a
          k = bq(ia)%qidx(kq)
          do jq = 1, np_b
             rb = blk(ib)%retained_idx(jq)
             v = g(k, rb)
             if (abs(v) < thr) then
                qp(i)%gba(jq, kq) = czero
             else if (.not. icgn_coupling_survives(v, blk(ia)%evals_full(k), &
                  blk(ib)%evals_full(rb), eps_u, elo_u, ehi_u)) then
                qp(i)%gba(jq, kq) = czero
             else
                qp(i)%gba(jq, kq) = conjg(v); nzp = nzp + 1
             end if
          end do
       end do
       ! W(Q_A, Q_B)
       if (need_w) then
          allocate(qp(i)%wab(nq_a, nq_b))
          do kq = 1, nq_b
             k = bq(ib)%qidx(kq)
             do jq = 1, nq_a
                ra = bq(ia)%qidx(jq)
                v = g(ra, k)
                if (abs(v) < thr) then
                   qp(i)%wab(jq, kq) = czero
                else if (.not. icgn_coupling_survives(v, blk(ia)%evals_full(ra), &
                     blk(ib)%evals_full(k), eps_u, elo_u, ehi_u)) then
                   qp(i)%wab(jq, kq) = czero
                else
                   qp(i)%wab(jq, kq) = v; nzq = nzq + 2
                end if
             end do
          end do
       end if
       if (own_g) deallocate(g)
    end do
    !$omp end parallel do
    if (.not. present(tb)) call icgn_restore_nesting(saved_levels)
    npq_nz = nzp; nqq_nz = nzq
  end subroutine icgn_build_neumann_blocks

  subroutine icgn_free_neumann_blocks(bq, qp, npair)
    type(ICGNQBlock), allocatable, intent(inout) :: bq(:)
    type(ICGNQPair), allocatable, intent(inout) :: qp(:)
    integer, intent(in) :: npair
    integer :: i
    if (allocated(qp)) then
       do i = 1, npair
          if (associated(qp(i)%gab)) deallocate(qp(i)%gab)
          if (associated(qp(i)%gba)) deallocate(qp(i)%gba)
          if (associated(qp(i)%wab)) deallocate(qp(i)%wab)
       end do
       deallocate(qp)
    end if
    if (allocated(bq)) then
       do i = 1, size(bq)
          if (associated(bq(i)%qidx)) deallocate(bq(i)%qidx)
          if (associated(bq(i)%res))  deallocate(bq(i)%res)
       end do
       deallocate(bq)
    end if
  end subroutine icgn_free_neumann_blocks

  ! Single cleanup path for failed icgn_prepare stages. Local allocatables are
  ! automatically released on procedure exit; these helpers release pointer-
  ! owned arrays and matrix components that Fortran will not deallocate for us.
  subroutine icgn_prepare_cleanup(upt, pairs, bq, qp, sig)
    type(OUPT), intent(inout) :: upt
    type(ICGNPair), allocatable, intent(inout) :: pairs(:)
    type(ICGNQBlock), allocatable, intent(inout) :: bq(:)
    type(ICGNQPair), allocatable, intent(inout) :: qp(:)
    type(ICGNSigRow), allocatable, intent(inout) :: sig(:)
    integer :: nqp

    if (allocated(pairs)) call destroy_pairs(pairs)
    nqp = 0
    if (allocated(qp)) nqp = size(qp)
    if (allocated(bq) .or. allocated(qp)) call icgn_free_neumann_blocks(bq, qp, nqp)
    if (allocated(sig)) call icgn_free_sigrows(sig)
    call icgn_free_block_eigensystems(upt)
    call icgn_clear(upt)
  end subroutine icgn_prepare_cleanup

  ! Block-sparse Neumann Sigma in P-P space.  Only upper block triangle
  ! (a2 >= a) is stored; lower follows from Hermiticity at CSR merge time.
  subroutine icgn_sigma_blocksparse(nb, order, maxdim, bq, qp, &
       nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, sig, tb)
    integer, intent(in) :: nb, order, maxdim
    integer, intent(in) :: nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, intent(in) :: nbr_isA(:)
    type(ICGNQBlock), intent(in) :: bq(:)
    type(ICGNQPair), intent(in) :: qp(:)
    type(ICGNSigRow), intent(inout) :: sig(:)
    integer, intent(in), optional :: tb
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)
    real(dp), parameter :: thr = 1.0e-14_dp
    integer :: a, b, b2, a2, e, f, ord, kq, nth, saved_levels, slot
    integer :: npa, nqb, nqb2, npa2
    type(ICGNMat), allocatable :: Acur(:), Anext(:)
    complex(dp), pointer :: Gm(:,:), Wm(:,:)
    character(len=1) :: opw
    real(dp) :: amax
    logical :: any_cur

    if (present(tb)) then
       nth = max(1, min(tb, nb))
    else
       nth = min(icgn_thread_budget(), max(1, nb))
       nth = icgn_outer_threads(nb, nth, maxdim)
       call icgn_flatten_nesting(saved_levels)
    end if
    !$omp parallel default(shared) num_threads(nth) if(nth > 1) &
    !$omp& private(a, b, b2, a2, e, f, ord, kq, npa, nqb, nqb2, npa2, Gm, Wm, opw, amax, &
    !$omp&         any_cur, Acur, Anext, slot)
    allocate(Acur(nb), Anext(nb))
    !$omp do schedule(dynamic,1)
    do a = 1, nb
       npa = bq(a)%np
       if (npa == 0) cycle
       ! order 0:  A_0(a,b) = G(P_a,Q_b) R_b
       do e = nbr_ptr(a), nbr_ptr(a+1)-1
          b = nbr_blk(e); nqb = bq(b)%nq
          if (nqb == 0) cycle
          Gm => icgn_qg(qp, nbr_pair, nbr_isA, e)
          allocate(Acur(b)%m(npa, nqb))
          do kq = 1, nqb
             Acur(b)%m(:, kq) = Gm(:, kq) * bq(b)%res(kq)
          end do
       end do

       do ord = 0, order
          if (ord > 0) then
             ! A_next(a,b2) = sum_b A(a,b) W(b,b2) R_b2
             do b = 1, nb
                if (.not. associated(Acur(b)%m)) cycle
                do f = nbr_ptr(b), nbr_ptr(b+1)-1
                   b2 = nbr_blk(f); nqb2 = bq(b2)%nq
                   if (nqb2 == 0) cycle
                   if (.not. associated(Anext(b2)%m)) then
                      allocate(Anext(b2)%m(npa, nqb2)); Anext(b2)%m = (0.0_dp, 0.0_dp)
                   end if
                   Wm => qp(nbr_pair(f))%wab
                   if (nbr_isA(f)) then; opw = 'N'; else; opw = 'C'; end if
                   call icgn_zgemm_gen('N', opw, cone, Acur(b)%m, Wm, cone, Anext(b2)%m)
                end do
             end do
             do b = 1, nb
                if (associated(Acur(b)%m)) deallocate(Acur(b)%m)
             end do
             amax = 0.0_dp; any_cur = .false.
             do b2 = 1, nb
                if (.not. associated(Anext(b2)%m)) cycle
                do kq = 1, bq(b2)%nq
                   Anext(b2)%m(:, kq) = Anext(b2)%m(:, kq) * bq(b2)%res(kq)
                   amax = max(amax, maxval(abs(Anext(b2)%m(:, kq))))
                end do
                Acur(b2)%m => Anext(b2)%m
                nullify(Anext(b2)%m)
                any_cur = .true.
             end do
             if (.not. any_cur .or. amax < thr) exit
          end if

          ! close the chain into block-sparse Sigma(P_a, P_a2), a2 >= a
          do b = 1, nb
             if (.not. associated(Acur(b)%m)) cycle
             do f = nbr_ptr(b), nbr_ptr(b+1)-1
                a2 = nbr_blk(f)
                if (a2 < a) cycle
                npa2 = bq(a2)%np
                if (npa2 == 0) cycle
                Gm => icgn_qg(qp, nbr_pair, nbr_isA, nbr_rev(f))     ! G(P_a2, Q_b)
                slot = icgn_sig_slot(sig, a, a2, npa, npa2)
                call icgn_zgemm_gen('N', 'C', cone, Acur(b)%m, Gm, cone, sig(a)%mat(slot)%m)
             end do
          end do
       end do

       do b = 1, nb
          if (associated(Acur(b)%m))  deallocate(Acur(b)%m)
          if (associated(Anext(b)%m)) deallocate(Anext(b)%m)
       end do
    end do
    !$omp end do
    deallocate(Acur, Anext)
    !$omp end parallel
    if (.not. present(tb)) call icgn_restore_nesting(saved_levels)
  end subroutine icgn_sigma_blocksparse

  ! Merge H_PP CSR with block-sparse Sigma into a new full CSR (sorted columns).
  ! Lower triangle of Sigma is filled by Hermiticity; diagonal blocks are
  ! Hermitianized in place.
  subroutine icgn_csr_add_sigma(ham, roff_blk, sig, nred, ierr)
    type(CSR), intent(inout) :: ham
    integer, intent(in) :: roff_blk(:), nred
    type(ICGNSigRow), intent(in) :: sig(:)
    integer, intent(out) :: ierr
    integer :: a, t, a2, i, j, r, c, k, nnew, pos
    integer(8) :: ne
    integer, allocatable :: ii(:), jj(:), ord(:)
    complex(dp), allocatable :: vv(:)
    integer, allocatable :: rowcount(:), next(:)
    type(CSR) :: compact_ham
    real(dp), parameter :: thr = 1.0e-14_dp
    complex(dp) :: val

    ierr = 0
    ! Upper bound of CSR entries after merge.  Diagonal blocks are Hermitianized
    ! before the threshold test, so count after combining (i,j)/(j,i).
    ne = ham%nnz
    do a = 1, size(sig)
       do t = 1, sig(a)%n
          a2 = sig(a)%blk(t)
          do i = 1, size(sig(a)%mat(t)%m, 1)
             do j = 1, size(sig(a)%mat(t)%m, 2)
                if (a2 == a) then
                   if (j < i) cycle
                   if (j > i) then
                      val = (sig(a)%mat(t)%m(i,j) + conjg(sig(a)%mat(t)%m(j,i))) / 2.0_dp
                      if (abs(val) <= thr) cycle
                      ne = ne + 2
                   else
                      if (abs(sig(a)%mat(t)%m(i,j)) <= thr) cycle
                      ne = ne + 1
                   end if
                else
                   if (abs(sig(a)%mat(t)%m(i,j)) <= thr) cycle
                   ne = ne + 2
                end if
             end do
          end do
       end do
    end do
    if (ne > int(huge(0), kind=8)) then
       ierr = 1
       return
    end if
    allocate(ii(max(1, int(ne))), jj(max(1, int(ne))), vv(max(1, int(ne))))
    nnew = 0
    ! Existing H_PP
    do r = 1, ham%nrow
       do k = ham%Mi(r), ham%Mi(r+1)-1
          nnew = nnew + 1
          ii(nnew) = r; jj(nnew) = ham%Mj(k); vv(nnew) = ham%M(k)
       end do
    end do
    ! Sigma blocks
    do a = 1, size(sig)
       do t = 1, sig(a)%n
          a2 = sig(a)%blk(t)
          do i = 1, size(sig(a)%mat(t)%m, 1)
             do j = 1, size(sig(a)%mat(t)%m, 2)
                val = sig(a)%mat(t)%m(i,j)
                if (a2 == a) then
                   ! Hermitianize diagonal block on the fly (keep upper + diag)
                   if (j < i) cycle
                   if (j > i) then
                      val = (val + conjg(sig(a)%mat(t)%m(j,i))) / 2.0_dp
                   else
                      val = cmplx(real(val,dp), 0.0_dp, kind=dp)
                   end if
                   if (abs(val) <= thr) cycle
                   nnew = nnew + 1
                   ii(nnew) = roff_blk(a)+i-1; jj(nnew) = roff_blk(a)+j-1; vv(nnew) = val
                   if (j > i) then
                      nnew = nnew + 1
                      ii(nnew) = roff_blk(a)+j-1; jj(nnew) = roff_blk(a)+i-1; vv(nnew) = conjg(val)
                   end if
                else
                   if (abs(val) <= thr) cycle
                   nnew = nnew + 1
                   ii(nnew) = roff_blk(a)+i-1; jj(nnew) = roff_blk(a2)+j-1; vv(nnew) = val
                   nnew = nnew + 1
                   ii(nnew) = roff_blk(a2)+j-1; jj(nnew) = roff_blk(a)+i-1; vv(nnew) = conjg(val)
                end if
             end do
          end do
       end do
    end do

    ! Counting sort by row, then insertion-sort by col within each row and merge.
    allocate(rowcount(nred), next(nred), ord(max(1,nnew)))
    rowcount = 0
    do k = 1, nnew
       rowcount(ii(k)) = rowcount(ii(k)) + 1
    end do
    next(1) = 1
    do r = 1, nred-1
       next(r+1) = next(r) + rowcount(r)
    end do
    do k = 1, nnew
       r = ii(k); ord(next(r)) = k; next(r) = next(r) + 1
    end do
    ! next(r) now past-the-end; start of row r is next(r)-rowcount(r)
    call destroy_matrix(ham)
    call create_matrix(ham, nred, nred, max(1, nnew))
    ham%sparse_fmt = 'F'
    ham%Mi(1) = 1
    pos = 0
    do r = 1, nred
       i = next(r) - rowcount(r)          ! first index in ord for this row
       j = next(r) - 1                    ! last index
       ! insertion sort ord(i:j) by column
       do k = i+1, j
          a = ord(k); t = k - 1
          do while (t >= i)
             if (jj(ord(t)) <= jj(a)) exit
             ord(t+1) = ord(t); t = t - 1
          end do
          ord(t+1) = a
       end do
       ! emit merged entries
       do k = i, j
          a = ord(k)
          if (pos >= ham%Mi(r)) then
             if (ham%Mj(pos) == jj(a)) then
                ham%M(pos) = ham%M(pos) + vv(a)
             else
                pos = pos + 1
                ham%Mj(pos) = jj(a)
                ham%M(pos)  = vv(a)
             end if
          else
             pos = pos + 1
             ham%Mj(pos) = jj(a)
             ham%M(pos)  = vv(a)
          end if
       end do
       ham%Mi(r+1) = pos + 1
    end do
    ham%nnz = pos

    ! The merge can collapse duplicate (row,column) entries, so the allocation
    ! made with nnew entries may have unused tail elements.  Some downstream
    ! routines use size(ham%M), not ham%nnz, as the stored-entry count; compact
    ! both value/index arrays to the actual CSR length before returning.
    if (size(ham%M) /= pos .or. size(ham%Mj) /= pos) then
       call create_matrix(compact_ham, nred, nred, pos)
       if (pos > 0) then
          compact_ham%M(1:pos)  = ham%M(1:pos)
          compact_ham%Mj(1:pos) = ham%Mj(1:pos)
       end if
       compact_ham%Mi = ham%Mi
       compact_ham%nnz = pos
       compact_ham%sparse_fmt = ham%sparse_fmt
       call destroy_matrix(ham)
       ham = compact_ham
    end if

    deallocate(ii, jj, vv, ord, rowcount, next)
  end subroutine icgn_csr_add_sigma

  ! Sigma = sum_{ord=0..order} H_PQ R (W R)^ord H_QP on the block structure.
  ! Parallel over the P blocks a (each owns the block row Sigma(P_a, :)); only
  ! the blocks a2 >= a are accumulated, the rest follows from Hermiticity.
  ! Locate or create the Sigma block (a, a2) inside sig(a).  a2 >= a.
  function icgn_sig_slot(sig, a, a2, npa, npa2) result(slot)
    type(ICGNSigRow), intent(inout) :: sig(:)
    integer, intent(in) :: a, a2, npa, npa2
    integer :: slot, t
    integer, pointer :: blk_new(:)
    type(ICGNMat), pointer :: mat_new(:)
    slot = 0
    do t = 1, sig(a)%n
       if (sig(a)%blk(t) == a2) then; slot = t; return; end if
    end do
    ! grow by one
    slot = sig(a)%n + 1
    allocate(blk_new(slot), mat_new(slot))
    if (sig(a)%n > 0) then
       blk_new(1:sig(a)%n) = sig(a)%blk(1:sig(a)%n)
       do t = 1, sig(a)%n
          mat_new(t)%m => sig(a)%mat(t)%m
          nullify(sig(a)%mat(t)%m)
       end do
       deallocate(sig(a)%blk, sig(a)%mat)
    end if
    sig(a)%blk => blk_new; sig(a)%mat => mat_new
    sig(a)%n = slot
    sig(a)%blk(slot) = a2
    allocate(sig(a)%mat(slot)%m(npa, npa2))
    sig(a)%mat(slot)%m = (0.0_dp, 0.0_dp)
  end function icgn_sig_slot

  subroutine icgn_free_sigrows(sig)
    type(ICGNSigRow), allocatable, intent(inout) :: sig(:)
    integer :: a, t
    if (.not. allocated(sig)) return
    do a = 1, size(sig)
       do t = 1, sig(a)%n
          if (associated(sig(a)%mat(t)%m)) deallocate(sig(a)%mat(t)%m)
       end do
       if (associated(sig(a)%blk)) deallocate(sig(a)%blk)
       if (associated(sig(a)%mat)) deallocate(sig(a)%mat)
       sig(a)%n = 0
    end do
    deallocate(sig)
  end subroutine icgn_free_sigrows

  ! G(P_a, Q_b) of the adjacency entry e (block a, neighbour b)
  function icgn_qg(qp, nbr_pair, nbr_isA, e) result(p)
    type(ICGNQPair), intent(in) :: qp(:)
    integer, intent(in) :: nbr_pair(:), e
    logical, intent(in) :: nbr_isA(:)
    complex(dp), pointer :: p(:,:)
    if (nbr_isA(e)) then
       p => qp(nbr_pair(e))%gab
    else
       p => qp(nbr_pair(e))%gba
    end if
  end function icgn_qg

  ! out_X = [R] * sum_Y (W_{Y,X})^T in_Y   (the operator W^T = conj(W) of the power iteration)
  !   subroutine icgn_blk_wt_apply(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, qoff, vin, vout, apply_res)
  !     integer, intent(in) :: nb, nbr_ptr(:), nbr_blk(:), nbr_pair(:), qoff(:)
  !     logical, intent(in) :: nbr_isA(:), apply_res
  !     type(ICGNQBlock), intent(in) :: bq(:)
  !     type(ICGNQPair), intent(in) :: qp(:)
  !     complex(dp), intent(in) :: vin(:)
  !     complex(dp), intent(out) :: vout(:)
  !     complex(dp), parameter :: cone = (1.0_dp, 0.0_dp), czero = (0.0_dp, 0.0_dp)
  !     integer :: X, Y, f, nqx, nqy, oX, oY, maxnq, nth, q
  !     complex(dp), allocatable :: tmp(:), yv(:)
  !     external :: zgemv
  ! 
  !     maxnq = 1
  !     do X = 1, nb
  !        maxnq = max(maxnq, bq(X)%nq)
  !     end do
  !     nth = min(icgn_thread_budget(), max(1, nb))
  !     !$omp parallel default(shared) num_threads(nth) if(nth > 1) private(X, Y, f, nqx, nqy, oX, oY, q, tmp, yv)
  !     allocate(tmp(maxnq), yv(maxnq))
  !     !$omp do schedule(dynamic,1)
  !     do X = 1, nb
  !        nqx = bq(X)%nq
  !        if (nqx == 0) cycle
  !        oX = qoff(X)
  !        vout(oX+1:oX+nqx) = czero
  !        do f = nbr_ptr(X), nbr_ptr(X+1)-1
  !           Y = nbr_blk(f); nqy = bq(Y)%nq
  !           if (nqy == 0) cycle
  !           oY = qoff(Y)
  !           if (nbr_isA(f)) then
  !              ! X = A, Y = B, wab is (nq_X x nq_Y):  out_X += conj(wab) in_Y = conj(wab * conj(in_Y))
  !              do q = 1, nqy
  !                 tmp(q) = conjg(vin(oY+q))
  !              end do
  !              call zgemv('N', nqx, nqy, cone, qp(nbr_pair(f))%wab, max(1,nqx), tmp, 1, czero, yv, 1)
  !              do q = 1, nqx
  !                 vout(oX+q) = vout(oX+q) + conjg(yv(q))
  !              end do
  !           else
  !              ! X = B, Y = A, wab is (nq_Y x nq_X):  out_X += wab^T in_Y
  !              do q = 1, nqy
  !                 tmp(q) = vin(oY+q)
  !              end do
  !              call zgemv('T', nqy, nqx, cone, qp(nbr_pair(f))%wab, max(1,nqy), tmp, 1, czero, yv, 1)
  !              do q = 1, nqx
  !                 vout(oX+q) = vout(oX+q) + yv(q)
  !              end do
  !           end if
  !        end do
  !        if (apply_res) then
  !           do q = 1, nqx
  !              vout(oX+q) = vout(oX+q) * bq(X)%res(q)
  !           end do
  !        end if
  !     end do
  !     !$omp end do
  !     deallocate(tmp, yv)
  !     !$omp end parallel
  !     ! blocks without Q states keep zeros
  !     do X = 1, nb
  !        if (bq(X)%nq == 0) cycle
  !     end do
  !   end subroutine icgn_blk_wt_apply

  ! DISABLED Neumann convergence checker; retained for future reconsideration.
  !   ! Power iteration for ||T||_2, T = R W, on the block structure (same iteration
  !   ! as the CSR matvec form: u = R W^T v, w = R u, v <- W^T w).
  !   subroutine icgn_power_iteration_blk(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, &
  !        maxiter, tol, sigma_out, converged_out)
  !     integer, intent(in) :: nb, nbr_ptr(:), nbr_blk(:), nbr_pair(:), maxiter
  !     logical, intent(in) :: nbr_isA(:)
  !     type(ICGNQBlock), intent(in) :: bq(:)
  !     type(ICGNQPair), intent(in) :: qp(:)
  !     real(dp), intent(in) :: tol
  !     real(dp), intent(out) :: sigma_out
  !     logical, intent(out) :: converged_out
  !     integer :: k, X, Y, j, ntot, q, f, nqx, nqy
  !     integer, allocatable :: qoff(:)
  !     real(dp) :: nrm, diff2, rr, ri, frob2
  !     complex(dp), allocatable :: v(:), u(:), w(:), v_next(:)
  !     real(dp), allocatable :: res(:)
  !
  !     sigma_out = 0.0_dp; converged_out = .true.
  !     allocate(qoff(nb+1)); qoff(1) = 0
  !     do X = 1, nb
  !        qoff(X+1) = qoff(X) + bq(X)%nq
  !     end do
  !     ntot = qoff(nb+1)
  !     if (ntot == 0) return
  !     allocate(v(ntot), u(ntot), w(ntot), v_next(ntot), res(ntot))
  !     do X = 1, nb
  !        do q = 1, bq(X)%nq
  !           res(qoff(X)+q) = bq(X)%res(q)
  !        end do
  !     end do
  !
  !     ! Power iteration yields a lower estimate of the spectral norm and can
  !     ! therefore miss non-convergence. Use ||R W||_F as a rigorous upper bound
  !     ! for ||R W||_2; each undirected Q-Q block contributes both block rows.
  !     frob2 = 0.0_dp
  !     do X = 1, nb
  !        nqx = bq(X)%nq
  !        if (nqx == 0) cycle
  !        do f = nbr_ptr(X), nbr_ptr(X+1)-1
  !           if (.not. nbr_isA(f)) cycle
  !           Y = nbr_blk(f); nqy = bq(Y)%nq
  !           if (nqy == 0) cycle
  !           if (.not. associated(qp(nbr_pair(f))%wab)) cycle
  !           do j = 1, nqy
  !              do q = 1, nqx
  !                 frob2 = frob2 + (res(qoff(X)+q)**2 + res(qoff(Y)+j)**2) * &
  !                      abs(qp(nbr_pair(f))%wab(q,j))**2
  !              end do
  !           end do
  !        end do
  !     end do
  !     if (frob2 >= 0.0_dp) then
  !        sigma_out = sqrt(frob2)
  !     else
  !        ! Non-finite arithmetic must never turn into a zero/false-safe estimate.
  !        sigma_out = huge(1.0_dp)
  !     end if
  !
  !     do j = 1, ntot
  !        call random_number(rr); call random_number(ri)
  !        v(j) = cmplx(rr - 0.5_dp, ri - 0.5_dp, kind=dp)
  !     end do
  !     nrm = sqrt(sum(real(v,dp)**2 + aimag(v)**2))
  !     if (nrm < 1.0e-14_dp) then
  !        converged_out = (sigma_out == 0.0_dp)
  !        deallocate(v, u, w, v_next, res, qoff); return
  !     end if
  !     v = v / nrm
  !     converged_out = .false.
  !     do k = 1, maxiter
  !        call icgn_blk_wt_apply(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, qoff, v, u, .true.)
  !        do j = 1, ntot
  !           w(j) = u(j) * res(j)
  !        end do
  !        call icgn_blk_wt_apply(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, qoff, w, v_next, .false.)
  !        nrm = sqrt(sum(real(v_next,dp)**2 + aimag(v_next)**2))
  !        if (nrm < 1.0e-14_dp) then
  !           ! Keep the Frobenius upper bound computed above; never replace it
  !           ! with a zero lower estimate from a degenerate iterate. A zero bound
  !           ! means W is exactly zero in this representation, so convergence is
  !           ! trivially certified in that special case.
  !           converged_out = (sigma_out == 0.0_dp)
  !           deallocate(v, u, w, v_next, res, qoff); return
  !        end if
  !        v_next = v_next / nrm
  !        diff2 = sum(real(v_next - v, dp)**2 + aimag(v_next - v)**2)
  !        v = v_next
  !        if (sqrt(diff2) < tol) then
  !           converged_out = .true.
  !           exit
  !        end if
  !     end do
  !     if (icgn_timing_on()) write(*,'(a,i0,a,l1)') '[ICGN-TIME] power iteration sweeps: ', &
  !          min(k, maxiter), ' converged=', converged_out
  !     ! Report the Frobenius upper bound instead of the iteration's lower
  !     ! estimate, so a value < 1 is a valid sufficient convergence test.
  !     deallocate(v, u, w, v_next, res, qoff)
  !   end subroutine icgn_power_iteration_blk
  !
  !   ! Warnings about the Neumann convergence bound are routed through the host logger.
  !   subroutine icgn_report_convergence(upt)
  !     type(OUPT), intent(in) :: upt
  !     character(len=512) :: message
  !     if (.not. upt%icgn_pi_converged) then
  !        write(message,'(a,i0,a)') &
  !             'WARNING: power iteration did not converge in ', &
  !             upt%icgn_pi_maxiter, ' sweeps; using the Frobenius upper bound for safety.'
  !        call icgn_log_important(upt, trim(message))
  !     end if
  !     if (upt%icgn_sigma_T2 >= 1.0_dp) then
  !        write(message,'(a,f10.6,a)') &
  !             'WARNING: Frobenius upper bound for ||T||_2 = ', upt%icgn_sigma_T2, &
  !             ' >= 1; convergence is not certified at this E0.'
  !        call icgn_log_important(upt, trim(message))
  !        if (upt%icgn_selfenergy_order == 0) call icgn_log_important(upt, &
  !             'NOTE: order=0 does not use the Q-Q coupling W; this warning is informational only.')
  !     else if (upt%verbose > 0) then
  !        write(message,'(a,f10.6,a,l1)') &
  !             'Frobenius upper bound for ||T||_2 = ', upt%icgn_sigma_T2, &
  !             ', power iteration converged: ', upt%icgn_pi_converged
  !        call icgn_log_important(upt, trim(message))
  !        if (upt%icgn_selfenergy_order == 0) call icgn_log_important(upt, &
  !             'NOTE: order=0 ignores W; the norm bound is for diagnostics only.')
  !     end if
  !   end subroutine icgn_report_convergence

! ===========================================================================
! SECTION D — Dense eigensolvers & BLAS helpers
! ===========================================================================

  subroutine dense_eigh(a, w, ierr)
    complex(dp), intent(inout) :: a(:,:)
    real(dp), intent(out) :: w(:)
    integer, intent(out) :: ierr
    integer :: n, info, lwork, lrwork, liwork
    complex(dp) :: workq(1)
    real(dp) :: rworkq(1)
    integer :: iworkq(1)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    n=size(w)
    if (n == 0) then; ierr=0; return; end if
    call zheevd('V','U',n,a,n,w,workq,-1,rworkq,-1,iworkq,-1,info)
    if (info /= 0) then; ierr=10; return; end if
    lwork=max(1,int(real(workq(1)))); lrwork=max(1,int(rworkq(1))); liwork=max(1,iworkq(1))
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    call zheevd('V','U',n,a,n,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    deallocate(work,rwork,iwork); ierr=info
  end subroutine dense_eigh

  ! Expand keep_mask at the energy-window edges when consecutive (sorted)
  ! eigenvalues differ by at most delta.  Only computed eigenvectors participate
  ! (zero S_full columns from the ZHEEVR fast path are ignored).
  subroutine icgn_expand_keep_boundary(blk, keep, elo, ehi, delta)
    type(ICGNBlock), intent(in) :: blk
    logical, intent(inout) :: keep(:)
    real(dp), intent(in) :: elo, ehi, delta
    integer :: nn, ncomp, j, k, p, tmp
    integer, allocatable :: ord(:)
    real(dp) :: ej

    if (.not. associated(blk%evals_full) .or. .not. associated(blk%S_full)) return
    nn = blk%nrow
    if (nn <= 0) return

    allocate(ord(nn))
    ncomp = 0
    do j = 1, nn
       if (j > size(keep)) exit
       if (maxval(abs(blk%S_full(:, j))) <= 0.0_dp) cycle
       ncomp = ncomp + 1
       ord(ncomp) = j
    end do
    if (ncomp <= 1) then
       deallocate(ord)
       return
    end if

    ! Insertion sort by eigenvalue (ncomp is modest per block).
    do k = 2, ncomp
       tmp = ord(k)
       ej = blk%evals_full(tmp)
       p = k - 1
       do while (p >= 1)
          if (blk%evals_full(ord(p)) <= ej) exit
          ord(p + 1) = ord(p)
          p = p - 1
       end do
       ord(p + 1) = tmp
    end do

    ! High edge: from last kept state with E <= ehi, walk upward.
    p = 0
    do k = 1, ncomp
       j = ord(k)
       if (keep(j) .and. blk%evals_full(j) <= ehi) p = k
    end do
    do while (p >= 1 .and. p < ncomp)
       if (abs(blk%evals_full(ord(p+1)) - blk%evals_full(ord(p))) > delta) exit
       keep(ord(p+1)) = .true.
       p = p + 1
    end do

    ! Low edge: from first kept state with E >= elo, walk downward.
    p = 0
    do k = ncomp, 1, -1
       j = ord(k)
       if (keep(j) .and. blk%evals_full(j) >= elo) p = k
    end do
    do while (p >= 2)
       if (abs(blk%evals_full(ord(p)) - blk%evals_full(ord(p-1))) > delta) exit
       keep(ord(p-1)) = .true.
       p = p - 1
    end do

    deallocate(ord)
  end subroutine icgn_expand_keep_boundary

  ! Partial Hermitian eigenproblem on energy window (VL, VU] via ZHEEVR (MRRR).
  ! ABSTOL = 0 selects the relative-accuracy path (EPS*|T|), matching the
  ! high-accuracy regime of ZHEEVD for the eigenvalues that are requested.
  ! On failure the caller should fall back to dense_eigh.
  ! On success: m eigenvalues in w(1:m), eigenvectors in z(:,1:m).
  subroutine dense_eigh_window(a, vl, vu, w, z, m, ierr)
    complex(dp), intent(inout) :: a(:,:)
    real(dp), intent(in) :: vl, vu
    real(dp), intent(out) :: w(:)
    complex(dp), intent(out) :: z(:,:)
    integer, intent(out) :: m, ierr
    integer :: n, info, lwork, lrwork, liwork
    complex(dp) :: workq(1)
    real(dp) :: rworkq(1)
    integer :: iworkq(1)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:), isuppz(:)
    real(dp), parameter :: abstol = 0.0_dp
    n = size(a, 1)
    m = 0; ierr = 0; w = 0.0_dp; z = (0.0_dp, 0.0_dp)
    if (n == 0) return
    if (.not. (vu > vl)) then; ierr = 11; return; end if
    allocate(isuppz(2 * max(1, n)))
    call zheevr('V', 'V', 'U', n, a, n, vl, vu, 1, n, abstol, m, w, z, n, &
         isuppz, workq, -1, rworkq, -1, iworkq, -1, info)
    if (info /= 0) then; deallocate(isuppz); ierr = 10; return; end if
    lwork  = max(1, int(real(workq(1), kind=dp)))
    lrwork = max(1, int(rworkq(1)))
    liwork = max(1, iworkq(1))
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevr('V', 'V', 'U', n, a, n, vl, vu, 1, n, abstol, m, w, z, n, &
         isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    deallocate(work, rwork, iwork, isuppz)
    ierr = info
  end subroutine dense_eigh_window

  ! Full spectrum (ZHEEVD). Used when outer / Neumann need discarded states.
  subroutine coarse_eigh(a, w, ierr)
    complex(dp), intent(inout) :: a(:,:)
    real(dp), intent(out) :: w(:)
    integer, intent(out) :: ierr
    call dense_eigh(a, w, ierr)
  end subroutine coarse_eigh

  logical function stored_entry(fmt, r, c)
    character(1), intent(in) :: fmt
    integer, intent(in) :: r, c
    select case (fmt)
    case ('U', 'u')
       stored_entry = (r <= c)
    case ('L', 'l')
       stored_entry = (r >= c)
    case default
       ! The Hamiltonian builder may emit a complete matrix even when the
       ! sparse format is labelled F.  Process one triangle here; the
       ! projected matrix is completed by its Hermitian conjugate below.
       stored_entry = (r <= c)
    end select
  end function stored_entry

  ! C = op(A) * B with BLAS zgemm (opa = 'N', 'T', or 'C').  Much faster than the
  ! Fortran intrinsic matmul and thread-safe (runs single-threaded inside our
  ! OpenMP regions because nesting is flattened, threaded when called serially).
  subroutine icgn_zgemm(opa, A, B, C)
    character(len=1), intent(in) :: opa
    complex(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp), intent(inout) :: C(:,:)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp), czero = (0.0_dp, 0.0_dp)
    integer :: m, n, kk
    external :: zgemm
    m = size(C,1); n = size(C,2)
    if (m == 0 .or. n == 0) return
    if (opa == 'N') then; kk = size(A,2); else; kk = size(A,1); end if
    if (kk == 0) then; C = czero; return; end if
    call zgemm(opa, 'N', m, n, kk, cone, A, max(1,size(A,1)), B, max(1,size(B,1)), &
         czero, C, max(1,size(C,1)))
  end subroutine icgn_zgemm

  ! C = alpha*op(A)*op(B) + beta*C with BLAS zgemm (op = 'N','T','C').
  subroutine icgn_zgemm_gen(opa, opb, alpha, A, B, beta, C)
    character(len=1), intent(in) :: opa, opb
    complex(dp), intent(in) :: alpha, beta
    complex(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp), intent(inout) :: C(:,:)
    integer :: m, n, kk
    external :: zgemm
    m = size(C,1); n = size(C,2)
    if (m == 0 .or. n == 0) return
    if (opa == 'N') then; kk = size(A,2); else; kk = size(A,1); end if
    if (kk == 0) return
    call zgemm(opa, opb, m, n, kk, alpha, A, max(1,size(A,1)), B, max(1,size(B,1)), &
         beta, C, max(1,size(C,1)))
  end subroutine icgn_zgemm_gen

! ===========================================================================
! SECTION E — Logging, timing, OpenMP thread control
! ===========================================================================

   ! Statistics/progress message: printed only when ICGN_LOG_STATS is .true.
   ! (disabled by default so that a normal run is not flooded with ICGN details).
   subroutine icgn_log_progress(upt, message)
      type(OUPT), intent(in) :: upt
      character(*), intent(in) :: message
      if (ICGN_LOG_STATS) call icgn_log_important(upt, message)
   end subroutine icgn_log_progress

   ! Message that is always forwarded to the host application log (through the
   ! C++ callback upt_log_message): summaries, fall-backs, diagnostics.
   subroutine icgn_log_important(upt, message)
      type(OUPT), intent(in) :: upt
      character(*), intent(in) :: message
      character(kind=c_char), allocatable :: c_message(:)
      integer :: i, message_length

      message_length = len_trim(message)
      allocate(c_message(max(1, message_length)))
      do i = 1, message_length
          c_message(i) = message(i:i)
      end do
      call upt_log_message(c_message, int(message_length, c_int))
      deallocate(c_message)
   end subroutine icgn_log_important

   subroutine icgn_log_block(upt, mode, block_number, block_dimension, phase)
      type(OUPT), intent(in) :: upt
      character(*), intent(in) :: mode, phase
      integer, intent(in) :: block_number, block_dimension
      character(len=256) :: message
      write(message,'(a,a,a,i0,a,i0)') 'mode=', trim(mode), ', block=', block_number, &
             ', dimension=', block_dimension
      message = trim(message)//', phase='//trim(phase)
      call icgn_log_progress(upt, message)
   end subroutine icgn_log_block

   subroutine icgn_log_timing(upt, phase, seconds)
      type(OUPT), intent(in) :: upt
      character(*), intent(in) :: phase
      real(dp), intent(in) :: seconds
      integer(kind=8) :: total_seconds, hours, minutes, secs
      character(len=256) :: message

      total_seconds = int(max(0.0_dp, seconds), kind=8)
      hours = total_seconds / 3600_8
      minutes = mod(total_seconds, 3600_8) / 60_8
      secs = mod(total_seconds, 60_8)
      write(message,'(a,a,a,i0,a,i2.2,a,i2.2,a)') 'mode=icgn timing: ', trim(phase), &
           ' = ', hours, 'h ', minutes, 'm ', secs, 's'
      call icgn_log_progress(upt, trim(message))
   end subroutine icgn_log_timing

   ! One-shot summary of a finished preparation: original -> reduced dimension,
   ! number of blocks, cut fraction and (ICGN) Neumann-series norm estimate.
   ! This is the only ICGN statistic shown by default.
   subroutine icgn_log_info(upt, mode, original_dim, reduced_dim, &
     nblocks, cut_fraction, sigma_t2, pi_converged, has_norm)
    type(OUPT), intent(in) :: upt
    character(*), intent(in) :: mode
    integer, intent(in) :: original_dim, reduced_dim, nblocks
    real(dp), intent(in) :: cut_fraction, sigma_t2
    logical, intent(in) :: pi_converged, has_norm
      character(len=512) :: line

      write(line,'(a,a,a,i0,a,i0,a,i0,a,f8.4)') 'mode=', trim(mode), &
        ', dimension=', original_dim, ' -> ', reduced_dim, ', blocks=', nblocks, &
        ', retained fraction=', real(reduced_dim,dp)/max(1.0_dp,real(original_dim,dp))
      call icgn_log_important(upt, trim(line))
      write(line,'(a,f8.4)') 'cut fraction=', cut_fraction
      call icgn_log_important(upt, trim(line))
    ! if (has_norm) then
    !      write(line,'(a,f12.6,a,l1)') 'Neumann ||T||_2 upper bound (Frobenius)=', sigma_t2, &
    !      ', power iteration converged=', pi_converged
    !      call icgn_log_important(upt, trim(line))
    ! end if
  end subroutine icgn_log_info

  ! Wall-clock timing (portable, works with or without OpenMP).
  function icgn_wall() result(t)
    real(dp) :: t
    integer(8) :: cnt, rate
    call system_clock(cnt, rate)
    t = real(cnt, dp) / real(max(1_8, rate), dp)
  end function icgn_wall

  ! Prints "[ICGN-TIME] phase: x s" when env UPT_CG_TIMING is set (any value).
  subroutine icgn_report_time(phase, t0)
    character(*), intent(in) :: phase
    real(dp), intent(in) :: t0
    character(len=8) :: s
    integer :: ln, ios
    call get_environment_variable('UPT_CG_TIMING', s, length=ln, status=ios)
    if (ios == 0 .or. ios == -1) then
       if (ln > 0) write(*,'(a,a,a,f10.3,a)') '[ICGN-TIME] ', phase, ': ', icgn_wall() - t0, ' s'
    end if
  end subroutine icgn_report_time

  ! How sparse are the block-block couplings?  Counts the rows/columns of each
  ! H_ab that contain at least one nonzero (timing diagnostics only).
  subroutine icgn_report_boundary(pairs, npair)
    type(ICGNPair), intent(in) :: pairs(:)
    integer, intent(in) :: npair
    integer :: i, nr, nc, max_nr, max_nc
    integer(8) :: sum_nr, sum_nc, sum_ra, sum_cb
    sum_nr = 0; sum_nc = 0; sum_ra = 0; sum_cb = 0; max_nr = 0; max_nc = 0
    do i = 1, npair
       if (.not. associated(pairs(i)%h)) cycle
       nr = count(any(pairs(i)%h /= (0.0_dp, 0.0_dp), dim=2))
       nc = count(any(pairs(i)%h /= (0.0_dp, 0.0_dp), dim=1))
       sum_nr = sum_nr + nr; sum_nc = sum_nc + nc
       sum_ra = sum_ra + size(pairs(i)%h, 1); sum_cb = sum_cb + size(pairs(i)%h, 2)
       max_nr = max(max_nr, nr); max_nc = max(max_nc, nc)
    end do
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '[ICGN-TIME] boundary: npair=', npair, &
         ' nonzero_rows=', sum_nr, '/', sum_ra, ' nonzero_cols=', sum_nc, '/', sum_cb, &
         ' max_rows_per_pair=', max_nr, ' max_cols_per_pair=', max_nc
  end subroutine icgn_report_boundary

  ! True when env UPT_CG_TIMING is set (any non-empty value).
  function icgn_timing_on() result(on)
    logical :: on
    character(len=8) :: s
    integer :: ln, ios
    call get_environment_variable('UPT_CG_TIMING', s, length=ln, status=ios)
    on = (ln > 0 .and. (ios == 0 .or. ios == -1))
  end function icgn_timing_on

  function icgn_thread_budget() result(nth)
    !$ use omp_lib
    integer :: nth
    integer :: v, ios, ln, req, avail
    logical, save :: warned = .false.
    character(len=32) :: s
    nth = 1; req = 1; avail = 1
    !$ req = omp_get_max_threads(); avail = max(1, omp_get_num_procs())
    !$ nth = min(req, avail)
    !$ if (omp_in_parallel()) nth = 1
    ! Safety guard: never use more threads than the CPUs this process may use
    ! (affinity mask / cgroup / SLURM allocation); warn the user once.
    !$ if (req > avail .and. .not. warned .and. .not. omp_in_parallel()) then
    !$    warned = .true.
    !$    call icgn_warn_threads(req, avail)
    !$ end if
    call get_environment_variable('UPT_CG_THREADS', s, length=ln, status=ios)
    if (ios == 0 .and. ln > 0) then
       read(s(1:ln), *, iostat=ios) v
       if (ios == 0 .and. v >= 1) nth = min(nth, v)
    end if
    nth = max(1, nth)
  end function icgn_thread_budget

  ! Number of OUTER threads for a loop over nwork independent items.
  ! default_nth is the automatic choice of the caller.  Env UPT_CG_OUTER:
  !   0 -> no outer parallelism (items one after another; threaded BLAS/LAPACK
  !        inside uses the whole thread budget)
  !   1 -> always parallel over items (BLAS single-threaded inside)
  !   unset/other -> default_nth, except that LARGE items (dim >= 512) are run
  !        one after another so the threaded BLAS/LAPACK works on each of them
  !        (measured: blocks of 1280 are faster this way, blocks of 160 the other way)
  function icgn_outer_threads(nwork, default_nth, dim) result(nth)
    integer, intent(in) :: nwork, default_nth
    integer, intent(in), optional :: dim
    integer :: nth
    character(len=8) :: s
    integer :: ln, ios
    integer, parameter :: big_dim = 512
    nth = default_nth
    call get_environment_variable('UPT_CG_OUTER', s, length=ln, status=ios)
    if (ios == 0 .and. ln > 0) then
       if (s(1:1) == '0') then
          nth = 1
       else if (s(1:1) == '1') then
          nth = min(icgn_thread_budget(), max(1, nwork))
       else if (present(dim)) then
          if (dim >= big_dim) nth = 1
       end if
    else if (present(dim)) then
       if (dim >= big_dim) nth = 1
    end if
    nth = max(1, nth)
  end function icgn_outer_threads

  ! Printed once (rank 0 only) when more threads are requested than the CPUs
  ! this process is allowed to use.
  subroutine icgn_warn_threads(req, avail)
    integer, intent(in) :: req, avail
    character(len=16) :: rk
    integer :: ln, ios, rank
    call get_environment_variable('OMPI_COMM_WORLD_RANK', rk, length=ln, status=ios)
    if (ios == 0 .and. ln > 0) then
       read(rk(1:ln), *, iostat=ios) rank
       if (ios == 0 .and. rank /= 0) return
    end if
    write(*,'(a,i0,a,i0,a)') '  (icgn) WARNING: ', req, &
         ' OpenMP threads requested (OMP_NUM_THREADS / --cpus-per-task) but this process may use only ', &
         avail, ' CPU(s).'
    write(*,'(a,i0,a)') '  (icgn) Falling back to ', avail, &
         ' thread(s) for coarse-graining: more threads than CPUs only slows the run down.'
    write(*,'(a)') '  (icgn) Choose OMP_NUM_THREADS (or --cpus-per-task) <= the CPUs available to each MPI rank.'
  end subroutine icgn_warn_threads

  ! Allow only ONE active parallel level while our outer loop runs, so that a
  ! threaded BLAS/LAPACK called from inside it runs single-threaded instead of
  ! spawning nproc threads per outer thread (oversubscription).
  subroutine icgn_flatten_nesting(saved)
    !$ use omp_lib
    integer, intent(out) :: saved
    saved = 1
    !$ saved = omp_get_max_active_levels()
    !$ call omp_set_max_active_levels(1)
  end subroutine icgn_flatten_nesting

  subroutine icgn_restore_nesting(saved)
    !$ use omp_lib
    integer, intent(in) :: saved
    !$ call omp_set_max_active_levels(saved)
  end subroutine icgn_restore_nesting

end module coarse_grain
