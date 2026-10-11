! Coarse-grained tight-binding projection after Liu et al. (2022).
!
! ======================================================================
! OVERVIEW (read this first)
! ======================================================================
! This module is a "mid-layer" between the assembly of the orbital
! Hamiltonian H (sparse_ham, stored in upt%ham as CSR) and the eigensolvers
! (Lanczos, Jacobi-Davidson, LAPACK, external SLEPc). It replaces H by a much
! smaller reduced Hamiltonian H_red = Q^H H Q so that the eigensolver has less
! work to do, then maps the eigenvectors back to the orbital basis.
!
! ICGN is the single coarse-graining implementation; neumann_order < 0 gives the plain ICG/ICGN special case:
!
!   ICGN   (mode 1) original Liu method. The atoms are partitioned into blocks
!                 (METIS, or a built-in graph-BFS fallback). Each block
!                 Hamiltonian is diagonalized; the eigenstates whose energy
!                 lies in [icgn_emin, icgn_emax] are retained (columns of Q_I).
!                 Inter-block couplings are projected, T = H_IJ Q_J and
!                 H_red(I,J) = Q_I^H T, and small entries are dropped with
!                 the threshold epsilon.
!   ICG  (mode 2) "improved" ICGN — this mode is now folded into ICGN.
!                 Set neumann_order < 0 to skip the Neumann correction and
!                 recover the plain ICG behaviour.
!   ICGN (mode 2) improved ICGN with optional Neumann-series self-energy
!                 correction:
!                 Sigma = H_PQ (E0 - D)^-1 [W (E0 - D)^-1]^order H_QP
!                 (P = retained states, Q = discarded states) added to the
!                 reduced Hamiltonian.
!                 neumann_order < 0 disables the correction (plain ICG);
!                 neumann_order = 0, 1, 2, ... applies the series to that order.
!
! Life cycle of one calculation (driven by uptight.F90):
!   1. UPT_hamiltonian builds H, then calls icgn_legacy_prepare / icgn_legacy_prepare,
!      which fill upt%icgn_ham / icgn_ham (CSR) and the per-block data
!      (upt%icgn_blocks / icgn_blocks, type ICGNBlock: rows, eigenvalues,
!      Q = retained eigenvectors).
!   2. icgn_set_original remembers the orbital operator, icgn_swap_in makes
!      upt%ham point at the reduced operator. From now on every solver sees
!      only the reduced matrix and does not know about coarse-graining.
!   3. After the solve, icgn_finalize_eigenvectors lifts the eigenvectors back
!      (physical = Q * reduced, see icgn_legacy_lift / icgn_legacy_lift) and
!      icgn_swap_out restores the original operator.
!
! Parallelism (OpenMP, icgn_prepare only): block diagonalization and the
! inter-block pair projections run in parallel over blocks / pairs.  Thread
! budget = omp_get_max_threads() (capped by env UPT_CG_THREADS if set); nesting
! is flattened inside the regions so threaded BLAS/LAPACK never oversubscribes.
!
! Sparse storage convention used everywhere: CSR with Mi = row pointer
! (size nrow+1), Mj = column index, M = values; sparse_fmt is 'F' (full),
! 'U' (upper triangle) or 'L' (lower triangle).
! ======================================================================
module coarse_grain
   use, intrinsic :: iso_c_binding, only : c_int, c_double, c_char
  use precision, only : dp
  use upt_param, only : OUPT, ICGNBlock
  use sparse_matrix, only : CSR, create_matrix, destroy_matrix, write_csr
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
     ! the acquaintance path so the Neumann correction need not project again.
     complex(dp), dimension(:,:), pointer :: g => null()
  end type ICGNPair

  ! ---- Block-dense Neumann correction (see icgn_sigma_blockdense) ----------
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

contains

  ! Called by the orchestrator right after the orbital operator has been built
  ! and before the reduced one is prepared.
  subroutine icgn_set_original(upt)
    type(OUPT), intent(inout) :: upt
    ! write(*,'(A,L1,A,L1,A,I0,A,I0)') '[ICGN-TRACE] set_original ham M/Mi=', associated(upt%ham%M), '/', &
         ! associated(upt%ham%Mi), ' nrow=', upt%ham%nrow, ' nnz=', upt%ham%nnz
    ! write(*,'(A,L1,A,L1,A,I0,A,I0)') '[ICGN-TRACE] set_original U M/Mi=', associated(upt%U%M), '/', &
         ! associated(upt%U%Mi), ' nrow=', upt%U%nrow, ' nnz=', upt%U%nnz
    icgn_original_ham = upt%ham    ! value copy of struct header (pointer addresses, not heap data)
    icgn_original_u   = upt%U
    icgn_original_valid = .true.
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] set_original snapshot ham M/Mi=', associated(icgn_original_ham%M), '/', &
         ! associated(icgn_original_ham%Mi)
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] set_original snapshot U M/Mi=', associated(icgn_original_u%M), '/', &
         ! associated(icgn_original_u%Mi)
  end subroutine icgn_set_original

  ! Common coupling filter shared by ICGN and ICGN. A coupling V between two
  ! block eigenstates with energies ea, eb survives only if
  !   |V|^2 > epsilon * |ea - eb|        (strict: equality is filtered out)
  ! with epsilon >= 0. This is the product form of |V|^2/|ea - eb| > epsilon,
  ! so no division by a (nearly) vanishing energy difference occurs: for
  ! degenerate states the right-hand side is 0 and any non-zero coupling
  ! survives, without overflow or NaN. For epsilon = 0 the test reduces to
  ! |V|^2 > 0, i.e. only exact zeros are dropped.
  pure logical function icgn_coupling_survives(v, ea, eb, epsilon)
    complex(dp), intent(in) :: v
    real(dp), intent(in) :: ea, eb, epsilon
    real(dp) :: v2, de
    v2 = real(v,dp)**2 + aimag(v)**2
    icgn_coupling_survives = .false.
    if (.not. (v2 > 0.0_dp)) return            ! exact zero (or NaN): no coupling
    de = abs(ea - eb)
    if (de <= 0.0_dp) then
       icgn_coupling_survives = .true.           ! degenerate: |V|^2 > 0 = epsilon*0
    else
       icgn_coupling_survives = (v2 > max(epsilon, 0.0_dp)*de)
    end if
  end function icgn_coupling_survives

  subroutine icgn_legacy_configure(upt, enabled, nblocks, emin, emax, epsilon, metis_imbalance)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled
    integer, intent(in) :: nblocks
    real(dp), intent(in) :: emin, emax, epsilon, metis_imbalance
    call icgn_legacy_clear(upt)
    upt%icgn_enabled = enabled
    upt%icgn_num_blocks = nblocks
    upt%icgn_emin = emin
    upt%icgn_emax = emax
    upt%icgn_epsilon = epsilon
    upt%icgn_metis_imbalance = metis_imbalance
  end subroutine icgn_legacy_configure

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
       call icgn_legacy_lift(upt, reduced, physical)
    else
       physical = reduced
    end if
  end subroutine icgn_lift_active

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

    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] swap_in snapshot ham M/Mi=', associated(icgn_original_ham%M), '/', &
         ! associated(icgn_original_ham%Mi)
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] swap_in reduced ham M/Mi=', associated(reduced_ham%M), '/', &
         ! associated(reduced_ham%Mi)
    icgn_saved_ham => icgn_original_ham
    icgn_saved_u   => icgn_original_u
    icgn_saved_shift_init = shift_init
    icgn_saved_shift_end = shift_end
    icgn_saved_shift_init_mi = shift_init_Mi(id)
    icgn_saved_shift_end_mi = shift_end_Mi(id)
    icgn_saved_n_spin = upt%n_spin

    upt%ham = reduced_ham
    upt%U   = reduced_u
    ! The reduced operator already carries the spin degeneracy of the states it
    ! retained, so a solver must not generate degenerate partners again.
    upt%n_spin = 1

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

    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] swap_out saved_ham/saved_u=', associated(icgn_saved_ham), '/', &
         ! associated(icgn_saved_u)
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] swap_out snapshot ham M/Mi=', associated(icgn_original_ham%M), '/', &
         ! associated(icgn_original_ham%Mi)
    if (icgn_original_valid) then
      upt%ham = icgn_original_ham
      upt%U   = icgn_original_u
    end if
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] swap_out restored ham M/Mi=', associated(upt%ham%M), '/', &
         ! associated(upt%ham%Mi)

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
    if (nfull <= nred) then
      ! No larger basis to lift into: restore the operator and leave the
      ! vectors alone rather than lift into something of the wrong size.
      write(*,'(a,i0,a,i0)') &
           '(icgn) finalize: no orbital basis recorded (nred=', nred, &
           ', nfull=', nfull, '), keeping reduced vectors'
      call icgn_swap_out(upt)
      return
    end if

    write(*,'(a,i0,a,i0,a,i0)') '(icgn) finalize: lifting nred=', nred, &
         ' -> nfull=', nfull, ', num_ev=', num_ev

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
  end subroutine icgn_forget

  ! Drops every pointer into the reduced operator. Only for teardown: icgn_legacy_clear()
  ! runs at the start of every icgn_legacy_prepare(), so it must not touch the record of
  ! the original operator.
  subroutine icgn_destroy_operator(upt)
    type(OUPT), intent(inout), target :: upt
    integer :: i
    call icgn_forget(upt)
    ! write(*,"(A,L1)") "[ICGN-TRACE] clear ham associated=", associated(upt%icgn_ham%M)
    if (associated(upt%icgn_ham%M)) then
       ! write(*,"(A)") "[ICGN-TRACE] clear destroy ham BEGIN"
       call destroy_matrix(upt%icgn_ham)
       ! write(*,"(A)") "[ICGN-TRACE] clear destroy ham END"
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
       ! write(*,'(A)') '[ICGN-TRACE] clear: deallocated icgn_blocks'
    end if
    upt%icgn_ready = .false.
    upt%icgn_original_dim = 0
    upt%icgn_reduced_dim = 0
    upt%icgn_cut_fraction = 0.0_dp
  end subroutine icgn_destroy_operator

  subroutine icgn_legacy_clear(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    if (associated(upt%icgn_ham%M)) call destroy_matrix(upt%icgn_ham)
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
       ! write(*,"(A)") "[ICGN-TRACE] clear blocks deallocated"
    end if
    upt%icgn_ready = .false.
    upt%icgn_original_dim = 0
    upt%icgn_reduced_dim = 0
    upt%icgn_cut_fraction = 0.0_dp
  end subroutine icgn_legacy_clear

  subroutine icgn_legacy_prepare(upt, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(out) :: ierr
    integer :: n, na, nb, i, j, k, p, q, status, nedge, maxedge
    integer :: r, c, br, bc, lr, lc, npair, total_ret, pos
    integer, allocatable :: atom_of(:), local_of(:), label(:), row_of(:)
    integer, allocatable :: counts(:), cursor(:), offsets(:), bsize(:)
    integer(c_int), allocatable :: xadj(:), adjncy(:), vwgt(:), adjwgt(:), part(:)
    real(dp), allocatable :: edge_weight(:)
    real(dp) :: max_weight, all_weight, cut_weight, workspace_mib
    integer(kind=8) :: timer_start, timer_end, timer_rate
    real(dp) :: timer_seconds
    character(len=256) :: line
    type(ICGNPair), allocatable :: pairs(:)

    ierr = 0
    ! write(*,'(A)') '[ICGN-TRACE] prepare ENTER'
    call icgn_legacy_clear(upt)
    if (.not.upt%icgn_enabled) return
    call icgn_log_progress(upt, 'mode=icgn preparation started')
    if (num_procs /= 1) then
       ierr = 1; write(*,*) '(icgn) MPI runs are not supported'; return
    end if
    na = upt%basis%n_basis
    n = upt%ham%nrow
    if (na < 1 .or. upt%icgn_num_blocks < 1 .or. upt%icgn_num_blocks > na) then
       ierr = 2; write(*,*) '(icgn) invalid number of blocks'; return
    end if
    if (upt%icgn_emin >= upt%icgn_emax .or. upt%icgn_metis_imbalance < 0.0_dp) then
       ierr = 3; write(*,*) '(icgn) invalid energy window or metis_imbalance'; return
    end if
    if (.not.associated(upt%ham%M)) then
       ierr = 4; write(*,*) '(icgn) Hamiltonian is not initialized'; return
    end if
    if (upt%icgn_num_blocks == 1 .and. upt%verbose > 0) then
       write(*,*) '(icgn) one block selected; energy window controls rank'
    end if

    ! write(*,'(A,I0,A,I0)') '[ICGN-TRACE] prepare before mapping allocate n=', n, ' na=', na
    allocate(atom_of(n), local_of(n), offsets(na+1), bsize(na))
    ! write(*,'(A)') '[ICGN-TRACE] prepare after mapping allocate'
    pos = 1
    ! write(*,'(A)') '[ICGN-TRACE] prepare before mapping loop'
    do i = 1, na
       offsets(i) = pos
       bsize(i) = upt%n_spin * upt%basis%n_st(i)
       do j = 1, bsize(i)
          atom_of(pos) = i; local_of(pos) = j; pos = pos + 1
       end do
    end do
    offsets(na+1) = pos
    if (pos-1 /= n) then
       ierr = 5; write(*,*) '(icgn) atom/orbital mapping is inconsistent'; return
    end if

    ! Parallel graph edges are intentional: their summed weights are the
    ! Frobenius norm squared of an atom-to-atom Hamiltonian block.
    maxedge = max(1, upt%ham%nnz)
    allocate(edge_weight(maxedge), counts(na), label(na))
    nedge = 0; max_weight = 0.0_dp; counts = 0
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          nedge = nedge + 1
          if (nedge > maxedge) then
             ierr = 6; return
          end if
          edge_weight(nedge) = abs(upt%ham%M(k))**2
          counts(atom_of(r)) = counts(atom_of(r)) + 1
          counts(atom_of(c)) = counts(atom_of(c)) + 1
          max_weight = max(max_weight, edge_weight(nedge))
       end do
    end do
    if (nedge == 0 .or. max_weight == 0.0_dp) then
       ierr = 7; write(*,*) '(icgn) atom graph has no couplings'; return
    end if
    write(line,'(a,i0,a,a1,a,i0,a,i0)') 'mode=icgn graph scan complete, orbitals=', n, &
         ', sparse_format=', upt%ham%sparse_fmt, ', stored_nnz=', upt%ham%nnz, &
         ', interblock_entries=', nedge
    call icgn_log_progress(upt, trim(line))
    allocate(xadj(na+1), cursor(na), adjncy(2*nedge), adjwgt(2*nedge), vwgt(na), part(na))
    xadj(1) = 0_c_int
    do i = 1, na
       xadj(i+1) = xadj(i) + int(counts(i), c_int)
       cursor(i) = int(xadj(i)) + 1
       vwgt(i) = int(bsize(i), c_int)
    end do
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          br = atom_of(r); bc = atom_of(c)
          p = max(1, nint(abs(upt%ham%M(k))**2 / max_weight * 1000000.0_dp))
          adjncy(cursor(br)) = int(bc-1, c_int); adjwgt(cursor(br)) = int(p, c_int); cursor(br)=cursor(br)+1
          adjncy(cursor(bc)) = int(br-1, c_int); adjwgt(cursor(bc)) = int(p, c_int); cursor(bc)=cursor(bc)+1
       end do
    end do
    status = icgn_metis_partition(int(na,c_int), xadj, adjncy, vwgt, adjwgt, &
         int(upt%icgn_num_blocks,c_int), int(nint(1000.0_dp*upt%icgn_metis_imbalance),c_int), 42_c_int, part)
    if (status /= 0) then       ! METIS not available — fallback to connectivity-aware partition
       call icgn_log_progress(upt, 'METIS unavailable, using built-in graph-BFS fallback partitioning')
       call icgn_graph_partition(na, upt%icgn_num_blocks, vwgt, xadj, adjncy, adjwgt, part)
    else
       call icgn_log_progress(upt, 'METIS partitioning done')
    end if
    do i = 1, na
       label(i) = int(part(i)) + 1
    end do
    ! write(*,'(A)') '[ICGN-TRACE] prepare: graph partition complete'
    call icgn_log_progress(upt, 'mode=icgn graph partition complete')

    all_weight = 0.0_dp; cut_weight = 0.0_dp
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_weight = all_weight + abs(upt%ham%M(k))**2
          if (label(atom_of(r)) /= label(atom_of(c))) cut_weight = cut_weight + abs(upt%ham%M(k))**2
       end do
    end do
    upt%icgn_cut_fraction = cut_weight / all_weight

    deallocate(counts)
    allocate(upt%icgn_blocks(upt%icgn_num_blocks), counts(upt%icgn_num_blocks))
    counts = 0
    do i = 1, na
       counts(label(i)) = counts(label(i)) + bsize(i)
    end do
    do i = 1, upt%icgn_num_blocks
       upt%icgn_blocks(i)%nrow = counts(i)
       allocate(upt%icgn_blocks(i)%rows(counts(i)))
    end do
    ! ZHEEVD needs the dense matrix plus work arrays.  This is deliberately
    ! only a prediction: the allocation itself remains inside diagonalize_block.
    workspace_mib = 16.0_dp * real(maxval(counts),dp)**2 / (1024.0_dp**2)
    if (upt%verbose > 0) write(*,'(a,i0,a,f10.2,a)') '(icgn) largest dense block ', &
         maxval(counts), ', matrix workspace at least ', workspace_mib, ' MiB'
    cursor = 0
    do i = 1, na
       br = label(i)
       do j = offsets(i), offsets(i+1)-1
          cursor(br) = cursor(br) + 1
          upt%icgn_blocks(br)%rows(cursor(br)) = j
       end do
    end do

    allocate(row_of(n)); row_of = 0
    call system_clock(timer_start, timer_rate)
    do i = 1, upt%icgn_num_blocks
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
       call diagonalize_block(upt, i, atom_of, row_of, ierr)
       if (ierr /= 0) return
       row_of(upt%icgn_blocks(i)%rows) = 0
    end do
    call system_clock(timer_end)
    timer_seconds = real(timer_end-timer_start,dp) / real(max(1_8,timer_rate),dp)
    call icgn_log_timing(upt, 'block diagonalization', timer_seconds)
    do i = 1, upt%icgn_num_blocks
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
    end do
    total_ret = 0
    do i = 1, upt%icgn_num_blocks
       total_ret = total_ret + upt%icgn_blocks(i)%nret
    end do
    if (total_ret == 0) then
       ierr = 9; write(*,*) '(icgn) no block state retained, please expand the energy window'; return
    end if
    ! No check needed - we will compute ALL eigenvalues of reduced matrix
    upt%icgn_original_dim = n; upt%icgn_reduced_dim = total_ret
    write(line,'(a,i0,a,i0,a,i0,a,f8.4)') 'mode=icgn block diagonalization phase complete, total_retained=', &
         total_ret, ', reduced/original=', total_ret, '/', n, ', cut_fraction=', upt%icgn_cut_fraction
    call icgn_log_progress(upt, trim(line))
    call icgn_log_progress(upt, 'mode=icgn reduced Hamiltonian build started')

    call build_reduced_hamiltonian(upt, atom_of, label, row_of, pairs, npair, ierr)
    if (ierr /= 0) return
    write(line,'(a,i0,a,a1,a,i0)') 'mode=icgn reduced Hamiltonian build complete, active_pairs=', npair, &
         ', sparse_format=', upt%icgn_ham%sparse_fmt, ', reduced_nnz=', upt%icgn_ham%nnz
    call icgn_log_progress(upt, trim(line))
    call destroy_pairs(pairs)
    upt%icgn_ready = .true.
       call icgn_log_info(upt, 'icgn', n, total_ret, &
          upt%icgn_num_blocks, upt%icgn_cut_fraction, -1.0_dp, .false., .false.)
  end subroutine icgn_legacy_prepare

  subroutine diagonalize_block(upt, ib, atom_of, local, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, atom_of(:), local(:)
    integer, intent(out) :: ierr
    integer :: n, i, k, r, c
    complex(dp), allocatable :: h(:,:)
    real(dp), allocatable :: w(:)
    ierr = 0; n = upt%icgn_blocks(ib)%nrow
    if (n == 0) then
       upt%icgn_blocks(ib)%nret = 0
       return
    end if
    allocate(h(n,n), w(n)); h = (0.0_dp,0.0_dp)
    do r = 1, upt%ham%nrow
       if (local(r) == 0) cycle
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (local(c) == 0) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          h(local(r),local(c)) = upt%ham%M(k)
          if (r /= c) h(local(c),local(r)) = conjg(upt%ham%M(k))
       end do
    end do
   call coarse_eigh(h, w, ierr)
    if (ierr /= 0) return
    ! Store the FULL eigensystem — projection happens in build_reduced_hamiltonian
    ! using all nrow columns of S before the energy window truncation.
    allocate(upt%icgn_blocks(ib)%evals_full(n), upt%icgn_blocks(ib)%S_full(n,n))
    upt%icgn_blocks(ib)%evals_full = w
    upt%icgn_blocks(ib)%S_full     = h   ! columns are eigenvectors after zheevd
    ! Count how many states fall in the energy window (needed for reduced dim)
    upt%icgn_blocks(ib)%nret = count(w >= upt%icgn_emin .and. w <= upt%icgn_emax)
    deallocate(h, w)
  end subroutine diagonalize_block

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

  ! Block diagonalization always uses LAPACK (ZHEEVD via dense_eigh).
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


  subroutine build_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(ICGNPair), allocatable, intent(out) :: pairs(:)
    integer, intent(out) :: npair, ierr
    integer :: i,j,k,r,c,a,b,ia,ib,nnz,pos,slot,nred,na
    integer(kind=8) :: nnz64, timer_start, timer_end, timer_rate
    real(dp) :: timer_seconds
    character(len=256) :: line
    character(1) :: output_fmt
    integer, allocatable :: roff(:), rowcount(:), next(:), pair_map(:,:)
    integer, allocatable :: win_a(:)   ! indices of retained states within S_full
    ierr=0; nred=upt%icgn_reduced_dim

    ! --- Step 1: Build q/eval for each block from S_full by applying window ---
    ! (This is done here so S_full is available; S_full is freed at the end.)
    do i = 1, upt%icgn_num_blocks
       na = upt%icgn_blocks(i)%nrow
       if (na == 0 .or. upt%icgn_blocks(i)%nret == 0) cycle
       ! collect indices of retained states
       allocate(win_a(upt%icgn_blocks(i)%nret))
       k = 0
       do j = 1, na
          if (upt%icgn_blocks(i)%evals_full(j) >= upt%icgn_emin .and. &
              upt%icgn_blocks(i)%evals_full(j) <= upt%icgn_emax) then
             k = k + 1; win_a(k) = j
          end if
       end do
       allocate(upt%icgn_blocks(i)%eval(upt%icgn_blocks(i)%nret))
       allocate(upt%icgn_blocks(i)%q(na, upt%icgn_blocks(i)%nret))
       do j = 1, upt%icgn_blocks(i)%nret
          upt%icgn_blocks(i)%eval(j) = upt%icgn_blocks(i)%evals_full(win_a(j))
          upt%icgn_blocks(i)%q(:,j)  = upt%icgn_blocks(i)%S_full(:, win_a(j))
       end do
       deallocate(win_a)
    end do

    allocate(roff(upt%icgn_num_blocks+1)); roff(1)=1
    do i=1,upt%icgn_num_blocks; roff(i+1)=roff(i)+upt%icgn_blocks(i)%nret; end do
    ! At most one pair slot per distinct block pair (a<b), not per Hamiltonian entry.
    allocate(pairs(max(1, upt%icgn_num_blocks*(upt%icgn_num_blocks-1)/2))); npair=0
    allocate(pair_map(upt%icgn_num_blocks, upt%icgn_num_blocks)); pair_map = 0

    ! --- Step 2: Project inter-block couplings directly from global CSR ---
    ! For each active pair (A,B), accumulate T = H_AB * Q_B while scanning
    ! the existing sparse Hamiltonian.  This avoids allocating the dense
    ! physical-space H_AB (nA*nB entries) and avoids a dense H_AB*Q_B GEMM.
    call system_clock(timer_start, timer_rate)
    do r=1,upt%ham%nrow
       do k=upt%ham%Mi(r),upt%ham%Mi(r+1)-1
          c=upt%ham%Mj(k); a=label(atom_of(r)); b=label(atom_of(c))
          if(a==b) cycle
          if(upt%icgn_blocks(a)%nret==0 .or. upt%icgn_blocks(b)%nret==0) cycle
          if(.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          ia=min(a,b); ib=max(a,b)
          slot=pair_slot_cg(pairs,npair,ia,ib,upt,pair_map)
          if(slot==0) then; ierr=11; return; end if
          ! Accumulate T = H_AB * Q_B directly from the CSR entries, without
          ! building the dense H_AB. Pairs are stored once with ia < ib, so
          ! T always has the rows of block ia and the columns of Q_ib.
          !   a < b : entry (r,c) already is an (ia,ib) element:
          !           T(local(r),:) += M * Q_b(local(c),:)
          !   a > b : the entry is the (ib,ia) element, so use its Hermitian
          !           conjugate, which lives at (local(c), local(r)):
          !           T(local(c),:) += conjg(M) * Q_a(local(r),:)
          !           (Q of block a = ib, the row's block; NOT block b.)
          if(a < b) then
             pairs(slot)%t(local(r),:) = pairs(slot)%t(local(r),:) + &
                  upt%ham%M(k) * upt%icgn_blocks(b)%q(local(c),:)
          else
             pairs(slot)%t(local(c),:) = pairs(slot)%t(local(c),:) + &
                  conjg(upt%ham%M(k)) * upt%icgn_blocks(a)%q(local(r),:)
          end if
       end do
    end do

    deallocate(pair_map)
    call system_clock(timer_end)
    timer_seconds = real(timer_end-timer_start,dp) / real(max(1_8,timer_rate),dp)
    call icgn_log_timing(upt, 'ICGN CSR projection assembly', timer_seconds)
    write(line,'(a,i0,a,a1,a,i0)') 'mode=icgn interblock CSR projection assembly complete, active_pairs=', npair, &
         ', sparse_format=', upt%ham%sparse_fmt, ', ham_nnz=', upt%ham%nnz
    call icgn_log_progress(upt, trim(line))

    ! Finish V_AB = Q_A^H T with one dense GEMM per active pair.
    call icgn_log_progress(upt, 'mode=icgn interblock projection started')
    call system_clock(timer_start, timer_rate)
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       call project_accumulated_pair_cg(pairs(i), upt%icgn_blocks(a)%q)

       ! ICGN keeps all retained states, but sparsifies the projected couplings
       ! with the common filter (see icgn_coupling_survives):
       !   |V_ij|^2 > epsilon * |E_i-E_j|   (strict)
       ! Couplings failing it are set to exact zero; the states themselves
       ! remain in the reduced Hamiltonian. epsilon = 0 only removes exact zeros.
       do j = 1, size(pairs(i)%v,1)
          do k = 1, size(pairs(i)%v,2)
             if (.not. icgn_coupling_survives(pairs(i)%v(j,k), upt%icgn_blocks(a)%eval(j), &
                 upt%icgn_blocks(b)%eval(k), upt%icgn_epsilon)) then
                pairs(i)%v(j,k) = (0.0_dp,0.0_dp)
             end if
          end do
       end do

       if (i <= 5 .or. mod(i,100) == 0 .or. i == npair) then
          write(line,'(a,i0,a,i0,a,i0,a,i0)') 'mode=icgn interblock projection progress, pair=', i, &
               '/', npair, ', block_a=', a, ', block_b=', b
          call icgn_log_progress(upt, trim(line))
       end if
    end do
    call system_clock(timer_end)
    timer_seconds = real(timer_end-timer_start,dp) / real(max(1_8,timer_rate),dp)
    call icgn_log_timing(upt, 'ICGN dense left projection + epsilon filtering', timer_seconds)
    write(line,'(a,es12.4)') 'mode=icgn coupling filter epsilon=', upt%icgn_epsilon
    call icgn_log_progress(upt, trim(line))
    call icgn_log_progress(upt, 'mode=icgn interblock projection + epsilon filtering complete')
    call icgn_log_progress(upt, 'mode=icgn reduced Hamiltonian: projection phase exited')

    ! --- Step 3: Free S_full (no longer needed) ---
    call system_clock(timer_start, timer_rate)
    do i = 1, upt%icgn_num_blocks
       if (associated(upt%icgn_blocks(i)%S_full)) then
          deallocate(upt%icgn_blocks(i)%S_full)
          nullify(upt%icgn_blocks(i)%S_full)
       end if
       if (associated(upt%icgn_blocks(i)%evals_full)) then
          deallocate(upt%icgn_blocks(i)%evals_full)
          nullify(upt%icgn_blocks(i)%evals_full)
       end if
    end do
    call system_clock(timer_end)
    timer_seconds = real(timer_end-timer_start,dp) / real(max(1_8,timer_rate),dp)
    call icgn_log_timing(upt, 'ICGN S_full/evals_full deallocation', timer_seconds)
    call icgn_log_progress(upt, 'mode=icgn reduced Hamiltonian: S_full/evals_full deallocated')

    ! --- Step 4: Build CSR reduced Hamiltonian ---
    ! Keep the reduced Hamiltonian in the same sparse format as the input.
    ! In particular, CUDA JD consumes the CSR entries directly; it does not
    ! reconstruct the missing Hermitian triangle for sparse_fmt='U'.
    ! After the ICGN coupling filter, the full representation fits safely in
    ! the current 32-bit CSR index space for the present test case.
    output_fmt = upt%ham%sparse_fmt

    call system_clock(timer_start, timer_rate)
    allocate(rowcount(nred),next(nred)); rowcount=1
    do i=1,npair
       a=pairs(i)%a; b=pairs(i)%b
       select case(output_fmt)
       case('L', 'l')
          do j = 1, upt%icgn_blocks(b)%nret
             rowcount(roff(b)+j-1) = rowcount(roff(b)+j-1) + &
                  count(abs(pairs(i)%v(:,j)) > 0.0_dp)
          end do
       case('U', 'u')
          do j = 1, upt%icgn_blocks(a)%nret
             rowcount(roff(a)+j-1) = rowcount(roff(a)+j-1) + &
                  count(abs(pairs(i)%v(j,:)) > 0.0_dp)
          end do
       case default
          ! F: emit both Hermitian triangles.
          do j = 1, upt%icgn_blocks(a)%nret
             rowcount(roff(a)+j-1) = rowcount(roff(a)+j-1) + &
                  count(abs(pairs(i)%v(j,:)) > 0.0_dp)
          end do
          do j = 1, upt%icgn_blocks(b)%nret
             rowcount(roff(b)+j-1) = rowcount(roff(b)+j-1) + &
                  count(abs(pairs(i)%v(:,j)) > 0.0_dp)
          end do
       end select
    end do

    nnz64 = sum(int(rowcount,kind=8))
    if (nnz64 > int(huge(0),kind=8)) then
       ierr = 15
       write(*,'(a,a1,a,i0)') '(icgn) reduced CSR still exceeds default INTEGER capacity, format=', output_fmt, &
            ', nnz=', nnz64
       return
    end if
    nnz = int(nnz64)
    write(line,'(a,i0,a,i0,a,a1)') 'mode=icgn reduced Hamiltonian: CSR sizing, nred=', nred, &
         ', nnz=', nnz64, ', format=', output_fmt
    call icgn_log_progress(upt, trim(line))
    call create_matrix(upt%icgn_ham,nred,nred,nnz)
    call icgn_log_progress(upt, 'mode=icgn reduced Hamiltonian: CSR allocation complete')
    upt%icgn_ham%sparse_fmt=output_fmt; upt%icgn_ham%Mi(1)=1
    do i=1,nred; upt%icgn_ham%Mi(i+1)=upt%icgn_ham%Mi(i)+rowcount(i); end do
    next=upt%icgn_ham%Mi(1:nred)
    do a=1,upt%icgn_num_blocks
       do i=1,upt%icgn_blocks(a)%nret
          pos=next(roff(a)+i-1); upt%icgn_ham%Mj(pos)=roff(a)+i-1; upt%icgn_ham%M(pos)=upt%icgn_blocks(a)%eval(i)
          next(roff(a)+i-1)=pos+1
       end do
    end do
    do i=1,npair
       a=pairs(i)%a; b=pairs(i)%b
       call emit_pair(upt%icgn_ham,pairs(i),roff(a),roff(b),output_fmt,next)
    end do
    call icgn_log_progress(upt, 'mode=icgn reduced Hamiltonian: pair emission complete')
    upt%icgn_ham%nnz=nnz

    ! Dump the exact reduced CSR passed to the eigenvalue solver for
    ! independent validation outside TIBERCAD.  The binary layout is:
    !   int32 nrow, int32 nnz, complex(dp) M, int32 Mj, int32 Mi.
    ! A small text sidecar records the sparse format used by the CSR.
    open(unit=91, file='icgn_reduced_hamiltonian.dat', access='stream', &
         form='unformatted', status='replace')
    write(91) nred, nnz
    write(91) upt%icgn_ham%M
    write(91) upt%icgn_ham%Mj
    write(91) upt%icgn_ham%Mi
    close(91)
    open(unit=91, file='icgn_reduced_hamiltonian.meta', form='formatted', status='replace')
    write(91,'(a,i0)') 'n=', nred
    write(91,'(a,i0)') 'nnz=', nnz
    write(91,'(a,a1)') 'format=', output_fmt
    close(91)
    call icgn_log_progress(upt, 'mode=icgn reduced Hamiltonian dump: icgn_reduced_hamiltonian.dat')

    ! write(*,'(A,I0)') '[ICGN-TRACE] prepare: create icgn_U nred=', nred
    call create_matrix(upt%icgn_U,nred,nred,nred)
    ! write(*,'(A)') '[ICGN-TRACE] prepare: create icgn_U returned'
    upt%icgn_U%sparse_fmt='F'; upt%icgn_U%Mi(1)=1
    do i=1,nred
       upt%icgn_U%Mi(i)=i; upt%icgn_U%Mj(i)=i; upt%icgn_U%M(i)=(1.0_dp,0.0_dp)
    end do
    upt%icgn_U%Mi(nred+1)=nred+1; upt%icgn_U%nnz=nred
    deallocate(roff,rowcount,next)
    call system_clock(timer_end)
    timer_seconds = real(timer_end-timer_start,dp) / real(max(1_8,timer_rate),dp)
    call icgn_log_timing(upt, 'ICGN reduced CSR sizing/allocation/emission', timer_seconds)
  end subroutine build_reduced_hamiltonian

  integer function pair_slot(pairs,npair,a,b,upt,pair_map)
    type(ICGNPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer,intent(in)::a,b
    type(OUPT),intent(in)::upt
    integer, intent(inout) :: pair_map(:,:)

    ! pair_map makes the repeated pair lookup O(1) instead of scanning
    ! all previously created block pairs for every physical Hamiltonian entry.
    if (pair_map(a,b) /= 0) then
       pair_slot = pair_map(a,b)
       return
    end if

    npair=npair+1
    if(npair>size(pairs)) then; pair_slot=0; return; end if
    pairs(npair)%a=a; pairs(npair)%b=b
    allocate(pairs(npair)%h(upt%icgn_blocks(a)%nrow,upt%icgn_blocks(b)%nrow)); pairs(npair)%h=(0.0_dp,0.0_dp)
    pair_map(a,b) = npair
    pair_slot=npair
  end function pair_slot

  ! ICGN-specific pair slot: allocate only T = H_AB * Q_B.  Unlike pair_slot,
  ! this never materializes the dense physical-space H_AB matrix.
  integer function pair_slot_cg(pairs,npair,a,b,upt,pair_map)
    type(ICGNPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(in) :: a,b
    type(OUPT), intent(in) :: upt
    integer, intent(inout) :: pair_map(:,:)

    if (pair_map(a,b) /= 0) then
       pair_slot_cg = pair_map(a,b)
       return
    end if

    npair=npair+1
    if(npair>size(pairs)) then; pair_slot_cg=0; return; end if
    pairs(npair)%a=a; pairs(npair)%b=b
    allocate(pairs(npair)%t(upt%icgn_blocks(a)%nrow, upt%icgn_blocks(b)%nret))
    pairs(npair)%t=(0.0_dp,0.0_dp)
    pair_map(a,b)=npair
    pair_slot_cg=npair
  end function pair_slot_cg

  subroutine project_accumulated_pair_cg(p, qa)
    type(ICGNPair), intent(inout) :: p
    complex(dp), intent(in) :: qa(:,:)

    ! The CSR scan has already computed T = H_AB * Q_B.  Only the dense
    ! left projection remains: V_AB = Q_A^H * T.
    allocate(p%v(size(qa,2), size(p%t,2)))
    p%v = matmul(conjg(transpose(qa)), p%t)
    deallocate(p%t)
  end subroutine project_accumulated_pair_cg

  subroutine project_pair(p, qa, qb, keep_h)
    type(ICGNPair), intent(inout) :: p
    complex(dp), intent(in) :: qa(:,:), qb(:,:)
    logical, intent(in), optional :: keep_h
    complex(dp), allocatable :: tmp(:,:)

    ! First H_AB * Q_B, then Q_A^H * (H_AB * Q_B).
    ! The expensive arithmetic is now done at block level instead of once per
    ! physical Hamiltonian entry.
    allocate(p%v(size(qa,2), size(qb,2)))
    call icgn_project_sparse(p%h, qa, qb, p%v)
    if (.not. present(keep_h) .or. .not. keep_h) deallocate(p%h)
  end subroutine project_pair

  ! Write the coupling block p%v of the block pair (a,b) into the reduced CSR h.
  !   oa, ob     : offsets of the retained states of blocks a and b in the reduced basis
  !   fmt        : 'F' writes the block and its Hermitian conjugate, 'U' only
  !                the block (upper), 'L' only the conjugate (lower)
  !   next(:)    : running insertion position of every reduced row (updated)
  !   keep_zeros : if .true. exact zeros are written too. Default: zeros are
  !                skipped. The caller must size the CSR with the same rule
  !                (ICGN and ICGN both count only the non-zero entries).
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
    ! Both ICGN and ICGN count only non-zero entries when sizing the CSR, so
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

  subroutine destroy_pairs(pairs)
    type(ICGNPair),allocatable,intent(inout)::pairs(:)
    integer::i,nv,nh,nt
    if(.not.allocated(pairs)) then
       ! write(*,'(A)') '[ICGN-TRACE] destroy_pairs: not allocated'
       return
    end if
    nv=0; nh=0; nt=0
    do i=1,size(pairs)
       if(associated(pairs(i)%v)) nv=nv+1
       if(associated(pairs(i)%h)) nh=nh+1
       if(associated(pairs(i)%t)) nt=nt+1
    end do
    ! write(*,'(A,I0,A,I0,A,I0,A,I0)') '[ICGN-TRACE] destroy_pairs count=', size(pairs), &
         ! ' v=', nv, ' h=', nh, ' t=', nt
    do i=1,size(pairs)
       if(associated(pairs(i)%v)) deallocate(pairs(i)%v)
       if(associated(pairs(i)%h)) deallocate(pairs(i)%h)
       if(associated(pairs(i)%t)) deallocate(pairs(i)%t)
       if(associated(pairs(i)%g)) deallocate(pairs(i)%g)
    end do
    deallocate(pairs)
    ! write(*,'(A)') '[ICGN-TRACE] destroy_pairs DONE'
  end subroutine destroy_pairs

  subroutine icgn_legacy_lift(upt, reduced, physical)
    type(OUPT),intent(in)::upt
    complex(dp),intent(in)::reduced(:,:)
    complex(dp),intent(out)::physical(:,:)
    integer::i,j,k,off
    physical=(0.0_dp,0.0_dp); off=1
    do i=1,size(upt%icgn_blocks)
       if(upt%icgn_blocks(i)%nret>0) then
          do j=1,size(upt%icgn_blocks(i)%rows)
             physical(upt%icgn_blocks(i)%rows(j),:)=matmul(upt%icgn_blocks(i)%q(j,:),reduced(off:off+upt%icgn_blocks(i)%nret-1,:))
          end do
       end if
       off=off+upt%icgn_blocks(i)%nret
    end do
  end subroutine icgn_legacy_lift

  ! ===========================================================================
  ! ==========================================================================
  ! ICGN: improved ICGN with optional Neumann self-energy correction (mode 2)
  !
  ! Strategy: icgn_legacy_prepare partitions atoms → diagonalizes blocks (full S) →
  ! selects retained states with a core+buffer criterion → projects inter-block
  ! couplings → builds icgn_ham CSR. When neumann_order >= 0, the Neumann
  ! self-energy  Sigma = H_PQ * (E0-D)^-1 * [W*(E0-D)^-1]^order * H_QP
  ! is extracted from 'pairs' and added into the dense form of icgn_ham before
  ! converting back to CSR. When neumann_order < 0 the correction is skipped.
  ! ==========================================================================

  subroutine icgn_configure(upt, enabled, nblocks, core_emin, core_emax, &
     top_buffer, bottom_buffer, epsilon, add_core_acquaintances, selfenergy_order, E0, metis_imbalance, &
       check_convergence, pi_maxiter, pi_tol)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled, check_convergence
    integer, intent(in) :: nblocks, add_core_acquaintances, selfenergy_order, pi_maxiter
   real(dp), intent(in) :: core_emin, core_emax, top_buffer, bottom_buffer, epsilon, E0, metis_imbalance, pi_tol
    upt%icgn_enabled = enabled
    upt%icgn_num_blocks = nblocks
    upt%icgn_core_emin = core_emin
    upt%icgn_core_emax = core_emax
   upt%icgn_top_buffer = top_buffer
   upt%icgn_bottom_buffer = bottom_buffer
    upt%icgn_epsilon = epsilon
    upt%icgn_add_core_acquaintances = max(0, add_core_acquaintances)
    upt%icgn_selfenergy_order = selfenergy_order
    upt%icgn_E0 = E0
    upt%icgn_metis_imbalance = metis_imbalance
    upt%icgn_check_convergence = check_convergence
    upt%icgn_pi_maxiter = pi_maxiter
    upt%icgn_pi_tol     = pi_tol
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

  subroutine icgn_clear(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    ! write(*,'(A)') '[ICGN-TRACE] clear ENTER'
    ! write(*,'(A,L1)') '[ICGN-TRACE] clear ham associated=', associated(upt%icgn_ham%M)
    if (associated(upt%icgn_ham%M)) then
       ! write(*,'(A)') '[ICGN-TRACE] clear destroy ham BEGIN'
       call destroy_matrix(upt%icgn_ham)
       ! write(*,'(A)') '[ICGN-TRACE] clear destroy ham END'
    end if
    ! write(*,'(A,L1)') '[ICGN-TRACE] clear U associated=', associated(upt%icgn_U%M)
    if (associated(upt%icgn_U%M)) then
       ! write(*,'(A)') '[ICGN-TRACE] clear destroy U BEGIN'
       call destroy_matrix(upt%icgn_U)
       ! write(*,'(A)') '[ICGN-TRACE] clear destroy U END'
    end if
    if (associated(upt%icgn_blocks)) then
       ! write(*,'(A,I0)') '[ICGN-TRACE] clear blocks count=', size(upt%icgn_blocks)
       do i = 1, size(upt%icgn_blocks)
          if (associated(upt%icgn_blocks(i)%rows))         deallocate(upt%icgn_blocks(i)%rows)
          if (associated(upt%icgn_blocks(i)%eval))         deallocate(upt%icgn_blocks(i)%eval)
          if (associated(upt%icgn_blocks(i)%q))            deallocate(upt%icgn_blocks(i)%q)
          if (associated(upt%icgn_blocks(i)%evals_full))   deallocate(upt%icgn_blocks(i)%evals_full)
          if (associated(upt%icgn_blocks(i)%S_full))       deallocate(upt%icgn_blocks(i)%S_full)
          if (associated(upt%icgn_blocks(i)%retained_idx)) deallocate(upt%icgn_blocks(i)%retained_idx)
       end do
       deallocate(upt%icgn_blocks)
       ! write(*,'(A)') '[ICGN-TRACE] clear blocks deallocated'
    end if
    upt%icgn_ready        = .false.
    upt%icgn_original_dim = 0
    upt%icgn_reduced_dim  = 0
    upt%icgn_cut_fraction = 0.0_dp
    ! write(*,'(A)') '[ICGN-TRACE] clear EXIT'
  end subroutine icgn_clear


  ! ---------------------------------------------------------------------------
  ! icgn_prepare: same P-space selection as icgn_prepare, but after building the
  ! reduced Hamiltonian we add a Neumann-series self-energy correction for the
  ! discarded Q states. 'pairs' holds the full (nrow_a x nrow_b) coupling
  ! matrices in the block eigenbasis; we reuse them to extract H_PQ and H_QQ.
  ! ---------------------------------------------------------------------------
  subroutine icgn_prepare(upt, ierr)
    character(len=*), parameter :: core_acq_warning = &
         'Including core-state acquaintances higher than level 1 is not yet implemented --> ' // &
         'Fall back to level-1 acquaintances only.'
    type(OUPT), intent(inout) :: upt
    integer, intent(out) :: ierr

    integer :: n, na, nb_atoms, i, j, k, p, status, nedge, maxedge
    integer :: r, c, br, bc, npair, total_ret, pos
    integer, allocatable :: atom_of(:), local_of(:), label(:), row_of(:)
    integer, allocatable :: counts(:), cursor(:), offsets(:), bsize(:)
    integer(c_int), allocatable :: xadj(:), adjncy(:), vwgt(:), adjwgt(:), part(:)
    real(dp), allocatable :: edge_weight(:)
    real(dp) :: max_weight, all_weight, cut_weight
    type(ICGNPair), allocatable :: pairs(:)

    logical, allocatable :: is_core(:,:), keep_mask(:,:)

    integer :: ia, ib, slot, nred, nnz
    integer, allocatable :: roff(:), rowcount(:), next(:), win_a(:), win_b(:), pair_map(:,:)
    ! g_full points either to pairs(i)%g (kept projection) or to a locally computed array
    complex(dp), pointer :: g_full(:,:)

    ! For Neumann self-energy
    integer :: nb, ord, ip, iq, gp_row, gp_col, q_idx, qb_idx
    integer, allocatable :: block_of_state(:), local_of_state(:), ret_offset(:)
    real(dp), allocatable :: evals_q(:)
    complex(dp), allocatable :: H_dense(:,:), Sigma(:,:)
    complex(dp), allocatable :: amp_vec(:), next_vec(:)
    real(dp) :: E0_used, resolvent_val
    complex(dp) :: amp_val, contrib
    ! pq_list: list of (global_p_in_red, global_q_flat, coupling_value)
    integer, allocatable :: pq_p(:), pq_q(:)
    complex(dp), allocatable :: pq_v(:)
    ! qq_list: Q-Q coupling edges in block eigenbasis
    integer, allocatable :: qq_i(:), qq_j(:)
    complex(dp), allocatable :: qq_v(:)
    integer :: npq, nqq, cnt_pq, cnt_qq
    integer :: nstates_total
    ! block→reduced-row offset
    integer, allocatable :: roff_blk(:)
    ! local state index → retained index (0 if discarded)
    integer, allocatable :: local_ret_idx(:,:)
    ! OpenMP work arrays (block diagonalization / pair projection)
    integer :: nth, saved_levels, mypq, myqq
    integer, allocatable :: blk_of(:), ierr_blk(:)
    integer, allocatable :: pq_off(:), pq_cnt(:), qq_off(:), qq_cnt(:)
    logical :: need_qq
    real(dp) :: t_all, t_stage, amax, t_mark
    complex(dp), allocatable :: tmp_g(:,:)
    ! Sparse Sigma work arrays
    real(dp), allocatable :: res_q(:)
    integer, allocatable :: rp_ptr(:), rp_ord(:), cq_ptr(:), cq_ord(:), qs_ptr(:), qs_ord(:)
    integer, allocatable :: supp(:), supp2(:)
    logical, allocatable :: mark(:), mark2(:)
    integer :: nsupp, nsupp2, tt, ee
    logical :: dummy_qq
    ! Block-dense Neumann correction
    logical :: use_blocked
    type(ICGNQBlock), allocatable :: bq(:)
    type(ICGNQPair), allocatable :: qp(:)
    integer, allocatable :: nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, allocatable :: nbr_isA(:)
    integer(8) :: npq_nz, nqq_nz
    ! Shared edge graphs for the convergence check and Sigma (built once)
    integer, allocatable :: p_ptr(:), p_col(:), c_ptr(:), c_row(:), s_ptr(:), s_col(:)
    complex(dp), allocatable :: p_val(:), c_val(:), s_val(:), s_wres(:)

    ierr = 0
    ! write(*,"(A)") "[ICGN-TRACE] prepare ENTER"
    ! write(*,"(A)") "[ICGN-TRACE] prepare before clear"
    call icgn_clear(upt)
    ! write(*,"(A)") "[ICGN-TRACE] prepare after clear"
    if (.not. upt%icgn_enabled) return
   call icgn_log_progress(upt, 'mode=icgn preparation started')
    if (num_procs /= 1) then
       ierr = 1; write(*,*) '(icgn) MPI not supported'; return
    end if

    na = upt%basis%n_basis
    n  = upt%ham%nrow
    if (na < 1 .or. upt%icgn_num_blocks < 1 .or. upt%icgn_num_blocks > na) then
       ierr = 2; write(*,*) '(icgn) invalid number of blocks'; return
    end if
   if (upt%icgn_core_emin >= upt%icgn_core_emax .or. upt%icgn_top_buffer < 0.0_dp .or. &
      upt%icgn_bottom_buffer < 0.0_dp) then
       ierr = 3; write(*,*) '(icgn) invalid core window or buffer'; return
    end if
    ! write(*,'(A,L1)') '[ICGN-TRACE] prepare: ham associated=', associated(upt%ham%M)
    if (.not. associated(upt%ham%M)) then
       ierr = 4; write(*,*) '(icgn) Hamiltonian not initialized'; return
    end if
    ! write(*,'(A,I0,A,I0,A,I0)') '[ICGN-TRACE] prepare ham nrow=', upt%ham%nrow, ' nnz=', upt%ham%nnz, ' n_basis=', upt%basis%n_basis
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] prepare ham Mi=', associated(upt%ham%Mi), ' Mj=', associated(upt%ham%Mj)

    ! ---- Build atom→orbital mapping (identical to icgn_prepare) ------------
    ! write(*,'(A,I0,A,I0)') '[ICGN-TRACE] prepare before mapping allocate n=', n, ' na=', na
    allocate(atom_of(n), local_of(n), offsets(na+1), bsize(na))
    ! write(*,'(A)') '[ICGN-TRACE] prepare after mapping allocate'
    pos = 1
    ! write(*,'(A)') '[ICGN-TRACE] prepare before mapping loop'
    do i = 1, na
       offsets(i) = pos
       bsize(i)   = upt%n_spin * upt%basis%n_st(i)
       do j = 1, bsize(i)
          atom_of(pos) = i; local_of(pos) = j; pos = pos + 1
       end do
    end do
    ! write(*,'(A,I0)') '[ICGN-TRACE] prepare after mapping loop pos=', pos
    call icgn_log_progress(upt, 'mode=icgn graph partition complete')
    offsets(na+1) = pos
    if (pos-1 /= n) then
       ierr = 5; write(*,*) '(icgn) atom/orbital mapping inconsistent'; return
    end if

    ! ---- METIS partition (same as icgn_prepare) ----------------------------
    maxedge = max(1, upt%ham%nnz)
    allocate(edge_weight(maxedge), counts(na), label(na))
    nedge = 0; max_weight = 0.0_dp; counts = 0
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          nedge = nedge + 1
          if (nedge > maxedge) then; ierr = 6; return; end if
          edge_weight(nedge) = abs(upt%ham%M(k))**2
          counts(atom_of(r)) = counts(atom_of(r)) + 1
          counts(atom_of(c)) = counts(atom_of(c)) + 1
          max_weight = max(max_weight, edge_weight(nedge))
       end do
    end do
    if (nedge == 0 .or. max_weight == 0.0_dp) then
       ierr = 7; write(*,*) '(icgn) atom graph has no couplings'; return
    end if
    allocate(xadj(na+1), cursor(na), adjncy(2*nedge), adjwgt(2*nedge), vwgt(na), part(na))
    xadj(1) = 0_c_int
    do i = 1, na
       xadj(i+1) = xadj(i) + int(counts(i), c_int)
       cursor(i)  = int(xadj(i)) + 1
       vwgt(i)    = int(bsize(i), c_int)
    end do
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          br = atom_of(r); bc = atom_of(c)
          p  = max(1, nint(abs(upt%ham%M(k))**2 / max_weight * 1000000.0_dp))
          adjncy(cursor(br)) = int(bc-1,c_int); adjwgt(cursor(br)) = int(p,c_int); cursor(br)=cursor(br)+1
          adjncy(cursor(bc)) = int(br-1,c_int); adjwgt(cursor(bc)) = int(p,c_int); cursor(bc)=cursor(bc)+1
       end do
    end do
    status = icgn_metis_partition(int(na,c_int), xadj, adjncy, vwgt, adjwgt, &
         int(upt%icgn_num_blocks,c_int), int(nint(1000.0_dp*upt%icgn_metis_imbalance),c_int), 42_c_int, part)
    if (status /= 0) then
       call icgn_log_progress(upt, 'METIS unavailable, using built-in graph-BFS fallback partitioning')
       call icgn_graph_partition(na, upt%icgn_num_blocks, vwgt, xadj, adjncy, adjwgt, part)
    else
       call icgn_log_progress(upt, 'METIS partitioning done')
    end if
    do i = 1, na; label(i) = int(part(i)) + 1; end do

    all_weight = 0.0_dp; cut_weight = 0.0_dp
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_weight = all_weight + abs(upt%ham%M(k))**2
          if (label(atom_of(r)) /= label(atom_of(c))) &
               cut_weight = cut_weight + abs(upt%ham%M(k))**2
       end do
    end do
    upt%icgn_cut_fraction = cut_weight / all_weight

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
          ierr = ierr_blk(i); return
       end if
    end do
    deallocate(ierr_blk, blk_of)

    ! ---- Improved keep-mask: core + buffer + acquaintance -------------------
    allocate(is_core(nb, maxval(counts)), keep_mask(nb, maxval(counts)))
    is_core = .false.; keep_mask = .false.

    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          if (.not. associated(upt%icgn_blocks(i)%evals_full)) cycle
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
    end do

    row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
    end do

    ! Level-1 acquaintance: collect physical block couplings, then project
    ! each pair once at matrix level.  The zero-acquaintance/no-Neumann case
    ! defers pair construction until q is available and uses the ICGN-style
    ! sparse projection path instead of materializing dense physical H_ab.
    if (upt%icgn_add_core_acquaintances > 0 .or. upt%icgn_selfenergy_order >= 0) then
       allocate(pairs(max(1, upt%ham%nnz)))
       t_stage = icgn_wall()
       call icgn_collect_pairs(upt, atom_of, label, row_of, pairs, npair, ierr)
       if (ierr /= 0) return
       call icgn_report_time('collect pair couplings H_ab', t_stage)
       if (icgn_timing_on()) call icgn_report_boundary(pairs, npair)
    else
       allocate(pairs(max(1, upt%icgn_num_blocks*(upt%icgn_num_blocks-1)/2))); npair = 0
    end if
    if (upt%icgn_add_core_acquaintances > 0) then
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
       call icgn_report_time('project acquaintance pairs', t_stage)
    end if

    ! Acquaintance selection is optional. Level 1 is the currently implemented
    ! core-state acquaintance level; higher requested levels fall back to it.
    t_stage = icgn_wall()
    if (upt%icgn_add_core_acquaintances > 1) then
       call icgn_log_important(upt, core_acq_warning)
       write(*,'(a)') core_acq_warning
    end if
    if (upt%icgn_add_core_acquaintances > 0) then
       do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       do j = 1, upt%icgn_blocks(ia)%nrow
          if (is_core(ia, j)) then
             do k = 1, upt%icgn_blocks(ib)%nrow
                if (.not. keep_mask(ib, k)) then
                   if (icgn_coupling_survives(pairs(i)%v(j,k), upt%icgn_blocks(ia)%evals_full(j), &
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
                   if (icgn_coupling_survives(pairs(i)%v(j,k), upt%icgn_blocks(ib)%evals_full(k), &
                       upt%icgn_blocks(ia)%evals_full(j), upt%icgn_epsilon)) &
                       keep_mask(ia, j) = .true.
                end if
             end do
          end if
       end do
    end do
       end if
    call icgn_report_time('acquaintance selection (keep_mask, serial)', t_stage)
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
       ierr = 9; write(*,*) '(icgn) no states retained'; return
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
         upt%icgn_add_core_acquaintances > 0, ierr)
    if (ierr /= 0) return
    call icgn_report_time('build reduced H (pair projection)', t_stage)

    ! ====================================================================
    ! Self-energy correction via Neumann series (skipped when neumann_order < 0)
    ! ====================================================================
    if (upt%icgn_selfenergy_order >= 0) then
    ! pairs(s)%v is now (nret_a x nret_b): H_PP sliced couplings.
    ! We need H_PQ: from the *unsliced* g_full before retained_idx slicing.
    ! Strategy: re-project directly from S_full to extract P-Q and Q-Q.

    t_mark = icgn_wall()
    nred = total_ret
    allocate(H_dense(nred, nred))
    H_dense = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Fill H_dense from icgn_ham CSR (may be triangular — symmetrize immediately)
    do r = 1, upt%icgn_ham%nrow
       do k = upt%icgn_ham%Mi(r), upt%icgn_ham%Mi(r+1)-1
          c = upt%icgn_ham%Mj(k)
          H_dense(r, c) = upt%icgn_ham%M(k)
          H_dense(c, r) = conjg(upt%icgn_ham%M(k))
       end do
    end do

    ! Build roff: global reduced-row offset per block
    allocate(roff_blk(nb+1)); roff_blk(1) = 1
    do i = 1, nb; roff_blk(i+1) = roff_blk(i) + upt%icgn_blocks(i)%nret; end do

    ! Build local_ret_idx(block, local_state) → retained index in block (0=discarded)
    allocate(local_ret_idx(nb, maxval([(upt%icgn_blocks(i)%nrow, i=1,nb)])))
    local_ret_idx = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nret
          local_ret_idx(i, upt%icgn_blocks(i)%retained_idx(j)) = j
       end do
    end do

    ! Count total block states for flat Q-space indexing
    nstates_total = sum([(upt%icgn_blocks(i)%nrow, i=1,nb)])
    allocate(block_of_state(nstates_total), local_of_state(nstates_total), &
         evals_q(nstates_total), ret_offset(nb))
    pos = 0
    do i = 1, nb
       ret_offset(i) = pos
       do j = 1, upt%icgn_blocks(i)%nrow
          pos = pos + 1
          block_of_state(pos) = i
          local_of_state(pos) = j
          evals_q(pos) = upt%icgn_blocks(i)%evals_full(j)
       end do
    end do

    ! Exact P-Q count per pair: nret_a * (nrow_b - nret_b) + nret_b * (nrow_a - nret_a)
    cnt_pq = 0
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       cnt_pq = cnt_pq + upt%icgn_blocks(ia)%nret * &
            (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
       cnt_pq = cnt_pq + upt%icgn_blocks(ib)%nret * &
            (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret)
    end do
    ! pairs%v now holds the *sliced* P-P block. We need the g_full from
    ! the projection step. Since we already freed g_full inside
    ! build_icgn_reduced_hamiltonian, we need to re-project.
    ! Re-project inter-block couplings to get the full g matrices.
    call icgn_report_time('Neumann setup (H_dense fill, index tables)', t_mark)
    use_blocked = icgn_neumann_blocked()
    if (use_blocked) then
       ! ---- Block-dense Neumann correction (dense BLAS3 on the block structure) ----
       E0_used = upt%icgn_E0
       if (E0_used == 0.0_dp) E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp
       need_qq = (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence)
       t_stage = icgn_wall()
       call icgn_build_neumann_blocks(upt, pairs, npair, local_ret_idx, E0_used, need_qq, &
            bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, npq_nz, nqq_nz)
       npq = int(npq_nz); nqq = int(nqq_nz)
       call icgn_report_time('Neumann: dense P-Q/Q-Q blocks', t_stage)
       t_mark = icgn_wall()
       ! Free S_full now that the couplings are projected
       do i = 1, nb
          if (associated(upt%icgn_blocks(i)%S_full))     deallocate(upt%icgn_blocks(i)%S_full)
          if (associated(upt%icgn_blocks(i)%evals_full)) deallocate(upt%icgn_blocks(i)%evals_full)
          nullify(upt%icgn_blocks(i)%S_full, upt%icgn_blocks(i)%evals_full)
       end do
       call destroy_pairs(pairs)
       ! dummies so the common clean-up below stays valid
       allocate(pq_p(1), pq_q(1), pq_v(1)); pq_p = 1; pq_q = 1; pq_v = (0.0_dp, 0.0_dp)

       upt%icgn_sigma_T2     = -1.0_dp
       upt%icgn_pi_converged = .true.
       if (upt%icgn_check_convergence .and. nqq > 0) then
          t_stage = icgn_wall()
          call icgn_power_iteration_blk(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, &
               upt%icgn_pi_maxiter, upt%icgn_pi_tol, upt%icgn_sigma_T2, upt%icgn_pi_converged)
          call icgn_report_convergence(upt)
          call icgn_report_time('power iteration (block matvec)', t_stage)
       end if
       call icgn_report_time('free S_full/pairs + convergence check', t_mark)

       allocate(Sigma(nred, nred))
       Sigma = cmplx(0.0_dp, 0.0_dp, kind=dp)
       t_stage = icgn_wall()
       if (npq > 0) call icgn_sigma_blockdense(nb, upt%icgn_selfenergy_order, &
            maxval(upt%icgn_blocks(1:nb)%nrow), roff_blk, bq, qp, &
            nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, Sigma)
       call icgn_report_time('Sigma (block-dense zgemm)', t_stage)
       call icgn_free_neumann_blocks(bq, qp, npair)
       deallocate(nbr_ptr, nbr_blk, nbr_pair, nbr_rev, nbr_isA)
    else
    allocate(pq_p(max(1,cnt_pq)), pq_q(max(1,cnt_pq)))
    allocate(pq_v(max(1,cnt_pq)))
    npq = 0

    ! Each pair owns a disjoint, exactly-bounded slice of the flat edge lists
    ! ([off(i)+1, off(i+1)]), so pairs can be projected concurrently.  The slices
    ! are compacted afterwards in pair order, which reproduces the serial edge
    ! ordering (hence the summation order of Sigma) exactly.
    need_qq = (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence)
    allocate(pq_off(npair+1), pq_cnt(npair), qq_off(npair+1), qq_cnt(npair))
    pq_off(1) = 0; qq_off(1) = 0; pq_cnt = 0; qq_cnt = 0
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       pq_off(i+1) = pq_off(i) &
            + upt%icgn_blocks(ia)%nret * (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret) &
            + upt%icgn_blocks(ib)%nret * (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret)
       qq_off(i+1) = qq_off(i)
       if (need_qq) qq_off(i+1) = qq_off(i) + 2 * &
            (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret) * &
            (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
    end do
    if (need_qq) then
       allocate(qq_i(max(1,qq_off(npair+1))), qq_j(max(1,qq_off(npair+1))))
       allocate(qq_v(max(1,qq_off(npair+1))))
    end if

    ! Re-project each coupling pair (a,b): g_full = S_a^H H_ab S_b is computed once
    ! and used for both the P-Q and (if requested) the Q-Q edges.
    nth = min(icgn_thread_budget(), max(1, npair))
    t_stage = icgn_wall()
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel do default(shared) private(i, ia, ib, j, k, g_full, tmp_g, mypq, myqq) &
    !$omp& num_threads(nth) schedule(dynamic,1) if(nth > 1)
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       if (associated(pairs(i)%g)) then
          ! Acquaintance path: the full projection S_a^H H_ab S_b was kept.
          g_full => pairs(i)%g
       else
          allocate(g_full(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
          allocate(tmp_g(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
          ! The physical block coupling was already collected in pairs(i)%h.
          call icgn_zgemm('N', pairs(i)%h, upt%icgn_blocks(ib)%S_full, tmp_g)
          call icgn_zgemm('C', upt%icgn_blocks(ia)%S_full, tmp_g, g_full)
          deallocate(tmp_g)
       end if
       mypq = pq_off(i)
       ! Extract P(ia)-Q(ib) entries
       do j = 1, upt%icgn_blocks(ia)%nret
          do k = 1, upt%icgn_blocks(ib)%nrow
             if (local_ret_idx(ib, k) /= 0) cycle ! k is in P
             if (abs(g_full(upt%icgn_blocks(ia)%retained_idx(j), k)) < 1.0e-14_dp) cycle
             if (.not. icgn_coupling_survives(g_full(upt%icgn_blocks(ia)%retained_idx(j), k), &
                 upt%icgn_blocks(ia)%evals_full(upt%icgn_blocks(ia)%retained_idx(j)), &
                 upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) cycle
             mypq = mypq + 1
             pq_p(mypq) = roff_blk(ia) + j - 1
             pq_q(mypq) = ret_offset(ib) + k
             pq_v(mypq) = g_full(upt%icgn_blocks(ia)%retained_idx(j), k)
          end do
       end do
       ! Extract P(ib)-Q(ia) entries (g_ba = g_ab^†)
       do j = 1, upt%icgn_blocks(ib)%nret
          do k = 1, upt%icgn_blocks(ia)%nrow
             if (local_ret_idx(ia, k) /= 0) cycle ! k is in P
             if (abs(g_full(k, upt%icgn_blocks(ib)%retained_idx(j))) < 1.0e-14_dp) cycle
             if (.not. icgn_coupling_survives(g_full(k, upt%icgn_blocks(ib)%retained_idx(j)), &
                 upt%icgn_blocks(ia)%evals_full(k), &
                 upt%icgn_blocks(ib)%evals_full(upt%icgn_blocks(ib)%retained_idx(j)), &
                 upt%icgn_epsilon)) cycle
             mypq = mypq + 1
             pq_p(mypq) = roff_blk(ib) + j - 1
             pq_q(mypq) = ret_offset(ia) + k
             pq_v(mypq) = conjg(g_full(k, upt%icgn_blocks(ib)%retained_idx(j)))
          end do
       end do
       pq_cnt(i) = mypq - pq_off(i)
       if (need_qq) then
          myqq = qq_off(i)
          ! Store BOTH directions (j->k and k->j with conj(v))
          do j = 1, upt%icgn_blocks(ia)%nrow
             if (local_ret_idx(ia, j) /= 0) cycle ! j in P, skip
             do k = 1, upt%icgn_blocks(ib)%nrow
                if (local_ret_idx(ib, k) /= 0) cycle ! k in P, skip
                if (abs(g_full(j, k)) < 1.0e-14_dp) cycle
                if (.not. icgn_coupling_survives(g_full(j, k), upt%icgn_blocks(ia)%evals_full(j), &
                    upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) cycle
                myqq = myqq + 1
                qq_i(myqq) = ret_offset(ia) + j
                qq_j(myqq) = ret_offset(ib) + k
                qq_v(myqq) = g_full(j, k)
                ! Reverse direction: k -> j, conj(v)
                myqq = myqq + 1
                qq_i(myqq) = ret_offset(ib) + k
                qq_j(myqq) = ret_offset(ia) + j
                qq_v(myqq) = conjg(g_full(j, k))
             end do
          end do
          qq_cnt(i) = myqq - qq_off(i)
       end if
       if (.not. associated(pairs(i)%g)) deallocate(g_full)
    end do
    !$omp end parallel do
    call icgn_restore_nesting(saved_levels)
    call icgn_report_time('Neumann P-Q/Q-Q projection', t_stage)

    ! Compact the per-pair slices into contiguous edge lists (pair order).
    ! Destination index never exceeds the source index, so the ascending copy is safe.
    npq = 0
    do i = 1, npair
       do k = 1, pq_cnt(i)
          npq = npq + 1
          pq_p(npq) = pq_p(pq_off(i)+k)
          pq_q(npq) = pq_q(pq_off(i)+k)
          pq_v(npq) = pq_v(pq_off(i)+k)
       end do
    end do
    nqq = 0
    if (need_qq) then
       do i = 1, npair
          do k = 1, qq_cnt(i)
             nqq = nqq + 1
             qq_i(nqq) = qq_i(qq_off(i)+k)
             qq_j(nqq) = qq_j(qq_off(i)+k)
             qq_v(nqq) = qq_v(qq_off(i)+k)
          end do
       end do
    end if
    deallocate(pq_off, pq_cnt, qq_off, qq_cnt)
    t_mark = icgn_wall()

    ! Free S_full now that projections are done
    do i = 1, nb
       if (associated(upt%icgn_blocks(i)%S_full))     deallocate(upt%icgn_blocks(i)%S_full)
       if (associated(upt%icgn_blocks(i)%evals_full)) deallocate(upt%icgn_blocks(i)%evals_full)
       nullify(upt%icgn_blocks(i)%S_full, upt%icgn_blocks(i)%evals_full)
    end do
    call destroy_pairs(pairs)

    ! ---- Convergence check: estimate ||T||_2 = ||(E0-D)^-1 W||_2 -----------
    ! Uses qq_i/qq_j/qq_v already built. Only meaningful when order >= 1.
    upt%icgn_sigma_T2     = -1.0_dp
    upt%icgn_pi_converged = .true.
    E0_used = upt%icgn_E0
    if (E0_used == 0.0_dp) E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp

    ! Edge graphs (contiguous CSR) shared by the convergence check and Sigma:
    ! built once here and reused, together with the resolvent R(q).
    t_stage = icgn_wall()
    allocate(res_q(nstates_total))
    do i = 1, nstates_total
       res_q(i) = 1.0_dp / (E0_used - evals_q(i))
    end do
    dummy_qq = .not. allocated(qq_i)
    if (dummy_qq) allocate(qq_i(1), qq_j(1), qq_v(1))   ! Q-Q lists exist only if needed
    call icgn_build_edge_csr(npq, pq_p, pq_q, pq_v, nqq, qq_i, qq_j, qq_v, nred, nstates_total, res_q, &
         nqq > 0 .and. (upt%icgn_check_convergence .or. upt%icgn_selfenergy_order >= 1), &
         p_ptr, p_col, p_val, c_ptr, c_row, c_val, s_ptr, s_col, s_val, s_wres)
    if (dummy_qq) deallocate(qq_i, qq_j, qq_v)
    call icgn_report_time('build edge CSR (shared by conv. check + Sigma)', t_stage)
    t_stage = icgn_wall()
    if (upt%icgn_check_convergence .and. nqq > 0) then
       call icgn_power_iteration_csr(evals_q, nstates_total, s_ptr, s_col, s_val, E0_used, &
            upt%icgn_pi_maxiter, upt%icgn_pi_tol, &
            upt%icgn_sigma_T2, upt%icgn_pi_converged)
       if (.not. upt%icgn_pi_converged) then
          write(*,'(a,i0,a)') &
               '  (icgn) WARNING: power iteration did not converge in ', &
               upt%icgn_pi_maxiter, ' sweeps — ||T||_2 estimate may be unreliable.'
       end if
       if (upt%icgn_sigma_T2 >= 1.0_dp) then
          write(*,'(a,f10.6,a)') &
               '  (icgn) WARNING: ||T||_2 = ', upt%icgn_sigma_T2, &
               ' >= 1 — Neumann series not guaranteed to converge at this E0.'
          if (upt%icgn_selfenergy_order == 0) write(*,'(a)') &
               '  (icgn) NOTE: order=0 does not use the Q-Q coupling W, so this' // &
               ' warning is informational only for this run.'
       else if (upt%verbose > 0) then
          write(*,'(a,f10.6,a,l1)') &
               '  (icgn) ||T||_2 = ', upt%icgn_sigma_T2, &
               ', power-iter converged: ', upt%icgn_pi_converged
          if (upt%icgn_selfenergy_order == 0) write(*,'(a)') &
               '  (icgn) NOTE: order=0 ignores W — norm above is for diagnostics only.'
       end if
    end if

    if (upt%icgn_check_convergence .and. nqq > 0) call icgn_report_time('power iteration (parallel)', t_stage)
    call icgn_report_time('free S_full/pairs + convergence check', t_mark)
    ! ---- Compute Sigma and add to H_dense ----------------------------------
    allocate(Sigma(nred, nred))
    Sigma = cmplx(0.0_dp, 0.0_dp, kind=dp)

    E0_used = upt%icgn_E0
    if (E0_used == 0.0_dp) E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp

    ! Sigma = sum_{ord=0..N} H_PQ R (W R)^ord H_QP on the shared edge graphs
    ! (built above, also used by the convergence check): row-blocked, threaded,
    ! upper triangle + Hermitian mirror.  Summation order differs from the old
    ! per-row scan only at rounding level; the result does not depend on the
    ! number of threads.
    t_stage = icgn_wall()
    if (npq > 0) call icgn_sigma_blocked(nred, nstates_total, upt%icgn_selfenergy_order, 16, &
         p_ptr, p_col, p_val, c_ptr, c_row, c_val, s_ptr, s_col, s_wres, Sigma)
    deallocate(res_q, p_ptr, p_col, p_val, c_ptr, c_row, c_val, s_ptr, s_col, s_val, s_wres)
    call icgn_report_time('Sigma (row-blocked, threaded)', t_stage)
    end if   ! use_blocked / legacy edge-list Neumann path
    if (icgn_timing_on()) then
       write(*,'(a,i0,a,i0,a,i0,a,i0,a,es10.3)') '[ICGN-TIME] stats: nred=', nred, ' nstates_total=', &
            nstates_total, ' npq=', npq, ' nqq=', nqq, ' max|Sigma|=', maxval(abs(Sigma))
    end if
    t_mark = icgn_wall()
    deallocate(pq_p, pq_q, pq_v)
    if (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence) then
       if (allocated(qq_i)) deallocate(qq_i, qq_j, qq_v)
    end if
    deallocate(block_of_state, local_of_state, evals_q, ret_offset)
    deallocate(roff_blk, local_ret_idx)

    ! Add Sigma into H_dense and re-Hermitize
    H_dense = H_dense + Sigma
    nth = min(icgn_thread_budget(), max(1, nred))
    !$omp parallel do default(shared) private(i, j) num_threads(nth) schedule(static) if(nth > 1)
    do i = 1, nred
       do j = i+1, nred
          H_dense(i,j) = (H_dense(i,j) + conjg(H_dense(j,i))) / 2.0_dp
          H_dense(j,i) = conjg(H_dense(i,j))
       end do
    end do
    !$omp end parallel do
    deallocate(Sigma)
    call icgn_report_time('free edge lists + H_dense+=Sigma, hermitize', t_mark)
    t_mark = icgn_wall()

    ! Rebuild icgn_ham from H_dense
    call destroy_matrix(upt%icgn_ham)
    nnz = count(abs(H_dense) > 1.0e-14_dp)
    ! write(*,'(A,I0,I0,I0)') '[ICGN-TRACE] prepare: create icgn_ham nred/nnz=', nred, nred, nnz
    call create_matrix(upt%icgn_ham, nred, nred, nnz)
    ! write(*,'(A)') '[ICGN-TRACE] prepare: create icgn_ham returned'
    upt%icgn_ham%sparse_fmt = 'F'  ! always full — H_dense is complete after Sigma+Hermitianize
    upt%icgn_ham%Mi(1) = 1
    k = 0
    do i = 1, nred
       do j = 1, nred
          if (abs(H_dense(i,j)) > 1.0e-14_dp) then
             k = k + 1
             upt%icgn_ham%Mj(k) = j
             upt%icgn_ham%M(k)  = H_dense(i,j)
          end if
       end do
       upt%icgn_ham%Mi(i+1) = k + 1
    end do
    upt%icgn_ham%nnz = k
    deallocate(H_dense)
    call icgn_report_time('rebuild icgn_ham from dense (serial)', t_mark)

    end if  ! icgn_selfenergy_order >= 0 (Neumann correction)

    call icgn_report_time('TOTAL from block diagonalization to end of Neumann/Sigma', t_all)

    ! When the Neumann branch is skipped (selfenergy_order < 0) nothing above has
    ! released the pair matrices (v/h/t) nor the full block eigensystems.
    ! Free them here; both calls are no-ops if the Neumann branch already did it.
    if (allocated(pairs)) call destroy_pairs(pairs)
    do i = 1, nb
       if (associated(upt%icgn_blocks(i)%S_full))     deallocate(upt%icgn_blocks(i)%S_full)
       if (associated(upt%icgn_blocks(i)%evals_full)) deallocate(upt%icgn_blocks(i)%evals_full)
       nullify(upt%icgn_blocks(i)%S_full, upt%icgn_blocks(i)%evals_full)
    end do

    ! Build icgn_U (identity)
    call create_matrix(upt%icgn_U, nred, nred, nred)
    upt%icgn_U%sparse_fmt = 'F'; upt%icgn_U%Mi(1) = 1
    do i = 1, nred
       upt%icgn_U%Mj(i) = i; upt%icgn_U%M(i) = (1.0_dp, 0.0_dp)
       upt%icgn_U%Mi(i+1) = i + 1
    end do
    upt%icgn_U%nnz = nred

    upt%icgn_ready = .true.
       call icgn_log_info(upt, 'icgn', n, total_ret, &
          upt%icgn_num_blocks, upt%icgn_cut_fraction, upt%icgn_sigma_T2, &
          upt%icgn_pi_converged, upt%icgn_check_convergence)

    ! Cleanup
    ! write(*,'(A,L1,A,L1,A,I0,A,I0)') '[ICGN-TRACE] prepare final ham M/Mi=', associated(upt%icgn_ham%M), '/', &
         ! associated(upt%icgn_ham%Mi), ' nrow=', upt%icgn_ham%nrow, ' nnz=', upt%icgn_ham%nnz
    ! write(*,'(A,L1,A,L1,A,I0,A,I0)') '[ICGN-TRACE] prepare final U M/Mi=', associated(upt%icgn_U%M), '/', &
         ! associated(upt%icgn_U%Mi), ' nrow=', upt%icgn_U%nrow, ' nnz=', upt%icgn_U%nnz
    ! write(*,'(A,L1,A,L1)') '[ICGN-TRACE] prepare final original snapshot ham M/Mi=', associated(icgn_original_ham%M), '/', &
         ! associated(icgn_original_ham%Mi)
    deallocate(atom_of, local_of, offsets, bsize, label, row_of)
    deallocate(counts, edge_weight, xadj, adjncy, adjwgt, vwgt, part, cursor)
    ! write(*,'(A)') '[ICGN-TRACE] prepare EXIT cleanup done'

  end subroutine icgn_prepare

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
   if (has_norm) then
         write(line,'(a,f12.6,a,l1)') 'Neumann norm=', sigma_t2, &
         ', power iteration converged=', pi_converged
         call icgn_log_important(upt, trim(line))
    end if
  end subroutine icgn_log_info

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
  ! Stable counting sort: groups the indices 1..n by key(i) in 1..nkeys.
  ! Entries of group g are ord(ptr(g):ptr(g+1)-1), in ascending index order.
  subroutine icgn_group_by(key, n, nkeys, ptr, ord)
    integer, intent(in) :: key(:), n, nkeys
    integer, allocatable, intent(out) :: ptr(:), ord(:)
    integer, allocatable :: fill(:)
    integer :: i, g
    allocate(ptr(nkeys+1), ord(max(1,n)), fill(nkeys))
    ptr = 0
    do i = 1, n
       ptr(key(i)+1) = ptr(key(i)+1) + 1
    end do
    ptr(1) = 1
    do g = 1, nkeys
       ptr(g+1) = ptr(g+1) + ptr(g)
    end do
    fill = ptr(1:nkeys)
    do i = 1, n
       g = key(i); ord(fill(g)) = i; fill(g) = fill(g) + 1
    end do
    deallocate(fill)
  end subroutine icgn_group_by

  ! C = op(A) * B with BLAS zgemm (opa = 'N' or 'C').  Much faster than the
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

  ! Collect the dense physical block couplings H_ab (nrow_a x nrow_b) of every
  ! pair of blocks.  Same pair order and same accumulation order as the old
  ! serial scan (so the result is identical), but the dense allocation/zeroing
  ! and scatter run in parallel over pairs: a serial pass only discovers the
  ! pairs and buckets the CSR entries per pair.
  subroutine icgn_collect_pairs(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(in) :: upt
    integer, intent(in) :: atom_of(:), label(:), local(:)
    type(ICGNPair), intent(inout) :: pairs(:)
    integer, intent(out) :: npair, ierr
    integer :: r, k, c, ia, ib, slot, ne, e, m, i, nth, saved_levels
    integer, allocatable :: pair_map(:,:), e_slot(:), e_row(:), e_k(:), e_ord(:), p_start(:), p_fill(:)
    ierr = 0; npair = 0
    allocate(pair_map(upt%icgn_num_blocks, upt%icgn_num_blocks)); pair_map = 0
    allocate(e_slot(max(1,upt%ham%nnz)), e_row(max(1,upt%ham%nnz)), e_k(max(1,upt%ham%nnz)))
    ne = 0
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          if (upt%icgn_blocks(ia)%nrow == 0 .or. upt%icgn_blocks(ib)%nrow == 0) cycle
          slot = pair_map(min(ia,ib), max(ia,ib))
          if (slot == 0) then
             npair = npair + 1
             if (npair > size(pairs)) then; ierr = 11; return; end if
             pairs(npair)%a = min(ia,ib); pairs(npair)%b = max(ia,ib)
             pair_map(min(ia,ib), max(ia,ib)) = npair; slot = npair
          end if
          ne = ne + 1
          if (ia < ib) then; e_slot(ne) = slot; else; e_slot(ne) = -slot; end if
          e_row(ne) = r; e_k(ne) = k
       end do
    end do
    deallocate(pair_map)
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

  ! True when env UPT_CG_TIMING is set (any non-empty value).
  function icgn_timing_on() result(on)
    logical :: on
    character(len=8) :: s
    integer :: ln, ios
    call get_environment_variable('UPT_CG_TIMING', s, length=ln, status=ios)
    on = (ln > 0 .and. (ios == 0 .or. ios == -1))
  end function icgn_timing_on

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

  ! ---------------------------------------------------------------------------
  ! Block-dense Neumann correction.
  !
  ! The P-Q and Q-Q couplings between blocks are nearly dense (the block
  ! couplings of a Hamiltonian in the eigenbasis of the blocks), so the Neumann
  ! series is evaluated with dense BLAS3 products on the block structure instead
  ! of with edge lists:
  !     A_0 = G R,   A_{k+1} = A_k W R,   Sigma = sum_k A_k G^H
  ! where G (P x Q), W (Q x Q) and the diagonal resolvent R are stored per block
  ! (pair), with the same filtering (|g| >= 1e-14 and icgn_coupling_survives) as
  ! the edge-list path.  Memory is proportional to the number of coupled block
  ! pairs.  The legacy edge-list path stays available with UPT_CG_NEUMANN=old.
  ! ---------------------------------------------------------------------------
  function icgn_neumann_blocked() result(on)
    logical :: on
    character(len=8) :: s
    integer :: ln, ios
    on = .true.
    call get_environment_variable('UPT_CG_NEUMANN', s, length=ln, status=ios)
    if (ios == 0 .and. ln >= 3) then
       if (s(1:3) == 'old') on = .false.
    end if
  end function icgn_neumann_blocked

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

  subroutine icgn_build_neumann_blocks(upt, pairs, npair, local_ret_idx, E0, need_w, &
       bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, npq_nz, nqq_nz)
    type(OUPT), intent(in) :: upt
    type(ICGNPair), intent(in) :: pairs(:)
    integer, intent(in) :: npair, local_ret_idx(:,:)
    real(dp), intent(in) :: E0
    logical, intent(in) :: need_w
    type(ICGNQBlock), allocatable, intent(out) :: bq(:)
    type(ICGNQPair), allocatable, intent(out) :: qp(:)
    integer, allocatable, intent(out) :: nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, allocatable, intent(out) :: nbr_isA(:)
    integer(8), intent(out) :: npq_nz, nqq_nz
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    real(dp), parameter :: thr = 1.0e-14_dp
    integer :: nb, b, i, j, k, ia, ib, jq, kq, ra, rb, e, np_a, np_b, nq_a, nq_b
    integer :: nth, saved_levels, nra, nrb
    integer, allocatable :: cnt(:), fillp(:), eA(:), eB(:)
    complex(dp) :: v
    complex(dp), pointer :: g(:,:)
    complex(dp), allocatable :: tmp(:,:)
    logical :: own_g
    integer(8) :: nzp, nzq

    nb = upt%icgn_num_blocks
    ! ---- per-block Q data ----
    allocate(bq(nb))
    do b = 1, nb
       bq(b)%np = upt%icgn_blocks(b)%nret
       bq(b)%nq = upt%icgn_blocks(b)%nrow - upt%icgn_blocks(b)%nret
       allocate(bq(b)%qidx(max(0,bq(b)%nq)), bq(b)%res(max(0,bq(b)%nq)))
       kq = 0
       do k = 1, upt%icgn_blocks(b)%nrow
          if (local_ret_idx(b, k) == 0) then
             kq = kq + 1
             bq(b)%qidx(kq) = k
             bq(b)%res(kq) = 1.0_dp / (E0 - upt%icgn_blocks(b)%evals_full(k))
          end if
       end do
    end do

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
    nth = min(icgn_thread_budget(), max(1, npair))
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel do default(shared) private(i, ia, ib, jq, kq, k, ra, rb, v, g, tmp, own_g, &
    !$omp&   np_a, np_b, nq_a, nq_b, nra, nrb) reduction(+:nzp, nzq) num_threads(nth) &
    !$omp&   schedule(dynamic,1) if(nth > 1)
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       np_a = bq(ia)%np; nq_a = bq(ia)%nq; np_b = bq(ib)%np; nq_b = bq(ib)%nq
       nra = upt%icgn_blocks(ia)%nrow; nrb = upt%icgn_blocks(ib)%nrow
       own_g = .false.
       if (associated(pairs(i)%g)) then
          g => pairs(i)%g                      ! full projection kept (acquaintance path)
       else
          allocate(g(nra, nrb)); own_g = .true.
          call icgn_project_sparse(pairs(i)%h, upt%icgn_blocks(ia)%S_full, &
               upt%icgn_blocks(ib)%S_full, g)
       end if

       ! G(P_A, Q_B)
       allocate(qp(i)%gab(np_a, nq_b))
       do kq = 1, nq_b
          k = bq(ib)%qidx(kq)
          do jq = 1, np_a
             ra = upt%icgn_blocks(ia)%retained_idx(jq)
             v = g(ra, k)
             if (abs(v) < thr) then
                qp(i)%gab(jq, kq) = czero
             else if (.not. icgn_coupling_survives(v, upt%icgn_blocks(ia)%evals_full(ra), &
                  upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) then
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
             rb = upt%icgn_blocks(ib)%retained_idx(jq)
             v = g(k, rb)
             if (abs(v) < thr) then
                qp(i)%gba(jq, kq) = czero
             else if (.not. icgn_coupling_survives(v, upt%icgn_blocks(ia)%evals_full(k), &
                  upt%icgn_blocks(ib)%evals_full(rb), upt%icgn_epsilon)) then
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
                else if (.not. icgn_coupling_survives(v, upt%icgn_blocks(ia)%evals_full(ra), &
                     upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) then
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
    call icgn_restore_nesting(saved_levels)
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

  ! Sigma = sum_{ord=0..order} H_PQ R (W R)^ord H_QP on the block structure.
  ! Parallel over the P blocks a (each owns the block row Sigma(P_a, :)); only
  ! the blocks a2 >= a are accumulated, the rest follows from Hermiticity.
  subroutine icgn_sigma_blockdense(nb, order, maxdim, roff_blk, bq, qp, &
       nbr_ptr, nbr_blk, nbr_pair, nbr_isA, nbr_rev, Sigma)
    integer, intent(in) :: nb, order, maxdim
    integer, intent(in) :: roff_blk(:), nbr_ptr(:), nbr_blk(:), nbr_pair(:), nbr_rev(:)
    logical, intent(in) :: nbr_isA(:)
    type(ICGNQBlock), intent(in) :: bq(:)
    type(ICGNQPair), intent(in) :: qp(:)
    complex(dp), intent(inout) :: Sigma(:,:)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)
    real(dp), parameter :: thr = 1.0e-14_dp
    integer :: a, b, b2, a2, e, f, ord, kq, nth, saved_levels, nred, i, gp
    integer :: r1, r2, c1, c2, npa, nqb, nqb2
    type(ICGNMat), allocatable :: Acur(:), Anext(:)
    complex(dp), pointer :: Gm(:,:), Wm(:,:)
    character(len=1) :: opw
    real(dp) :: amax
    logical :: any_cur

    nred = size(Sigma, 1)
    nth = min(icgn_thread_budget(), max(1, nb))
    nth = icgn_outer_threads(nb, nth, maxdim)
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel default(shared) num_threads(nth) if(nth > 1) &
    !$omp& private(a, b, b2, a2, e, f, ord, kq, r1, r2, c1, c2, npa, nqb, nqb2, Gm, Wm, opw, amax, &
    !$omp&         any_cur, Acur, Anext)
    allocate(Acur(nb), Anext(nb))
    !$omp do schedule(dynamic,1)
    do a = 1, nb
       npa = bq(a)%np
       if (npa == 0) cycle
       r1 = roff_blk(a); r2 = r1 + npa - 1
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

          ! close the chain:  Sigma(P_a, P_a2) += A(a,b) G(P_a2, Q_b)^H,  a2 >= a
          do b = 1, nb
             if (.not. associated(Acur(b)%m)) cycle
             do f = nbr_ptr(b), nbr_ptr(b+1)-1
                a2 = nbr_blk(f)
                if (a2 < a) cycle
                if (bq(a2)%np == 0) cycle
                Gm => icgn_qg(qp, nbr_pair, nbr_isA, nbr_rev(f))     ! G(P_a2, Q_b)
                c1 = roff_blk(a2); c2 = c1 + bq(a2)%np - 1
                call icgn_zgemm_gen('N', 'C', cone, Acur(b)%m, Gm, cone, Sigma(r1:r2, c1:c2))
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
    call icgn_restore_nesting(saved_levels)

    ! lower triangle from Hermiticity
    nth = min(icgn_thread_budget(), max(1, nred))
    !$omp parallel do default(shared) private(i, gp) num_threads(nth) schedule(static) if(nth > 1)
    do gp = 1, nred
       do i = gp+1, nred
          Sigma(i, gp) = conjg(Sigma(gp, i))
       end do
    end do
    !$omp end parallel do
  end subroutine icgn_sigma_blockdense

  ! out_X = [R] * sum_Y (W_{Y,X})^T in_Y   (the operator W^T = conj(W) of the power iteration)
  subroutine icgn_blk_wt_apply(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, qoff, vin, vout, apply_res)
    integer, intent(in) :: nb, nbr_ptr(:), nbr_blk(:), nbr_pair(:), qoff(:)
    logical, intent(in) :: nbr_isA(:), apply_res
    type(ICGNQBlock), intent(in) :: bq(:)
    type(ICGNQPair), intent(in) :: qp(:)
    complex(dp), intent(in) :: vin(:)
    complex(dp), intent(out) :: vout(:)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp), czero = (0.0_dp, 0.0_dp)
    integer :: X, Y, f, nqx, nqy, oX, oY, maxnq, nth, q
    complex(dp), allocatable :: tmp(:), yv(:)
    external :: zgemv

    maxnq = 1
    do X = 1, nb
       maxnq = max(maxnq, bq(X)%nq)
    end do
    nth = min(icgn_thread_budget(), max(1, nb))
    !$omp parallel default(shared) num_threads(nth) if(nth > 1) private(X, Y, f, nqx, nqy, oX, oY, q, tmp, yv)
    allocate(tmp(maxnq), yv(maxnq))
    !$omp do schedule(dynamic,1)
    do X = 1, nb
       nqx = bq(X)%nq
       if (nqx == 0) cycle
       oX = qoff(X)
       vout(oX+1:oX+nqx) = czero
       do f = nbr_ptr(X), nbr_ptr(X+1)-1
          Y = nbr_blk(f); nqy = bq(Y)%nq
          if (nqy == 0) cycle
          oY = qoff(Y)
          if (nbr_isA(f)) then
             ! X = A, Y = B, wab is (nq_X x nq_Y):  out_X += conj(wab) in_Y = conj(wab * conj(in_Y))
             do q = 1, nqy
                tmp(q) = conjg(vin(oY+q))
             end do
             call zgemv('N', nqx, nqy, cone, qp(nbr_pair(f))%wab, max(1,nqx), tmp, 1, czero, yv, 1)
             do q = 1, nqx
                vout(oX+q) = vout(oX+q) + conjg(yv(q))
             end do
          else
             ! X = B, Y = A, wab is (nq_Y x nq_X):  out_X += wab^T in_Y
             do q = 1, nqy
                tmp(q) = vin(oY+q)
             end do
             call zgemv('T', nqy, nqx, cone, qp(nbr_pair(f))%wab, max(1,nqy), tmp, 1, czero, yv, 1)
             do q = 1, nqx
                vout(oX+q) = vout(oX+q) + yv(q)
             end do
          end if
       end do
       if (apply_res) then
          do q = 1, nqx
             vout(oX+q) = vout(oX+q) * bq(X)%res(q)
          end do
       end if
    end do
    !$omp end do
    deallocate(tmp, yv)
    !$omp end parallel
    ! blocks without Q states keep zeros
    do X = 1, nb
       if (bq(X)%nq == 0) cycle
    end do
  end subroutine icgn_blk_wt_apply

  ! Power iteration for ||T||_2, T = R W, on the block structure (same iteration
  ! as icgn_power_iteration_csr: u = R W^T v, w = R u, v <- W^T w).
  subroutine icgn_power_iteration_blk(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, &
       maxiter, tol, sigma_out, converged_out)
    integer, intent(in) :: nb, nbr_ptr(:), nbr_blk(:), nbr_pair(:), maxiter
    logical, intent(in) :: nbr_isA(:)
    type(ICGNQBlock), intent(in) :: bq(:)
    type(ICGNQPair), intent(in) :: qp(:)
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: sigma_out
    logical, intent(out) :: converged_out
    integer :: k, X, j, ntot, q
    integer, allocatable :: qoff(:)
    real(dp) :: nrm, sigma_est, diff2, rr, ri
    complex(dp), allocatable :: v(:), u(:), w(:), v_next(:)
    real(dp), allocatable :: res(:)

    sigma_out = 0.0_dp; converged_out = .true.
    allocate(qoff(nb+1)); qoff(1) = 0
    do X = 1, nb
       qoff(X+1) = qoff(X) + bq(X)%nq
    end do
    ntot = qoff(nb+1)
    if (ntot == 0) return
    allocate(v(ntot), u(ntot), w(ntot), v_next(ntot), res(ntot))
    do X = 1, nb
       do q = 1, bq(X)%nq
          res(qoff(X)+q) = bq(X)%res(q)
       end do
    end do
    do j = 1, ntot
       call random_number(rr); call random_number(ri)
       v(j) = cmplx(rr - 0.5_dp, ri - 0.5_dp, kind=dp)
    end do
    nrm = sqrt(sum(real(v,dp)**2 + aimag(v)**2))
    if (nrm < 1.0e-14_dp) then
       deallocate(v, u, w, v_next, res, qoff); return
    end if
    v = v / nrm
    converged_out = .false.
    sigma_est = 0.0_dp
    do k = 1, maxiter
       call icgn_blk_wt_apply(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, qoff, v, u, .true.)
       sigma_est = sqrt(sum(real(u,dp)**2 + aimag(u)**2))
       do j = 1, ntot
          w(j) = u(j) * res(j)
       end do
       call icgn_blk_wt_apply(nb, bq, qp, nbr_ptr, nbr_blk, nbr_pair, nbr_isA, qoff, w, v_next, .false.)
       nrm = sqrt(sum(real(v_next,dp)**2 + aimag(v_next)**2))
       if (nrm < 1.0e-14_dp) then
          sigma_out = 0.0_dp; converged_out = .true.
          deallocate(v, u, w, v_next, res, qoff); return
       end if
       v_next = v_next / nrm
       diff2 = sum(real(v_next - v, dp)**2 + aimag(v_next - v)**2)
       v = v_next
       if (sqrt(diff2) < tol) then
          converged_out = .true.
          exit
       end if
    end do
    if (icgn_timing_on()) write(*,'(a,i0,a,l1)') '[ICGN-TIME] power iteration sweeps: ', &
         min(k, maxiter), ' converged=', converged_out
    sigma_out = sigma_est
    deallocate(v, u, w, v_next, res, qoff)
  end subroutine icgn_power_iteration_blk

  ! Warnings about the Neumann convergence estimate (same messages as the edge-list path)
  subroutine icgn_report_convergence(upt)
    type(OUPT), intent(in) :: upt
    if (.not. upt%icgn_pi_converged) then
       write(*,'(a,i0,a)') &
            '  (icgn) WARNING: power iteration did not converge in ', &
            upt%icgn_pi_maxiter, ' sweeps — ||T||_2 estimate may be unreliable.'
    end if
    if (upt%icgn_sigma_T2 >= 1.0_dp) then
       write(*,'(a,f10.6,a)') &
            '  (icgn) WARNING: ||T||_2 = ', upt%icgn_sigma_T2, &
            ' >= 1 — Neumann series not guaranteed to converge at this E0.'
       if (upt%icgn_selfenergy_order == 0) write(*,'(a)') &
            '  (icgn) NOTE: order=0 does not use the Q-Q coupling W, so this' // &
            ' warning is informational only for this run.'
    else if (upt%verbose > 0) then
       write(*,'(a,f10.6,a,l1)') &
            '  (icgn) ||T||_2 = ', upt%icgn_sigma_T2, &
            ', power-iter converged: ', upt%icgn_pi_converged
       if (upt%icgn_selfenergy_order == 0) write(*,'(a)') &
            '  (icgn) NOTE: order=0 ignores W — norm above is for diagnostics only.'
    end if
  end subroutine icgn_report_convergence

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
    !$    ! also cut the thread count seen by threaded BLAS/LAPACK (MKL follows it)
    !$    call omp_set_num_threads(avail)
    !$ end if
    call get_environment_variable('UPT_CG_THREADS', s, length=ln, status=ios)
    if (ios == 0 .and. ln > 0) then
       read(s(1:ln), *, iostat=ios) v
       if (ios == 0 .and. v >= 1) nth = min(nth, v)
    end if
    nth = max(1, nth)
  end function icgn_thread_budget

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

  ! Diagonalize block ib of icgn_blocks.  Thread-safe for distinct ib: it reads
  ! upt%ham and writes only upt%icgn_blocks(ib).  local(r) = position of row r
  ! inside its own block, blk_of(r) = block that owns row r (both cover all rows).
  subroutine icgn_diagonalize_block(upt, ib, local, blk_of, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, local(:), blk_of(:)
    integer, intent(out) :: ierr
    integer :: nn, k, r, c, jr
    complex(dp), allocatable :: h(:,:)
    real(dp), allocatable :: w(:)
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
   call coarse_eigh(h, w, ierr)
    if (ierr /= 0) return
    allocate(upt%icgn_blocks(ib)%evals_full(nn), upt%icgn_blocks(ib)%S_full(nn,nn))
    upt%icgn_blocks(ib)%evals_full = w
    upt%icgn_blocks(ib)%S_full     = h
    deallocate(h, w)
  end subroutine icgn_diagonalize_block

  ! pair_slot variant for icgn: allocates v(nrow_a, nrow_b) using icgn_blocks.
  integer function pair_slot_icgn(pairs, npair, a, b, upt, pair_map)
    type(ICGNPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(in) :: a, b
    type(OUPT), intent(in) :: upt
    integer, intent(inout) :: pair_map(:,:)
    if (pair_map(a,b) /= 0) then
       pair_slot_icgn = pair_map(a,b); return
    end if
    npair = npair + 1
    if (npair > size(pairs)) then; pair_slot_icgn = 0; return; end if
    pairs(npair)%a = a; pairs(npair)%b = b
    allocate(pairs(npair)%h(upt%icgn_blocks(a)%nrow, upt%icgn_blocks(b)%nrow))
    pairs(npair)%h = (0.0_dp, 0.0_dp)
    pair_map(a,b) = npair
    pair_slot_icgn = npair
  end function pair_slot_icgn

  integer function pair_slot_icgn_fast(pairs,npair,a,b,upt,pair_map)
    type(ICGNPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(in) :: a,b
    type(OUPT), intent(in) :: upt
    integer, intent(inout) :: pair_map(:,:)

    if (pair_map(a,b) /= 0) then
       pair_slot_icgn_fast = pair_map(a,b)
       return
    end if

    npair = npair + 1
    if (npair > size(pairs)) then; pair_slot_icgn_fast = 0; return; end if
    pairs(npair)%a = a; pairs(npair)%b = b
    allocate(pairs(npair)%t(upt%icgn_blocks(a)%nrow, upt%icgn_blocks(b)%nret))
    pairs(npair)%t = (0.0_dp, 0.0_dp)
    pair_map(a,b) = npair
    pair_slot_icgn_fast = npair
  end function pair_slot_icgn_fast

  subroutine project_accumulated_pair_icgn(p, qa)
    type(ICGNPair), intent(inout) :: p
    complex(dp), intent(in) :: qa(:,:)

    allocate(p%v(size(qa,2), size(p%t,2)))
    call icgn_zgemm('C', qa, p%t, p%v)
    deallocate(p%t)
  end subroutine project_accumulated_pair_icgn

  ! Build reduced Hamiltonian for ICGN.
  subroutine build_icgn_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, reuse_full, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(ICGNPair), allocatable, intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    ! .true.: pairs already hold the full projected coupling v (and h) from
    ! icgn_prepare (acquaintance path); only slicing to retained states is left.
    logical, intent(in) :: reuse_full
    integer, intent(out) :: ierr
    integer :: i, j, k, r, c, a, b, ia, ib, nnz, pos, slot, nred, nth, saved_levels
    integer :: ne, m, e
    integer, allocatable :: e_slot(:), e_row(:), e_k(:), e_ord(:), p_start(:), p_fill(:)
    real(dp) :: t_stage
    integer, allocatable :: roff(:), rowcount(:), next(:), pair_map(:,:)
    complex(dp), allocatable :: g_full(:,:)
    ierr = 0; nred = upt%icgn_reduced_dim
    allocate(roff(upt%icgn_num_blocks+1)); roff(1) = 1
    do i = 1, upt%icgn_num_blocks; roff(i+1) = roff(i) + upt%icgn_blocks(i)%nret; end do

    if (upt%icgn_add_core_acquaintances == 0 .and. upt%icgn_selfenergy_order < 0) then
       ! Fast path: accumulate T = H_AB * Q_B directly from the global CSR,
       ! exactly like ICGN.  No dense physical H_AB is needed because Neumann
       ! corrections are disabled and the pairs are used only to build H_PP.
       npair = 0
       allocate(pair_map(upt%icgn_num_blocks, upt%icgn_num_blocks)); pair_map = 0
       t_stage = icgn_wall()
       ! Phase A (serial, integer work only): discover the pairs (same order as
       ! the old scan) and record every off-block entry as (row r, CSR position k)
       ! with a signed slot: +slot if label(r) < label(c), -slot otherwise.
       allocate(e_slot(max(1,upt%ham%nnz)), e_row(max(1,upt%ham%nnz)), e_k(max(1,upt%ham%nnz)))
       ne = 0
       do r = 1, upt%ham%nrow
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
             if (a == b) cycle
             if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             ia = min(a,b); ib = max(a,b)
             slot = pair_map(ia, ib)
             if (slot == 0) then
                npair = npair + 1
                if (npair > size(pairs)) then; ierr = 11; return; end if
                pairs(npair)%a = ia; pairs(npair)%b = ib
                pair_map(ia, ib) = npair; slot = npair
             end if
             ne = ne + 1
             if (a < b) then; e_slot(ne) = slot; else; e_slot(ne) = -slot; end if
             e_row(ne) = r; e_k(ne) = k
          end do
       end do
       deallocate(pair_map)
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
                    upt%icgn_epsilon)) pairs(i)%v(j,k) = (0.0_dp, 0.0_dp)
             end do
          end do
       end do
       !$omp end parallel do
       call icgn_restore_nesting(saved_levels)
       call icgn_report_time('  fast path: accumulate + project (parallel)', t_stage)
       deallocate(e_slot, e_row, e_k, e_ord, p_start, p_fill)
    else if (upt%icgn_add_core_acquaintances == 0) then
       ! Zero-acquaintance path with Neumann enabled: keep the full physical
       ! pair matrices because the correction needs P-Q and Q-Q couplings.
       nth = min(icgn_thread_budget(), max(1, npair))
       call icgn_flatten_nesting(saved_levels)
       !$omp parallel do default(shared) private(i, a, b, j, k) num_threads(nth) schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
          call project_pair(pairs(i), upt%icgn_blocks(a)%q, &
               upt%icgn_blocks(b)%q, .true.)
          do j = 1, upt%icgn_blocks(a)%nret
             do k = 1, upt%icgn_blocks(b)%nret
                if (.not. icgn_coupling_survives(pairs(i)%v(j,k), &
                    upt%icgn_blocks(a)%eval(j), upt%icgn_blocks(b)%eval(k), &
                    upt%icgn_epsilon)) pairs(i)%v(j,k) = (0.0_dp, 0.0_dp)
             end do
          end do
       end do
       !$omp end parallel do
       call icgn_restore_nesting(saved_levels)
    else
       if (.not. reuse_full) then
          if (allocated(pairs)) call destroy_pairs(pairs)
          allocate(pairs(max(1, upt%ham%nnz)))
          ! Collect physical H_ab first; project each pair once at matrix level.
          call icgn_collect_pairs(upt, atom_of, label, local, pairs, npair, ierr)
          if (ierr /= 0) return
       end if
       ! Per pair: project H_ab with the full eigenbasis (unless already done in
       ! icgn_prepare), then slice to the retained states.  Pairs are
       ! independent; g_full is a per-thread temporary.
       nth = min(icgn_thread_budget(), max(1, npair))
       call icgn_flatten_nesting(saved_levels)
       !$omp parallel do default(shared) private(i, a, b, j, k, g_full) num_threads(nth) &
       !$omp& schedule(dynamic,1) if(nth > 1)
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
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
                       upt%icgn_blocks(b)%eval(j), upt%icgn_epsilon)) &
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
                       upt%icgn_blocks(b)%eval(j), upt%icgn_epsilon)) &
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
    end if

    t_stage = icgn_wall()
    ! Build CSR.  Count exactly the entries that emit_pair will write:
    ! zero-valued reduced couplings are omitted to keep the sparse matrix small.
    allocate(rowcount(nred), next(nred)); rowcount = 1
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       do j = 1, size(pairs(i)%v,1)
          do k = 1, size(pairs(i)%v,2)
             if (abs(pairs(i)%v(j,k)) == 0.0_dp) cycle
             if (upt%ham%sparse_fmt /= 'L') rowcount(roff(a)+j-1) = rowcount(roff(a)+j-1) + 1
             if (upt%ham%sparse_fmt == 'F' .or. upt%ham%sparse_fmt == 'L') &
                  rowcount(roff(b)+k-1) = rowcount(roff(b)+k-1) + 1
          end do
       end do
    end do
    nnz = sum(rowcount); call create_matrix(upt%icgn_ham, nred, nred, nnz)
    upt%icgn_ham%sparse_fmt = upt%ham%sparse_fmt; upt%icgn_ham%Mi(1) = 1
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
       a = pairs(i)%a; b = pairs(i)%b
       call emit_pair(upt%icgn_ham, pairs(i), roff(a), roff(b), upt%ham%sparse_fmt, next)
    end do
    upt%icgn_ham%nnz = nnz
    deallocate(roff, rowcount, next)
    call icgn_report_time('  assemble reduced CSR (serial)', t_stage)
  end subroutine build_icgn_reduced_hamiltonian

  ! Lift ICGN eigenvectors from reduced basis back to physical space.
  subroutine icgn_lift(upt, reduced, physical)
    type(OUPT), intent(in) :: upt
    complex(dp), intent(in) :: reduced(:,:)
    complex(dp), intent(out) :: physical(:,:)
    integer :: i, j, off
    physical = (0.0_dp, 0.0_dp); off = 1
    do i = 1, size(upt%icgn_blocks)
       if (upt%icgn_blocks(i)%nret > 0) then
          do j = 1, size(upt%icgn_blocks(i)%rows)
             physical(upt%icgn_blocks(i)%rows(j),:) = &
                  matmul(upt%icgn_blocks(i)%q(j,:), reduced(off:off+upt%icgn_blocks(i)%nret-1,:))
          end do
       end if
       off = off + upt%icgn_blocks(i)%nret
    end do
  end subroutine icgn_lift

  ! ---------------------------------------------------------------------------
  ! Power iteration to estimate ||T||_2, T = (E0-D)^{-1} W
  ! where D = diag(evals_q restricted to Q) and W = Q-Q coupling matrix
  ! stored as the edge list (qq_i, qq_j, qq_v, nqq), BOTH directions.
  !
  ! Algorithm: standard power method on T^dagger T (same as Julia prototype).
  !   v_0   = random unit vector on Q states with at least one Q-Q edge
  !   u     = T * v_k         (apply_T)
  !   w     = T^dagger * u    (apply_T_adjoint)
  !   v_{k+1} = w / ||w||
  !   sigma_est = ||u|| = ||T v_k||   (Rayleigh quotient)
  !   stop when ||v_{k+1} - v_k|| < tol or k == maxiter
  !
  ! Uses dense arrays indexed by flat state index (faster than Dict).
  ! Only states appearing in the Q-Q graph need nonzero entries.
  subroutine icgn_power_iteration(evals_q, nstates_total, &
       qq_i_arr, qq_j_arr, qq_v_arr, nqq, E0, maxiter, tol, &
       sigma_out, converged_out)
    real(dp),     intent(in)  :: evals_q(:)
    integer,      intent(in)  :: nstates_total
    integer,      intent(in)  :: qq_i_arr(:), qq_j_arr(:)
    complex(dp),  intent(in)  :: qq_v_arr(:)
    integer,      intent(in)  :: nqq, maxiter
    real(dp),     intent(in)  :: E0, tol
    real(dp),     intent(out) :: sigma_out
    logical,      intent(out) :: converged_out

    integer  :: k, e, qi, qj
    real(dp) :: nrm, sigma_est, diff2, res_qi
    complex(dp) :: amp
    complex(dp), allocatable :: v(:), u(:), w(:), v_next(:)

    sigma_out     = 0.0_dp
    converged_out = .true.

    if (nqq == 0) return

    allocate(v(nstates_total), u(nstates_total), w(nstates_total), v_next(nstates_total))

    ! ---- initialise v as random unit vector on Q-Q support -----------------
    call random_number_cx(v, nstates_total, qq_i_arr, nqq)
    nrm = sqrt(real(dot_product(v, v), kind=dp))
    if (nrm < 1.0e-14_dp) then; sigma_out = 0.0_dp; converged_out = .true.; return; end if
    v = v / nrm

    converged_out = .false.
    sigma_est     = 0.0_dp

    do k = 1, maxiter

       ! u = T * v = (E0 - D)^{-1} W v
       ! Step 1: w_tmp = W * v  (sparse matvec over qq edges, both directions stored)
       u = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do e = 1, nqq
          qi = qq_i_arr(e); qj = qq_j_arr(e)
          u(qj) = u(qj) + qq_v_arr(e) * v(qi)
       end do
       ! Step 2: scale by resolvent  u(q) = u(q) / (E0 - evals_q(q))
       do qi = 1, nstates_total
          if (abs(u(qi)) < 1.0e-300_dp) cycle
          res_qi = 1.0_dp / (E0 - evals_q(qi))
          u(qi) = u(qi) * res_qi
       end do

       ! sigma_est = ||u|| = ||T v||
       sigma_est = sqrt(real(dot_product(u, u), kind=dp))

       ! w = T^dagger * u = W^dagger * (E0-D)^{-1} * u
       ! Step 1: scale by resolvent (E0-D is real diagonal, self-adjoint)
       w = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do qi = 1, nstates_total
          if (abs(u(qi)) < 1.0e-300_dp) cycle
          w(qi) = u(qi) / (E0 - evals_q(qi))
       end do
       ! Step 2: apply W^dagger via conjugated edge walk
       v_next = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do e = 1, nqq
          qi = qq_i_arr(e); qj = qq_j_arr(e)
          v_next(qi) = v_next(qi) + conjg(qq_v_arr(e)) * w(qj)
       end do

       ! normalize v_{k+1}
       nrm = sqrt(real(dot_product(v_next, v_next), kind=dp))
       if (nrm < 1.0e-14_dp) then
          sigma_out = 0.0_dp; converged_out = .true.; return
       end if
       v_next = v_next / nrm

       ! convergence: ||v_{k+1} - v_k||
       diff2 = real(dot_product(v_next - v, v_next - v), kind=dp)
       v = v_next

       if (sqrt(diff2) < tol) then
          converged_out = .true.
          exit
       end if
    end do

    sigma_out = sigma_est
    deallocate(v, u, w, v_next)
  end subroutine icgn_power_iteration

  ! Fill v with random complex values on the support of qq edges, zeros elsewhere
  subroutine random_number_cx(v, n, qi_arr, nqq)
    complex(dp), intent(out) :: v(:)
    integer,     intent(in)  :: n, nqq, qi_arr(:)
    real(dp) :: rr, ri
    integer  :: e
    v = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do e = 1, nqq
       call random_number(rr); call random_number(ri)
       v(qi_arr(e)) = cmplx(rr - 0.5_dp, ri - 0.5_dp, kind=dp)
    end do
  end subroutine random_number_cx

  ! ---------------------------------------------------------------------------
  ! Edge graphs for the Neumann/Sigma stage, built ONCE and shared by the
  ! convergence check and the Sigma kernel (contiguous CSR-ordered arrays, so
  ! the kernels stream through memory instead of gathering through index maps).
  !   P-grouped : for every P row p, its Q neighbours q and  g(p,q)*R(q)
  !   Q-grouped : for every Q state q, the P rows p (ascending) and conj(g(p,q))
  !   S (W)     : for every source Q state q, its targets q', raw coupling w and
  !               w*R(q')  (only if build_w)
  ! ---------------------------------------------------------------------------
  subroutine icgn_build_edge_csr(npq, pq_p, pq_q, pq_v, nqq, qq_i, qq_j, qq_v, &
       nred, nstates, res_q, build_w, &
       p_ptr, p_col, p_val, c_ptr, c_row, c_val, s_ptr, s_col, s_val, s_wres)
    integer, intent(in) :: npq, nqq, nred, nstates
    integer, intent(in) :: pq_p(:), pq_q(:), qq_i(:), qq_j(:)
    complex(dp), intent(in) :: pq_v(:), qq_v(:)
    real(dp), intent(in) :: res_q(:)
    logical, intent(in) :: build_w
    integer, allocatable, intent(out) :: p_ptr(:), p_col(:), c_ptr(:), c_row(:), s_ptr(:), s_col(:)
    complex(dp), allocatable, intent(out) :: p_val(:), c_val(:), s_val(:), s_wres(:)
    integer, allocatable :: ord(:), fill(:)
    integer :: e, i, p, q

    ! ---- P-grouped (stable: entries of a row keep their original order) ----
    call icgn_group_by(pq_p, npq, nred, p_ptr, ord)
    allocate(p_col(max(1,npq)), p_val(max(1,npq)))
    do e = 1, npq
       i = ord(e)
       p_col(e) = pq_q(i)
       p_val(e) = pq_v(i) * res_q(pq_q(i))
    end do

    ! ---- Q-grouped, rows p ascending inside each q (needed for the triangle cut) ----
    allocate(c_ptr(nstates+1), c_row(max(1,npq)), c_val(max(1,npq)), fill(nstates))
    c_ptr = 0
    do e = 1, npq
       c_ptr(pq_q(e)+1) = c_ptr(pq_q(e)+1) + 1
    end do
    c_ptr(1) = 1
    do q = 1, nstates
       c_ptr(q+1) = c_ptr(q+1) + c_ptr(q)
    end do
    fill = c_ptr(1:nstates)
    do p = 1, nred
       do e = p_ptr(p), p_ptr(p+1)-1
          i = ord(e); q = pq_q(i)
          c_row(fill(q)) = p
          c_val(fill(q)) = conjg(pq_v(i))
          fill(q) = fill(q) + 1
       end do
    end do
    deallocate(ord, fill)

    ! ---- W grouped by source state ----
    if (build_w) then
       call icgn_group_by(qq_i, nqq, nstates, s_ptr, ord)
       allocate(s_col(max(1,nqq)), s_val(max(1,nqq)), s_wres(max(1,nqq)))
       !$omp parallel do default(shared) private(e, i) schedule(static)
       do e = 1, nqq
          i = ord(e)
          s_col(e)  = qq_j(i)
          s_val(e)  = qq_v(i)
          s_wres(e) = qq_v(i) * res_q(qq_j(i))
       end do
       !$omp end parallel do
       deallocate(ord)
    else
       allocate(s_ptr(nstates+1), s_col(1), s_val(1), s_wres(1))
       s_ptr = 1; s_col = 1; s_val = (0.0_dp, 0.0_dp); s_wres = (0.0_dp, 0.0_dp)
    end if
  end subroutine icgn_build_edge_csr

  ! ---------------------------------------------------------------------------
  ! Sigma = sum_{ord=0..order} H_PQ R (W R)^ord H_QP, row-blocked and threaded.
  !
  ! The rows of Sigma are processed in blocks of rb rows at once: the amplitudes
  ! of the rb rows are stored side by side (A(r,q), r contiguous), so every
  ! coupling read from memory is used rb times.  The sparse support of the
  ! amplitudes is tracked (union over the block), so the work stays proportional
  ! to the paths actually followed.  Only the upper triangle (gp >= first row of
  ! the block) is accumulated; the lower triangle is filled from Hermiticity.
  ! Thread-private scratch only; every block owns its rows of Sigma.
  ! ---------------------------------------------------------------------------
  subroutine icgn_sigma_blocked(nred, nstates, order, rb_req, &
       p_ptr, p_col, p_val, c_ptr, c_row, c_val, s_ptr, s_col, s_wres, Sigma)
    integer, intent(in) :: nred, nstates, order, rb_req
    integer, intent(in) :: p_ptr(:), p_col(:), c_ptr(:), c_row(:), s_ptr(:), s_col(:)
    complex(dp), intent(in) :: p_val(:), c_val(:), s_wres(:)
    complex(dp), intent(inout) :: Sigma(:,:)
    real(dp), parameter :: thr2 = 1.0e-28_dp          ! (1e-14)^2
    real(dp), parameter :: scratch_bytes = 2.56e8_dp   ! per-thread budget for A and Nx
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    integer :: rb, nblk, nth, saved_levels
    integer :: ib, ip0, nr, r, ord, tt, ee, q, qb, nsupp, nsupp2, lo, hi, mid, gp, i
    real(dp) :: amax2, m2
    logical :: anynz
    complex(dp), allocatable :: A(:,:), Nx(:,:), a_loc(:)
    integer, allocatable :: supp(:), supp2(:)
    logical, allocatable :: mark(:), mark2(:)

    if (nred <= 0) return
    rb = max(1, min(rb_req, nred, int(scratch_bytes / (32.0_dp * real(max(1,nstates), dp)))))
    nblk = (nred + rb - 1) / rb
    nth = min(icgn_thread_budget(), max(1, nblk))
    call icgn_flatten_nesting(saved_levels)
    !$omp parallel default(shared) num_threads(nth) if(nth > 1) &
    !$omp& private(ib, ip0, nr, r, ord, tt, ee, q, qb, nsupp, nsupp2, lo, hi, mid, gp, i, &
    !$omp&         amax2, m2, anynz, A, Nx, a_loc, supp, supp2, mark, mark2)
    allocate(A(rb, nstates), Nx(rb, nstates), a_loc(rb))
    allocate(supp(nstates), supp2(nstates), mark(nstates), mark2(nstates))
    A = czero; Nx = czero; mark = .false.; mark2 = .false.
    !$omp do schedule(dynamic,1)
    do ib = 1, nblk
       ip0 = (ib-1)*rb + 1
       nr  = min(rb, nred - ip0 + 1)

       ! Order-0 amplitudes of the rows of this block: A(r,q) = g(ip,q) * R(q)
       nsupp = 0
       do r = 1, nr
          do ee = p_ptr(ip0+r-1), p_ptr(ip0+r)-1
             q = p_col(ee)
             A(r,q) = A(r,q) + p_val(ee)
             if (.not. mark(q)) then
                mark(q) = .true.; nsupp = nsupp + 1; supp(nsupp) = q
             end if
          end do
       end do

       do ord = 0, order
          if (ord > 0) then
             ! A <- (W R) A over the active support
             nsupp2 = 0
             do tt = 1, nsupp
                q = supp(tt)
                anynz = .false.
                do r = 1, nr
                   m2 = real(A(r,q), dp)**2 + aimag(A(r,q))**2
                   if (m2 < thr2) then
                      a_loc(r) = czero
                   else
                      a_loc(r) = A(r,q); anynz = .true.
                   end if
                end do
                if (.not. anynz) cycle
                do ee = s_ptr(q), s_ptr(q+1)-1
                   qb = s_col(ee)
                   do r = 1, nr
                      Nx(r,qb) = Nx(r,qb) + a_loc(r) * s_wres(ee)
                   end do
                   if (.not. mark2(qb)) then
                      mark2(qb) = .true.; nsupp2 = nsupp2 + 1; supp2(nsupp2) = qb
                   end if
                end do
             end do
             ! clear the old amplitude, move Nx -> A
             do tt = 1, nsupp
                q = supp(tt)
                do r = 1, nr
                   A(r,q) = czero
                end do
                mark(q) = .false.
             end do
             amax2 = 0.0_dp
             do tt = 1, nsupp2
                q = supp2(tt)
                do r = 1, nr
                   A(r,q) = Nx(r,q); Nx(r,q) = czero
                   amax2 = max(amax2, real(A(r,q), dp)**2 + aimag(A(r,q))**2)
                end do
                mark(q) = .true.; mark2(q) = .false.
                supp(tt) = q
             end do
             nsupp = nsupp2
             if (amax2 < thr2) exit
          end if

          ! Close the chain: Sigma(ip,gp) += A(ip,q) * conj(g(gp,q)), gp >= ip0 only
          do tt = 1, nsupp
             q = supp(tt)
             anynz = .false.
             do r = 1, nr
                m2 = real(A(r,q), dp)**2 + aimag(A(r,q))**2
                if (m2 < thr2) then
                   a_loc(r) = czero
                else
                   a_loc(r) = A(r,q); anynz = .true.
                end if
             end do
             if (.not. anynz) cycle
             lo = c_ptr(q); hi = c_ptr(q+1)
             do while (lo < hi)                 ! first entry with c_row >= ip0
                mid = (lo + hi) / 2
                if (c_row(mid) < ip0) then
                   lo = mid + 1
                else
                   hi = mid
                end if
             end do
             do ee = lo, c_ptr(q+1)-1
                gp = c_row(ee)
                do r = 1, nr
                   Sigma(ip0+r-1, gp) = Sigma(ip0+r-1, gp) + a_loc(r) * c_val(ee)
                end do
             end do
          end do
       end do

       ! reset the scratch for the next block
       do tt = 1, nsupp
          q = supp(tt)
          do r = 1, nr
             A(r,q) = czero
          end do
          mark(q) = .false.
       end do
    end do
    !$omp end do
    deallocate(A, Nx, a_loc, supp, supp2, mark, mark2)
    !$omp end parallel
    call icgn_restore_nesting(saved_levels)

    ! Lower triangle from Hermiticity: Sigma(i,gp) = conj(Sigma(gp,i)), i > gp
    nth = min(icgn_thread_budget(), max(1, nred))
    !$omp parallel do default(shared) private(i, gp) num_threads(nth) schedule(static) if(nth > 1)
    do gp = 1, nred
       do i = gp+1, nred
          Sigma(i, gp) = conjg(Sigma(gp, i))
       end do
    end do
    !$omp end parallel do
  end subroutine icgn_sigma_blocked

  ! Power iteration for ||T||_2, T = (E0-D)^-1 W, on the shared S graph.
  ! Same method as icgn_power_iteration but both W and W^dagger products are
  ! row-parallel matvecs on the CSR-by-source arrays (W is stored in both
  ! directions with conjugate values, so W^dagger(r,c) = conj(w(r,c)) too).
  subroutine icgn_power_iteration_csr(evals_q, nstates, s_ptr, s_col, s_val, E0, maxiter, tol, &
       sigma_out, converged_out)
    real(dp),    intent(in)  :: evals_q(:)
    integer,     intent(in)  :: nstates, s_ptr(:), s_col(:), maxiter
    complex(dp), intent(in)  :: s_val(:)
    real(dp),    intent(in)  :: E0, tol
    real(dp),    intent(out) :: sigma_out
    logical,     intent(out) :: converged_out
    integer  :: k, j, ee, nth
    real(dp) :: nrm, sigma_est, diff2, rr, ri
    complex(dp) :: acc
    complex(dp), allocatable :: v(:), u(:), w(:), v_next(:)
    real(dp), allocatable :: res(:)

    sigma_out = 0.0_dp; converged_out = .true.
    allocate(v(nstates), u(nstates), w(nstates), v_next(nstates), res(nstates))
    do j = 1, nstates
       res(j) = 1.0_dp / (E0 - evals_q(j))
    end do
    ! random start vector on the states that have at least one Q-Q edge
    v = (0.0_dp, 0.0_dp)
    do j = 1, nstates
       if (s_ptr(j+1) > s_ptr(j)) then
          call random_number(rr); call random_number(ri)
          v(j) = cmplx(rr - 0.5_dp, ri - 0.5_dp, kind=dp)
       end if
    end do
    nrm = sqrt(sum(real(v,dp)**2 + aimag(v)**2))
    if (nrm < 1.0e-14_dp) then
       deallocate(v, u, w, v_next, res); return
    end if
    v = v / nrm
    converged_out = .false.
    sigma_est = 0.0_dp
    nth = icgn_thread_budget()

    do k = 1, maxiter
       ! u = (E0-D)^-1 W v      (row j: sum_i conj(w_ji) v_i = sum_i w_ij v_i)
       !$omp parallel do default(shared) private(j, ee, acc) num_threads(nth) &
       !$omp& schedule(dynamic,256) if(nth > 1)
       do j = 1, nstates
          acc = (0.0_dp, 0.0_dp)
          do ee = s_ptr(j), s_ptr(j+1)-1
             acc = acc + conjg(s_val(ee)) * v(s_col(ee))
          end do
          u(j) = acc * res(j)
       end do
       !$omp end parallel do
       sigma_est = sqrt(sum(real(u,dp)**2 + aimag(u)**2))

       ! w = (E0-D)^-1 u ; v_next = W^dagger w  (row i: sum_j conj(w_ij) w_j)
       do j = 1, nstates
          w(j) = u(j) * res(j)
       end do
       !$omp parallel do default(shared) private(j, ee, acc) num_threads(nth) &
       !$omp& schedule(dynamic,256) if(nth > 1)
       do j = 1, nstates
          acc = (0.0_dp, 0.0_dp)
          do ee = s_ptr(j), s_ptr(j+1)-1
             acc = acc + conjg(s_val(ee)) * w(s_col(ee))
          end do
          v_next(j) = acc
       end do
       !$omp end parallel do

       nrm = sqrt(sum(real(v_next,dp)**2 + aimag(v_next)**2))
       if (nrm < 1.0e-14_dp) then
          sigma_out = 0.0_dp; converged_out = .true.
          deallocate(v, u, w, v_next, res); return
       end if
       v_next = v_next / nrm
       diff2 = sum(real(v_next - v, dp)**2 + aimag(v_next - v)**2)
       v = v_next
       if (sqrt(diff2) < tol) then
          converged_out = .true.
          exit
       end if
    end do
    sigma_out = sigma_est
    deallocate(v, u, w, v_next, res)
  end subroutine icgn_power_iteration_csr

  ! Helper: grow integer array
  subroutine grow_int_array(arr, new_size)
    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    integer, allocatable :: tmp(:)
    allocate(tmp(new_size))
    tmp(1:size(arr)) = arr
    call move_alloc(tmp, arr)
  end subroutine grow_int_array

  ! Helper: grow complex array
  subroutine grow_cx_array(arr, new_size)
    complex(dp), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    complex(dp), allocatable :: tmp(:)
    allocate(tmp(new_size))
    tmp(1:size(arr)) = arr
    call move_alloc(tmp, arr)
  end subroutine grow_cx_array

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

end module coarse_grain
