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
  end type ICGNPair

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
    allocate(tmp(size(qa,1), size(qb,2)))
    tmp = matmul(p%h, qb)
    allocate(p%v(size(qa,2), size(qb,2)))
    p%v = matmul(conjg(transpose(qa)), tmp)
    deallocate(tmp)
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
    complex(dp), allocatable :: g_full(:,:)

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
    allocate(row_of(n)); row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
      call icgn_log_block(upt, 'icgn', i, upt%icgn_blocks(i)%nrow, 'processing')
      call icgn_diagonalize_block(upt, i, row_of, ierr)
       if (ierr /= 0) return
       row_of(upt%icgn_blocks(i)%rows) = 0
    end do

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
       allocate(pairs(max(1, upt%ham%nnz))); npair = 0
       allocate(pair_map(upt%icgn_num_blocks, upt%icgn_num_blocks)); pair_map = 0
       do r = 1, upt%ham%nrow
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k)
             ia = label(atom_of(r)); ib = label(atom_of(c))
             if (ia == ib) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             if (upt%icgn_blocks(ia)%nrow == 0 .or. upt%icgn_blocks(ib)%nrow == 0) cycle
             slot = pair_slot_icgn(pairs, npair, min(ia,ib), max(ia,ib), upt, pair_map)
             if (slot == 0) then; ierr = 11; return; end if
             if (ia < ib) then
                pairs(slot)%h(row_of(r),row_of(c)) = pairs(slot)%h(row_of(r),row_of(c)) + upt%ham%M(k)
             else
                pairs(slot)%h(row_of(c),row_of(r)) = pairs(slot)%h(row_of(c),row_of(r)) + conjg(upt%ham%M(k))
             end if
          end do
       end do
       deallocate(pair_map)
    else
       allocate(pairs(max(1, upt%icgn_num_blocks*(upt%icgn_num_blocks-1)/2))); npair = 0
    end if
    if (upt%icgn_add_core_acquaintances > 0) then
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          call project_pair(pairs(i), upt%icgn_blocks(ia)%S_full, upt%icgn_blocks(ib)%S_full)
       end do
    end if

    ! Acquaintance selection is optional. Level 1 is the currently implemented
    ! core-state acquaintance level; higher requested levels fall back to it.
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
    ! Keep physical pair matrices for the zero-acquaintance fast path.
    if (upt%icgn_add_core_acquaintances > 0) call destroy_pairs(pairs)

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

    ! ---- Build reduced Hamiltonian (keeping pairs for Sigma) ----------------
    call build_icgn_reduced_hamiltonian(upt, atom_of, label, row_of, pairs, npair, ierr)
    if (ierr /= 0) return

    ! ====================================================================
    ! Self-energy correction via Neumann series (skipped when neumann_order < 0)
    ! ====================================================================
    if (upt%icgn_selfenergy_order >= 0) then
    ! pairs(s)%v is now (nret_a x nret_b): H_PP sliced couplings.
    ! We need H_PQ: from the *unsliced* g_full before retained_idx slicing.
    ! Strategy: re-project directly from S_full to extract P-Q and Q-Q.

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

    ! Build P-Q edge list: (global_red_row_p, flat_q_idx, coupling g(p,q))
    ! These come from pairs before slicing. Re-project with S_full.
    ! For each inter-block physical CSR entry (r,c) where r in P, c in Q:
    cnt_pq = 0
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          ! We count all P-Q edges across the full g matrix
          cnt_pq = cnt_pq + upt%icgn_blocks(ia)%nret * &
               (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
          cnt_pq = cnt_pq + upt%icgn_blocks(ib)%nret * &
               (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret)
       end do
    end do
    ! Upper bound — just pre-allocate generously using npair info
    ! Exact P-Q: for each pair, nret_a * (nrow_b - nret_b) + nret_b * (nrow_a - nret_a)
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
    allocate(pq_p(max(1,cnt_pq)), pq_q(max(1,cnt_pq)))
    allocate(pq_v(max(1,cnt_pq)))
    npq = 0

    ! Re-project: for each coupling pair (a,b), build g_full(nrow_a, nrow_b) fresh.
    ! Then extract P-Q entries.
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       allocate(g_full(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
       ! The physical block coupling was already collected in pairs(i)%h.
       ! Reuse it directly; do not rescan the full sparse Hamiltonian.
       block
          complex(dp), allocatable :: tmp_g(:,:)
          allocate(tmp_g(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
          tmp_g = matmul(pairs(i)%h, upt%icgn_blocks(ib)%S_full)
          g_full = matmul(conjg(transpose(upt%icgn_blocks(ia)%S_full)), tmp_g)
          deallocate(tmp_g)
       end block
       ! Extract P(ia)-Q(ib) entries
       do j = 1, upt%icgn_blocks(ia)%nret
          do k = 1, upt%icgn_blocks(ib)%nrow
             if (local_ret_idx(ib, k) /= 0) cycle ! k is in P
             if (abs(g_full(upt%icgn_blocks(ia)%retained_idx(j), k)) < 1.0e-14_dp) cycle
             if (.not. icgn_coupling_survives(g_full(upt%icgn_blocks(ia)%retained_idx(j), k), &
                 upt%icgn_blocks(ia)%evals_full(upt%icgn_blocks(ia)%retained_idx(j)), &
                 upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) cycle
             npq = npq + 1
             if (npq > size(pq_p)) then
                call grow_int_array(pq_p, 2*size(pq_p))
                call grow_int_array(pq_q, 2*size(pq_q))
                call grow_cx_array(pq_v, 2*size(pq_v))
             end if
             pq_p(npq) = roff_blk(ia) + j - 1
             pq_q(npq) = ret_offset(ib) + k
             pq_v(npq) = g_full(upt%icgn_blocks(ia)%retained_idx(j), k)
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
             npq = npq + 1
             if (npq > size(pq_p)) then
                call grow_int_array(pq_p, 2*size(pq_p))
                call grow_int_array(pq_q, 2*size(pq_q))
                call grow_cx_array(pq_v, 2*size(pq_v))
             end if
             pq_p(npq) = roff_blk(ib) + j - 1
             pq_q(npq) = ret_offset(ia) + k
             pq_v(npq) = conjg(g_full(k, upt%icgn_blocks(ib)%retained_idx(j)))
          end do
       end do
       deallocate(g_full)
    end do

    ! Build Q-Q edge list if order >= 1 or convergence check is requested
    if (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence) then
       ! Both directions: cnt_qq * 2 (upper bound)
       cnt_qq = 0
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          cnt_qq = cnt_qq + 2 * (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret) * &
               (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
       end do
       allocate(qq_i(max(1,cnt_qq)), qq_j(max(1,cnt_qq)))
       allocate(qq_v(max(1,cnt_qq)))
       nqq = 0
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          allocate(g_full(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
          ! Reuse the already-collected physical coupling for this pair.
          block
             complex(dp), allocatable :: tmp_g(:,:)
             allocate(tmp_g(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
             tmp_g = matmul(pairs(i)%h, upt%icgn_blocks(ib)%S_full)
             g_full = matmul(conjg(transpose(upt%icgn_blocks(ia)%S_full)), tmp_g)
             deallocate(tmp_g)
          end block
          ! Store BOTH directions (j->k and k->j with conj(v))
          do j = 1, upt%icgn_blocks(ia)%nrow
             if (local_ret_idx(ia, j) /= 0) cycle ! j in P, skip
             do k = 1, upt%icgn_blocks(ib)%nrow
                if (local_ret_idx(ib, k) /= 0) cycle ! k in P, skip
                if (abs(g_full(j, k)) < 1.0e-14_dp) cycle
                if (.not. icgn_coupling_survives(g_full(j, k), upt%icgn_blocks(ia)%evals_full(j), &
                    upt%icgn_blocks(ib)%evals_full(k), upt%icgn_epsilon)) cycle
                nqq = nqq + 1
                if (nqq > size(qq_i)) then
                   call grow_int_array(qq_i, 2*size(qq_i))
                   call grow_int_array(qq_j, 2*size(qq_j))
                   call grow_cx_array(qq_v, 2*size(qq_v))
                end if
                qq_i(nqq) = ret_offset(ia) + j
                qq_j(nqq) = ret_offset(ib) + k
                qq_v(nqq) = g_full(j, k)
                ! Reverse direction: k -> j, conj(v)
                nqq = nqq + 1
                if (nqq > size(qq_i)) then
                   call grow_int_array(qq_i, 2*size(qq_i))
                   call grow_int_array(qq_j, 2*size(qq_j))
                   call grow_cx_array(qq_v, 2*size(qq_v))
                end if
                qq_i(nqq) = ret_offset(ib) + k
                qq_j(nqq) = ret_offset(ia) + j
                qq_v(nqq) = conjg(g_full(j, k))
             end do
          end do
          deallocate(g_full)
       end do
    else
       nqq = 0
    end if

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
    if (upt%icgn_check_convergence .and. nqq > 0) then
       call icgn_power_iteration(evals_q, nstates_total, &
            qq_i, qq_j, qq_v, nqq, E0_used, &
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

    ! ---- Compute Sigma and add to H_dense ----------------------------------
    allocate(Sigma(nred, nred))
    Sigma = cmplx(0.0_dp, 0.0_dp, kind=dp)
    allocate(amp_vec(nstates_total), next_vec(nstates_total))

    E0_used = upt%icgn_E0
    if (E0_used == 0.0_dp) E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp

    ! For each unique P-state, propagate into Q-space and close back
    do ip = 1, nred
       amp_vec = cmplx(0.0_dp, 0.0_dp, kind=dp)
       ! Order-0 propagation: amp_vec(q) = g_{ip,q} * resolvent(q)
       do i = 1, npq
          if (pq_p(i) /= ip) cycle
          q_idx = pq_q(i)
          resolvent_val = 1.0_dp / (E0_used - evals_q(q_idx))
          amp_vec(q_idx) = amp_vec(q_idx) + pq_v(i) * resolvent_val
       end do

      ! Neumann truncation includes the zeroth term and every term through
      ! the requested order: Sigma = sum_{ord=0..N} H_PQ R (W R)^ord H_QP.
      do ord = 0, upt%icgn_selfenergy_order
          if (ord > 0) then
             next_vec = cmplx(0.0_dp, 0.0_dp, kind=dp)
             do i = 1, nqq
                q_idx  = qq_i(i)
                amp_val = amp_vec(q_idx)
                if (abs(amp_val) < 1.0e-14_dp) cycle
                qb_idx = qq_j(i)
                resolvent_val = 1.0_dp / (E0_used - evals_q(qb_idx))
                next_vec(qb_idx) = next_vec(qb_idx) + amp_val * qq_v(i) * resolvent_val
             end do
             amp_vec = next_vec
             if (maxval(abs(amp_vec)) < 1.0e-14_dp) exit
          end if

          ! Close chain: for each Q state with nonzero amp, sum over P' connected to Q
          do i = 1, npq
             q_idx = pq_q(i)
             amp_val = amp_vec(q_idx)
             if (abs(amp_val) < 1.0e-14_dp) cycle
             gp_col = pq_p(i)
             contrib = amp_val * conjg(pq_v(i))
             Sigma(ip, gp_col) = Sigma(ip, gp_col) + contrib
          end do
       end do
    end do

    deallocate(amp_vec, next_vec)
    deallocate(pq_p, pq_q, pq_v)
    if (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence) then
       if (allocated(qq_i)) deallocate(qq_i, qq_j, qq_v)
    end if
    deallocate(block_of_state, local_of_state, evals_q, ret_offset)
    deallocate(roff_blk, local_ret_idx)

    ! Add Sigma into H_dense and re-Hermitize
    H_dense = H_dense + Sigma
    do i = 1, nred
       do j = i+1, nred
          H_dense(i,j) = (H_dense(i,j) + conjg(H_dense(j,i))) / 2.0_dp
          H_dense(j,i) = conjg(H_dense(i,j))
       end do
    end do
    deallocate(Sigma)

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

    end if  ! icgn_selfenergy_order >= 0 (Neumann correction)

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

  ! Diagonalize block ib of icgn_blocks (same logic as icgn_diagonalize_block).
  subroutine icgn_diagonalize_block(upt, ib, local, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, local(:)
    integer, intent(out) :: ierr
    integer :: nn, k, r, c
    complex(dp), allocatable :: h(:,:)
    real(dp), allocatable :: w(:)
    ierr = 0; nn = upt%icgn_blocks(ib)%nrow
    if (nn == 0) then; upt%icgn_blocks(ib)%nret = 0; return; end if
    allocate(h(nn,nn), w(nn)); h = (0.0_dp, 0.0_dp)
    do r = 1, upt%ham%nrow
       if (local(r) == 0) cycle
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (local(c) == 0) cycle
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
    p%v = matmul(conjg(transpose(qa)), p%t)
    deallocate(p%t)
  end subroutine project_accumulated_pair_icgn

  ! Build reduced Hamiltonian for ICGN.
  subroutine build_icgn_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(ICGNPair), allocatable, intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(out) :: ierr
    integer :: i, j, k, r, c, a, b, ia, ib, nnz, pos, slot, nred
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
       do r = 1, upt%ham%nrow
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
             if (a == b) cycle
             if (upt%icgn_blocks(a)%nret == 0 .or. upt%icgn_blocks(b)%nret == 0) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             ia = min(a,b); ib = max(a,b)
             slot = pair_slot_icgn_fast(pairs, npair, ia, ib, upt, pair_map)
             if (slot == 0) then; ierr = 11; return; end if
             if (a < b) then
                pairs(slot)%t(local(r),:) = pairs(slot)%t(local(r),:) + &
                     upt%ham%M(k) * upt%icgn_blocks(b)%q(local(c),:)
             else
                pairs(slot)%t(local(c),:) = pairs(slot)%t(local(c),:) + &
                     conjg(upt%ham%M(k)) * upt%icgn_blocks(a)%q(local(r),:)
             end if
          end do
       end do
       deallocate(pair_map)
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
          call project_accumulated_pair_icgn(pairs(i), upt%icgn_blocks(a)%q)
          do j = 1, upt%icgn_blocks(a)%nret
             do k = 1, upt%icgn_blocks(b)%nret
                if (.not. icgn_coupling_survives(pairs(i)%v(j,k), &
                    upt%icgn_blocks(a)%eval(j), upt%icgn_blocks(b)%eval(k), &
                    upt%icgn_epsilon)) pairs(i)%v(j,k) = (0.0_dp, 0.0_dp)
             end do
          end do
       end do
    else if (upt%icgn_add_core_acquaintances == 0) then
       ! Zero-acquaintance path with Neumann enabled: keep the full physical
       ! pair matrices because the correction needs P-Q and Q-Q couplings.
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
    else
       if (allocated(pairs)) call destroy_pairs(pairs)
       allocate(pairs(max(1, upt%ham%nnz))); npair = 0
       allocate(pair_map(upt%icgn_num_blocks, upt%icgn_num_blocks)); pair_map = 0

       ! Collect physical H_ab first; project each pair once at matrix level.
       do r = 1, upt%ham%nrow
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
             if (a == b) cycle
             if (upt%icgn_blocks(a)%nrow == 0 .or. upt%icgn_blocks(b)%nrow == 0) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             ia = min(a,b); ib = max(a,b)
             slot = pair_slot_icgn(pairs, npair, ia, ib, upt, pair_map)
             if (slot == 0) then; ierr = 11; return; end if
             if (a < b) then
                pairs(slot)%h(local(r),local(c)) = pairs(slot)%h(local(r),local(c)) + upt%ham%M(k)
             else
                pairs(slot)%h(local(c),local(r)) = pairs(slot)%h(local(c),local(r)) + conjg(upt%ham%M(k))
             end if
          end do
       end do
       deallocate(pair_map)
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          call project_pair(pairs(i), upt%icgn_blocks(ia)%S_full, &
               upt%icgn_blocks(ib)%S_full, .true.)
       end do

       ! Slice to retained states using retained_idx.
       do i = 1, npair
          a = pairs(i)%a; b = pairs(i)%b
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
       end do
    end if

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
