! Test driver for 4 solve modes:
!   (1) Standard full diagonalization
!   (2) Original coarse-graining (Liu et al. 2022)
!   (3) Improved coarse-graining (core + buffer + level-1 acquaintance)
!   (4) Improved CG + Neumann self-energy correction (ICGN)
!
! TIMING: The sparse Hamiltonian is built once and reused.  For each mode
! the clock starts when that mode receives the already-built H and stops
! after eigenvectors have been lifted back to the original orbital basis.
! For modes 2/3/4 this therefore includes: reduced-H preparation
! (cg_prepare / icg_prepare / icgn_prepare) + eigensolver + lift.
! For mode 1 it includes only: eigensolver (no lift needed).
!
! AAD: eigenvalues from modes 2/3/4 are compared to the reference from
! mode 1.  Both sets are sorted ascending before comparison.
!   n_up:   the n_up smallest *positive* eigenvalues (closest to 0 from above)
!   n_down: the n_down largest *negative* eigenvalues (closest to 0 from below)
!
! Config file format (one value per line):
!   1:  structure file (.upg)
!   2:  relativistic (.true./.false.)
!   3:  Harrison scaling (.true./.false.)
!   4:  c-axis (3 floats)
!   5:  solver (LK / JD / LO)
!   6:  nVB  (standard mode)
!   7:  nCB  (standard mode)
!   8:  lambda_vb  (eV)
!   9:  lambda_cb  (eV)
!   10: n_blocks   (CG & ICG & ICGN)
!   11: cg_emin    (eV)
!   12: cg_emax    (eV)
!   13: imbalance  (METIS)
!   14: icg_core_emin  (eV)
!   15: icg_core_emax  (eV)
!   16: icg_e_buffer   (eV)
!   17: icg_epsilon    (threshold factor)
!   18: icgn_selfenergy_order (0,1,2,...)
!   19: icgn_E0     (eV, 0.0 = auto = core window midpoint)
!   20: n_up   (n smallest positive eigenvalues for AAD, 0 to skip)
!   21: n_down (n largest  negative eigenvalues for AAD, 0 to skip)
!   22: icgn_check_convergence (.true./.false.)
!   23: icgn_pi_tol  (power-iteration tolerance, e.g. 1e-6)
program test_supercell

  USE precision
  USE globals,             only : MST, LST
  USE mpi_globals,         only : upt_mpi_init, upt_mpi_end
  USE upt_param,           only : OUPT
  USE struct_building,     only : init_structure, make_basis, init_basis, subs_dg_ions
  USE neighbours,          only : refine_neighbours_map, check_input_nn_list, &
                                  check_nn_map, write_neighbours_map
  USE input_data,          only : read_data
  USE type_defs,           only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   sort_states, set_max_order
  USE alloys,              only : init_mat_ion
  USE uptight,             only : UPT_configure_coarse_graining,  &
                                  UPT_get_coarse_graining_info,    &
                                  UPT_configure_improved_cg,       &
                                  UPT_get_improved_cg_info,        &
                                  UPT_configure_icgn,              &
                                  UPT_get_icgn_info,               &
                                  upt_hamiltonian
  USE lapack_driver,       only : lapack, lapack_icg, lapack_icgn
  USE JD_driver,           only : jd
  USE lanczos_driver,      only : lanczos
  USE sparse_matrix,       only : destroy_matrix
  USE clock,               only : set_clock, get_sclock

  IMPLICIT NONE

  TYPE(OUPT), TARGET  :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER        :: i, err, n_ham, num_ev
  CHARACTER(LST) :: config_file
  CHARACTER(MST) :: solver_choice
  INTEGER        :: n_blocks, nVB, nCB
  REAL(dp)       :: cg_emin, cg_emax, imbalance
  REAL(dp)       :: icg_core_emin, icg_core_emax, icg_e_buffer, icg_epsilon
  INTEGER        :: icgn_selfenergy_order
  REAL(dp)       :: icgn_E0
  LOGICAL        :: icgn_check_conv
  REAL(dp)       :: icgn_pi_tol
  INTEGER        :: icgn_pi_maxiter
  INTEGER        :: n_up, n_down
  REAL(sp)       :: solve_time
  LOGICAL        :: cg_ready, icg_ready_flag, icgn_ready_flag
  INTEGER        :: orig_dim, red_dim, nb_out
  REAL(dp)       :: cut_frac
  ! Reference eigenvalues from mode 1 (sorted ascending, allocated after mode 1)
  REAL(dp), ALLOCATABLE :: ref_evals(:)
  ! ICGN convergence results (from icgn_get_info)
  REAL(dp)   :: sigma_T2
  LOGICAL    :: pi_converged

  call upt_mpi_init(0)
  pupt => upt

  ! ---- default paths --------------------------------------------------------
  upt%database_path = './'
  upt%work_path     = './'
  upt%out_path      = './'
  upt%gen_out       = 'out.gen'
  upt%state_file    = 'states.data'
  upt%sparse_format = 'U'
  upt%verbose       = 10

  ! ---- read config ----------------------------------------------------------
  config_file = 'config'
  open(10, file=trim(config_file), status='old', action='read', iostat=err)
  if (err /= 0) then
     write(*,*) 'ERROR: cannot open config file: ', trim(config_file); stop 1
  end if
  read(10,*) upt%gen_filename
  read(10,*) upt%relat
  read(10,*) upt%scaling
  read(10,*) upt%c_axis(:)
  read(10,*) solver_choice
  read(10,*) nVB
  read(10,*) nCB
  read(10,*) upt%lambda_vb
  read(10,*) upt%lambda_cb
  read(10,*) n_blocks
  read(10,*) cg_emin
  read(10,*) cg_emax
  read(10,*) imbalance
  read(10,*) icg_core_emin
  read(10,*) icg_core_emax
  read(10,*) icg_e_buffer
  read(10,*) icg_epsilon
  read(10,*) icgn_selfenergy_order
  read(10,*) icgn_E0
  read(10,*) n_up
  read(10,*) n_down
  read(10,*) icgn_check_conv
  read(10,*) icgn_pi_tol
  close(10)

  icgn_pi_maxiter = 2000   ! default for power iteration

  write(*,'(a)') '========================================'
  write(*,'(a,a)')     ' Structure:    ', trim(upt%gen_filename)
  write(*,'(a,l1)')    ' Relativistic: ', upt%relat
  write(*,'(a,l1)')    ' Scaling:      ', upt%scaling
  write(*,'(a,a)')     ' Solver:       ', trim(solver_choice)
  write(*,'(a,i0)')    ' n_blocks:     ', n_blocks
  write(*,'(a,2f8.3)') '  CG window:   ', cg_emin, cg_emax
  write(*,'(a,2f8.3)') '  ICG core:    ', icg_core_emin, icg_core_emax
  write(*,'(a,f8.3)')  '  ICG buffer:  ', icg_e_buffer
  write(*,'(a,es10.2)')'  ICG epsilon: ', icg_epsilon
  write(*,'(a,i0)')    '  ICGN order:  ', icgn_selfenergy_order
  write(*,'(a,f8.3)')  '  ICGN E0:     ', icgn_E0
  write(*,'(a,l1)')    '  ICGN check convergence: ', icgn_check_conv
  if (icgn_check_conv) then
     write(*,'(a,es10.2,a,i0)') '  ICGN PI tol: ', icgn_pi_tol, '  maxiter: ', icgn_pi_maxiter
  end if
  write(*,'(a,i0,a,i0)') '  AAD bands: +', n_up, ' / -', n_down
  write(*,'(a)') '========================================'

  ! ---- misc params ----------------------------------------------------------
  upt%structure%gen_filename = upt%gen_filename
  upt%d_onsite_shift_flag = .true.
  upt%potential_flag      = .false.
  upt%syst_rotated        = .false.
  upt%ioutput_flag        = .false.
  upt%optmat              = .false.
  upt%poldir  = 3
  upt%d_H     = 0.1d0
  upt%E_H     = -200.0d0
  upt%estimate_factor = 1.0
  upt%check_bondmap   = .false.
  upt%n_spin = merge(2, 1, upt%relat)

  ! ---- solver bookkeeping ---------------------------------------------------
  upt%num_vb   = nVB;  upt%num_cb   = nCB
  upt%start_vb = 1;    upt%start_cb = 1
  upt%min_iter = 2;    upt%long_iter = 30;  upt%max_iter = 100000
  upt%fast_tol = 1.0d-1; upt%long_tol = 1.0d-10; upt%ort_tol = 1.0d-5
  upt%solver_flag = 0;  upt%dynamic = .true.
  upt%seed_flag   = .false.; upt%bitoff = 0.1_dp
  upt%k_point = (/ 0.0d0, 0.0d0, 0.0d0 /)

  ! ---- build structure (done once) -----------------------------------------
  write(*,'(a)') ' Building structure...'
  call set_machine_acc
  call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat, &
                      upt%interfaces, upt%nr_int)
  call check_input_nn_list(upt%structure)
  do i = 1, upt%nr_mat
     call read_data(upt%materials(i), upt%work_path, upt%database_path)
  end do
  call make_basis(upt%verbose, upt%structure, upt%basis)
  call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, &
                                upt%ref_couplings, upt%n_ref_cpl)
  do i = 1, upt%nr_mat
     call sort_states(upt%materials(i), upt%ref_states, upt%ref_couplings)
  end do
  call set_max_order(upt%materials)
  do i = 1, upt%nr_mat; call init_mat_ion(upt%materials(i)); end do
  call write_materials(upt%nr_mat, upt%materials)
  call init_basis(upt%basis, upt%materials)
  call write_basis(upt%basis, 0)
  call refine_neighbours_map(upt%structure, upt%basis, upt%materials, upt%nn_map)
  call write_neighbours_map(upt%nn_map)
  call check_nn_map(upt%nn_map)
  call subs_dg_ions(upt%basis, upt%materials, upt%nn_map)
  call init_n_st(upt%basis, upt%materials)

  ! ---- build the sparse Hamiltonian once (no CG); record n_ham -------------
  upt%cg_enabled   = .false.
  upt%icg_enabled  = .false.
  upt%icgn_enabled = .false.
  upt%verbose = 0
  call upt_hamiltonian(pupt)
  n_ham = upt%ham%nrow
  write(*,'(a,i0)') ' Full Hamiltonian dimension: ', n_ham
  ! Keep upt%ham alive — all modes reuse it.

  ! ==========================================================================
  ! MODE 1: Standard full diagonalization
  !   Clock starts: eigensolver receives H (no prepare step)
  !   Clock stops:  eigenvectors are in original orbital basis (no lift)
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 1: Standard full diagonalization'
  write(*,'(a)') '========================================'

  upt%cg_enabled   = .false.
  upt%icg_enabled  = .false.
  upt%icgn_enabled = .false.
  num_ev = nVB + nCB
  allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
           upt%particles(num_ev), stat=err)
  upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0

  call set_clock()          ! --- start timing ---
  select case (trim(solver_choice))
  case ('LK'); call lapack(upt)
  case ('JD'); call jd(upt)
  case ('LO'); call lanczos(upt)
  case default; write(*,*) 'Unknown solver: ', trim(solver_choice); stop 1
  end select
  solve_time = get_sclock() ! --- stop timing ---

  write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
  write(*,'(a,f10.3)')  ' Total time:   ', solve_time
  write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)

  ! Save sorted reference eigenvalues for AAD
  allocate(ref_evals(size(upt%eigen_values)))
  ref_evals = upt%eigen_values
  call sort_ascending(ref_evals)

  call write_eigenvalues('eigenvalues_standard.dat', upt%eigen_values, solve_time, &
       'STANDARD', n_ham, n_ham, 0.0_dp, ref_evals, 0, 0)
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)

  ! ==========================================================================
  ! MODE 2: Original coarse-graining (Liu et al.)
  !   Clock starts: cg_prepare receives H  (inside upt_hamiltonian)
  !   Clock stops:  eigenvectors lifted to original orbital basis
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 2: Original coarse-graining (Liu et al.)'
  write(*,'(a)') '========================================'

  ! Re-build with only CG enabled so prepare is included in the timed region
  call destroy_matrix(upt%ham)
  call UPT_configure_coarse_graining(upt, .true.,  n_blocks, cg_emin, cg_emax, imbalance)
  call UPT_configure_improved_cg    (upt, .false., n_blocks, icg_core_emin, icg_core_emax, &
       icg_e_buffer, icg_epsilon, imbalance)
  call UPT_configure_icgn           (upt, .false., n_blocks, icg_core_emin, icg_core_emax, &
       icg_e_buffer, icg_epsilon, icgn_selfenergy_order, icgn_E0, imbalance, &
       icgn_check_conv, icgn_pi_maxiter, icgn_pi_tol)

  call set_clock()          ! --- start timing (includes cg_prepare) ---
  call upt_hamiltonian(pupt)

  call UPT_get_coarse_graining_info(upt, cg_ready, orig_dim, red_dim, nb_out, cut_frac)
  if (.not. cg_ready) then
     write(*,*) ' WARNING: CG not ready, skipping mode 2'
     call destroy_matrix(upt%ham); goto 300
  end if
  write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
       '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(orig_dim,dp)), '% reduction)'
  num_ev = red_dim
  allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
           upt%particles(num_ev), stat=err)
  upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0
  upt%cg_enabled = .true.
  select case (trim(solver_choice))
  case ('LK'); call lapack(upt)   ! lapack calls cg_lift internally
  case ('JD'); call jd(upt)
  case ('LO'); call lanczos(upt)
  end select
  solve_time = get_sclock() ! --- stop timing (after lift) ---

  write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
  write(*,'(a,f10.3)')  ' Total time:   ', solve_time
  write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
  call write_eigenvalues('eigenvalues_cg.dat', upt%eigen_values, solve_time, &
       'CG', orig_dim, red_dim, cut_frac, ref_evals, n_up, n_down)
  call destroy_matrix(upt%ham)
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
300 continue

  ! ==========================================================================
  ! MODE 3: Improved coarse-graining
  !   Clock starts: icg_prepare receives H
  !   Clock stops:  eigenvectors lifted to original orbital basis
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 3: Improved coarse-graining'
  write(*,'(a)') '========================================'

  call UPT_configure_coarse_graining(upt, .false., n_blocks, cg_emin, cg_emax, imbalance)
  call UPT_configure_improved_cg    (upt, .true.,  n_blocks, icg_core_emin, icg_core_emax, &
       icg_e_buffer, icg_epsilon, imbalance)
  call UPT_configure_icgn           (upt, .false., n_blocks, icg_core_emin, icg_core_emax, &
       icg_e_buffer, icg_epsilon, icgn_selfenergy_order, icgn_E0, imbalance, &
       icgn_check_conv, icgn_pi_maxiter, icgn_pi_tol)

  call set_clock()          ! --- start timing (includes icg_prepare) ---
  call upt_hamiltonian(pupt)

  call UPT_get_improved_cg_info(upt, icg_ready_flag, orig_dim, red_dim, nb_out, cut_frac)
  if (.not. icg_ready_flag) then
     write(*,*) ' WARNING: ICG not ready, skipping mode 3'
     call destroy_matrix(upt%ham); goto 400
  end if
  write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
       '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(orig_dim,dp)), '% reduction)'
  num_ev = red_dim
  allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
           upt%particles(num_ev), stat=err)
  upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0
  upt%icg_enabled = .true.
  select case (trim(solver_choice))
  case ('LK'); call lapack_icg(upt)   ! includes icg_lift
  case default
     write(*,*) ' ICG currently supports LK solver only'
     call destroy_matrix(upt%ham)
     deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles); goto 400
  end select
  solve_time = get_sclock() ! --- stop timing (after lift) ---

  write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
  write(*,'(a,f10.3)')  ' Total time:   ', solve_time
  write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
  call write_eigenvalues('eigenvalues_icg.dat', upt%eigen_values, solve_time, &
       'ICG', orig_dim, red_dim, cut_frac, ref_evals, n_up, n_down)
  call destroy_matrix(upt%ham)
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
400 continue

  ! ==========================================================================
  ! MODE 4: Improved CG + Neumann self-energy correction (ICGN)
  !   Clock starts: icgn_prepare receives H (includes self-energy build)
  !   Clock stops:  eigenvectors lifted to original orbital basis
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 4: ICGN (Neumann self-energy)'
  write(*,'(a)') '========================================'

  call UPT_configure_coarse_graining(upt, .false., n_blocks, cg_emin, cg_emax, imbalance)
  call UPT_configure_improved_cg    (upt, .false., n_blocks, icg_core_emin, icg_core_emax, &
       icg_e_buffer, icg_epsilon, imbalance)
  call UPT_configure_icgn           (upt, .true.,  n_blocks, icg_core_emin, icg_core_emax, &
       icg_e_buffer, icg_epsilon, icgn_selfenergy_order, icgn_E0, imbalance, &
       icgn_check_conv, icgn_pi_maxiter, icgn_pi_tol)

  call set_clock()          ! --- start timing (includes icgn_prepare) ---
  call upt_hamiltonian(pupt)

  call UPT_get_icgn_info(upt, icgn_ready_flag, orig_dim, red_dim, nb_out, cut_frac, &
       sigma_T2, pi_converged)
  if (.not. icgn_ready_flag) then
     write(*,*) ' WARNING: ICGN not ready, skipping mode 4'
     call destroy_matrix(upt%ham); goto 500
  end if
  write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
       '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(orig_dim,dp)), '% reduction)'
  if (icgn_check_conv) then
     write(*,'(a,l1)') '  Power iteration converged: ', pi_converged
     if (sigma_T2 < 0.0_dp) then
        write(*,'(a)') '  Spectral norm sigma_T2: N/A (no Q-Q edges found)'
     else
        if (sigma_T2 >= 1.0_dp) then
           write(*,'(a,es12.4,a)') '  Spectral norm ||T||_2 = ', sigma_T2, &
                ' >= 1 (WARNING: full Neumann series not guaranteed to converge)'
        else
           write(*,'(a,es12.4,a)') '  Spectral norm ||T||_2 = ', sigma_T2, &
                ' < 1  (Neumann series on solid footing)'
        end if
     end if
  end if
  num_ev = red_dim
  allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
           upt%particles(num_ev), stat=err)
  upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0
  upt%icgn_enabled = .true.
  select case (trim(solver_choice))
  case ('LK'); call lapack_icgn(upt)  ! includes icgn_lift
  case default
     write(*,*) ' ICGN currently supports LK solver only'
     call destroy_matrix(upt%ham)
     deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles); goto 500
  end select
  solve_time = get_sclock() ! --- stop timing (after lift) ---

  write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
  write(*,'(a,f10.3)')  ' Total time:   ', solve_time
  write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
  call write_eigenvalues('eigenvalues_icgn.dat', upt%eigen_values, solve_time, &
       'ICGN', orig_dim, red_dim, cut_frac, ref_evals, n_up, n_down, &
       sigma_T2=sigma_T2, pi_conv=pi_converged)
  call destroy_matrix(upt%ham)
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
500 continue

  ! ---- cleanup --------------------------------------------------------------
  if (allocated(ref_evals)) deallocate(ref_evals)
  call upt_mpi_end

contains

  ! ---------------------------------------------------------------------------
  ! In-place ascending sort (insertion sort)
  subroutine sort_ascending(a)
    real(dp), intent(inout) :: a(:)
    integer  :: i, j
    real(dp) :: tmp
    do i = 2, size(a)
       tmp = a(i); j = i - 1
       do while (j >= 1 .and. a(j) > tmp)
          a(j+1) = a(j); j = j - 1
       end do
       a(j+1) = tmp
    end do
  end subroutine sort_ascending

  ! ---------------------------------------------------------------------------
  ! Compute AAD for:
  !   n_up   smallest positive eigenvalues (closest to 0 from above)
  !   n_down largest  negative eigenvalues (closest to 0 from below)
  ! Both `evals` and `ref` must be sorted ascending on entry.
  ! Returns -1.0 when there are not enough matching bands.
  subroutine compute_aad(evals, ref, n_up, n_down, aad_up, aad_down)
    real(dp), intent(in)  :: evals(:), ref(:)
    integer,  intent(in)  :: n_up, n_down
    real(dp), intent(out) :: aad_up, aad_down

    integer  :: n_ev, n_ref, i
    integer  :: first_pos_ev, first_pos_ref   ! first index ≥ 0 in each array
    integer  :: last_neg_ev,  last_neg_ref    ! last  index <  0 in each array

    n_ev  = size(evals)
    n_ref = size(ref)

    ! ---- locate sign boundaries ----
    first_pos_ev  = n_ev  + 1
    first_pos_ref = n_ref + 1
    do i = 1, n_ev;  if (evals(i) >= 0.0_dp) then; first_pos_ev  = i; exit; end if; end do
    do i = 1, n_ref; if (ref(i)   >= 0.0_dp) then; first_pos_ref = i; exit; end if; end do

    last_neg_ev  = first_pos_ev  - 1
    last_neg_ref = first_pos_ref - 1

    ! ---- AAD for n_up smallest positives ----
    aad_up = -1.0_dp
    if (n_up > 0) then
       ! Need at least n_up positive values in both arrays
       if ((n_ev - first_pos_ev + 1 >= n_up) .and. &
           (n_ref - first_pos_ref + 1 >= n_up)) then
          aad_up = 0.0_dp
          do i = 0, n_up - 1
             aad_up = aad_up + abs(evals(first_pos_ev + i) - ref(first_pos_ref + i))
          end do
          aad_up = aad_up / real(n_up, dp)
       end if
    end if

    ! ---- AAD for n_down largest negatives ----
    aad_down = -1.0_dp
    if (n_down > 0) then
       ! Need at least n_down negative values in both arrays
       if ((last_neg_ev >= n_down) .and. (last_neg_ref >= n_down)) then
          aad_down = 0.0_dp
          do i = 0, n_down - 1
             aad_down = aad_down + abs(evals(last_neg_ev - i) - ref(last_neg_ref - i))
          end do
          aad_down = aad_down / real(n_down, dp)
       end if
    end if
  end subroutine compute_aad

  ! ---------------------------------------------------------------------------
  subroutine write_eigenvalues(fname, evals, t, mode, ndim, nred, cut, ref, nu, nd, &
       sigma_T2, pi_conv)
    character(*), intent(in) :: fname, mode
    real(dp),     intent(in) :: evals(:), cut, ref(:)
    real(sp),     intent(in) :: t
    integer,      intent(in) :: ndim, nred, nu, nd
    real(dp),     intent(in), optional :: sigma_T2
    logical,      intent(in), optional :: pi_conv

    integer  :: n, fu, ii
    real(dp), allocatable :: se(:)
    real(dp) :: aad_up, aad_down

    n = size(evals)
    allocate(se(n))
    se = evals
    call sort_ascending(se)

    aad_up = -1.0_dp; aad_down = -1.0_dp
    if ((nu > 0 .or. nd > 0) .and. size(ref) > 0) then
       ! ref is already sorted ascending (saved that way from mode 1)
       call compute_aad(se, ref, nu, nd, aad_up, aad_down)
    end if

    open(newunit=fu, file=trim(fname), status='replace', action='write')
    write(fu,'(a,a)')     '# Mode: ', trim(mode)
    write(fu,'(a,i0)')    '# Full dimension:    ', ndim
    write(fu,'(a,i0)')    '# Reduced dimension: ', nred
    write(fu,'(a,f8.2,a)') '# Rank reduction: ', &
         100.0_dp*(1.0_dp - real(nred,dp)/real(ndim,dp)), ' %'
    write(fu,'(a,f10.4)') '# Cut fraction:      ', cut
    write(fu,'(a,f12.6)') '# Total time (s):    ', t
    write(fu,'(a,i0)')    '# Total bands:       ', n
    ! ICGN convergence metadata
    if (present(sigma_T2)) then
       if (sigma_T2 < 0.0_dp) then
          write(fu,'(a)') '# Neumann ||T||_2: N/A (no Q-Q edges)'
       else
          if (sigma_T2 >= 1.0_dp) then
             write(fu,'(a,es14.6,a)') '# Neumann ||T||_2 = ', sigma_T2, &
                  ' (WARNING: >= 1, convergence not guaranteed)'
          else
             write(fu,'(a,es14.6,a)') '# Neumann ||T||_2 = ', sigma_T2, &
                  ' (< 1, series on solid footing)'
          end if
       end if
       if (present(pi_conv)) write(fu,'(a,l1)') '# Power iter converged: ', pi_conv
    end if
    if (nu > 0) then
       if (aad_up >= 0.0_dp) then
          write(fu,'(a,i0,a,es14.6,a)') &
               '# AAD ', nu, ' smallest positive evals (eV): ', aad_up, ''
       else
          write(fu,'(a,i0,a)') '# AAD ', nu, ' smallest positive evals (eV): N/A'
       end if
    end if
    if (nd > 0) then
       if (aad_down >= 0.0_dp) then
          write(fu,'(a,i0,a,es14.6,a)') &
               '# AAD ', nd, ' largest  negative evals (eV): ', aad_down, ''
       else
          write(fu,'(a,i0,a)') '# AAD ', nd, ' largest  negative evals (eV): N/A'
       end if
    end if
    write(fu,'(a)') '#'
    write(fu,'(a)') '# Index    Energy(eV)'
    do ii = 1, n
       write(fu,'(i6,2x,f16.8)') ii, se(ii)
    end do
    close(fu)

    write(*,'(a,a)') ' Output: ', trim(fname)
    if (nu > 0 .and. aad_up   >= 0.0_dp) &
         write(*,'(a,i0,a,es12.4)') '  AAD ', nu, ' smallest pos: ', aad_up
    if (nd > 0 .and. aad_down >= 0.0_dp) &
         write(*,'(a,i0,a,es12.4)') '  AAD ', nd, ' largest  neg: ', aad_down
    deallocate(se)
  end subroutine write_eigenvalues

end program test_supercell
