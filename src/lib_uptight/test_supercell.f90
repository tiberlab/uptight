! Test driver for 3 solve modes:
!   (1) Standard full diagonalization
!   (2) Original coarse-graining (Liu et al. 2022)
!   (3) Improved coarse-graining (core + buffer + level-1 acquaintance)
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
!   10: n_blocks   (CG & ICG)
!   11: cg_emin    (eV)
!   12: cg_emax    (eV)
!   13: imbalance  (METIS)
!   14: icg_core_emin  (eV)
!   15: icg_core_emax  (eV)
!   16: icg_e_buffer   (eV)
!   17: icg_epsilon    (threshold factor)
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
                                  upt_hamiltonian
  USE lapack_driver,       only : lapack, lapack_icg_solve
  USE JD_driver,           only : jd
  USE lanczos_driver,      only : lanczos
  USE sparse_matrix,       only : destroy_matrix
  USE clock,               only : set_clock, get_sclock
  USE coarse_grain,        only : icg_lift

  IMPLICIT NONE

  TYPE(OUPT), TARGET  :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER        :: i, j, err, n_ham, num_ev
  CHARACTER(LST) :: config_file
  CHARACTER(MST) :: solver_choice
  INTEGER        :: n_blocks, nVB, nCB
  REAL(dp)       :: cg_emin, cg_emax, imbalance
  REAL(dp)       :: icg_core_emin, icg_core_emax, icg_e_buffer, icg_epsilon
  REAL(sp)       :: solve_time
  INTEGER        :: file_out
  LOGICAL        :: cg_ready, icg_ready_flag
  INTEGER        :: orig_dim, red_dim, nb_out
  REAL(dp)       :: cut_frac
  REAL(dp), ALLOCATABLE :: sorted_e(:)
  INTEGER,  ALLOCATABLE :: sidx(:)
  REAL(dp)       :: tmp_swap

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
  close(10)

  write(*,'(a)') '========================================'
  write(*,'(a,a)')   ' Structure:      ', trim(upt%gen_filename)
  write(*,'(a,l1)')  ' Relativistic:   ', upt%relat
  write(*,'(a,l1)')  ' Scaling:        ', upt%scaling
  write(*,'(a,a)')   ' Solver:         ', trim(solver_choice)
  write(*,'(a,i0)')  ' n_blocks:       ', n_blocks
  write(*,'(a,2f8.3)') ' CG window:    ', cg_emin, cg_emax
  write(*,'(a,2f8.3)') ' ICG core:     ', icg_core_emin, icg_core_emax
  write(*,'(a,f8.3)')  ' ICG buffer:   ', icg_e_buffer
  write(*,'(a,es10.2)')' ICG epsilon:  ', icg_epsilon
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

  ! ---- configure both CG modes (always enabled) ----------------------------
  call UPT_configure_coarse_graining(upt, .true., n_blocks, &
       cg_emin, cg_emax, imbalance)
  call UPT_configure_improved_cg(upt, .true., n_blocks, &
       icg_core_emin, icg_core_emax, icg_e_buffer, icg_epsilon, imbalance)

  ! ---- solver bookkeeping ---------------------------------------------------
  upt%num_vb   = nVB;  upt%num_cb   = nCB
  upt%start_vb = 1;    upt%start_cb = 1
  upt%min_iter = 2;    upt%long_iter = 30;  upt%max_iter = 100000
  upt%fast_tol = 1.0d-1; upt%long_tol = 1.0d-10; upt%ort_tol = 1.0d-5
  upt%solver_flag = 0;  upt%dynamic = .true.
  upt%seed_flag   = .false.; upt%bitoff = 0.1_dp
  upt%k_point = (/ 0.0d0, 0.0d0, 0.0d0 /)

  ! ---- build structure ------------------------------------------------------
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

  ! ---- build Hamiltonian (triggers cg_prepare + icg_prepare if enabled) ----
  write(*,'(a)') '========================================'
  write(*,'(a)') ' Building Hamiltonian...'
  call set_clock()
  call upt_hamiltonian(pupt)
  n_ham = upt%ham%nrow
  write(*,'(a,i0)') ' Full Hamiltonian dimension: ', n_ham

  ! query CG info
  call UPT_get_coarse_graining_info(upt, cg_ready, orig_dim, red_dim, nb_out, cut_frac)
  call UPT_get_improved_cg_info(upt, icg_ready_flag, orig_dim, red_dim, nb_out, cut_frac)

  ! ==========================================================================
  ! MODE 1: Standard full diagonalization
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 1: Standard full diagonalization'
  write(*,'(a)') '========================================'

  ! Temporarily disable CG for the standard solve
  upt%cg_enabled  = .false.
  upt%icg_enabled = .false.
  upt%num_vb = nVB; upt%num_cb = nCB

  num_ev = nVB + nCB
  allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
           upt%particles(num_ev), stat=err)
  upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0
  upt%verbose = 0
  call set_clock()
  select case (trim(solver_choice))
  case ('LK'); call lapack(upt)
  case ('JD'); call jd(upt)
  case ('LO'); call lanczos(upt)
  case default; write(*,*) 'Unknown solver: ', trim(solver_choice); stop 1
  end select
  solve_time = get_sclock()
  write(*,'(a,i0)')   ' Bands found:  ', size(upt%eigen_values)
  write(*,'(a,f10.3)') ' Solve time:   ', solve_time
  write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
  call write_eigenvalues('eigenvalues_standard.dat', upt%eigen_values, solve_time, &
       'STANDARD', n_ham, n_ham, 0.0_dp)
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)

  ! ==========================================================================
  ! MODE 2: Original coarse-graining (Liu et al.)
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 2: Original coarse-graining (Liu et al.)'
  write(*,'(a)') '========================================'
  call UPT_get_coarse_graining_info(upt, cg_ready, orig_dim, red_dim, nb_out, cut_frac)
  upt%cg_enabled  = .true.
  upt%icg_enabled = .false.
     if (cg_ready) then
        write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
             '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(orig_dim,dp)), '% reduction)'
        write(*,'(a,f6.4)') ' Cut fraction: ', cut_frac
        num_ev = red_dim
     else
        write(*,*) ' WARNING: CG not ready, skipping mode 2'
        goto 300
     end if
     allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
              upt%particles(num_ev), stat=err)
     upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0
     upt%verbose = 0
     call set_clock()
     select case (trim(solver_choice))
     case ('LK'); call lapack(upt)
     case ('JD'); call jd(upt)
     case ('LO'); call lanczos(upt)
     end select
     solve_time = get_sclock()
     write(*,'(a,i0)')   ' Bands found:  ', size(upt%eigen_values)
     write(*,'(a,f10.3)') ' Solve time:   ', solve_time
     write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
     call write_eigenvalues('eigenvalues_cg.dat', upt%eigen_values, solve_time, &
          'CG', orig_dim, red_dim, cut_frac)
     deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
300 continue

  ! ==========================================================================
  ! MODE 3: Improved coarse-graining
  ! ==========================================================================
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 3: Improved coarse-graining'
  write(*,'(a)') '========================================'
  call UPT_get_improved_cg_info(upt, icg_ready_flag, orig_dim, red_dim, nb_out, cut_frac)
  upt%cg_enabled  = .false.
  upt%icg_enabled = .true.
     if (icg_ready_flag) then
        write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
             '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(orig_dim,dp)), '% reduction)'
        write(*,'(a,f6.4)') ' Cut fraction: ', cut_frac
        num_ev = red_dim
     else
        write(*,*) ' WARNING: ICG not ready, skipping mode 3'
        goto 400
     end if
     allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
              upt%particles(num_ev), stat=err)
     upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0
     upt%verbose = 0
     call set_clock()
     select case (trim(solver_choice))
     case ('LK'); call lapack_icg(upt)
     case default
        write(*,*) ' ICG currently supports LK solver only'
        deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
        goto 400
     end select
     solve_time = get_sclock()
     write(*,'(a,i0)')   ' Bands found:  ', size(upt%eigen_values)
     write(*,'(a,f10.3)') ' Solve time:   ', solve_time
     write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
     call write_eigenvalues('eigenvalues_icg.dat', upt%eigen_values, solve_time, &
          'ICG', orig_dim, red_dim, cut_frac)
     deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
400 continue

  ! ---- cleanup --------------------------------------------------------------
  call destroy_matrix(upt%ham)
  call upt_mpi_end

contains

  ! ---------------------------------------------------------------------------
  subroutine write_eigenvalues(fname, evals, t, mode, ndim, nred, cut)
    character(*), intent(in) :: fname, mode
    real(dp), intent(in)     :: evals(:), cut
    real(sp), intent(in)     :: t
    integer, intent(in)      :: ndim, nred
    integer :: n, fu, ii, jj, itmp
    real(dp), allocatable :: se(:)
    integer,  allocatable :: si(:)
    real(dp) :: stmp
    n = size(evals)
    allocate(se(n), si(n))
    se = evals
    do ii = 1, n; si(ii) = ii; end do
    do ii = 1, n-1          ! bubble sort descending
       do jj = ii+1, n
          if (se(jj) > se(ii)) then
             itmp=si(ii); si(ii)=si(jj); si(jj)=itmp
             stmp=se(ii); se(ii)=se(jj); se(jj)=stmp
          end if
       end do
    end do
    open(newunit=fu, file=trim(fname), status='replace', action='write')
    write(fu,'(a,a)')   '# Mode: ', trim(mode)
    write(fu,'(a,i0)')  '# Full dimension:    ', ndim
    write(fu,'(a,i0)')  '# Reduced dimension: ', nred
    write(fu,'(a,f8.2,a)') '# Rank reduction: ', &
         100.0_dp*(1.0_dp - real(nred,dp)/real(ndim,dp)), ' %'
    write(fu,'(a,f10.4)') '# Cut fraction:   ', cut
    write(fu,'(a,f12.6)') '# Solve time (s): ', t
    write(fu,'(a,i0)')  '# Total bands:    ', n
    write(fu,'(a)') '#'
    write(fu,'(a)') '# Index    Energy(eV)'
    do ii = 1, n
       write(fu,'(i6,2x,f16.8)') ii, se(ii)
    end do
    close(fu)
    write(*,'(a,a)') ' Output: ', trim(fname)
    deallocate(se, si)
  end subroutine write_eigenvalues

  ! ---------------------------------------------------------------------------
  ! ICG LAPACK solver: assembles dense icg_ham, diagonalizes, lifts to physical.
  subroutine lapack_icg(upt)
    use lapack_driver, only : lapack_icg_solve
    use coarse_grain,  only : icg_lift
    use precision,     only : dp
    use upt_param,     only : OUPT
    type(OUPT), intent(inout) :: upt
    integer :: nred, nfull, ii
    complex(dp), allocatable :: h(:,:)
    real(dp),    allocatable :: eval(:)
    nred  = upt%icg_ham%nrow
    nfull = upt%ham%nrow
    call lapack_icg_solve(upt, h, eval)
    if (associated(upt%eigen_values))  deallocate(upt%eigen_values)
    if (associated(upt%eigen_vectors)) deallocate(upt%eigen_vectors)
    if (associated(upt%particles))     deallocate(upt%particles)
    allocate(upt%eigen_values(nred), upt%eigen_vectors(nfull,nred), upt%particles(nred))
    do ii = 1, nred
       upt%eigen_values(ii) = eval(ii)
       upt%particles(ii)    = 0
    end do
    call icg_lift(upt, h, upt%eigen_vectors)
    deallocate(h, eval)
  end subroutine lapack_icg

end program test_supercell
