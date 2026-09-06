! Standalone test for coarse-graining feature
! Reads config from TEST/supercell_coarse-grain/config
! Solves at Gamma point only
program test_supercell

  USE precision
  USE globals, only : MST, LST
  USE mpi_globals, only : upt_mpi_init, upt_mpi_end
  USE upt_param, only : OUPT
  USE struct_building, only : init_structure, make_basis, init_basis, subs_dg_ions
  USE neighbours, only : refine_neighbours_map, check_input_nn_list, &
                         check_nn_map, write_neighbours_map
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   sort_states, set_max_order
  USE alloys, only : init_mat_ion
  USE uptight, only : UPT_configure_coarse_graining, UPT_get_coarse_graining_info, &
                      upt_hamiltonian
  USE lapack_driver, only : lapack
  USE JD_driver, only : jd
  USE lanczos_driver, only : lanczos
  USE sparse_matrix, only : destroy_matrix
  USE clock, only : set_clock, get_sclock

  IMPLICIT NONE

  TYPE(OUPT), TARGET :: upt
  TYPE(OUPT), POINTER :: pupt
  
  INTEGER :: i, j, err, n_ham, num_ev
  CHARACTER(LST) :: config_file, output_file
  CHARACTER(MST) :: solver_choice
  INTEGER :: n_blocks, nVB, nCB
  LOGICAL :: use_coarse_grain
  REAL(dp) :: Emin, Emax, imbalance
  REAL(sp) :: solve_time
  INTEGER :: file_out
  LOGICAL :: cg_ready
  INTEGER :: original_dim, reduced_dim, nb_out
  REAL(dp) :: cut_fraction
  REAL(dp), ALLOCATABLE :: sorted_energies(:)
  INTEGER, ALLOCATABLE :: sort_idx(:)
  REAL(dp) :: temp_swap

  ! Initialize MPI
  call upt_mpi_init(0)

  ! Point to the upt structure
  pupt => upt

  ! Set default paths
  upt%database_path = './'
  upt%work_path = './'
  upt%out_path = './'
  upt%gen_out = 'out.gen'
  upt%state_file = 'states.data'
  upt%sparse_format = 'U'
  upt%verbose = 10

  ! Read configuration file
  config_file = 'config'
  open(10, file=trim(config_file), status='old', action='read', iostat=err)
  if (err /= 0) then
     write(*,*) 'ERROR: Cannot open config file: ', trim(config_file)
     stop 1
  end if

  write(*,*) '========================================='
  write(*,*) 'Reading configuration from: ', trim(config_file)
  write(*,*) '========================================='

  read(10,*) upt%gen_filename      ! Line 1: alloy.upg file
  read(10,*) upt%relat             ! Line 2: relativistic flag
  read(10,*) upt%scaling           ! Line 3: Harrison scaling
  read(10,*) upt%c_axis(:)         ! Line 4: c-axis direction
  read(10,*) solver_choice         ! Line 5: solver (LK=LAPACK, JD=Jacobi-Davidson, LO=Lanczos)
  read(10,*) nVB                   ! Line 6: number of valence bands
  read(10,*) nCB                   ! Line 7: number of conduction bands
  read(10,*) upt%lambda_vb         ! Line 8: valence band guess (eV)
  read(10,*) upt%lambda_cb         ! Line 9: conduction band guess (eV)
  read(10,*) n_blocks              ! Line 10: number of blocks
  read(10,*) use_coarse_grain      ! Line 11: enable coarse-graining (.true./.false.)
  read(10,*) Emin                  ! Line 12: energy window minimum (eV)
  read(10,*) Emax                  ! Line 13: energy window maximum (eV)
  read(10,*) imbalance             ! Line 14: METIS imbalance tolerance
  close(10)

  write(*,*) 'Structure file:       ', trim(upt%gen_filename)
  write(*,*) 'Relativistic:         ', upt%relat
  write(*,*) 'Harrison scaling:     ', upt%scaling
  write(*,*) 'C-axis:               ', upt%c_axis
  write(*,*) 'Solver:               ', trim(solver_choice)
  write(*,*) 'Valence bands:        ', nVB
  write(*,*) 'Conduction bands:     ', nCB
  write(*,*) 'VB guess (eV):        ', upt%lambda_vb
  write(*,*) 'CB guess (eV):        ', upt%lambda_cb
  write(*,*) 'Number of blocks:     ', n_blocks
  write(*,*) 'Use coarse-graining:  ', use_coarse_grain
  write(*,*) 'Energy window (eV):   ', Emin, ' to ', Emax
  write(*,*) 'METIS imbalance:      ', imbalance
  write(*,*) '========================================='

  ! Set other parameters
  upt%structure%gen_filename = upt%gen_filename
  upt%d_onsite_shift_flag = .true.
  upt%potential_flag = .false.
  upt%syst_rotated = .false.
  upt%ioutput_flag = .false.
  upt%optmat = .false.
  upt%poldir = 3
  upt%d_H = 0.1d0
  upt%E_H = -200.0d0
  upt%estimate_factor = 1.0
  upt%check_bondmap = .false.

  ! Set spin based on relativistic flag
  if (upt%relat) then
     upt%n_spin = 2
     write(*,*) 'Relativistic calculation (n_spin=2)'
  else
     upt%n_spin = 1
     write(*,*) 'Non-relativistic calculation (n_spin=1)'
  end if

  ! Configure coarse-graining
  if (use_coarse_grain) then
     write(*,*) '========================================='
     write(*,*) 'COARSE-GRAINING ENABLED (n_blocks=', n_blocks, ')'
     write(*,*) '========================================='
     call UPT_configure_coarse_graining(upt, .true., n_blocks, Emin, Emax, imbalance)
     upt%num_vb = 1
     upt%num_cb = 1
     upt%start_vb = 1
     upt%start_cb = 1
  else
     write(*,*) '========================================='
     write(*,*) 'STANDARD MODE (no coarse-graining)'
     write(*,*) '========================================='
     call UPT_configure_coarse_graining(upt, .false., 1, Emin, Emax, imbalance)
     upt%num_vb = nVB
     upt%num_cb = nCB
     upt%start_vb = 1
     upt%start_cb = 1
  end if

  ! Solver parameters
  upt%min_iter = 2
  upt%long_iter = 30
  upt%max_iter = 100000
  upt%fast_tol = 1.0d-1
  upt%long_tol = 1.0d-10
  upt%ort_tol = 1.0d-5
  upt%solver_flag = 0
  upt%dynamic = .true.
  upt%seed_flag = .false.
  upt%bitoff = 0.1_dp

  ! Set k-point to Gamma
  upt%k_point = (/ 0.0d0, 0.0d0, 0.0d0 /)
  write(*,*) 'K-point: Gamma (0, 0, 0)'

  ! Build structure
  write(*,*) '========================================='
  write(*,*) 'BUILDING STRUCTURE'
  write(*,*) '========================================='
  
  call set_machine_acc

  write(*,*) 'Reading and initializing structure...'
  call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat, &
                      upt%interfaces, upt%nr_int)

  write(*,*) 'Checking input nearest neighbor list...'
  call check_input_nn_list(upt%structure)

  write(*,*) 'Reading material database...'
  do i = 1, upt%nr_mat
     call read_data(upt%materials(i), upt%work_path, upt%database_path)
  end do

  write(*,*) 'Making basis...'
  call make_basis(upt%verbose, upt%structure, upt%basis)

  write(*,*) 'Setting up states and couplings...'
  call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, &
                                upt%ref_couplings, upt%n_ref_cpl)

  write(*,*) 'Sorting states...'
  do i = 1, upt%nr_mat
     call sort_states(upt%materials(i), upt%ref_states, upt%ref_couplings)
  end do
  do i = 1, size(upt%materials)
     do j = 1, size(upt%materials)
        if (associated(upt%interfaces(i,j)%nr_parents)) then
           call sort_states(upt%interfaces(i,j), upt%ref_states, upt%ref_couplings)
        end if
     end do
  end do

  write(*,*) 'Setting max order...'
  call set_max_order(upt%materials)

  write(*,*) 'Initializing ions...'
  do i = 1, upt%nr_mat
     call init_mat_ion(upt%materials(i))
  end do

  call write_materials(upt%nr_mat, upt%materials)

  write(*,*) 'Initializing ion valences...'
  call init_basis(upt%basis, upt%materials)

  call write_basis(upt%basis, 0)

  write(*,*) 'Building nearest neighbor table...'
  call refine_neighbours_map(upt%structure, upt%basis, upt%materials, upt%nn_map)

  call write_neighbours_map(upt%nn_map)
  call check_nn_map(upt%nn_map)

  write(*,*) 'Treating dangling bonds...'
  call subs_dg_ions(upt%basis, upt%materials, upt%nn_map)

  write(*,*) 'Initializing n_st...'
  call init_n_st(upt%basis, upt%materials)

  ! Build Hamiltonian (this triggers coarse-graining if enabled)
  write(*,*) '========================================='
  write(*,*) 'BUILDING HAMILTONIAN'
  write(*,*) '========================================='
  
  ! Start timing from here for fair comparison (includes coarse-graining preparation)
  call set_clock()
  
  call upt_hamiltonian(pupt)
  
  n_ham = upt%ham%nrow
  write(*,*) 'Hamiltonian dimension: ', n_ham

  ! Query coarse-graining status AFTER sparse_ham which triggers coarse-graining
  call UPT_get_coarse_graining_info(upt, cg_ready, original_dim, reduced_dim, &
                                    nb_out, cut_fraction)
  
  if (cg_ready) then
     write(*,*) '----------------------------------------'
     write(*,*) 'COARSE-GRAINING STATISTICS:'
     write(*,*) '  Original dimension:  ', original_dim
     write(*,*) '  Reduced dimension:   ', reduced_dim
     write(*,*) '  Rank reduction:      ', 100.0_dp*(1.0_dp-real(reduced_dim,dp)/real(original_dim,dp)), '%'
     write(*,*) '  Number of blocks:    ', nb_out
     write(*,*) '  Cut edge fraction:   ', cut_fraction
     write(*,*) '----------------------------------------'
     ! For coarse-graining, we will get ALL reduced_dim eigenvalues
     num_ev = reduced_dim
  else
     ! For standard mode, use requested VB and CB counts
     num_ev = upt%num_vb + upt%num_cb
  end if

  ! Allocate eigensystem arrays
  allocate(upt%eigen_values(num_ev), stat=err)
  allocate(upt%eigen_vectors(n_ham, num_ev), stat=err)
  allocate(upt%particles(num_ev), stat=err)
  if (err /= 0) stop 'Allocation error for eigensystem'
  
  upt%eigen_values = 0.0d0
  upt%eigen_vectors = (0.0d0, 0.0d0)
  upt%particles = 0

  ! Solve eigensystem
  write(*,*) '========================================='
  write(*,*) 'SOLVING EIGENSYSTEM'
  write(*,*) '========================================='
  write(*,*) 'Solver: ', trim(solver_choice)
  write(*,*) 'Number of eigenvalues to compute: ', num_ev

  upt%verbose = 0  ! Reduce verbosity during solve

  select case (trim(solver_choice))
  case ('LK')
     write(*,*) 'Using LAPACK dense solver...'
     call lapack(upt)
  case ('JD')
     write(*,*) 'Using Jacobi-Davidson iterative solver...'
     call jd(upt)
  case ('LO')
     write(*,*) 'Using Lanczos iterative solver...'
     call lanczos(upt)
  case default
     write(*,*) 'ERROR: Unknown solver: ', trim(solver_choice)
     stop 1
  end select

  solve_time = get_sclock()

  ! Report results
  write(*,*) '========================================='
  write(*,*) 'RESULTS'
  write(*,*) '========================================='
  write(*,*) 'Total bands found:    ', size(upt%eigen_values)
  write(*,*) 'Solve time (s):       ', solve_time
  
  if (size(upt%eigen_values) > 0) then
     write(*,*) 'Energy range (eV):    ', minval(upt%eigen_values), ' to ', maxval(upt%eigen_values)
  end if
  
  if (cg_ready) then
     write(*,*) 'Config energy window: ', Emin, ' to ', Emax
  end if

  ! Sort eigenvalues in descending order
  num_ev = size(upt%eigen_values)
  allocate(sorted_energies(num_ev), stat=err)
  allocate(sort_idx(num_ev), stat=err)
  if (err /= 0) stop 'Allocation error for sorting'

  sorted_energies = upt%eigen_values
  
  ! Simple bubble sort (descending)
  do i = 1, num_ev
     sort_idx(i) = i
  end do
  
  do i = 1, num_ev - 1
     do j = i + 1, num_ev
        if (sorted_energies(j) > sorted_energies(i)) then
           ! Swap indices
           err = sort_idx(i)
           sort_idx(i) = sort_idx(j)
           sort_idx(j) = err
           
           ! Swap energies
           temp_swap = sorted_energies(i)
           sorted_energies(i) = sorted_energies(j)
           sorted_energies(j) = temp_swap
        end if
     end do
  end do

  ! Write energies to file
  output_file = 'eigenvalues.dat'
  open(file_out, file=trim(output_file), status='replace', action='write', iostat=err)
  if (err /= 0) then
     write(*,*) 'ERROR: Cannot open output file: ', trim(output_file)
     stop 1
  end if

  write(file_out, '(A)') '# Eigenvalues (eV) - sorted descending'
  write(file_out, '(A, I0)') '# Total number of bands: ', size(upt%eigen_values)
  write(file_out, '(A, F12.6)') '# Solve time (s): ', solve_time
  if (cg_ready) then
     write(file_out, '(A)') '# Coarse-graining: ENABLED'
     write(file_out, '(A, I0, A, I0)') '# Dimension: ', original_dim, ' -> ', reduced_dim
     write(file_out, '(A, F8.2, A)') '# Rank reduction: ', &
          100.0_dp*(1.0_dp-real(reduced_dim,dp)/real(original_dim,dp)), '%'
     write(file_out, '(A, F8.3, A, F8.3, A)') '# Config energy window: [', Emin, ', ', Emax, '] eV'
  else
     write(file_out, '(A)') '# Coarse-graining: DISABLED'
  end if
  write(file_out, '(A)') '#'
  write(file_out, '(A)') '# Index    Energy(eV)'
  
  do i = 1, num_ev
     write(file_out, '(I6, 2X, F16.8)') i, sorted_energies(i)
  end do

  close(file_out)

  write(*,*) '========================================='
  write(*,*) 'Output written to: ', trim(output_file)
  write(*,*) '========================================='

  ! Cleanup
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles, stat=err)
  deallocate(sorted_energies, sort_idx, stat=err)
  call destroy_matrix(upt%ham)

  call upt_mpi_end

end program test_supercell
