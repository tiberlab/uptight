program test_magnetic
  ! Test program for Peierls substitution implementation

  USE precision
  USE globals, only : MST, LST, GAUGE_NONE, GAUGE_LANDAU_Z, GAUGE_SYMMETRIC_Z, &
                      magnetic_field_vector, gauge_choice ! Import B-field globals
  USE mpi_globals, only : upt_mpi_init, upt_mpi_end, id0, num_procs, upt_comm, ierr
  USE upt_param, only : OUPT
  USE struct_building, only : init_structure, make_basis, init_basis, subs_dg_ions
  USE neighbours, only : refine_neighbours_map, check_input_nn_list, check_nn_map
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials, destroy_basis, destroy_material, &
                        destroy_structure, nullify_basis, nullify_structure
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   sort_states, set_max_order
  USE alloys, only : init_mat_ion
  USE TB_ham, only : sparse_ham, check_if_hermitian
  USE lapack_driver, only : lapack
  USE lanczos_driver, only : lanczos ! Added Lanczos as an option
  USE sparse_matrix, only : destroy_matrix
  USE uptight, only : upt_nullify_all, upt_destruct, upt_version, upt_set_defaults, &
                      upt_init, upt_hamiltonian, upt_lapack, upt_lanczos, upt_alloc_eigv
  USE clock

  IMPLICIT NONE

  TYPE(OUPT), TARGET :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER :: i, j, k, n_ham, num_ev, err
  CHARACTER(2) :: calc_method = 'LK' ! Default to Lapack

  ! --- Initialize MPI ---
  call upt_mpi_init(0)

  ! --- Initialize UPT object ---
  pupt => upt
  call upt_nullify_all(pupt)
  call upt_set_defaults(pupt) ! Use defaults from upt_param

  if (id0) call upt_version(pupt)

  ! --- Read Basic Parameters from Standard Input ---
  if (id0) write(*,*) 'Reading input parameters...'
  read(*,*) upt%database_path ! line 1: database path (e.g., ../../parameters/Jancu)
  upt%work_path = "./"
  read(*,*) upt%gen_filename  ! line 2: name of .upg file (e.g., ../../src/tests/GaN_column/temp_c.gen)
  read(*,*) upt%relat         ! line 3: relativistic flag (.TRUE. or .FALSE.)
  read(*,*) upt%scaling       ! line 4: harrison scaling (.TRUE. or .FALSE.)
  read(*,*) upt%c_axis(:)     ! line 5: c axis (e.g., 0.0 0.0 1.0)

  ! --- Read Magnetic Field Parameters ---
  read(*,*) upt%use_magnetic_field ! line 6: Use magnetic field? (.TRUE. or .FALSE.)
  if (upt%use_magnetic_field) then
    read(*,*) upt%magnetic_field_vector(:) ! line 7: Bx, By, Bz (Tesla)
    read(*,*) upt%gauge_choice             ! line 8: Gauge choice (1=Landau Z, 2=Symmetric Z)
    if (id0) then
       write(*,*) '(main) Magnetic field enabled:'
       write(*,'(A,3F10.4)') ' B = ', upt%magnetic_field_vector
       write(*,'(A,I2)')     ' Gauge = ', upt%gauge_choice
    endif
  else
    if (id0) write(*,*) '(main) Magnetic field disabled.'
  endif

  ! --- Read Solver Parameters ---
  read(*,*) calc_method       ! line 9: Solver ('LK' for Lapack, 'LO' for Lanczos)
  read(*,*) upt%num_vb        ! line 10: number of valence bands (below guess)
  read(*,*) upt%num_cb        ! line 11: number of conduction bands (above guess)
  read(*,*) upt%lambda_vb     ! line 12: valence guess (for Lanczos)
  read(*,*) upt%lambda_cb     ! line 13: conduction guess (for Lanczos)

  ! --- Set some defaults suitable for testing ---
  upt%verbose = 2             ! Moderate verbosity
  upt%sparse_format = "F"     ! Use Full format for easier debugging
  upt%k_point = (/0.0d0, 0.0d0, 0.0d0/) ! Gamma point for finite systems
  upt%check_bondmap = .true.
  upt%hybrid_passivation = .true. ! Often needed for finite systems

  ! --- Initialize Uptight Library ---
  if (id0) write(*,*) '(main) Initializing Uptight...'
  call upt_init(pupt)

  ! --- Build Hamiltonian ---
  if (id0) write(*,*) '(main) Building Hamiltonian...'
  call upt_hamiltonian(pupt)

  if (upt%sparse_format.eq.'F') then
     if (id0) write(*,*) '(main) Checking Hamiltonian Hermiticity...'
     call check_if_hermitian(upt%ham)
  end if

  ! --- Allocate Eigenvalue/vector Storage ---
  n_ham = upt%ham%nrow
  num_ev = upt%num_vb + upt%num_cb
  if (num_ev .le. 0) then
     if (id0) write(*,*) 'ERROR: num_vb + num_cb must be > 0'
     call upt_mpi_end
     stop 1
  endif
  if (id0) write(*,*) '(main) Hamiltonian dimension:', n_ham
  if (id0) write(*,*) '(main) Number of eigenvalues requested:', num_ev

  call upt_alloc_eigv(pupt) ! Use library allocation routine

  ! --- Diagonalize ---
  if (id0) write(*,*) '(main) Diagonalizing Hamiltonian using ', TRIM(calc_method), '...'
  call set_clock()
  select case ( TRIM(calc_method) )
     case ('LK')
        !if (num_ev .ne. n_ham) then
        !   if (id0) write(*,*) 'WARNING: Lapack (LK) computes all eigenvalues. num_vb/num_cb ignored.'
        !   upt%num_vb = n_ham ! Adjust for Lapack output
        !   upt%num_cb = 0
        !endif
        call upt_lapack(pupt)
     case ('LO')
        call upt_lanczos(pupt) ! Use default Lanczos settings from upt_param
     case default
        if (id0) write(*,*) 'ERROR: Unknown solver method: ', TRIM(calc_method)
        call upt_mpi_end
        stop 1
  end select
  if (id0) call write_clock()

  ! --- Print Results (Eigenvalues) ---
  if (id0) then
     write(*,*) '----------------------------------------'
     write(*,*) 'Computed Eigenvalues (eV):'
     write(*,*) '----------------------------------------'
     if (TRIM(calc_method) .eq. 'LK') then
        ! Lapack computes all eigenvalues
        do i = 1, num_ev
           write(*,'(I5, F15.8)') i, upt%eigen_values(i)
        enddo
     else
        ! Lanczos computes requested number
        do i = 1, num_ev
           write(*,'(I5, F15.8)') i, upt%eigen_values(i)
        enddo
     endif
     write(*,*) '----------------------------------------'
  endif

  ! --- Clean up ---
  if (id0) write(*,*) '(main) Cleaning up...'
  call upt_destruct(pupt)

  call upt_mpi_end

end program test_magnetic