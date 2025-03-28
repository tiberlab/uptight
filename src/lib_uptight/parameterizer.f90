program parameterizer

  USE precision
  USE globals, only : MST, LST
  USE mpi_globals, only : upt_mpi_init, upt_mpi_end
  USE upt_param, only : OUPT
  USE struct_building, only : init_structure, make_basis, init_basis, &
                              write_gen, read_gen, subs_dg_ions
  USE neighbours, only : make_neighbours_map, refine_neighbours_map, &
                         check_input_nn_list, check_nn_map, write_neighbours_map
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   sort_states, set_max_order
  USE alloys, only : init_mat_ion
  USE TB_ham, only : sparse_ham, check_if_hermitian, check_if_antihermitian, hermitianize
  USE lanczos_driver, only : lanczos
  USE savemofile, only : writemofile, write_eigenstates, write_cube
  USE lapack_driver, only : lapack
  USE sparse_matrix, only : destroy_matrix,  write_sprs_to_coo
  !USE arpack, only : arpack_dr
  USE uptight, only : upt_get_mat_el
  USE input_output, only : open_file

  IMPLICIT NONE


  TYPE(OUPT), TARGET :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER :: i, j, k, file_num, n_ham, num_ev, err, at, bg, st, ao, nk, n_col, n_ao_type
  CHARACTER(2) :: calc
  COMPLEX(dp) :: matel
  CHARACTER(LST) :: kfile, file_name
  REAL(dp), DIMENSION(3) :: k_xyz, k_old, strain
  REAL(dp) :: Lk, Lk_min, weight, pop, norm
  CHARACTER(2) :: id
  real(dp), dimension(16) :: coef2

  complex(dp), dimension(:,:), allocatable :: temp_eivec
  real(dp), dimension(:,:), allocatable :: eigensol
  !real(dp), dimension(:,:), allocatable :: dos(:,:)
  !real(dp), dimension(:), allocatable :: energy_grid(:)

  LOGICAL :: exist, coef2_cal

  integer :: comm
  integer :: upt_group

  call upt_mpi_init(0)

  read(*,*) upt%database_path ! line 1 of `parameterizer.in`: database path
  upt%work_path = "./"
  read(*,*) upt%gen_filename  ! line 2 of `parameterizer.in`: name of `.upg` file
  upt%gen_out = "out.gen"
  upt%state_file = "states.data"  
  upt%sparse_format = "U"
  upt%verbose = 10
  
  upt%structure%gen_filename = upt%gen_filename
  read(*,*) upt%relat  ! line 3 of `parameterizer.in`: relativistic flag
  upt%out_path = "./"
  read(*,*) upt%scaling  ! line 4 of `parameterizer.in`: harrison scaling

  upt%d_onsite_shift_flag=.true.
  upt%potential_flag=.false.
  upt%syst_rotated=.false.   ! if .true. means the structure is already
                             ! rotated along the (111) direction

  !upt%c_axis= (/ 0.d0, 0.d0, 1.d0 /)
  !upt%c_axis= (/ 1.d0, 0.d0, 0.d0 /)
  read(*,*) upt%c_axis(:)  ! line 5 of `parameterizer.in`: c axis (not relevant for Tan but important for Jancu wurtzite and 2D MoS2)

  upt%ioutput_flag=.false.

  upt%optmat=.false.
  upt%poldir = 3    

  upt%d_H = 0.1d0 
  upt%E_H = -200.0d0

  upt%estimate_factor = 1.0

  if(upt%relat) then 
     upt%n_spin = 2
     write(*,*) "(main) relativistic calculation"
     write(*,*) "(main) n. spins:",upt%n_spin
  else
     upt%n_spin = 1
     write(*,*) "(main) non relativistic calculation"
     write(*,*) "(main) n. spins:",upt%n_spin
  end if

  call set_machine_acc


  write(*,*) '(main) read and init structure'
  call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat, upt%interfaces, upt%nr_int)
  !call make_struct(upt%verbose, upt%structure, upt%basis, upt%materials, upt%nr_mat)
  
  write(*,*) '(main) check input nn_list'
  call check_input_nn_list(upt%structure)

  write(*,*) '(main) read material database:'
  do i = 1, upt%nr_mat
     call read_data(upt%materials(i), upt%work_path, upt%database_path)
  end do

  !==================================================================================

  write(*,*) '(main) make basis'
  call make_basis(upt%verbose, upt%structure, upt%basis)


  write(*,*) '(main) setup states and couplings:'
  call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, & 
       upt%ref_couplings, upt%n_ref_cpl)
  
  if(upt%verbose.gt.0) write(*,*) '(main) sorting states according to reference:'
  do i = 1, upt%nr_mat
      call sort_states(upt%materials(i), upt%ref_states, upt%ref_couplings)
  end do
                                    

  write(*,*) '(main) set max order in each material:'
  call set_max_order(upt%materials)

  write(*,*) '(main) init ions'
  do i=1,upt%nr_mat
     call init_mat_ion(upt%materials(i))
  enddo

  write(*,*) '(main) materials info:'  
  call write_materials(upt%nr_mat,upt%materials)

  write(*,*) '(main) init ion valences'
  call init_basis(upt%basis, upt%materials)


  write(*,*) '(main) structure info:'
  call write_basis(upt%basis,0)


  write(*,*) '(main) Build Nearest Neighbour Table'
  call refine_neighbours_map(upt%structure, upt%basis, upt%materials, upt%nn_map)


  write(*,*) '(main) write neighbour list:'
  call write_neighbours_map(upt%nn_map)
  
  write(*,*) '(main) check neighbour list'
  call check_nn_map(upt%nn_map)

  write(*,*) '(main) treat dangling bonds'
  call subs_dg_ions(upt%basis,upt%materials,upt%nn_map)

  !write(*,*) '(main) write gen'
  !call write_gen(upt%gen_out,upt%basis,upt%nn_map, upt%syst_rotated, .false.)


  write(*,*) '(main) init n_st'
  call init_n_st(upt%basis,upt%materials)

  !call write_gen(upt_in%gen_out, upt%basis, upt%nn_map)  

  write(*,*) '(main) read solver options'


  calc = 'DF'
  read(*,*) calc  ! line 6 of `parameterizer.in`: eigensolver (e.g. LK = lapack)
  
  upt%start_vb = 1  
  upt%start_cb = 1 

  read(*,*) upt%num_vb   ! line 7 of `parameterizer.in`: number of valence bands (actually, number of bands below valence guess)
  read(*,*) upt%num_cb   ! line 8 of `parameterizer.in`: number of conduction bands (actually, number of bands above conduction guess)

  read(*,*) upt%lambda_vb  ! line 9 of `parameterizer.in`: valence guess
  read(*,*) upt%lambda_cb  ! line 10 of `parameterizer.in`: conduction guess

  upt%min_iter = 2
  upt%long_iter = 30
  upt%max_iter = 100000

  upt%fast_tol = 1.0d-1
  upt%long_tol = 1.0d-10
  upt%ort_tol  = 1.0d-5
     
  upt%solver_flag = 0
  upt%dynamic        = .true. 
  upt%seed_flag      = .false.

  read(*,*) kfile   ! line 11 of `parameterizer.in`: name of file containing list of k-points in fractions coordinates


  read(*,*) strain(1:3)  ! line 12 of `parameterizer.in`: strain

  upt%basis%strain(1,:) = strain(1)
  upt%basis%strain(2,:) = strain(2)
  upt%basis%strain(3,:) = strain(3)

  read(*,*) coef2_cal   ! line 13 of `parameterizer.in`: whether calculating coef2 or not, T or F

  if (coef2_cal) then
    read(*,*) n_ao_type   ! line 14 of `parameterizer.in`: number of AO type per atom: 1 (s), 2 (s, p), 3 (s, p, d), 4 (s, p, d, f)
  end if

  ! read potential values
  !read(*,*) upt%potential_flag   ! line 13 of `parameterizer.in`: potential flag
!
  !if (upt%potential_flag) then
  !  upt%potential_file = "pot_on_atoms.dat"
  !  allocate(upt%pot_data(upt%basis%n_basis), STAT = err)
  !  write(*,*) "(main) pot_data size ", upt%basis%n_basis
  !  inquire(FILE=upt%potential_file,EXIST=exist)  
  !  if (exist) then
  !    write(*,*) '(main) Read Potential File: ',trim(upt%potential_file)
  !    open(51,file=upt%potential_file)
!
  !    do i= 1, upt%basis%n_basis
  !      read(51,*) upt%pot_data(i)
  !    end do
!
  !    close(51)
  !  else
  !    upt%pot_data = 0.d0     
  !  endif
  !endif 

  open(15,file=trim(kfile), STATUS='OLD')

  file_name = trim(upt%out_path)//'etb.csv'
  
  CALL open_file( file_name, file_num, operation = "write", &
       format_flag = .TRUE., replace_flag = .TRUE. )
  close(file_num)

  CALL open_file( file_name, file_num, operation = "write", &
       format_flag = .TRUE., replace_flag = .FALSE. )

  !write( file_num, '(a36)', advance='NO' ) '# k-path, Energies, Coefficients^2'
  !do i =1, upt%num_vb
  !  write(id,'(i2.2)') i
  !  write( file_num, '(a4,a2)', advance='NO' ) '  Ev',id
  !enddo
  !do i =1, upt%num_cb
  !  write(id,'(i2.2)') i
  !  write( file_num, '(a4,a2)', advance='NO' ) '  Ec',id
  !enddo
  !write( file_num, *)

  Lk = 0.d0
  k_xyz = 0.d0
  k_old = 0.d0

  upt%verbose=0

  num_ev = upt%num_vb + upt%num_cb 
  
  !import kfile and count how many k-points
  nk = 0
  do k = 1, 100000
     read(15,fmt=*,end=80) upt%k_point   !`kpoints` file needs to be in the format: kx,ky,kz
     nk = nk + 1
  enddo
80  rewind(15) !rewind the file 15 to reader once again from the beginning

  write(*,*) '(main) nk = ', nk

  n_col = 2
  if (coef2_cal) then
    n_col = n_col + n_ao_type*upt%basis%n_basis
  end if

  allocate(eigensol(num_ev*nk, n_col), STAT = err)
  ! 2 + 3*upt%basis%n_basis: 2 columns for k-path and eigenvalues, 3 columns for s,p,d coefficients for each atom
  IF (err.NE.0) stop 'Allocation error'
  eigensol = 0.d0

  do k = 1, 100000
     read(15,fmt=*,end=90) upt%k_point   !`kpoints` file needs to be in the format: kx,ky,kz
 
     write(*,'(a,i3, 3(f8.4))') '(k-loop) ', k, (upt%k_point(i), i=1, 3)

     call sparse_ham(upt)
     
     !if (upt%sparse_format.eq.'F') then
        !write(*,*) '(main) make full matrix Hermitian'
     !   call hermitianize(upt%ham)
        
        !write(*,*) '(main) check if hermitian'
     !   call check_if_hermitian(upt%ham)
     !end if

     !write(*,*) 'saving H_COO.DAT...'
     !call  write_sprs_to_coo(upt%ham%M,upt%ham%Mij)

     n_ham = upt%ham%nrow

     allocate(upt%eigen_values(num_ev), STAT = err)
     allocate(upt%eigen_vectors(n_ham,num_ev), STAT = err)
     IF (err.NE.0) stop 'Allocation error'
     upt%eigen_values = 0.d0
     upt%eigen_vectors = (0.d0,0.d0)

     select case ( calc ) 
     case ('LO')
        !write(*,*) '(main) lanczos diagonalization'
        call lanczos(upt)
     case ('LK')
        !write(*,*) '(main) lapack diagonalization'
        call lapack(upt)
     case ('AR')
        !write(*,*) '(main) arpack diagonalization'
        !call arpack_dr(upt)
     case ('DF')
        stop 'invalid calculation mode'
     end select     

     k_xyz = MATMUL(upt%basis%rec_latt, upt%k_point)

     !write(*,'(a, 3(f8.4))') '(k-loop) ', (k_xyz(i), i= 1, 3)

     Lk = Lk + sqrt((k_xyz(1)-k_old(1))**2 + (k_xyz(2)-k_old(2))**2 + &
                                             (k_xyz(3)-k_old(3))**2 )

     if (k.EQ.1) Lk_min = Lk
     Lk = Lk - Lk_min !shift k-path to 0.0 at the beginning of the path

     k_old = k_xyz

     allocate(temp_eivec(num_ev, n_ham), STAT = err)
     IF (err.NE.0) stop 'Allocation error'

     temp_eivec = transpose(upt%eigen_vectors)
     !write( file_num,'(ES19.8)', advance='NO' ) Lk

     do i = 1, num_ev

       eigensol(k+(i-1)*nk, 1) = Lk
       eigensol(k+(i-1)*nk, 2) = upt%eigen_values(i)

         !WRITE( file_num, '(( x1 f13.8 ))', advance='NO' ) &
         !     upt%eigen_values(i)      
        if (coef2_cal) then
          call AO_decomposition(temp_eivec(i,:), coef2, upt%basis%n_basis, upt%relat, n_ao_type)

            !write( file_coef2,'(ES15.8)', advance='NO' ) Lk
            !write( file_coef2,'(( x1 f13.8 ))', advance='NO' ) upt%eigen_values(i)
          do ao = 1, n_ao_type*upt%basis%n_basis !types of orbitals: s, p, d, ... for each atom
            !write( file_coef2, '(( x1 f8.5 ))', advance='NO' ) coef2(ao)
            eigensol(k+(i-1)*nk, 2+ao) = coef2(ao)
          end do
           !write(file_coef2,*)
        end if
        
     enddo

     !write(file_coef2,*)
     !write(file_num,*) 

     deallocate(upt%eigen_values, upt%eigen_vectors, temp_eivec, STAT = err)
     IF (err.NE.0) stop 'Deallocation error'
     !call write_eigenstates(upt)  
     
     !write(*,*) '(main) write mo file'
     
     !call writemofile(upt)
     !call write_cube(upt)

  end do

  call upt_mpi_end

  !write `eigensol` to file_num
90  do i = 1, nk*num_ev
    write(file_num,'(f8.6,a)', advance='NO') eigensol(i, 1)/Lk, ','
    do j = 2, size(eigensol, 2)
      write(file_num, '(f12.6,a)', advance='NO') eigensol(i, j), ','
    end do
    write(file_num, *)
    if (mod(i, nk) == 0) then
      write(file_num, *)
      write(file_num, *)
    endif
  end do
  close(file_num)

  deallocate(eigensol, STAT = err)
  IF (err.NE.0) stop 'Deallocation error'


contains

  subroutine AO_decomposition(eivec, coef2, n_atom, spin_flag, n_ao_type)

    complex(dp), dimension(:), intent(in) :: eivec
    LOGICAL, intent(in) :: spin_flag
    integer, intent(in) :: n_atom, n_ao_type
    integer :: atom, idx
    real(dp), dimension(3*n_atom), intent(out) :: coef2 !3 is the number of orbital types: s, p, d

    coef2 = 0.d0

    if (spin_flag) then !if there is spin degeneracy
      do atom = 1, n_atom
        idx = (atom-1)*10*2 !10 is number of orbitals per atom per spin, 2 spins
        coef2((atom-1)*n_ao_type + 1) = abs(eivec(idx+1))**2 + abs(eivec(idx+5))**2 +& 
                                abs(eivec(idx+11))**2 + abs(eivec(idx+15))**2
        coef2((atom-1)*n_ao_type + 2) = real(dot_product(eivec(idx+2:idx+4), eivec(idx+2:idx+4))) +&
                                real(dot_product(eivec(idx+12:idx+14), eivec(idx+12:idx+14)))
        coef2((atom-1)*n_ao_type + 3) = real(dot_product(eivec(idx+6:idx+10), eivec(idx+6:idx+10))) +&  
                                real(dot_product(eivec(idx+16:idx+20), eivec(idx+16:idx+20)))
      end do
    else !if there is no spin degeneracy
      do atom = 1, n_atom
        idx = (atom-1)*10 !10 is number of orbitals per atom per spin
        coef2((atom-1)*n_ao_type + 1) = abs(eivec(idx+1))**2 + abs(eivec(idx+5))**2
        coef2((atom-1)*n_ao_type + 2) = real(dot_product(eivec(idx+2:idx+4), eivec(idx+2:idx+4)))
        coef2((atom-1)*n_ao_type + 3) = real(dot_product(eivec(idx+6:idx+10), eivec(idx+6:idx+10)))
      end do
    end if

  end subroutine AO_decomposition

end program parameterizer
