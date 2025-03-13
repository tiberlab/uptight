! This file is part of uptight.
!
! uptight is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! uptight is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with uptight. If not, see <https://www.gnu.org/licenses/>.
!
program test7

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

  INTEGER :: i, j, k, file_num, file_coef2, filepdos, n_ham, num_ev, err, at, bg, st, ao
  CHARACTER(2) :: calc
  COMPLEX(dp) :: matel
  CHARACTER(LST) :: kfile, file_name, file_coef2_name
  REAL(dp), DIMENSION(3) :: k_xyz, k_old, strain
  REAL(dp) :: Lk, weight, pop, norm
  CHARACTER(2) :: id
  real(dp), dimension(16) :: coef2
  complex(dp), dimension(:,:), allocatable :: temp_eivec

  integer :: comm
  integer :: upt_group

  call upt_mpi_init(0)

  read(*,*) upt%database_path ! line 1 of `dr_bands`: database path
  upt%work_path = "./"
  read(*,*) upt%gen_filename  ! line 2 of `dr_bands`: name of `.upg` file
  upt%gen_out = "out.gen"
  upt%state_file = "states.data"  
  upt%sparse_format = "U"
  upt%verbose = 10
  
  upt%structure%gen_filename = upt%gen_filename
  read(*,*) upt%relat  ! line 3 of `dr_bands`: relativistic flag
  upt%out_path = "./"
  read(*,*) upt%scaling  ! line 4 of `dr_bands`: harrison scaling

  upt%d_onsite_shift_flag=.true.
  upt%potential_flag=.false.
  upt%syst_rotated=.false.   ! if .true. means the structure is already
                             ! rotated along the (111) direction

  !upt%c_axis= (/ 0.d0, 0.d0, 1.d0 /)
  !upt%c_axis= (/ 1.d0, 0.d0, 0.d0 /)
  read(*,*) upt%c_axis(:)  ! line 5 of `dr_bands`: c axis (not important for Tan but important for Jancu wurtzite and 2D MoS2)

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
<<<<<<< .working
  !call sort_states(upt%materials, upt%ref_states, upt%ref_couplings)
  do i=1,size(upt%materials)
    call sort_states(upt%materials(i), upt%ref_states, upt%ref_couplings)
    do j=1,size(upt%materials)
      if (associated(upt%interfaces(i,j)%nr_parents)) then
        call sort_states(upt%interfaces(i,j), upt%ref_states, upt%ref_couplings)
      end if
    end do
  end do

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
  read(*,*) calc  ! line 6 of `dr_bands`: eigensolver (LK = lapack)
  
  upt%start_vb = 1  
  upt%start_cb = 1 

  read(*,*) upt%num_vb   ! line 7 of `dr_bands`: number of valence bands (actually, number of bands below valence guess)
  read(*,*) upt%num_cb   ! line 8 of `dr_bands`: number of conduction bands (actually, number of bands above conduction guess)

  read(*,*) upt%lambda_vb  ! line 9 of `dr_bands`: valence guess
  read(*,*) upt%lambda_cb  ! line 10 of `dr_bands`: conduction guess

  upt%min_iter = 2
  upt%long_iter = 30
  upt%max_iter = 100000

  upt%fast_tol = 1.0d-1
  upt%long_tol = 1.0d-10
  upt%ort_tol  = 1.0d-5
     
  upt%solver_flag = 0
  upt%dynamic        = .true. 
  upt%seed_flag      = .false.

  read(*,*) kfile   ! line 11 of `dr_bands`: name of file containing list of k-points in fractions coordinates


  read(*,*) strain(1:3)  ! line 12 of `dr_bands`: strain

  upt%basis%strain(1,:) = strain(1)
  upt%basis%strain(2,:) = strain(2)
  upt%basis%strain(3,:) = strain(3)

  open(15,file=trim(kfile), STATUS='OLD')

  file_name = trim(upt%out_path)//'eigv.dat'
  file_coef2_name = trim(upt%out_path)//'coef2.dat'
  
  CALL open_file( file_name, file_num, operation = "write", &
       format_flag = .TRUE., replace_flag = .TRUE. )

  close(file_num)

  CALL open_file( file_coef2_name, file_coef2, operation = "write", &
       format_flag = .TRUE., replace_flag = .TRUE. )

  close(file_coef2)

  CALL open_file( file_name, file_num, operation = "write", &
       format_flag = .TRUE., replace_flag = .FALSE. )

  CALL open_file( file_coef2_name, file_coef2, operation = "write", &
       format_flag = .TRUE., replace_flag = .FALSE. )

  write( file_num, '(a4)', advance='NO' ) '# N  '
  do i =1, upt%num_vb
    write(id,'(i2.2)') i
    write( file_num, '(a4,a2)', advance='NO' ) '  Ev',id
  enddo
  do i =1, upt%num_cb
    write(id,'(i2.2)') i
    write( file_num, '(a4,a2)', advance='NO' ) '  Ec',id
  enddo
  write( file_num, *)

  write( file_coef2, '(a60)', advance='NO' ) '# k-path, eigenvalues, coefficients(s,p,s*,d|Ga1,Ga2,N3,N4)'
  write( file_coef2, *)

  Lk = 0.d0
  k_xyz = 0.d0
  k_old = 0.d0

  upt%verbose=0

  !filepdos = 400
  !OPEN(filepdos, file='pdos.dat')

  do k = 1, 100000

     read(15,fmt=*,end=100) upt%k_point!, weight
 
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
     num_ev = upt%num_vb + upt%num_cb 

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

     k_old = k_xyz

     allocate(temp_eivec(num_ev, n_ham), STAT = err)
     IF (err.NE.0) stop 'Allocation error'
     temp_eivec = transpose(upt%eigen_vectors)
     write( file_num,'(ES19.8)', advance='NO' ) Lk

     do i = 1, upt%num_vb + upt%num_cb 
        
        WRITE( file_num, '(( x1 f13.8 ))', advance='NO' ) &
             upt%eigen_values(i)      

        call AO_decomposition(temp_eivec(i,:), coef2)
        write( file_coef2,'(ES15.8)', advance='NO' ) Lk
        write( file_coef2,'(( x1 f13.8 ))', advance='NO' ) upt%eigen_values(i)
        do ao = 1, 16 !16 types of orbitals: s, p, s*, d for 4 atoms: Ga1, Ga2, N3, N4
          write( file_coef2, '(( x1 f8.5 ))', advance='NO' ) coef2(ao)
        end do
        write(file_coef2,*)
        
     enddo

     write(file_coef2,*)
     write(file_num,*) 

     deallocate(upt%eigen_values, upt%eigen_vectors, temp_eivec, STAT = err)
     IF (err.NE.0) stop 'Deallocation error'
     !call write_eigenstates(upt)  
     
     !write(*,*) '(main) write mo file'
     
     !call writemofile(upt)
     !call write_cube(upt)

  end do

  call upt_mpi_end

100 CLOSE(file_num)
close(file_coef2)
!CLOSE(filepdos)


contains

  subroutine pdos(nfile,upt,at,k,weight)
    integer, intent(in) :: nfile, at
    type(OUPT) :: upt
    integer :: k
    real(dp) :: weight
     
    integer bg, st
    real(dp) :: pop 
    !print*,'write PDOS'   

     bg= 20*(at-1) +1
     st= 20*(at-1) +20
     write(nfile,*) 'KPT', k, 'SPIN', 1, 'KWEIGHT', weight 
     do i = upt%num_vb, 1, -1
        pop=real(dot_product(upt%eigen_vectors(bg:st,i), &
                            upt%eigen_vectors(bg:st,i)))
        write(nfile,*) upt%eigen_values(i),  pop
     enddo 
     do i = upt%num_vb+1, upt%num_vb+upt%num_cb 
        pop=real(dot_product(upt%eigen_vectors(bg:st,i), &
                            upt%eigen_vectors(bg:st,i)))
        write(nfile,*) upt%eigen_values(i),  pop
     enddo 

  end subroutine pdos


  subroutine AO_decomposition(eivec, coef2)

    complex(dp), dimension(:), intent(in) :: eivec
    real(dp), dimension(16), intent(out) :: coef2

    coef2 = 0.d0
    !Ga1: s, p, s*, d
    coef2(1) = abs(eivec(1))**2 + abs(eivec(11))**2
    coef2(2) = abs(eivec(2))**2 + abs(eivec(3))**2 + abs(eivec(4))**2 +&
               abs(eivec(12))**2 + abs(eivec(13))**2 + abs(eivec(14))**2
    coef2(3) = abs(eivec(5))**2 + abs(eivec(15))**2
    coef2(4) = abs(eivec(6))**2 + abs(eivec(7))**2 + abs(eivec(8))**2 + abs(eivec(9))**2 + abs(eivec(10))**2 +&
               abs(eivec(16))**2 + abs(eivec(17))**2 + abs(eivec(18))**2 + abs(eivec(19))**2 + abs(eivec(20))**2

    !Ga2: s, p, s*, d
    coef2(5) = abs(eivec(21))**2 + abs(eivec(31))**2
    coef2(6) = abs(eivec(22))**2 + abs(eivec(23))**2 + abs(eivec(24))**2 +&
               abs(eivec(32))**2 + abs(eivec(33))**2 + abs(eivec(34))**2
    coef2(7) = abs(eivec(25))**2 + abs(eivec(35))**2
    coef2(8) = abs(eivec(26))**2 + abs(eivec(27))**2 + abs(eivec(28))**2 + abs(eivec(29))**2 + abs(eivec(30))**2 +&
               abs(eivec(36))**2 + abs(eivec(37))**2 + abs(eivec(38))**2 + abs(eivec(39))**2 + abs(eivec(40))**2

    !N3: s, p, s*, d
    coef2(9) = abs(eivec(41))**2 + abs(eivec(51))**2
    coef2(10) = abs(eivec(42))**2 + abs(eivec(43))**2 + abs(eivec(44))**2 +&
                abs(eivec(52))**2 + abs(eivec(53))**2 + abs(eivec(54))**2
    coef2(11) = abs(eivec(45))**2 + abs(eivec(55))**2
    coef2(12) = abs(eivec(46))**2 + abs(eivec(47))**2 + abs(eivec(48))**2 + abs(eivec(49))**2 + abs(eivec(50))**2 +&
                abs(eivec(56))**2 + abs(eivec(57))**2 + abs(eivec(58))**2 + abs(eivec(59))**2 + abs(eivec(60))**2

    !N4: s, p, s*, d
    coef2(13) = abs(eivec(61))**2 + abs(eivec(71))**2
    coef2(14) = abs(eivec(62))**2 + abs(eivec(63))**2 + abs(eivec(64))**2 +&
                abs(eivec(72))**2 + abs(eivec(73))**2 + abs(eivec(74))**2
    coef2(15) = abs(eivec(65))**2 + abs(eivec(75))**2
    coef2(16) = abs(eivec(66))**2 + abs(eivec(67))**2 + abs(eivec(68))**2 + abs(eivec(69))**2 + abs(eivec(70))**2 +&
                abs(eivec(76))**2 + abs(eivec(77))**2 + abs(eivec(78))**2 + abs(eivec(79))**2 + abs(eivec(80))**2

  end subroutine AO_decomposition

end program test7
