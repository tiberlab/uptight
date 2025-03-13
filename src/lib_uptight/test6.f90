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
program test6

  USE precision
  USE globals, only : LST
  USE upt_param, only : OUPT
  USE struct_building, only : init_structure, make_basis, init_basis, &
                              write_gen, read_gen, subs_dg_ions
  USE neighbours, only : make_neighbours_map, refine_neighbours_map, &
                         check_input_nn_list, check_nn_map, write_neighbours_map, &
                         write_am_coo
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   sort_states, set_max_order
  USE alloys, only : init_mat_ion
  USE TB_ham, only : sparse_ham, check_if_hermitian, check_if_real, hermitianize
  USE lanczos_driver, only : lanczos
  !USE feast_driver, only : feast
  USE savemofile, only : writemofile, write_eigenstates, write_cube, write_jvxl, &
                         read_eigenvalues, read_eigenstates
  USE lapack_driver, only : lapack
  USE sparse_matrix, only : destroy_matrix,  write_sprs_to_coo, &
                            write_csr_to_coo, write_csr, analyze
  !USE arpack, only : arpack_dr
  USE uptight, only : upt_get_mat_el, upt_nullify_all, upt_destruct, upt_version
  USE clock

  IMPLICIT NONE

  TYPE(OUPT), TARGET :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER :: i, j, k, err, n_ham, num_ev, nev, nec
  CHARACTER(2) :: calc
  CHARACTER(LST) :: filename
  COMPLEX(dp) :: matel
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: P
  LOGICAL :: exist

  pupt => upt
  call upt_nullify_all(pupt)

  call upt_version(pupt)

  read(*,*) upt%database_path 
  upt%work_path = "./"
  read(*,*) upt%gen_filename
  upt%gen_out = "out.gen"
  upt%state_file = "states.data"
  upt%sparse_format = "F" !WARNNG WITH F some TEST FAIL
  upt%verbose = 10
  
  upt%structure%gen_filename = upt%gen_filename
  upt%relat = .true.
  upt%out_path = "./"
  read(*,*) upt%scaling

  upt%d_onsite_shift_flag=.false.
  !if (upt%scaling) upt%d_onsite_shift_flag=.true.

  upt%potential_flag=.true.
  upt%potential_file="pot_on_atoms.dat"
  upt%syst_rotated=.false.   ! if .true. means the c-axis is already
                             ! rotated along the (111) direction

  !upt%c_axis= (/ 0.d0, 0.d0, 1.d0 /)
  !upt%c_axis= (/ 1.d0, 0.d0, 0.d0 /)
  read(*,*) upt%c_axis(:)

  upt%k_point= (/0.0d0, 0.0d0, 0.0d0/)  

  upt%ioutput_flag=.false.

  upt%optmat=.false.
  upt%poldir = 3    

  upt%hybrid_passivation=.true.
  upt%d_H = 100.d0  
  upt%E_H = -200.d0 

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
  call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat)
  !call make_struct(upt%verbose, upt%structure, upt%basis, upt%materials, upt%nr_mat)

  write(*,*) '(main) check input nn_list'
  call check_input_nn_list(upt%structure)

  write(*,*) '(main) read material database:'
  do i = 1, upt%nr_mat
     call read_data(upt%materials(i), upt%work_path, upt%database_path)
  end do

  

  !=================================================================================  
  ! can be || from here onward 
  write(*,*) '(main) make basis'
  call make_basis(upt%verbose, upt%structure, upt%basis)

  write(*,*) '(main) setup states and couplings:'
  call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, & 
       upt%ref_couplings, upt%n_ref_cpl)
  
  if(upt%verbose.gt.0) write(*,*) '(main) sorting states according to reference:'
  call sort_states(upt%materials, upt%ref_states, upt%ref_couplings)
                                    
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

  !write(*,*) '(main) write neighbour list:'
  !call write_am_coo(upt%basis,upt%nn_map)
  !call write_neighbour_list(upt%nn_list)
  
  !write(*,*) '(main) check neighbour list'
  !call check_nn_list(upt%nn_list)

  write(*,*) '(main) treat dangling bonds'
  call subs_dg_ions(upt%basis,upt%materials,upt%nn_map)

  write(*,*) '(main) write gen'
  call write_gen(upt%gen_out,upt%basis,upt%nn_map, upt%syst_rotated, .false.)

  write(*,*) '(main) init n_st'
  call init_n_st(upt%basis,upt%materials)


  if (upt%potential_flag) then
     allocate(upt%pot_data(upt%basis%n_basis), STAT = err)
     write(*,*) "(main) pot_data size ", upt%basis%n_basis
     inquire(FILE=upt%potential_file,EXIST=exist)  
     if (exist) then
       write(*,*) '(main) Read Potential File: ',trim(upt%potential_file)
       open(51,file=upt%potential_file)

       do i= 1, upt%basis%n_basis
         read(51,*) upt%pot_data(i)
       end do

       close(51)
     else
        upt%pot_data = 0.d0     
     endif
  endif      

  write(*,*) '(main) read input...'
  calc = 'DF'
  read(*,*) calc
  read(*,*) upt%k_point

  !if (calc.eq.'LK') then 
  !endif     

  if (calc.ne.'LD') then
     write(*,*) '(main) compute H'
     call sparse_ham(upt)

     !write(*,*) '(main) check if real H'
     !call check_if_real(upt%ham)

     if (upt%sparse_format.eq.'F') then
        !write(*,*) '(main) make full matrix Hermitian'
        !call hermitianize(upt%ham)
        
        write(*,*) '(main) check if Hermitian'
        call check_if_hermitian(upt%ham)
        
     end if
     !call analyze(upt%ham)
     !write(*,*) 'saving H_COO.DAT...'
     !call  write_csr_to_coo(upt%ham, upt%work_path)

  else
    upt%ham%nrow = 10*upt%n_spin*upt%basis%n_basis + upt%basis%n_dg_bond     
  end if


  
  upt%start_vb = 1  
  upt%start_cb = 1 

  read(*,*) upt%num_vb  
  read(*,*) upt%num_cb

  read(*,*) upt%lambda_vb 
  read(*,*) upt%lambda_cb 

  upt%min_iter = 30     ! 2 
  upt%long_iter = 32    ! 30
  upt%max_iter = 100000

  upt%fast_tol = 1.0d-1
  upt%long_tol = 1.0d-9 !
  upt%ort_tol  = 1.0d-6 ! 1d-5
     
  upt%solver_flag = 0
  upt%dynamic        = .true.
  upt%seed_flag      = .false.

  n_ham = upt%ham%nrow
  num_ev = upt%num_vb + upt%num_cb 
  
  if(num_ev.eq.0) goto 200  

  allocate(upt%eigen_values(num_ev), STAT = err)
  allocate(upt%eigen_vectors(n_ham,num_ev), STAT = err)
  IF (err.NE.0) stop 'Allocation error'
  upt%eigen_values = 0.d0
  upt%eigen_vectors = (0.d0,0.d0)

  select case ( calc ) 
  case ('LO')
     write(*,*) '(main) lanczos diagonalization'
     read(*,*) upt%max_iter
     read(*,*) upt%long_tol
     call set_clock()
     call lanczos(upt)
     call write_clock()
  case ('LC')
     write(*,*) '(main) lanczos with cuda'
     upt%solver_flag = 1
     read(*,*) upt%max_iter 
     read(*,*) upt%long_tol
     call lanczos(upt)
  case ('LR')
     write(*,*) '(main) lanczos with cuda - split'
     upt%solver_flag = 2
     read(*,*) upt%max_iter 
     read(*,*) upt%long_tol
     call lanczos(upt)
  case ('LK')
     write(*,*) '(main) lapack diagonalization'
     call lapack(upt)
  case ('FT')
     write(*,*) '(main) FEAST solver'
     !call feast(upt)
  case ('AR')
     write(*,*) '(main) arpack diagonalization'
     !call arpack_dr(upt)
  case ('LD')
     call read_eigenvalues(upt,filename,nev,nec)
     write(*,*) '(main) load from file',n_ham,'rows'
     filename = 'states.upt'
     call read_eigenstates(pupt,filename,nev,nec)
  case ('LL')
     call read_eigenvalues(upt,filename,nev,nec)
     write(*,*) '(main) load from file',n_ham,'rows'
     filename = 'states.upt'
     call read_eigenstates(pupt,filename,nev,nec)
     write(*,*) '(main) lanczos diagonalization'
     read(*,*) upt%max_iter
     call lanczos(upt)
  case ('DF')
     stop 'invalid calculation mode'
  end select
  
  upt%grid_step = 0.5d0

  call write_eigenstates(upt)  

  write(*,*) '(main) write cube files'
  call write_cube(upt)

  write(*,*) '(main) write jvxl files'
  call write_jvxl(upt)


  write(*,*) '(main) compute optical matrix elements'

  if(upt%num_cb.eq.0 .or. upt%num_vb.eq.0) goto 100

  upt%optmat=.true.
  upt%sparse_format = "F"

  allocate(P(upt%num_cb+5,upt%num_vb+5))

  do k= 1, 3

     write(*,*) 'polarization',k
     upt%poldir = k

     call destroy_matrix(upt%ham)
     call sparse_ham(upt)

     !call hermitianize(upt%ham)
     !call check_if_hermitian(upt%ham)

     do i=1, upt%num_cb
        do j= 1, upt%num_vb
           
           call upt_get_mat_el(pupt,upt%num_vb+i,j,matel)

           write(*,'(a2,i2,i2,a2,f15.8)') 'M(',i,j,')=',ABS(matel)*ABS(matel)
           P(i,j) = ABS(matel)*ABS(matel);

        end do
     end do

     write(*,*) 'C-HH:',P(1,1)+P(2,1)+P(1,2)+P(2,2)
     write(*,*) 'C-LH:',P(1,3)+P(2,3)+P(1,4)+P(2,4)
     write(*,*) 'C-SO:',P(1,5)+P(2,5)+P(1,6)+P(2,6)     


  end do
 
  deallocate(P)

100 continue

  deallocate(upt%eigen_values, upt%eigen_vectors)

200 continue

   write(*,*) 'freeing memory...'

   call upt_destruct(pupt)

end program test6
