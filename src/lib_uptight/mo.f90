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
program buildmo
  
  USE precision
  USE globals, only : LST
  USE upt_param, only : OUPT
  !USE struct_building, only : make_struct, init_basis, &
  !                            write_gen, read_gen, subs_dg_ions
  USE struct_building, only : init_structure, make_basis, init_basis, &
                              write_gen, read_gen, subs_dg_ions
  USE neighbours, only : make_neighbours_map, refine_neighbours_map, nearest_neighbours, &
                         check_input_nn_list, check_nn_map, write_neighbours_map, &
                         write_am_coo
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, set_max_order
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

  IMPLICIT NONE

  TYPE(OUPT), TARGET :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER :: i, j, n_atoms, ia, ib, k, offa, offb,  err
  TYPE(nearest_neighbours), POINTER :: current_near
  REAL(dp), dimension(:,:), allocatable :: ss
  REAL(dp) :: psum
  CHARACTER(2) :: calc
  CHARACTER(LST) :: filename
  logical :: exist

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

  !-----------------------------------------------------------
  ! Set globals
  !-----------------------------------------------------------
  !index_file = "index.dat"
  !ham_file = "ham.dat"
  !sparse_format = "full"
  !verbose = 10
  
  upt%structure%gen_filename = upt%gen_filename
  upt%relat = .true.
  upt%out_path = "./" 
  read(*,*) upt%scaling

  upt%d_onsite_shift_flag=.false.
  upt%potential_flag=.false.
  upt%potential_file="pot_on_atoms.dat"
  upt%syst_rotated=.false.   ! if .true. means the structure is already
                             ! rotated along the (111) direction

  read(*,*) upt%c_axis

  upt%k_point= (/0.0d0, 0.0d0, 0.0d0/)  
  upt%ioutput_flag=.false.

  upt%optmat=.false.
  upt%poldir = 3    

  upt%hybrid_passivation=.true.
  upt%d_H = 100.0d0 
  upt%E_H = -200.0d0

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
  call read_data(upt%nr_mat, upt%materials, upt%work_path, upt%database_path)

  !=================================================================================  
  ! can be || from here onward 
  write(*,*) '(main) make basis'
  call make_basis(upt%verbose, upt%structure, upt%basis)

  write(*,*) '(main) setup states and couplings:'
  call ref_states_and_couplings(upt%materials,upt%nr_mat, & 
                                    upt%ref_couplings, upt%n_ref_st, upt%n_ref_cpl)

  write(*,*) '(main) set max order in each material:'
  call set_max_order(upt%materials)

  write(*,*) '(main) init ions'
  do i=1,upt%nr_mat
     call init_mat_ion(upt%materials,i)
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

  !write(*,*) '(main) write gen'
  !call write_gen(upt%gen_out,upt%basis,upt%nn_map, upt%syst_rotated, .false.)


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


  calc = 'DF'
  read(*,*) calc
  read(*,*) upt%k_point
  
  read(*,*) upt%num_vb  
  read(*,*) upt%num_cb

  read(*,*) upt%lambda_vb 
  read(*,*) upt%lambda_cb 



  upt%ham%nrow = SUM( upt%basis%n_st(1:upt%basis%n_basis) ) * upt%n_spin 
  write(*,*) 'n_ham',upt%ham%nrow

  allocate( upt%eigen_values(upt%num_vb+upt%num_cb) )
  allocate( upt%eigen_vectors(upt%ham%nrow, upt%num_vb+upt%num_cb) ) 

  filename = "eigvec"
  write(*,*) '(main) read eigenvectors'
  call read_eigenstates(upt, filename, upt%num_vb, upt%num_cb)


  !write(*,*) '(main) write eigenvectors'
  !call write_eigenstates(upt)  


  !write(*,*) '(main) write cube file'
  !call write_cube(upt)

  write(*,*) 'OVERLAPS:'
  n_atoms = upt%basis%n_basis
  allocate( ss( upt%num_vb+upt%num_cb, upt%num_vb+upt%num_cb))
   
  do i = 1, upt%num_vb + upt%num_cb
    do j = 1, upt%num_vb + upt%num_cb

      psum = 0.0
      do ia = 1, n_atoms
        current_near => upt%nn_map(ia)
        offa = SUM( upt%basis%n_st(1:ia-1) ) * upt%n_spin 
        do k = 1, 5 
           ib = current_near%ind(k)  
           if (ib.gt.n_atoms .or. ib .eq. ia) cycle 
            
           offb = SUM( upt%basis%n_st(1:ib-1) ) * upt%n_spin 

           psum = psum + abs( dot_product(upt%eigen_vectors(offa+1:offa+20,i), upt%eigen_vectors(offb+1:offb+20,j)) )
              
        end do
      end do

      ss(i,j)= psum

   end do
 end do  

 do i = 1, upt%num_vb 
   do j = upt%num_vb+1, upt%num_vb + upt%num_cb
      write(*,*) i, j, ss(i,j)/sqrt(ss(i,i)*ss(j,j)) !dot_product(upt%eigen_vectors(:,i),upt%eigen_vectors(:,j))
   end do
 end do

end program buildmo
