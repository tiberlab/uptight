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
program test8
  USE clock
  USE precision
  USE mpi_globals
  USE globals, only : LST, SST
  USE upt_param, only : OUPT
  USE struct_building, only : init_structure, write_gen, read_gen, &
       subs_dg_ions, init_basis, make_basis
  USE neighbours, only : make_neighbours_map, refine_neighbours_map, &
                         check_input_nn_list, check_nn_map, write_neighbours_map
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   set_max_order, sort_states
  USE alloys, only : init_mat_ion
  USE TB_ham, only : sparse_ham, check_if_hermitian, check_if_real, hermitianize
  USE lanczos_driver, only : lanczos
  USE savemofile, only : writemofile, write_eigenstates, write_cube, write_jvxl, &
                         read_eigenvalues, read_eigenstates
  USE lapack_driver, only : lapack
  USE sparse_matrix, only : destroy_matrix,  write_sprs_to_coo, &
       write_csr_to_coo, analyze, write_csr
  USE uptight, only : upt_get_mat_el, upt_nullify_all, upt_destruct, upt_version
  USE jd_driver, only : jd

  IMPLICIT NONE

  TYPE(OUPT) :: upt

  INTEGER :: i, j, k, err, n_ham, num_ev, nev, nec
  CHARACTER(2) :: calc
  CHARACTER(SST) :: calc_str
  CHARACTER(LST) :: states,filename
  COMPLEX(dp) :: matel
  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: P
  LOGICAL :: exist

  !MPI STUFF
  integer :: ptype,num_cores,num_threads


  print*,'Hello!'

#ifdef __CUDA 
#ifdef UPT_MPI
    call setdevicebeforeinit()
#else
    call setdeviceinit()
#endif
#endif

  call upt_mpi_init(MPI_COMM_WORLD)
  
  upt%mpi_comm = MPI_COMM_WORLD
 
  OPEN (UNIT=12, FILE='dr')

  call upt_nullify_all(upt)

  call upt_version(upt)

  read(12,*) upt%database_path 
  read(12,*) upt%gen_filename
  read(12,*) upt%scaling
  read(12,*) upt%c_axis(:)  ! rotated along the (111) direction

  upt%work_path = "./"

  upt%gen_out = "out.gen"
  upt%state_file = "states.data"
  upt%sparse_format = "F" !WARNNG WITH F some TEST FAIL
  if (id0) then
     upt%verbose = 10
  else
     upt%verbose = 0
  endif

  upt%structure%gen_filename = upt%gen_filename
  upt%relat = .true.
  upt%out_path = "./"

  upt%d_onsite_shift_flag=.false.
  !if (upt%scaling) upt%d_onsite_shift_flag=.true.

  upt%potential_flag=.true.
  upt%potential_file="pot_on_atoms.dat"
  upt%syst_rotated=.false.   ! if .true. means the c-axis is already
 
  upt%k_point= 0.d0

  upt%ioutput_flag=.false.

  upt%optmat=.false.
  upt%poldir = 3    

  upt%hybrid_passivation=.true.
  upt%d_H = 100.0d0 
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


  call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat)

  if (id0) write(*,*) '(main) check input nn_list'
  call check_input_nn_list(upt%structure)


  if (id0) write(*,*) '(main) read material database:'
  do i = 1, upt%nr_mat
     call read_data(upt%materials(i), upt%work_path, upt%database_path)
  end do

  !==================================================================================

  if (id0) write(*,*) '(main) make basis'
  call make_basis(upt%verbose, upt%structure, upt%basis)

  if (id0) write(*,*) '(main) setup states and couplings:'
  call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, & 
       upt%ref_couplings, upt%n_ref_cpl)
  
  if (id0) write(*,*) '(main) sorting states according to reference:'
  call sort_states(upt%materials, upt%ref_states, upt%ref_couplings)

  if (id0) write(*,*) '(main) set max order in each material:',id
  call set_max_order(upt%materials)

  if (id0) write(*,*) '(main) init ions'
  do i=1,upt%nr_mat
     call init_mat_ion(upt%materials(i))
  enddo

  if (id0) write(*,*) '(main) materials info:'  
  call write_materials(upt%nr_mat,upt%materials)

  if (id0) write(*,*) '(main) init ion valences'
  call init_basis(upt%basis, upt%materials)

  if (id0) write(*,*) '(main) structure info:'
  call write_basis(upt%basis,0)

#ifdef UPT_MPI
  call mpi_barrier(upt_comm,err)
#endif
  !==================================================================================

  if (id0) write(*,*) '(main) Build Nearest Neighbour Table', id
  call refine_neighbours_map(upt%structure, upt%basis, upt%materials, upt%nn_map)


  if (id0) write(*,*) '(main) treat dangling bonds'
  call subs_dg_ions(upt%basis,upt%materials,upt%nn_map)


  !if (id0) write(*,*) '(main) write gen'
  !call write_gen(upt%gen_out,upt%basis,upt%nn_map, upt%syst_rotated,.true.)


  if (id0) write(*,*) '(main) init n_st'
  call init_n_st(upt%basis,upt%materials)


  if (upt%potential_flag) then
     allocate(upt%pot_data(upt%basis%n_basis), STAT = err)
     if (id0) write(*,*) "(main) pot_data size ",upt%basis%n_basis+upt%basis%n_dg_bond
     inquire(FILE=trim(upt%potential_file),EXIST=exist)  
     if (exist) then
        
        open(51,file=upt%potential_file)
        do i= 1, upt%basis%n_basis
           read(51,*) upt%pot_data(i)
        end do
        close(51)

        write(*,*) '(main) Read Potential File ',trim(upt%potential_file),' CPU',id
     else
        upt%pot_data = 0.d0     
     endif
  endif

#ifdef UPT_MPI
  call mpi_barrier(upt_comm,err)
#endif

  calc = 'DF'
  upt%solver_flag = 0 
  read(12,'(A20)') calc_str
  read(calc_str, '(A2,I5)') calc, upt%solver_flag
  if (id0) write(*,*) '(main) Solver: ',calc,' with flag:', upt%solver_flag  
  read(12,*) upt%k_point

  if (calc.ne.'LD') then

     INQUIRE(FILE="H.dat", EXIST=exist) 
     if (exist) then
        if (id0) write(*,*) "(main) Reading Hamiltonian from H.dat..."
        !call set_clock()  
        !call read_csr(upt%ham)
        !call write_clock()  
     else

        if(id0) write(*,*) '(main) Computing Hamiltonian, CPU', id
        if(id0) call set_clock()  
        call sparse_ham(upt)
        if(id0) call write_clock(">>> H done. Time for H: ")  
        if(id0) write(*,*)

#ifdef UPT_MPI
        call mpi_barrier(upt_comm,err)
#endif
     endif

  else
     upt%ham%nrow = 10*upt%n_spin*upt%basis%n_basis + upt%basis%n_dg_bond     
  end if !LD

  n_ham = upt%ham%nrow


  upt%start_vb = 1  
  upt%start_cb = 1 

  read(12,*) upt%num_vb  
  read(12,*) upt%num_cb

  read(12,*) upt%lambda_vb 
  read(12,*) upt%lambda_cb 

  upt%min_iter = 30     ! 2
  upt%long_iter = 70    ! 30
  upt%max_iter = 10000

  upt%fast_tol = 1.0d-1
  upt%long_tol = 1.0d-7 !
  upt%ort_tol  = 1.0d-5 ! 1d-5

  upt%seed_flag      = .false.
  upt%dynamic = .true.

  if (id0) write(*,*) "num valence:", upt%num_vb
  if (id0) write(*,*) "num conduction:", upt%num_cb
  
  num_ev = upt%num_vb + upt%num_cb 

  if (num_ev.eq.0) goto 200  

  !allocate(upt%eigen_values(num_ev), STAT = err)
  !allocate(upt%eigen_vectors(n_ham,num_ev), STAT = err)
  !IF (err.NE.0) stop 'Allocation error'
  !upt%eigen_values = 0.d0
  !upt%eigen_vectors = (0.d0,0.d0)


  select case ( calc ) 

  case ('LO')

     if (id0) write(*,*) '(main) lanczos diagonalization'
     read(12,*) upt%max_iter
     read(12,*) upt%long_tol
     call lanczos(upt) 

  case ('JD')

     if (id0) write(*,*) '(main) Jacobi-Davidson diagonalization'
     read(12,*) upt%max_iter
     read(12,*) upt%long_tol
     call jd(upt)

  case ('FT')

     if (id0) write(*,*) '(main) FEAST solver'
     !call feast(upt)

  case ('LC')

     write(*,*) '(main) lanczos cascade'
     read(12,*) upt%long_iter 
     read(12,*) upt%long_tol
     upt%max_iter  = 100000

     call lanczos(upt)
     if (id0) write(*,*) '(main) CASCADE HAS NOT BEEN IMPLEMENTED YET'
     stop

  case ('LK')

     if (id0) write(*,*) '(main) lapack diagonalization'
     call lapack(upt)

  case ('AR')

     if (id0) write(*,*) '(main) arpack diagonalization'
     !call arpack_dr(upt)

  case ('LD')

     if (id0) then
        write(*,*) '(main) load from file',num_ev,'eigenvalues'
     endif

     filename = 'eigv.dat'
     call read_eigenvalues(upt,filename,nev,nec)
     write(*,*) '(main) load from file',n_ham,'rows'
     filename = 'eigvec'
     call read_eigenstates(upt,filename,nev,nec)

  case ('LL')


     if (id0) write(*,*) '(main) load from file',num_ev,'eigenvalues'
     filename = 'eigv.dat'
     call read_eigenvalues(upt,filename,nev,nec)
     if (id0) write(*,*) '(main) load from file',n_ham,'rows'
     filename = 'eigvec'
     call read_eigenstates(upt,filename,nev,nec)
     if (id0) write(*,*) '(main) lanczos diagonalization'
     read(12,*) upt%max_iter
     call lanczos(upt)


  case ('DF')

     stop 'invalid calculation mode'

  end select


  upt%grid_step = 0.5d0

  if (id0) then
     call write_eigenstates(upt)  

     write(*,*) '(main) write cube files'
     call write_cube(upt)

  endif

goto 100

  if (id0) write(*,*) '(main) compute optical matrix elements'

  if(upt%num_cb.eq.0 .or. upt%num_vb.eq.0) goto 100

  upt%optmat=.true.

  allocate(P(upt%num_cb+5,upt%num_vb+5))

  do k= 1, 3

     if (id0) write(*,*) 'polarization',k
     upt%poldir = k

     call destroy_matrix(upt%ham)

     call sparse_ham(upt)

     do i=1, upt%num_cb
        do j= 1, upt%num_vb

           call upt_get_mat_el(upt,upt%num_vb+i,j,matel)

           if (id0) write(*,'(a2,i2,i2,a2,f15.8)') 'M(',i,j,')=',ABS(matel)*ABS(matel)
           P(i,j) = ABS(matel)*ABS(matel);

        end do
     end do

     if (id0) write(*,*) 'C-HH:',P(1,1)+P(2,1)+P(1,2)+P(2,2)
     if (id0) write(*,*) 'C-LH:',P(1,3)+P(2,3)+P(1,4)+P(2,4)
     if (id0) write(*,*) 'C-SO:',P(1,5)+P(2,5)+P(1,6)+P(2,6)     


  end do

  deallocate(P)

100 continue

  deallocate(upt%eigen_values, upt%eigen_vectors)

200 continue

  if (id0) write(*,*) 'freeing memory...'

  call upt_destruct(upt)

300 continue


#ifdef UPT_MPI
  write(*,*) "(main) Finishing UPT PARALLEL Version@rank0", id
  call upt_mpi_end()
#endif

end program test8
