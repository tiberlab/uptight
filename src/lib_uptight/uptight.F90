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
!!* This module contains the high level routines for creating and diagonalising
!!* a Uptight Hamiltonian.
module UPTIGHT

  USE precision
  USE globals, only : LST
  USE mpi_globals
  USE upt_param, only : OUPT, set_defaults
  USE neighbours, only : make_neighbours_map, destroy_neighbours_map, refine_neighbours_map, &
                         check_input_nn_list, check_nn_map, write_neighbours_map
  USE struct_building, only : init_structure, make_basis, init_basis, &
                              write_gen, read_gen, subs_dg_ions
  USE input_data, only : read_data
  USE type_defs, only : write_basis, write_materials, destroy_basis, destroy_material, &
                        destroy_structure, nullify_basis, nullify_structure
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, set_max_order, &
                                   sort_states
  USE alloys, only : init_mat_ion
  USE TB_ham, only : sparse_ham, hermitianize, check_if_hermitian
  USE coarse_grain, only : cg_prepare, cg_clear, cg_configure
  USE lanczos_driver, only : lanczos
  USE lapack_driver, only : lapack
  USE jd_driver, only : jd
  !USE arpack, only : arpack_dr
  USE savemofile, only : writemofile, write_eigenstates, write_cube, write_jvxl
  USE sparse_numrec, only : sprs_ax
  USE sparse_matrix, only : destroy_matrix, write_csr_to_coo
  !USE feast_driver, only : feast
  USE errors
  USE omp_lib 

  implicit none
  private

  PUBLIC :: upt_init, upt_destruct, upt_hamiltonian, upt_lanczos, upt_feast
  PUBLIC :: upt_jd, upt_lapack
  public :: upt_write_eigenvectors, upt_get_mat_el, upt_nullify_all
  public :: upt_set_defaults, upt_get_hamil, upt_version, upt_alloc_eigv
  public :: upt_configure_coarse_graining, upt_get_coarse_graining_info
  public :: upt_writehamiltonian

  !interface init
  !   module procedure UPT_init
  !end interface
  
  !interface destruct
  !   module procedure UPT_destruct
  !end interface
  
contains
  
  subroutine UPT_nullify_all(upt)

    TYPE(OUPT) :: upt  

    NULLIFY(upt%materials)
    NULLIFY(upt%pot_data)
    NULLIFY(upt%eigen_values)
    NULLIFY(upt%eigen_vectors)
    NULLIFY(upt%particles)
    NULLIFY(upt%nn_map)
    NULLIFY(upt%ham%M)
    NULLIFY(upt%ham%Mj)
    NULLIFY(upt%ham%Mi)
    NULLIFY(upt%cg_ham%M)
    NULLIFY(upt%cg_ham%Mj)
    NULLIFY(upt%cg_ham%Mi)
    NULLIFY(upt%cg_U%M)
    NULLIFY(upt%cg_U%Mj)
    NULLIFY(upt%cg_U%Mi)
    NULLIFY(upt%ref_states)
    NULLIFY(upt%ref_couplings)

    call nullify_basis(upt%basis)

    call nullify_structure(upt%structure)

    ! initialize n_basis for following checks
    upt%basis%n_basis = 0

  end subroutine UPT_nullify_all
  !---------------------------------------------------------------------

  subroutine UPT_set_defaults(upt)
    TYPE(OUPT) :: upt    

    call set_defaults(upt)

  end subroutine UPT_set_defaults 

  !-------------------------------------------------------------------------------
  ! FUNDAMENTAL INITIALIZING ROUTINE
  ! 
  ! 1. upt%n_spin 
  ! 2. call init_structure() : Read .upg file and init upt%str, upt%materials
  ! 3. call make_basis()          : Init upt%basis
  ! 4. call check_input_nn_list() : Check symmetry of nn_list
  ! 5. call read_data()           : Read material_data() from databases 
  ! 6. call ref_states_and_couplings() : setup states and coupl lists 
  ! 7. call set_max_order()       : Set mat_data()%max_order for each material
  ! 8. call init_mat_ion()        : Initialize ions in each alloy material
  ! 9. call init_basis()          : Set ion valences
  !10. call refine_neighbour_map() : Build real nn_map depending on TB order
  !11. call check_nn_map()       : check again bondmap consistency
  !12. call subs_dg_ions()        : substitute H with missing ion (for binaries)
  !13. call init_n_st()           : Initialize upt%n_st(i) (s. states per ion)
  !-------------------------------------------------------------------------------
  subroutine UPT_init(upt)
    
    TYPE(OUPT) :: upt    
    
    integer :: i, j, num_threads

    call upt_mpi_init(upt%mpi_comm)

    if (.not.id0) upt%verbose = 0

    if (upt%verbose.gt.0) THEN
       write(*,*) 
       write(*,*) 'MPI parallelized on ',num_procs,'processes'
#ifdef __OMP
       num_threads = omp_get_max_threads()/num_procs
       call omp_set_num_threads(num_threads)
       write(*,*) 'OpenMP parallelized on',omp_get_max_threads(),'threads'
#endif
       ! !$OMP PARALLEL
       !  write(*,*) 'Thread ',omp_get_thread_num(),'alive'
       ! !$OMP END PARALLEL
       write(*,*)
    end if

    if(upt%relat) then 
       upt%n_spin = 2
       if(upt%verbose.gt.0) write(*,*) "(init upt) relativistic calculation"
    else
       upt%n_spin = 1
       if(upt%verbose.gt.0) write(*,*) "(init upt) non relativistic calculation"
    end if

    call set_machine_acc
    
    upt%structure%gen_filename = upt%gen_filename

    if(upt%verbose.gt.0) write(*,*) '(init upt) init_structure'
    call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat, &
                                                    upt%interfaces, upt%nr_int)
    
    if(upt%verbose.gt.0) write(*,*) '(init upt) check input nn_map'
    call check_input_nn_list(upt%structure)

    if(upt%verbose.gt.0) write(*,*) '(init upt) read material database:'
    do i = 1, size(upt%materials)
       call read_data(upt%materials(i), upt%work_path, upt%database_path)
    end do

    if(upt%verbose.gt.0) write(*,*) '(init upt) read interface database:'
    do i = 1, upt%nr_mat
      do j = 1, upt%nr_mat
         if (associated(upt%interfaces(i,j)%nr_parents)) then
           call read_data(upt%interfaces(i,j), upt%work_path, upt%database_path)
         end if
      end do
    end do

    !================================================================================

    if(upt%verbose.gt.0) write(*,*) '(init upt) init_basis'
    call make_basis(upt%verbose, upt%structure, upt%basis)

    
    if(upt%verbose.gt.0) write(*,*) '(main) setup states and couplings:'
    call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, & 
                                    upt%ref_couplings, upt%n_ref_cpl)

    if(upt%verbose.gt.0) write(*,*) '(main) sorting states according to reference:'
    !call sort_states(upt%materials, upt%ref_states, upt%ref_couplings)
    do i=1,size(upt%materials)
      call sort_states(upt%materials(i), upt%ref_states, upt%ref_couplings)
      do j=1,size(upt%materials)
        if (associated(upt%interfaces(i,j)%nr_parents)) then
          call sort_states(upt%interfaces(i,j), upt%ref_states, upt%ref_couplings)
        end if
      end do
    end do
                                    

    if(upt%verbose.gt.0) write(*,*) '(init upt) set max order in each material'
    call set_max_order(upt%materials)
    ! do i=1,size(upt%materials)
    !   call set_max_order(upt%interfaces(i,:))
    ! enddo

    if(upt%verbose.gt.0) write(*,*) '(init upt) init ions in each material'
    do i=1,size(upt%materials)
      if (size(upt%materials(i)%nr_parents) .gt. 0) then
       call init_mat_ion(upt%materials(i))
      end if
      do j=1,size(upt%materials)
        if (associated(upt%interfaces(i,j)%nr_parents)) then
          if (size(upt%interfaces(i,j)%nr_parents) .gt. 0) then
            call init_mat_ion(upt%interfaces(i,j))
          end if
        end if
      enddo
    enddo

    
    if(upt%verbose.gt.2) then 
       write(*,*) '(init upt) materials info:'  
       call write_materials(upt%nr_mat,upt%materials)
    endif

    if(upt%verbose.gt.0) write(*,*) '(init upt) init ion valences'
    call init_basis(upt%basis, upt%materials)

    if(upt%verbose.gt.2) then
        write(*,*) '(init upt) structure info:'
        call write_basis(upt%basis,0)
    endif

    if(upt%verbose.gt.2) write(*,*) '(init upt) Build Nearest Neighbour Table'
    call refine_neighbours_map(upt%structure, upt%basis, upt%materials, upt%nn_map)

    if(upt%check_bondmap) then
       if(upt%verbose.gt.2) write(*,*) '(init upt) check neighbour list'
       call check_nn_map(upt%nn_map)
    end if
  
    if(upt%verbose.gt.2) write(*,*) '(init upt) treat dangling bonds'
    call subs_dg_ions(upt%basis,upt%materials,upt%nn_map)
    

    if(upt%verbose.gt.2) write(*,*) '(init upt) init n_st'
    call init_n_st(upt%basis,upt%materials)
    
        
  end subroutine UPT_init

  !---------------------------------------------------------------------

  subroutine UPT_destruct(upt)
  
    TYPE(OUPT) :: upt     
    INTEGER :: i

    call destroy_matrix(upt%ham)
    call cg_clear(upt)

    !write(*,*) '(debug) deallocate basis'
    call destroy_basis(upt%basis)

    !write(*,*) '(debug) deallocate structure'
    call destroy_structure(upt%structure)

    !write(*,*) '(debug) deallocate ref_couplings'  
    if ( associated(upt%ref_states) ) deallocate(upt%ref_states)
    if ( associated(upt%ref_couplings) ) deallocate(upt%ref_couplings)

    !write(*,*) '(debug) deallocate potential'  
    if ( associated(upt%pot_data) ) deallocate(upt%pot_data)
    
    !write(*,*) '(debug) deallocate eigvals'        
    if ( associated(upt%eigen_values) ) deallocate(upt%eigen_values)
 
    !write(*,*) '(debug) deallocate eigvects'     
    if ( associated(upt%eigen_vectors) ) deallocate(upt%eigen_vectors)

    !write(*,*) '(debug) deallocate particles'     
    if ( associated(upt%particles) ) deallocate(upt%particles)
        
    if(associated(upt%nn_map)) then
       !write(*,*) '(debug) deallocate nn_map'     
       call destroy_neighbours_map(upt%nn_map)
    end if

    if(associated(upt%materials)) then
       DO i = 1, upt%nr_mat
          !write(*,*) '(debug) destroy material',i,'/',SIZE(upt%materials)
          call destroy_material( upt%materials(i) ) 
       END DO
       DEALLOCATE(upt%materials)
    end if

    call kill_shifts()

  end subroutine UPT_destruct

  !---------------------------------------------------------------------
  subroutine UPT_version(upt)  
   
    TYPE(OUPT) :: upt
    character(3), parameter :: UPTVER= __REVISION
    character(10),parameter :: DATE= __COMPDATE
    character(3),parameter :: MODIF= __MODIFIED

    write(*,'(a19,a3,a2,2x,a10)') 'UPTIGHT version: ',TRIM(UPTVER), &
                                           TRIM(MODIF), TRIM(DATE)
    write(*,*)

  end subroutine UPT_version
  !---------------------------------------------------------------------
  subroutine UPT_hamiltonian(upt)

    TYPE(OUPT), pointer :: upt 
    integer :: ierr

    if(upt%verbose.gt.0) write(*,*) '(uptight) clear any prev. matrix'    
    call destroy_matrix(upt%ham)

    if(upt%verbose.gt.0) write(*,*) '(uptight) compute new matrix'
    call sparse_ham(upt)
    if (upt%cg_enabled) then
       call cg_prepare(upt, ierr)
       if (ierr /= 0) then
          write(*,*) '(uptight) coarse-graining preparation failed'
          stop 1
       end if
    end if

    ! VERY IMPORTANT NOTE: The full H is non-hermitian when there are 
    ! interfaces Alloy/Alloy or Alloy/Pure since there is no symmetry 
    ! between the interaction parameters a-b/b-a because they are seen 
    ! as pure in one case and alloy in the other. 
    ! This problem has been cured in check_interface_ion() therefore
    ! hermitianize() has been commented out.

    ! if (upt%sparse_format.eq.'F') then
    !   if(upt%verbose.gt.0) write(*,*) '(uptight) make full matrix Hermitian'
    !   call hermitianize(upt%ham)

    !   if(upt%verbose.gt.0) write(*,*) '(uptight) check if Hermitian'
    !   call check_if_hermitian(upt%ham)
    !end if

  end subroutine UPT_hamiltonian

  subroutine UPT_configure_coarse_graining(upt, enabled, nblocks, emin, emax, imbalance)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled
    integer, intent(in) :: nblocks
    real(dp), intent(in) :: emin, emax, imbalance
    call cg_configure(upt, enabled, nblocks, emin, emax, imbalance)
  end subroutine UPT_configure_coarse_graining

  subroutine UPT_get_coarse_graining_info(upt, ready, original_dim, reduced_dim, nblocks, cut_fraction)
    use coarse_grain, only : cg_get_info
    type(OUPT), intent(in) :: upt
    logical, intent(out) :: ready
    integer, intent(out) :: original_dim, reduced_dim, nblocks
    real(dp), intent(out) :: cut_fraction
    call cg_get_info(upt, ready, original_dim, reduced_dim, nblocks, cut_fraction)
  end subroutine UPT_get_coarse_graining_info

  !---------------------------------------------------------------------
  subroutine UPT_writehamiltonian(upt)
    TYPE(OUPT), pointer :: upt 

    if(upt%verbose.gt.0) write(*,*) '(uptight) write hamiltonian'    
    call write_csr_to_coo(upt%ham,upt%work_path)

  end subroutine UPT_writehamiltonian

  !---------------------------------------------------------------------
  
  subroutine UPT_feast(upt)

    TYPE(OUPT), pointer :: upt 
    
    if(upt%verbose.gt.0) write(*,*) '(uptight) Feast eigensolver'
    write(*,*) 'FEAST SOLVER DISABLED'    
    !call feast(upt)

  end subroutine UPT_feast
  
  !---------------------------------------------------------------------

  subroutine UPT_lanczos(upt)

    TYPE(OUPT), pointer :: upt 
    
    if(upt%verbose.gt.0) write(*,*) '(uptight) Lanczos diagonalization'
#ifdef UPT_MPI
    call mpi_barrier(upt_comm,ierr)
#endif
    call lanczos(upt)

  end subroutine UPT_lanczos
  
  !---------------------------------------------------------------------
  subroutine UPT_jd(upt)

    TYPE(OUPT), pointer :: upt 
    
    if(upt%verbose.gt.0) write(*,*) '(uptight) Jacobi-Davidson diagonalization'
    call jd(upt)

  end subroutine UPT_jd
  
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  subroutine UPT_lapack(upt)

    TYPE(OUPT), pointer :: upt

    if(upt%verbose.gt.0) write(*,*) '(uptight) LAPACK diagonalization'
    call lapack(upt)

  end subroutine UPT_lapack

  !---------------------------------------------------------------------

  subroutine UPT_alloc_eigv(upt)

    TYPE(OUPT), pointer :: upt 

    integer err, num_ev, nrow

    if (upt%ham%nrow.eq.0) then
       write(*,*) '(uptight) Error: Hamiltonian not initialized'
       return
    endif

    if (upt%num_vb.eq.0 .and. upt%num_cb.eq.0) then
       write(*,*) '(uptight) Error: number of states not initialized'
       return
    endif

    nrow = upt%ham%nrow
    num_ev = UPT%num_vb + UPT%num_cb

    if (.not.associated(upt%eigen_values)) then  
       allocate(upt%eigen_values(num_ev), STAT = err)
       IF (err.NE.0) CALL alloc_error('uptight','main','eigen_values')
    endif  
    
    if (.not.associated(upt%eigen_vectors)) then  
       allocate(upt%eigen_vectors(nrow, num_ev), STAT = err)
       IF (err.NE.0) CALL alloc_error('uptight','main','eigen_vectors')
    endif
   
    if (.not.associated(upt%particles)) then  
       allocate(upt%particles(num_ev), STAT = err)
       IF (err.NE.0) CALL alloc_error('uptight','main','particles')
    endif
   
    
  end subroutine UPT_alloc_eigv

  !---------------------------------------------------------------------

  subroutine UPT_write_eigenvectors(upt)

    TYPE(OUPT), pointer :: upt 
    
    if (id0 .eqv. .true.) then
      !if(upt%verbose.gt.0) write(*,*) '(uptight) Write out the eigenvectors'
      !call write_eigenstates(upt)  

      if(upt%verbose.gt.0) write(*,*) '(uptight) write mo file'

      if(trim(upt%out_format).eq.'cube') then
         call write_cube(upt)
      elseif(trim(upt%out_format).eq.'jvxl') then
        if (id0) call write_jvxl(upt)
      endif
    endif

  end subroutine UPT_write_eigenvectors

  !---------------------------------------------------------------------

  subroutine UPT_get_mat_el(upt,i,j,matel)
    
    TYPE(OUPT), target :: upt
    integer :: i, j
    complex(dp) :: matel
    
    complex(dp), DIMENSION( : ), POINTER  :: p_i, p_j 
    complex(dp), DIMENSION( : ), POINTER  :: p_aux
    integer :: num_ev
    complex(dp) :: matel_local
   
    if(.not.associated(upt%eigen_vectors)) then
       write(*,*) "ERROR: void eigenvectors"
       return
    endif


    num_ev = size(upt%eigen_vectors,2)

    p_i => upt%eigen_vectors(:,i)
    p_j => upt%eigen_vectors(:,j)

    allocate(p_aux(size(upt%eigen_vectors,1)))

    call sprs_ax(upt%ham%M,upt%ham%Mj,upt%ham%Mi, upt%ham%sparse_fmt, p_i, p_aux)
    
    matel_local = dot_product( p_j(shift_init:shift_end), p_aux(shift_init:shift_end) )

#ifdef UPT_MPI
    call mpi_barrier(upt_comm,ierr) 
    call mpi_allreduce(matel_local, matel, 1 , MPI_double_complex, MPI_SUM, upt_comm,ierr)
#else
    matel = matel_local
#endif

    deallocate(p_aux)
    
  end subroutine UPT_get_mat_el

  !---------------------------------------------------------------------
  subroutine UPT_get_hdim(upt,hdim)

    TYPE(OUPT), pointer :: upt
    integer hdim

    hdim = upt%ham%nrow

  end subroutine UPT_get_hdim

  !---------------------------------------------------------------------
  subroutine UPT_get_hamil(upt,nrow,ncol,fmt,M,Mj,Mi)

    TYPE(OUPT), pointer :: upt
    INTEGER :: nrow, ncol
    CHARACTER(1) :: fmt
    COMPLEX ( dp ), DIMENSION( : ), POINTER :: M
    INTEGER,        DIMENSION( : ), POINTER :: Mj
    INTEGER,        DIMENSION( : ), POINTER :: Mi

    if(.not.associated(upt%ham%M)) then
       write(*,*) "ERROR: matrix not created yet"
       return
    endif 
    
    nrow = upt%ham%nrow
    ncol = upt%ham%ncol
    fmt = upt%ham%sparse_fmt     ! L(ower),U(pper),F(ull)
    M => upt%ham%M 
    Mj => upt%ham%Mj
    Mi => upt%ham%Mi
     
  end subroutine UPT_get_hamil



end module UPTIGHT

!---------------------------------------------------------------------------
! Note:
! 
! - La struttura e' ri-ordinata da struct_building con gli H in fondo.
!
! - Attenzione:
!   basis%n_basis = n_atoms - n_dg_bonds, 
!   size(basis%coord) = n_atoms. (ci sono tutti gli atomi)
!
!   Questo e' sfruttato in TB_ham e probabilmente puo' generare confusione.
!
! - basis%dg_coord non e' usato da nessuna parte.
!
! - basis%n_dg(i): e' definito in neighbours. Non molto bello. : ((
!                  da spostare in struct_building.
!
! - Forse c'e' un bug in TB_ham per H periodici ? Da verificare
!
!
 
