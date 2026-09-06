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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Fortran 77 style subroutines for communication with UPTIGHT LIB
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The order you should invoke the subroutines during initialisation is the
!! following:
!!
!!   upt_initsession
!!   upt_setmpicomm          (after upt_initsession)
!!   upt_fillbasicparameters (after upt_initsession)
!!   upt_inituptight         (after upt_fillbasicparameters)
!!
!! After initialisation you can call the following methods:
!!
!!   upt_createhamiltonian    
!!   upt_lanczos              (after upt_createhamiltonian)
!!   upt_getzcsrhamiltonian   (after upt_createhamiltonian)
!!
!! In order to destroy a UPT instance call
!!
!!   upt_destructuptight
!!   upt_destructsession      (after upt_destructuptight)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!* Returns the size of the handler array
!!* @param  handlerSize  Contains the size of the handler array on exit.
subroutine upt_gethandlersize(handlerSize)
  use uptightAPICommon  ! if:mod:use
  implicit none
  integer :: handlerSize  ! if:var:out

  handlerSize = DAC_handlerSize

end subroutine upt_gethandlersize

!!* Initialises a new UPTIGHT instance
!!* @param  handler  Contains the handler for the new instance on return
subroutine upt_initsession(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_nullify_all, upt_set_defaults
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:out
  
  type(UPTPointers) :: pUPTs
  integer :: err
  
  IF( size(transfer(pUPTs, handler)) > size(handler) ) stop 'Handler size mismatch'
 
  NULLIFY(pUPTs%pUPT)
  ALLOCATE(pUPTs%pUPT, stat = err)
  if (err .ne. 0) then
     write(*,*) "UPT init session error"
     stop
  endif

  handler(:) = 0
  handler = transfer(pUPTs, handler, size(handler))

  call upt_nullify_all(pUPTs%pUPT)

  call upt_set_defaults(pUPTs%pUPT)

end subroutine upt_initsession

!!* Destroys a certain UPTIGHT instance
!!* @param  handler  Handler for the instance to destroy
subroutine upt_destructsession(handler)
  use uptightAPICommon  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout

  type(UPTPointers) :: pUPTs
  integer :: err

  pUPTs = transfer(handler, pUPTs)

  !deallocate(pUPTs%pUPTIn, stat= err)
  deallocate(pUPTs%pUPT, stat= err)

  if (err.ne.0) write(*,*) '(upt_destructsession) Deallocation error'

end subroutine upt_destructsession


!!* Initialises a new UPTIGHT instance
!!* @param  handler  Contains the handler for the new instance on return
subroutine upt_getversion(handler)
  use uptightAPICommon                 !if:mod:use
  use uptight, only : upt_version      !if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTS = transfer(handler, pUPTs)

  call upt_version(pUPTS%pUPT)

end subroutine upt_getversion

!!* Initialises the MPI communicator
!!* @param  handler  Contains the handler 
!!* @param  comm     contains the communicator
subroutine upt_setmpicomm(handler, comm)
  use uptightAPICommon                 !if:mod:use
  use uptight, only : upt_version      !if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: comm                      ! if:var:in

  type(UPTPointers) :: pUPTs
  
  pUPTS = transfer(handler, pUPTs)

  pUPTS%pUPT%mpi_comm = comm

end subroutine upt_setmpicomm


!!* Intialises a paths for UPTIGHT .
!!* @param handler                : Handler of the instance.
!!* @param verbose_lev            : verbose level 
!!* @param database_path (UPT_LC) : where to find the TB material database
!!* @param work_path     (UPT_LC) : current working directory
subroutine upt_set_paths(handler, databasePath, workPath, outPath)
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)   ! if:var:in
  character(LST) :: databasePath(1)     ! if:var:in
  character(LST) :: workPath(1)         ! if:var:in
  character(LST) :: outPath(1)          ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
    
  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  pUPT%database_path = trim( databasePath(1) ) // '/'
  pUPT%work_path = trim( workPath(1) ) // '/' 
  pUPT%out_path = trim( outPath(1) ) // '/'

end subroutine upt_set_paths

!!* Intialises a specific container for UPTIGHT input parameters.
!!* @param handler                : Handler of the instance.
!!* @param verbose_lev            : verbose level 
!!* @param gen_filename  (UPT_MC) : filename of input structure to load
!!* @param gen_out       (UPT_MC) : output filename
!!* @param sparse_fmt    (UPT_SC) : sparse format
!!* @param max_n_n                : to be initialized (usually 1 for n.n. tb)
!!* @param harrison_flag          : harrison scaling of TB parameters 1 = yes 
!!* @param relat_flag             : relativistic calculations         1 = yes
!!* @param potential_flag         : include potential shift           1 = yes
!!* @param optmat_flag            : computes optical matrix elements  1 = yes
!!* @param poldir                 : field polarization direction
!!* @param check_bondmap          : check bondmap internally
subroutine upt_fillbasicparameters(handler, verbose_lev, gen_filename, gen_outname, sparse_fmt, max_n_n, &
                                  harrison_flag, relat_flag, potential_flag, optmat_flag, & 
                                  poldir, c_axis_x, c_axis_y, c_axis_z, check_bondmap, &
                                  dg_coupl_scale, dg_onsite, hybrid_passivation)
  use precision         ! if:mod:use
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  use uptight           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)   ! if:var:in
  integer :: verbose_lev                ! if:var:in 
  character(MST) :: gen_filename(1)     ! if:var:in
  character(MST) :: gen_outname(1)      ! if:var:in
  character(MST) :: sparse_fmt(1)       ! if:var:in
  integer :: max_n_n                    ! if:var:in
  integer :: harrison_flag              ! if:var:in
  integer :: relat_flag                 ! if:var:in
  integer :: potential_flag             ! if:var:in
  integer :: optmat_flag                ! if:var:in
  integer :: poldir                     ! if:var:in
  integer :: check_bondmap              ! if:var:in
  integer :: hybrid_passivation         ! if:var:in
  real(dp) :: c_axis_x                  ! if:var:in  
  real(dp) :: c_axis_y                  ! if:var:in
  real(dp) :: c_axis_z                  ! if:var:in
  real(dp) :: dg_coupl_scale            ! if:var:in  
  real(dp) :: dg_onsite                 ! if:var:in  

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
    
  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  !! Fill in parameters
  pUPT%gen_filename = trim(gen_filename(1))
  pUPT%gen_out =  trim(gen_outname(1))   
  pUPT%state_file  = "states.upt"
  select case(trim(sparse_fmt(1)))
  case('upper')
     pUPT%sparse_format = 'U'
  case('lower')
     pUPT%sparse_format = 'L'
  case('full')
     pUPT%sparse_format = 'F'
  end select
  pUPT%verbose = verbose_lev

  !! Initialize globals

  !write(*,*) '(upt-debug) db_path: ', trim(database_path)
  !write(*,*) '(upt-debug) wk_path: ', trim(work_path)
  !write(*,*) '(upt-debug) gen: ', trim(pUPTIn%gen_filename)

  !! Fill in UPT parameters

  pUPT%potential_flag = .false.
  pUPT%relat = .false.
  pUPT%n_spin = 2
  pUPT%scaling = .false.
  pUPT%optmat = .false.
  pUPT%check_bondmap = .false.
  pUPT%hybrid_passivation = .false.
  if(relat_flag.eq.1) pUPT%relat = .true.
  if(relat_flag.eq.1) pUPT%n_spin = 1
  if(harrison_flag.eq.1) pUPT%scaling = .true.
  if(optmat_flag.eq.1)  pUPT%optmat = .true.
  if(potential_flag.eq.1)  pUPT%potential_flag = .true.
  if(check_bondmap.eq.1)  pUPT%check_bondmap = .true.
  if(hybrid_passivation.eq.1)  pUPT%hybrid_passivation = .true.

  !write(*,*) '(upt-debug) Harrison scaling: ', pUPT%fuzzy
  !write(*,*) '(upt-debug) Optical Matrix: ', pUPT%optmat

  pUPT%poldir = poldir    
  !pUPT%nn_map%max_near = max_n_n

  pUPT%c_axis = (/ c_axis_x, c_axis_y, c_axis_z /)

  !write(*,*) '(upt-debug) c-axis: ', pUPT%c_axis

  pUPT%d_H = dg_coupl_scale !1.d0/sqrt(dg_coupl_scale)
  pUPT%E_H = dg_onsite

  ! nullify some pointers

  NULLIFY(pUPT%materials)
  NULLIFY(pUPT%pot_data)

  NULLIFY(pUPT%eigen_values)
  NULLIFY(pUPT%eigen_vectors) 

end subroutine upt_fillbasicparameters


subroutine upt_setworkpath(handler, workPath)
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  implicit none  
  integer :: handler(DAC_handlerSize)     ! if:var:in
  character(LST) :: workPath(1)     ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
  pUPTs = transfer(handler, pUPTs)
 
  pUPT => pUPTs%pUPT

  pUPT%work_path = trim( workPath(1) ) // '/'

end subroutine upt_setworkpath


subroutine upt_setoutpath(handler, outPath)
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  implicit none  
  integer :: handler(DAC_handlerSize)     ! if:var:in
  character(LST) :: outPath(1)            ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
  
  pUPT%out_path = trim( outPath(1) ) // '/'

end subroutine upt_setoutpath

subroutine upt_setloadpath(handler, loadPath)
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  implicit none  
  integer :: handler(DAC_handlerSize)     ! if:var:in
  character(LST) :: loadPath(1)     ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
    
  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  pUPT%load_path = trim( loadPath(1) ) // '/'
  !print*,'(UPT) load_path: ',trim(pUPT%load_path)

end subroutine upt_setloadpath

!!* Set the states file to read from/write to
!!* @param handler
!!* @param filename
subroutine upt_setstatefile(handler, filename)
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)     ! if:var:in
  character(LST) :: filename(1)     ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  pUPT%state_file = trim( filename(1) )
  print*,'(UPT) states file: ',trim(pUPT%state_file)

end subroutine upt_setstatefile



!!* Set parameters for representation of output states
!!* @param handler
!!* @param format
!!* @param step
subroutine upt_setoutput(handler, out_form, step)
  use precision         ! if:mod:use
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  use uptight           ! if:mod:use
  implicit none  
  integer :: handler(DAC_handlerSize)     ! if:var:in
  integer :: out_form                     ! if:var:in
  real(dp) :: step                        ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  pUPT%grid_step = step
  if (out_form.eq.1) pUPT%out_format = "jvxl"
  if (out_form.eq.2) pUPT%out_format = "cube"

end subroutine upt_setoutput

!!* Set the library to assemble the optical matrix
!!* @param handler                : Handler of the instance.
!!* @param optmat_flag            : computes optical matrix elements  1 = yes
!!* @param poldir                 : field polarization direction
subroutine upt_setpmatrix(handler, optmat_flag, poldir) 
  use uptightAPICommon  ! if:mod:use
  use upt_param         ! if:mod:use
  use uptight           ! if:mod:use
  implicit none 
  integer :: handler(DAC_handlerSize)   ! if:var:in  
  integer :: optmat_flag                ! if:var:in
  integer :: poldir                     ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  pUPT%optmat = .false.
  if(optmat_flag.eq.1)  pUPT%optmat = .true.

  pUPT%poldir = poldir    

end subroutine upt_setpmatrix



!!* Initializes a given UPTIGHT instance. 
!!* Read structure file and database materials and initializes all of the
!!* UPTIGHT data containers.
!!* @param handler Number for the UPTIGHT instance to initialise and the number
!!* for the input data container, which should be used for the initialisation.
subroutine upt_inituptight(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_init
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs
 
  pUPTs = transfer(handler, pUPTs)
  call upt_init(pUPTs%pUPT)

end subroutine upt_inituptight



!!* Destructs a given UPTIGHT instance.
!!* @param handler Number for the UPTIGHT instance to destroy.
subroutine upt_destructuptight(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_destruct
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)
  call upt_destruct(pUPTs%pUPT)

end subroutine upt_destructuptight



!!* Add an external potential on the atoms (external, piezo, piro, etc.) .
!!* @param handler       : handler for the UPTIGHT instance.
!!* @param nAtoms        : number of atoms in the structure
!!* @param potential     : atom-projected potential vector 
subroutine upt_addpotential(handler,nAtoms,potential)
  use uptightAPICommon !if:mod:use
  use precision, only : dp
  use upt_param        
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in
  integer  :: nAtoms                    ! if:var:in
  real(dp) :: potential(nAtoms)         ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  if ( pUPT%basis%n_basis .eq. 0 ) then
     write(*,*) 'ERROR: upt_inituptight must be called first'
     return
  end if

  if ( pUPT%basis%n_basis .gt. nAtoms ) then
     write(*,*) 'ERROR: number of atoms mismatch'
     return
  end if

  if (.not.associated(pUPT%pot_data)) then
     allocate(pUPT%pot_data(nAtoms))
     pUPT%pot_data = 0.d0
  end if
  ! this will be deallocated in upt_destructuptight

  pUPT%pot_data(1:natoms) =  pUPT%pot_data(1:natoms) &
                            + potential(1:natoms)

  pUPT%potential_flag = .true. 

end subroutine upt_addpotential


subroutine upt_erasepotential(handler)
  use uptightAPICommon !if:mod:use
  use upt_param       
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  if(associated(pUPT%pot_data)) then
     pUPT%pot_data=0.d0
  end if

end subroutine upt_erasepotential


!!* Computes the system hamiltonian
!!* @param handler Number for the UPTIGHT instance.
!!* @param sparse_fmt             : sparse format of upper matrix
subroutine upt_createhamiltonian(handler,sparse_fmt)
  use uptightAPICommon  ! if:mod:use
  use upt_param       
  use globals           ! if:mod:use
  use uptight, only : upt_hamiltonian
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  character(MST) :: sparse_fmt(1)       ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  select case(trim(sparse_fmt(1)))
  case('upper')
     pUPT%sparse_format = 'U'
  case('lower')
     pUPT%sparse_format = 'L'
  case('full')
     pUPT%sparse_format = 'F'
  end select

  call upt_hamiltonian(pUPTs%pUPT)  

end subroutine upt_createhamiltonian

!!* Configure optional block-basis coarse graining before creating the Hamiltonian.
subroutine upt_set_coarse_graining(handler, enabled, n_blocks, emin, emax, imbalance)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_configure_coarse_graining ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize) ! if:var:in
  integer :: enabled ! if:var:in
  integer :: n_blocks ! if:var:in
  real(dp) :: emin ! if:var:in
  real(dp) :: emax ! if:var:in
  real(dp) :: imbalance ! if:var:in
  type(UPTPointers) :: pUPTs
  pUPTs = transfer(handler, pUPTs)
  call upt_configure_coarse_graining(pUPTs%pUPT, enabled /= 0, n_blocks, emin, emax, imbalance)
end subroutine upt_set_coarse_graining

!!* Return coarse-graining reduction statistics after Hamiltonian creation.
subroutine upt_get_coarse_graining_info(handler, ready, original_dim, reduced_dim, n_blocks, cut_fraction)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_get_coarse_graining_info_f => upt_get_coarse_graining_info ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize) ! if:var:in
  integer :: ready ! if:var:out
  integer :: original_dim ! if:var:out
  integer :: reduced_dim ! if:var:out
  integer :: n_blocks ! if:var:out
  real(dp) :: cut_fraction ! if:var:out
  logical :: is_ready
  type(UPTPointers) :: pUPTs
  pUPTs = transfer(handler, pUPTs)
  call upt_get_coarse_graining_info_f(pUPTs%pUPT,is_ready,original_dim,reduced_dim,n_blocks,cut_fraction)
  ready = merge(1,0,is_ready)
end subroutine upt_get_coarse_graining_info

!!* Computes the system hamiltonian
!!* @param handler Number for the UPTIGHT instance.
subroutine upt_printhamiltonian(handler)
  use uptightAPICommon  ! if:mod:use
  use upt_param       
  use globals           ! if:mod:use
  use uptight, only : upt_writehamiltonian
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)
  
  call upt_writehamiltonian(pUPTs%pUPT)  

end subroutine upt_printhamiltonian

!!* Add k-points 
!!* @param handler       : handler for the UPTIGHT instance.
!!* @param k_vec(3)      : k-vector in fractional coord.  
subroutine upt_setkpoint(handler,k_vec)
  use uptightAPICommon !if:mod:use
  use precision, only : dp
  use upt_param        
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in
  real(dp) :: k_vec(3)                  ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  ! this will be deallocated in upt_destructuptight
  pUPT%k_point(:) = k_vec(:)

end subroutine upt_setkpoint


!!* Diagonalize using lanczos iterative methods.
!!* @param handler                : Handler of the instance.
!!* @param n_vb                   : Number of valence eigenstates
!!* @param n_cb                   : Number of conduction eigenstates       
!!* @param min_iter               : minimum number of iterations (fast) (~2)
!!* @param long_iter              : number of lanczos iterations (~30) 
!!* @param max_iter               : maximum number of iterations (~3000)
!!* @param guess_vb               : valence bottom guess 
!!* @param guess_cb               : conduction bottom
!!* @param fast_tol               : fast tolerance  (1e-7) 
!!* @param long_tol               : long tolerance  (1e-10)
!!* @param ort_tol                : orthogonality tolerance (1e-6)
!!* @param dynamic                : shift towards converged eigenvalues
!!* @param bitoff                 : tiny offset away from the last eigenvalue
subroutine upt_lanczosdiag(handler, st_vb, st_cb, n_vb, n_cb, guess_vb, guess_cb, &
                           min_iter, long_iter, max_iter, fast_tol, long_tol, ort_tol, &
                           dynamic, bitoff)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use uptight, only : upt_lanczos
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: st_vb                     ! if:var:in 
  integer :: st_cb                     ! if:var:in  
  integer :: n_vb                      ! if:var:in 
  integer :: n_cb                      ! if:var:in
  integer :: min_iter                  ! if:var:in
  integer :: long_iter                 ! if:var:in
  integer :: max_iter                  ! if:var:in
  real(dp) :: guess_vb                 ! if:var:in
  real(dp) :: guess_cb                 ! if:var:in 
  real(dp) :: fast_tol                 ! if:var:in
  real(dp) :: long_tol                 ! if:var:in
  real(dp) :: ort_tol                  ! if:var:in
  integer :: dynamic                   ! if:var:in
  real(dp) :: bitoff                   ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: upt

  pUPTs = transfer(handler, pUPTs)

  upt => pUPTs%pUPT

  upt%start_vb = st_vb  
  upt%start_cb = st_cb

  upt%num_vb = n_vb  
  upt%num_cb = n_cb

  upt%lambda_vb = guess_vb
  upt%lambda_cb = guess_cb

  upt%min_iter = min_iter
  upt%long_iter = long_iter
  upt%max_iter = max_iter

  upt%fast_tol = fast_tol
  upt%long_tol = long_tol
  upt%ort_tol  = ort_tol
     

  upt%dynamic= .false.
  if(dynamic.eq.1) upt%dynamic  = .true. 
  upt%seed_flag      = .false.

  upt%bitoff = bitoff
  !..................................................................
  ! Checks consistence on the number of iterations
  !..................................................................
  if (min_iter.gt.long_iter) upt%long_iter = upt%min_iter

  if (min_iter.gt.max_iter .OR. long_iter.gt.max_iter) upt%max_iter = upt%long_iter
  
  call upt_lanczos(upt)

end subroutine upt_lanczosdiag

subroutine upt_solver_flag(handler, flag)
  use uptightAPICommon  ! if:mod:use
  use upt_param
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in  
  integer :: flag  ! if:var:in 

  type(UPTPointers) :: pUPTs
  pUPTs = transfer(handler, pUPTs)

  pUPTs%pUPT%solver_flag = flag

end subroutine upt_solver_flag


!!* Diagonalize using JD methods.
!!* @param handler                : Handler of the instance.
!!* @param n_vb                   : Number of valence eigenstates
!!* @param n_cb                   : Number of conduction eigenstates       
!!* @param guess_vb               : valence bottom guess 
!!* @param guess_cb               : conduction bottom
!!* @param long_tol               : long tolerance  (1e-10)
subroutine upt_jd_diag(handler, st_vb, st_cb, n_vb, n_cb, guess_vb, guess_cb, &
                       long_tol)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use uptight, only : upt_jd
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: st_vb                     ! if:var:in 
  integer :: st_cb                     ! if:var:in  
  integer :: n_vb                      ! if:var:in 
  integer :: n_cb                      ! if:var:in
  real(dp) :: guess_vb                 ! if:var:in
  real(dp) :: guess_cb                 ! if:var:in 
  real(dp) :: long_tol                 ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: upt

  pUPTs = transfer(handler, pUPTs)

  upt => pUPTs%pUPT

  upt%start_vb = st_vb  
  upt%start_cb = st_cb

  upt%num_vb = n_vb  
  upt%num_cb = n_cb

  upt%lambda_vb = guess_vb
  upt%lambda_cb = guess_cb

  upt%long_tol = long_tol

  !upt%twice_cb_flag = .false. 
  !if(twice_cb.eq.1) upt%twice_cb_flag  = .true. 
  !upt%twice_vb_flag = .false. 
  !if(twice_vb.eq.1) upt%twice_vb_flag  = .true. 
  !upt%dynamic= .false.
  !if(dynamic.eq.1) upt%dynamic  = .true. 
  upt%seed_flag      = .false.
  
  call upt_jd(upt)

end subroutine upt_jd_diag




SUBROUTINE upt_feastsolver(handler, emin, emax, m0)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use uptight, only : upt_feast
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in  
  integer :: m0                      ! if:var:in 
  real(dp) :: emin                 ! if:var:in
  real(dp) :: emax                 ! if:var:in 



  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: upt

  pUPTs = transfer(handler, pUPTs)

  upt => pUPTs%pUPT

 upt%num_vb = m0
 upt%num_cb = m0
 upt%lambda_vb = emin
 upt%lambda_cb = emax

 call upt_feast(upt)

END SUBROUTINE upt_feastsolver



subroutine upt_lapacksolver(handler, n_vb, n_cb, guess_vb, guess_cb)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use uptight, only : upt_lapack
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: n_vb                      ! if:var:in
  integer :: n_cb                      ! if:var:in
  real(dp) :: guess_vb                 ! if:var:in
  real(dp) :: guess_cb                 ! if:var:in



  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: upt

  pUPTs = transfer(handler, pUPTs)

  upt => pUPTs%pUPT


  upt%num_vb = n_vb
  upt%num_cb = n_cb
  upt%lambda_vb = guess_vb
  upt%lambda_cb = guess_cb

  call upt_lapack(upt)

end subroutine upt_lapacksolver



!!* Create containers for eigenvectors
!!* @param handler Number for the UPTIGHT instance.
subroutine upt_alloc_states(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_alloc_eigv
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in  

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  call upt_alloc_eigv(pUPTs%pUPT)

end subroutine upt_alloc_states

!!* Store eigenvectors on files
!!* @param handler Number for the UPTIGHT instance.
subroutine upt_write_states(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_write_eigenvectors
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)
  call upt_write_eigenvectors(pUPTs%pUPT)

end subroutine upt_write_states




!!* Set the number of states (used for initialize reading)
!!* @param handler Number for the UPTIGHT instance.
subroutine upt_set_num_states(handler, n_vb, n_cb)
  use uptightAPICommon                 ! if:mod:use
  implicit none 
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: n_vb                      ! if:var:in
  integer :: n_cb                      ! if:var:in
  
  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs)  
  
  pUPTs%pUPT%num_vb = n_vb  
  pUPTs%pUPT%num_cb = n_cb  
  
end subroutine upt_set_num_states




!!* Read eigenvectors from files
!!* @param handler Number for the UPTIGHT instance.
!!* @param nev Number of valence read
!!* @param nec Number of conduction read
subroutine upt_read_states(handler, filename, nev, nec)
  use uptightAPICommon  ! if:mod:use
  use globals, only : LST
  use precision, only : dp
  use uptight, only : upt_get_mat_el
  use savemofile, only : read_eigenstates !, read_eigenvalues
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: nev                       ! if:var:inout
  integer :: nec                       ! if:var:inout
  character(LST) :: filename(1)        ! if:var:in

  integer :: i
  complex(dp) :: matel

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  !filename = trim(pUPTs%pUPT%load_path)//"eigv.dat"

  !call read_eigenvalues(pUPTs%pUPT,filename,nev,nec) 

  !filename = trim(pUPTs%pUPT%load_path)//trim(pUPTs%pUPT%state_file)

  call read_eigenstates(pUPTs%pUPT,trim(filename(1)),nev,nec)

  ! do a sanity check
  do i = 1, nev+nec
    call upt_get_mat_el(pUPTs%pUPT, i, i, matel)
    write(*,*) i, ' : eigenvalue read ', pUPTs%pUPT%eigen_values(i)
    write(*,*)    '         <i|H|i> = ',  matel
    if (abs(matel - pUPTs%pUPT%eigen_values(i)) > 1d-6*abs(matel)) then
      write(*,*) 'Eigenstate ', i, ' from file ', filename, &
                 ' is inconsistent with the given Hamiltonian.'
      exit
      end if
  end do

end subroutine upt_read_states


!!* Get the Hamiltonian dimension
!!* @param handler     : Number for the UPTIGHT instance.
!!* @param hdim        : Hamiltonian number of rows
subroutine upt_get_hamildim(handler, hdim)
  use uptightAPICommon  ! if:mod:use
  use sparse_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: hdim                      ! if:var:out

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  hdim = pUPTs%pUPT%ham%nrow

  !call get_nrow(pUPTs%pUPT%ham, hdim)

end subroutine upt_get_hamildim

!!* Get the Hamiltonian non zero values
!!* @param handler     : Number for the UPTIGHT instance.
!!* @param hdim        : Number of non zero values
subroutine upt_get_hamilnnz(handler, nnz)
  use uptightAPICommon  ! if:mod:use
  use sparse_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: nnz                       ! if:var:out

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  nnz = pUPTs%pUPT%ham%nnz
  !call get_nnz(pUPTs%pUPT%ham, nnz)

end subroutine upt_get_hamilnnz

!!* Get the Hamiltonian row size 
!!* @param handler     : Number for the UPTIGHT instance.
!!* @param row         : row number
!!* @param sz          : size of a row
subroutine upt_get_hamil_rowsize(handler, row, sz)
  use uptightAPICommon   
  use sparse_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  !if:var:in
  integer :: row                      !if:var:in   
  integer :: sz                       !if:var:out

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  sz = pUPTs%pUPT%ham%Mi(row+1) - pUPTs%pUPT%ham%Mi(row)

end subroutine upt_get_hamil_rowsize


!!* Get a Hamiltonian row
!!* @param handler     : Number for the UPTIGHT instance.
!!* @param row         : row number
!!* @param colind      : column indeces 
!!* @param vals        : column values
subroutine upt_get_hamil_row(handler, row, colind, vals)
  use precision, only : dp
  use uptightAPICommon  
  use sparse_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  !if:var:in
  integer :: row                      !if:var:in  
  integer :: colind(*)                !if:var:inout
  complex(dp) :: vals(*)              !if:var:inout

  integer :: i, j
  type(UPTPointers) :: pUPTs
  type(CSR), pointer :: ham

  pUPTs = transfer(handler, pUPTs)

  ham => pUPTs%pUPT%ham

  j = 1
  do i = ham%Mi(row), ham%Mi(row+1)-1
     colind(j) = ham%Mj(i)
     vals(j) = ham%M(i)
     j = j + 1
  enddo

end subroutine upt_get_hamil_row

!!* Get the Hamiltonian in sparse format
subroutine upt_get_csr_hamiltonian(handler, nrow, fmt, A, JA, IA)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use uptight, only : upt_get_hamil
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: nrow     ! if:var:out
  character(1) :: fmt    ! if:var:out
  complex(dp) :: A(1) ! if:var:out 
  integer :: JA(1)    ! if:var:out
  integer :: IA(1)    ! if:var:out

  type(UPTPointers) :: pUPTs
  COMPLEX ( dp ), DIMENSION( : ), POINTER :: M
  INTEGER,        DIMENSION( : ), POINTER :: Mj
  INTEGER,        DIMENSION( : ), POINTER :: Mi
  integer :: ncol, i, nnz

  pUPTs = transfer(handler, pUPTs)

  call upt_get_hamil(pUPTs%pUPT,nrow,ncol,fmt,M,Mj,Mi)
  
  do i = 1, nrow+1
     IA(i) = Mi(i)
  enddo
  nnz = Mi(nrow+1)-1
  
  do i = 1, nnz
     A(i) = M(i)
  enddo
  
  do i = 1, nnz
     JA(i) = Mj(i)
  enddo

end subroutine upt_get_csr_hamiltonian

!!* Get the Hamiltonian in sparse format
subroutine upt_set_csr_hamiltonian(handler, nrow, fmt, A, JA, IA)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use sparse_matrix, only : create_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: nrow     ! if:var:out
  character(1) :: fmt    ! if:var:out
  complex(dp) :: A(*) ! if:var:out 
  integer :: JA(*)    ! if:var:out
  integer :: IA(*)    ! if:var:out

  type(UPTPointers) :: pUPTs
  integer :: i,nnz, err

  pUPTs = transfer(handler, pUPTs)

  nnz = IA(nrow+1)-1
  
  !pUPTs%pUPT%ham%M => A(1:nnz)      DOES NOT WORK :( 
  !pUPTs%pUPT%ham%Mj => JA(1:nnz)    WE MAKE A COPY
  !pUPTs%pUPT%ham%Mi => IA(1:nrow+1) 

  call create_matrix(pUPTs%pUPT%ham,nrow,nrow,nnz)

  pUPTs%pUPT%ham%sparse_fmt = fmt 

  do i = 1, nnz
     pUPTs%pUPT%ham%M(i) = A(i)
  enddo
  do i = 1, nnz
     pUPTs%pUPT%ham%Mj(i)=JA(i)
  enddo
  do i = 1, nrow+1
     pUPTs%pUPT%ham%Mi(i)=IA(i)
  enddo

end subroutine upt_set_csr_hamiltonian        

!!* Store eigenvectors on files
!!* @param handler    : number for the UPTIGHT instance.
!!* @param num_ev     : number of states
!!* @param hdim       : dimensione di 
!!* @param eigenvals(num_ev)
!!* @param eigenstates(hdim,num_ev) 
subroutine upt_get_states(handler, num_ev, hdim, eigenvals, eigenstates, types)
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use precision, only : dp
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: num_ev                        ! if:var:in
  integer :: hdim                          ! if:var:in
  real(dp) :: eigenvals(num_ev)             ! if:var:out
  complex(dp) :: eigenstates(hdim,num_ev)   ! if:var:out
  integer :: types(num_ev)                  ! if:var:out

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
  integer :: k

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  eigenvals(1:num_ev) = pUPT%eigen_values(1:num_ev)
  types(1:num_ev) = pUPT%particles(1:num_ev)

  if (hdim .ne. size(pUPT%eigen_vectors, 1)) then 
          write(*,*) 'ERROR: eigenvector dimension mismatch'
          return
  end if

  eigenstates = pUPT%eigen_vectors
  !do k = 1,num_ev
  !   print*,'norm ',k, dot_product( eigenstates(:,k),eigenstates(:,k) )
  !enddo
    !eigenstates_im(1:hdim,1:num_ev) = imag(pUPT%eigen_vectors(1:hdim,1:num_ev))  

end subroutine upt_get_states

!!* Store eigenvectors on files
!!* @param handler    : number for the UPTIGHT instance.
!!* @param id         : the index of the state
!!* @param hdim       : dimensione di 
!!* @param eigenval
!!* @param eigenstate(hdim)
subroutine upt_set_state(handler, num_ev, id, hdim, eigenval, eigenstate, particle)
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use precision, only : dp
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: id                        ! if:var:in
  integer :: num_ev                    ! if:var:in
  integer :: hdim                      ! if:var:in
  real(dp) :: eigenval                 ! if:var:in
  complex(dp) :: eigenstate(hdim)      ! if:var:in
  integer :: particle                  ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
  integer :: k

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  if (.not.associated(pUPT%eigen_vectors)) then
     write(*,*) "(upt) allocating eigenvectors:",pUPT%ham%nrow,"x",num_ev
     allocate(pUPT%eigen_vectors(pUPT%ham%nrow,num_ev))
     pUPT%eigen_vectors = (0.d0,0.d0)
  end if

  if (.not.associated(pUPT%eigen_values)) then
     write(*,*) "(upt) allocating eigenvalue array of size",num_ev
     allocate(pUPT%eigen_values(num_ev))
     pUPT%eigen_values = (0.d0,0.d0)
  end if

  if (.not.associated(pUPT%particles)) then
     !write(*,*) "(upt) allocating particles type array of size",num_ev
     allocate(pUPT%particles(num_ev))
     pUPT%particles = 0
  end if


  pUPT%eigen_values(id) = eigenval
  pUPT%particles(id) = particle

  if (hdim .ne. size(pUPT%eigen_vectors, 1)) then 
    write(*,*) 'ERROR: eigenvector dimension mismatch'
    return
  end if


  pUPT%eigen_vectors(:,id) = eigenstate


  
  id = 0
  do k = 1,num_ev
    if (pUPT%particles(k) .eq. -1) then
      id = id + 1
    end if
  end do

  pUPT%num_vb = id
  pUPT%num_cb = num_ev - id


end subroutine upt_set_state



!!* Store eigenvectors on files
!!* @param handler    : number for the UPTIGHT instance.
!!* @param i          : left state index 
!!* @param j          : right state index
!!* @param matel      : complex matrix element
subroutine upt_get_matel(handler, i, j, matel)
  use uptightAPICommon  ! if:mod:use
  use upt_param
  use precision, only : dp
  use uptight, only : upt_get_mat_el
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: i                        ! if:var:in
  integer :: j                          ! if:var:in  
  complex(dp) :: matel                 ! if:var:out  


  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
  
  call upt_get_mat_el(pUPT, i, j, matel)

end subroutine upt_get_matel
  

subroutine upt_project_pot(handler, i, potential, average)
  use uptightAPICommon  ! if:mod:use  use upt_param
  use precision, only : dp
  use upt_param
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: i                        ! if:var:in
  real(dp) :: potential(1)            ! if:var:in
  real(dp) :: average                 ! if:var:out

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT 
  integer :: j

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  if (.not.associated(pUPT%eigen_vectors)) then
     average = 1d+38
     return 
  endif
 
  do j = 1,pUPT%Ham%nrow
     ! to be completed
  end do
  
end subroutine upt_project_pot


subroutine upt_get_ion_numorbitals(handler,ion_block_vector)
  use uptightAPICommon  ! if:mod:use  use upt_param
  use precision, only : dp
  use TB_ham, only : get_ion_block_size
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: ion_block_vector(*)      ! if:var:inout
  
  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs)  
  
  call get_ion_block_size(pUPTs%pUPT, ion_block_vector)
  
end subroutine upt_get_ion_numorbitals


subroutine upt_get_ion_orbitals(handler, ion, orbitals)
  use uptightAPICommon  ! if:mod:use  use upt_param
  !use precision, only : dp
  use TB_ham, only : get_ion_orbitals
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: ion              ! if:var:inout
  integer :: orbitals(*)      ! if:var:inout
  
  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs)  
  
  call get_ion_orbitals(pUPTs%pUPT, ion, orbitals)
  
end subroutine upt_get_ion_orbitals

subroutine upt_set_verbosity(handler,verbose_lev)
  use uptightAPICommon  ! if:mod:use  use upt_param  
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: verbose_lev               ! if:var:in

  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs) 
   
  pUPTs%pUPT%verbose = verbose_lev

end subroutine upt_set_verbosity

subroutine real_test(re)
  use precision, only : dp
  real(dp) :: re     ! if:var:out  

  re = 3.14159265358979323844_dp
end subroutine real_test


subroutine complex_test(re,im,zz)
  use precision, only : dp
  implicit none
  real(dp) :: re     ! if:var:out
  real(dp) :: im     ! if:var:out
  complex(dp) :: zz  ! if:var:out

  re=4.d0; im=1.d0
  zz=re+(0.d0,1.d0)*im
 
end subroutine complex_test

!!* Add an external potential on the atoms (external, piezo, piro, etc.) .
!!* @param handler       : handler for the UPTIGHT instance.
!!* @param nAtoms        : number of atoms in the structure
!!* @param strain_xx     : atom-projected potential vector 
!!* @param strain_yy     : atom-projected potential vector 
!!* @param strain_zz     : atom-projected potential vector 
subroutine upt_setstrain(handler,nAtoms,strain_xx,strain_yy,strain_zz)
  use uptightAPICommon              !if:mod:use
  use precision, only : dp
  use upt_param        
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in
  integer  :: nAtoms                    ! if:var:in
  real(dp) :: strain_xx(nAtoms)         ! if:var:in
  real(dp) :: strain_yy(nAtoms)         ! if:var:in
  real(dp) :: strain_zz(nAtoms)         ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
  integer :: err

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  if ( pUPT%basis%n_basis .eq. 0 ) then
     write(*,*) 'ERROR: upt_inituptight must be called first'
     return
  end if

  if ( pUPT%basis%n_basis .ne. nAtoms ) then
     write(*,*) 'ERROR: number of atoms mismatch: '
     write(*,*) 'basis: ', pUPT%basis%n_basis, 'nAtoms: ',nAtoms
     return
  end if

  if (associated(pUPT%basis%strain)) deallocate(pUPT%basis%strain)       

  if (.not.associated(pUPT%basis%strain)) then
     allocate(pUPT%basis%strain(3,nAtoms), stat=err)
     if(err.ne.0) then
        write(*,*) "ALLOCATION ERROR: could not allocate strain"
        return
     endif
  endif
 
  pUPT%basis%strain(1,1:nAtoms) = strain_xx(1:nAtoms)
  pUPT%basis%strain(2,1:nAtoms) = strain_yy(1:nAtoms)
  pUPT%basis%strain(3,1:nAtoms) = strain_zz(1:nAtoms)

  pUPT%d_onsite_shift_flag = .true.

end subroutine upt_setstrain
