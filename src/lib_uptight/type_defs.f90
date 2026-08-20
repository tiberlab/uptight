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
module type_defs
  
  use precision
  use globals, only : MST, n_ref_states
  use mpi_globals
  use exceptions
  
  implicit none
  private

  public :: ion_basis, ion_orbit, pair_coupling, TStructure
  public :: parent_data, material_data, create_basis, write_basis
  public :: write_materials, allocate_mat_ion, allocate_parent_ion
  public :: allocate_ion, allocate_parent_pair, allocate_pair
  public :: allocate_ons_corr, allocate_so_corr, allocate_intra
  public :: destroy_material, destroy_basis, nullify_basis, nullify_ion
  public :: destroy_structure, nullify_structure

  integer, public, parameter :: ATOMLEN = 2
  integer, public, parameter :: TYPELEN = 2 ! TEMPORARY: should be 6
  integer, public, parameter :: CPLLEN = 5
  integer, public, parameter :: STATELEN = 5
  integer, public, parameter :: MAX_NUM_NEIGHBOURS = 24

  character(2), public :: hydrosat = "Hx"

  !===========================================================================
  !
  ! Definition of types.
  !
  !___________________________________________________________________________
  !
  ! => cplxs
  ! => ints
  ! => names
  ! => pair_names
  ! => unit_lattice
  ! => ion
  ! => ion_basis
  ! => ion_orbit
  ! => pair_coupling
  !
  !===========================================================================
  ! defines the super basis
    
  TYPE ion_basis
     
     INTEGER                                            :: n_basis
     INTEGER                                            :: n_dg_bond
     CHARACTER( TYPELEN ),   DIMENSION( : ),    POINTER :: atomtypes => null()
     REAL ( dp ),            DIMENSION( :, : ), POINTER :: coord => null()
     REAL ( dp ),            DIMENSION( :, : ), POINTER :: dg_coord => null()
     REAL ( dp ),            DIMENSION( :, : ), POINTER :: strain => null()
     INTEGER,                DIMENSION( : ),    POINTER :: mat => null()
     INTEGER,                DIMENSION( : ),    POINTER :: ion => null()
     INTEGER,                DIMENSION( : ),    POINTER :: type => null()
     INTEGER,                DIMENSION( : ),    POINTER :: n_dg => null()
     INTEGER,                DIMENSION( : ),    POINTER :: n_st => null()
     INTEGER,                DIMENSION( : ),    POINTER :: valence => null()
     LOGICAL                                            :: periodic
     CHARACTER ( LEN = 3 )                              :: periodic_BC
     REAL ( dp ),            DIMENSION( : ),    POINTER :: start => null()
     REAL ( dp ),            DIMENSION( :, : ), POINTER :: prim => null()
     REAL ( dp ),            DIMENSION( :, : ), POINTER :: rec_latt => null()
     
  END TYPE ion_basis
  
  !___________________________________________________________________________
    
  ! defines the onsite ETB parameters

  TYPE ion_orbit

     REAL ( dp ),            DIMENSION( : ), POINTER :: energy => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: so_energy => null()
     CHARACTER( STATELEN ),  DIMENSION( : ), POINTER :: state => null()
     CHARACTER ( LEN = 1 ),  DIMENSION( : ), POINTER :: so_state => null()
     INTEGER,                DIMENSION( : ), POINTER :: deg => null()
     INTEGER,                DIMENSION( : ), POINTER :: ind_ref => null()
     INTEGER,                DIMENSION( : ), POINTER :: ind_ref_inv => null()
     REAL ( dp )                                     :: offset
     REAL ( dp )                                     :: b_d
     CHARACTER ( TYPELEN )                           :: name
     INTEGER                                         :: valence
     REAL ( dp )                                     :: content

  END TYPE ion_orbit
  
  !___________________________________________________________________________
  
  ! defines the non - diagonal ETB elements ( couplings between neighbours )

  TYPE pair_coupling
     
     CHARACTER ( CPLLEN ),   DIMENSION( : ), POINTER :: cpl => null()
     CHARACTER ( CPLLEN ),   DIMENSION( : ), POINTER :: ons_corr => null()
     CHARACTER ( CPLLEN ),   DIMENSION( : ), POINTER :: so_corr => null()
     CHARACTER ( CPLLEN ),   DIMENSION( : ), POINTER :: intra => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: energy => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: scaling => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: fac_P => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: fac_S => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: fac_Q => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: fac_I => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: l_I => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: fac_so => null()
     REAL ( dp ),            DIMENSION( : ), POINTER :: fac_C => null()
     REAL ( dp )                                     :: fac_O
     REAL ( dp )                                     :: l_O
     REAL ( dp )                                     :: dist_ref
     REAL ( dp )                                     :: delta_d
     CHARACTER ( ATOMLEN ),  DIMENSION( 2 )          :: name
     INTEGER                                         :: order
          
  END TYPE pair_coupling
  
  !___________________________________________________________________________

  ! describes a parent material completely

  TYPE parent_data
     
     INTEGER                                         :: n_ion
     TYPE (ion_orbit),       DIMENSION( : ), POINTER :: ion    => null() 
     INTEGER                                         :: n_pair 
     TYPE (pair_coupling),   DIMENSION( : ), POINTER :: pair => null()
     REAL ( dp )                                     :: content
     REAL ( dp ),            DIMENSION( 3 )          :: lat_par
     REAL ( dp ),            DIMENSION( 3 )          :: ref_axis
     REAL ( dp )                                     :: u_par
     REAL ( dp )                                     :: int_str_par
     REAL ( dp )                                     :: e_v
     CHARACTER ( LEN = 100 )                         :: data_file
     CHARACTER ( LEN = 20 )                          :: name
     CHARACTER ( LEN = 20 )                          :: type
     CHARACTER ( LEN = 20 )                          :: cry     
     CHARACTER ( LEN = 10 )                          :: scheme

  END TYPE parent_data

  !___________________________________________________________________________

  ! describes completely one material

  TYPE material_data
     
     INTEGER,                DIMENSION( : ), POINTER :: nr_parents => null()
     TYPE (parent_data),     DIMENSION( : ), POINTER :: parent => null()
     TYPE (ion_orbit),       DIMENSION( : ), POINTER :: ion => null()  ! Alloy averaged   
     REAL ( dp ),            DIMENSION( : ), POINTER :: bowing => null()
     REAL ( dp ),            DIMENSION( 2 )          :: gap_corr
     REAL ( dp ),            DIMENSION( 3 )          :: lat_par
     REAL ( dp ),            DIMENSION( 3 )          :: ref_axis
     REAL ( dp )                                     :: u_par
     REAL ( dp )                                     :: int_str_par   
     INTEGER                                         :: max_order
     CHARACTER ( LEN = 20 )                          :: cry
     CHARACTER ( LEN = 20 )                          :: type
     CHARACTER ( LEN = 20 )                          :: name
     CHARACTER ( LEN = 10 )                          :: scheme
     LOGICAL                                         :: alloy
     LOGICAL                                         :: alloy_random
     LOGICAL                                         :: extra_cpl
     LOGICAL                                         :: average_cat 

  END TYPE material_data

  !___________________________________________________________________________
  
  ! defines the atomic structure read for the gen file

  TYPE TStructure
     
     CHARACTER(MST)                       :: gen_filename
     INTEGER                              :: n_atoms
     CHARACTER(1)                         :: C_S     
     INTEGER                              :: nr_atomtypes
     CHARACTER(TYPELEN), DIMENSION(:), POINTER  :: atomtypes => null()
     INTEGER                              :: n_danglings
     INTEGER                              :: hydro_type
     REAL ( dp ), DIMENSION(:,:), POINTER :: coord   => null()  
     INTEGER, DIMENSION(:,:), POINTER     :: nn_list => null()
     INTEGER, DIMENSION(:), POINTER       :: mat    => null()           
     INTEGER, DIMENSION(:), POINTER       :: atomtype => null()
     INTEGER, DIMENSION(:), POINTER       :: indexa => null()
     INTEGER, DIMENSION(:), POINTER       :: inv_indexa => null()
     REAL ( dp ), DIMENSION(3)            :: start
     REAL ( dp ), DIMENSION(3,3)          :: prim

  END TYPE TStructure

contains

  !===========================================================================
  !
  ! Subroutine allocate_ion( mat_data, n_ion )
  !
  !===========================================================================
  !
  ! * input arguments :
  !
  !   --> n_ion : number of ions in the list
  !
  !---------------------------------------------------------------------------
  !
  ! * input - output arguments :
  !
  !   --> mat_data : data list to be allocated
  !
  !===========================================================================

  SUBROUTINE allocate_mat_ion( mat_data, n_ion )

    !=========================================================================
    ! Input - output  arguments
    !=========================================================================

    TYPE(material_data) :: mat_data

    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: n_ion

    !=========================================================================
    ! Local variables
    !=========================================================================

    INTEGER :: i_ion, err

    !=========================================================================

    ALLOCATE( mat_data%ion( n_ion ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR)
       
    ! The pointer elements of the structure must be initialized
    DO i_ion = 1, n_ion
       call nullify_ion( mat_data%ion(i_ion) )    
    END DO

  END SUBROUTINE allocate_mat_ion
  !=========================================================================

  SUBROUTINE allocate_parent_ion( parent, n_ion )
    
    TYPE(parent_data) :: parent
    INTEGER, INTENT( IN ) :: n_ion    

    !=========================================================================
    ! Local variables
    !=========================================================================
    INTEGER :: i_ion, err


    ALLOCATE( parent%ion(n_ion), STAT = err ) 
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR)

    ! The pointer elements of the structure must be initialized
    DO i_ion = 1, n_ion
       call nullify_ion( parent%ion(i_ion) )    
    END DO


  END SUBROUTINE allocate_parent_ion

  !=========================================================================
  subroutine nullify_ion(ion)

    TYPE(ion_orbit) :: ion
   
    NULLIFY(ion%energy, ion%so_energy, ion%state, ion%so_state, ion%ind_ref, ion%ind_ref_inv)

  end subroutine nullify_ion
  
  !========================================================================= 
  SUBROUTINE allocate_ion(ion, n_state, n_so)

    TYPE(ion_orbit) :: ion
    INTEGER, INTENT(IN) :: n_state, n_so

    ! LOCALS:
    INTEGER :: err

    ALLOCATE( ion%energy( n_state ), STAT = err )

    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( ion%state( n_state ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( ion%ind_ref( n_state ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( ion%ind_ref_inv( n_ref_states ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

    ion%energy = 0.0D0
    ion%state = ' '
    ion%ind_ref = 0

    ALLOCATE( ion%so_state(n_so), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( ion%so_energy(n_so), STAT = err )

    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 
              
    ion%so_state = ' '
    ion%so_energy = 0.0D0

  END SUBROUTINE allocate_ion

  !=========================================================================
  SUBROUTINE allocate_parent_pair(parent, n_pair)

    TYPE(parent_data) :: parent
    INTEGER, INTENT(IN) :: n_pair

    INTEGER :: i_pair, err

    ALLOCATE(parent%pair( n_pair ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 
          
    ! The pointer elements of the structure must be initialized
          
    DO i_pair = 1, n_pair
       
       NULLIFY(parent%pair(i_pair)%energy, &
            parent%pair(i_pair)%cpl,&
            parent%pair(i_pair)%scaling )
       
    END DO
         

  END SUBROUTINE allocate_parent_pair
  !=========================================================================

  SUBROUTINE allocate_pair(pair, n_cpl)

    TYPE(pair_coupling) :: pair
    INTEGER, INTENT(IN) :: n_cpl

    INTEGER :: err


    ALLOCATE( pair%cpl( n_cpl ), STAT = err )

    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    pair%cpl = ' '          
        
    !......................................................................

    ALLOCATE( pair%fac_P( n_cpl ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%fac_P = 0.0D0
           
    !......................................................................

    ALLOCATE( pair%fac_S( n_cpl ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%fac_S = 0.0D0
           
    !......................................................................

    ALLOCATE( pair%fac_Q( n_cpl ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%fac_Q = 0.0D0
              
    !......................................................................

    ALLOCATE( pair%energy( n_cpl ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%energy = 0.0D0
       
    !......................................................................
    

    ALLOCATE( pair%scaling( n_cpl ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%scaling = 0.0D0
   
  END SUBROUTINE allocate_pair

  !=========================================================================
  SUBROUTINE allocate_ons_corr(pair, n_I)

    TYPE(pair_coupling) :: pair
    INTEGER, INTENT(IN) :: n_I

    INTEGER :: err


    ALLOCATE( pair%ons_corr( n_I ), STAT = err )

    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    pair%ons_corr = ' '          
        
    !......................................................................

    ALLOCATE( pair%fac_I( n_I ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%fac_I = 0.0D0
           
    !......................................................................

    ALLOCATE( pair%l_I( n_I ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%l_I = 0.0D0
           
  END SUBROUTINE allocate_ons_corr

  !=========================================================================
  SUBROUTINE allocate_so_corr(pair, n_so)

    TYPE(pair_coupling) :: pair
    INTEGER, INTENT(IN) :: n_so

    INTEGER :: err

    ALLOCATE( pair%so_corr( n_so ), STAT = err )

    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    pair%so_corr = ' '          
        
    !......................................................................

    ALLOCATE( pair%fac_so( n_so ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%fac_so = 0.0D0
           
  END SUBROUTINE allocate_so_corr

  !=========================================================================
  SUBROUTINE allocate_intra(pair, n_intra)

    TYPE(pair_coupling) :: pair
    INTEGER, INTENT(IN) :: n_intra

    INTEGER :: err

    ALLOCATE( pair%intra( n_intra ), STAT = err )

    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    pair%intra = ' '          
        
    !......................................................................

    ALLOCATE( pair%fac_C( n_intra ), STAT = err )

    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
       
    pair%fac_C = 0.0D0
           
  END SUBROUTINE allocate_intra

  !=========================================================================
  SUBROUTINE nullify_mat(mat_data)
    type(material_data) :: mat_data

    NULLIFY(mat_data%nr_parents)
    NULLIFY(mat_data%bowing)

  END SUBROUTINE NULLIFY_MAT


  !=========================================================================

  subroutine create_basis(basis,n_basis)
    
    type(ion_basis) :: basis
    integer :: n_basis
    integer :: err
    
    basis%n_basis = n_basis
    basis%n_dg_bond = 0

    !ALLOCATE( basis%name( n_basis ), STAT = err )
    !IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%coord( 3, n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%strain( 3, n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%mat( n_basis ), STAT = err )
    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%ion( n_basis ), STAT = err )
    IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%type( n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%dg_coord( 3, n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%n_dg( n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%n_st( n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%valence( n_basis ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%start( 3 ), STAT = err )
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%prim( 3, 3 ), STAT = err )    
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

    ALLOCATE( basis%rec_latt( 3, 3 ), STAT = err )    
    IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 
    
    !basis%name = ''
    basis%coord = 0.d0
    basis%strain = 0.d0
    basis%mat = 0
    basis%ion = 0
    basis%type = 0
    basis%dg_coord = 0.d0
    basis%n_dg = 0
    basis%n_st = 0
    basis%start = 0.d0
    basis%prim = 0.d0
    basis%rec_latt = 0.d0    

  end subroutine create_basis
  ! -----------------------------------------------------------------
  subroutine nullify_basis(basis)
    
     type(ion_basis) :: basis

    !NULLIFY( basis%name )
    NULLIFY( basis%coord )
    NULLIFY( basis%strain )
    NULLIFY( basis%mat )
    NULLIFY( basis%ion )
    NULLIFY( basis%type )
    NULLIFY( basis%dg_coord )
    NULLIFY( basis%n_dg )
    NULLIFY( basis%n_st )
    NULLIFY( basis%valence )

    NULLIFY( basis%start )
    NULLIFY( basis%prim )      
    NULLIFY( basis%rec_latt )      

  end subroutine nullify_basis


  ! -----------------------------------------------------------------

  subroutine write_basis(basis,level)
 
     type(ion_basis) :: basis    
     integer :: level

     integer :: i

     write(*,*) 'n_basis=',basis%n_basis
     write(*,*) 'n_dg=',basis%n_dg_bond

     if (level.gt.0) then
        write(*,*) 'ion     material      valence'
        do i=1,basis%n_basis + basis%n_dg_bond
        write(*,*) basis%atomtypes(basis%type(i)),basis%mat(i),basis%valence(i)
        enddo
     end if

     write(*,*) 'periodic= ',basis%periodic
     write(*,*) 'directions= ',basis%periodic_BC
     write(*,*) 'start=',basis%start
     write(*,*) 'primitive vectors='
     write(*,*) basis%prim(1,1:3)
     write(*,*) basis%prim(2,1:3)
     write(*,*) basis%prim(3,1:3)
     write(*,*) 'reciprocal vectors='
     write(*,*) basis%rec_latt(1,1:3)
     write(*,*) basis%rec_latt(2,1:3)
     write(*,*) basis%rec_latt(3,1:3)

   end subroutine write_basis
   ! ------------------------------------------------------
   
   subroutine destroy_basis(basis)
    
     type(ion_basis) :: basis
     integer :: err
     
     err = 0


     if(associated(basis%atomtypes)) DEALLOCATE( basis%atomtypes, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%coord)) DEALLOCATE( basis%coord, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%strain)) DEALLOCATE( basis%strain, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%mat)) DEALLOCATE( basis%mat, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%ion)) DEALLOCATE( basis%ion, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%type)) DEALLOCATE( basis%type, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%dg_coord)) DEALLOCATE( basis%dg_coord, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%n_dg)) DEALLOCATE( basis%n_dg, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%n_st)) DEALLOCATE( basis%n_st, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%start)) DEALLOCATE( basis%start, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%prim)) DEALLOCATE( basis%prim, STAT = err )
     IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%rec_latt)) DEALLOCATE( basis%rec_latt, STAT = err )
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 

     if(associated(basis%valence)) DEALLOCATE( basis%valence, STAT = err )     
     IF ( err .NE. 0 )  call throw_init_exception(ERR_ALLOC_ERR) 
     
 
   end subroutine destroy_basis
   ! ------------------------------------------------------

   subroutine write_materials(nr_materials,mat)
     
     integer :: nr_materials
     type(material_data), DIMENSION(:), POINTER :: mat    
     
     integer :: i,k,l
     
     write(*,*) '=========================================='
     do i = 1,nr_materials
        
        write(*,*) 'Material',i
        write(*,*) 'Material Name= ',mat(i)%name
        write(*,*) 'Material type= ',mat(i)%type
        write(*,*) 'Material cystal= ',mat(i)%cry
        write(*,*) 'Is random alloy=',mat(i)%alloy_random
        write(*,*) 'n_parents=',sum(mat(i)%nr_parents)      
        
        do k=1, sum(mat(i)%nr_parents)

           write(*,*) '======================================='  
           write(*,*) 'parent',k,'name= ',mat(i)%parent(k)%name         
           write(*,*) 'parent',k,'type= ',mat(i)%parent(k)%type
           write(*,*) 'parent',k,'crystal= ',mat(i)%parent(k)%cry
           write(*,*) 'database=',trim(mat(i)%parent(k)%data_file)  
           write(*,*) 'n. ions=',mat(i)%parent(k)%n_ion
           write(*,*) 'n. pairs=',mat(i)%parent(k)%n_pair
           write(*,*) 'content=',mat(i)%parent(k)%content
           
           !call write_ion_orbits(mat(i)%parent(k)%ion)
           
           !call write_ion_couplings(mat(i)%parent(k)%pair)

        enddo

        write(*,*) 'max n.n. order=',mat(i)%max_order        
        write(*,*) '======================================='  
        
        if(.NOT.mat(i)%alloy_random.AND.mat(i)%alloy) then
           
           !write(*,*) 'averaged ion properties:'
           
           !call write_ion_orbits(mat(i)%ion)
           
        end if
        
        
     end do
      
   end subroutine write_materials
   ! ---------------------------------------------------------- 
   subroutine write_ion_orbits(ion)
     
     TYPE(ion_orbit), DIMENSION(:), POINTER :: ion

     INTEGER :: i,k

     
     do i=1,SIZE(ion)

        write(*,*) '----------------------------------------'
        write(*,*) 'ion= ',ion(i)%name
        write(*,*) 'N. of orbitals=',SIZE(ion(i)%state)
        
        do k=1,SIZE(ion(i)%state)
           
           write(*,'(i3,2x,a6,f10.6)') ion(i)%ind_ref(k),&
                                 ion(i)%state(k),ion(i)%energy(k)

        enddo

        write(*,*) 'spin-orbit couplings'

        do k=1,SIZE(ion(i)%so_state)
           
           write(*,*) ion(i)%so_state(k),ion(i)%so_energy(k)

        enddo        

        write(*,*) 'offset=',ion(i)%offset

        write(*,*) 'strain band coupling=',ion(i)%b_d

     end do

   end subroutine write_ion_orbits
    ! ---------------------------------------------------------- 
   subroutine write_ion_couplings(pair)

     TYPE(pair_coupling), DIMENSION(:), POINTER :: pair

     INTEGER :: i,k

     !CHARACTER ( LEN = 5 ),  DIMENSION( : ), POINTER :: cpl
     !REAL ( dp ),            DIMENSION( : ), POINTER :: energy
     !REAL ( dp ),            DIMENSION( : ), POINTER :: scaling
     !REAL ( dp )                                     :: dist_ref
     !CHARACTER ( LEN = 2 ),  DIMENSION( 2 )          :: name
     !INTEGER,                DIMENSION( 2 )          :: mat
     !INTEGER                                         :: order
     
     do i=1,SIZE(pair)

        write(*,*) '........................................'
        write(*,*) 'pair name=',pair(i)%name
        write(*,*) 'pair order=',pair(i)%order
        write(*,*) 'reference distance=',pair(i)%dist_ref
        write(*,*) 'N. of pair couplings=',SIZE(pair(i)%cpl)

        
        do k=1,SIZE(pair(i)%cpl)
           
           write(*,'(i3,2x,a6,f10.6, f10.6)') k, &
                pair(i)%cpl(k),pair(i)%energy(k),pair(i)%scaling(k)

        enddo

     end do       

   end subroutine write_ion_couplings
   ! ---------------------------------------------------------- 

   subroutine destroy_ion_orbit(ion)

     TYPE(ion_orbit) :: ion

     !write(*,*) 'loc:', %LOC(ion%energy)

     if(associated(ion%energy))    DEALLOCATE(ion%energy)
     if(associated(ion%so_energy)) DEALLOCATE(ion%so_energy)
     if(associated(ion%state))     DEALLOCATE(ion%state)
     if(associated(ion%so_state))  DEALLOCATE(ion%so_state)
     if(associated(ion%ind_ref))   DEALLOCATE(ion%ind_ref)

   end subroutine destroy_ion_orbit
   ! ---------------------------------------------------------- 
   
   subroutine destroy_pair_coupling(pair)
   
     TYPE(pair_coupling) :: pair

     if(associated(pair%cpl))     DEALLOCATE( pair%cpl )
     if(associated(pair%energy))  DEALLOCATE( pair%energy )
     if(associated(pair%scaling)) DEALLOCATE( pair%scaling )
     
     ! Deallocate Tan scheme arrays (memory leak fix)
     if(associated(pair%fac_I))   DEALLOCATE( pair%fac_I )
     if(associated(pair%l_I))     DEALLOCATE( pair%l_I )
     if(associated(pair%fac_so))  DEALLOCATE( pair%fac_so )
     if(associated(pair%ons_corr)) DEALLOCATE( pair%ons_corr )
     if(associated(pair%so_corr)) DEALLOCATE( pair%so_corr )
     if(associated(pair%intra))   DEALLOCATE( pair%intra )
     if(associated(pair%fac_P))   DEALLOCATE( pair%fac_P )
     if(associated(pair%fac_S))   DEALLOCATE( pair%fac_S )
     if(associated(pair%fac_Q))   DEALLOCATE( pair%fac_Q )
     if(associated(pair%fac_C))   DEALLOCATE( pair%fac_C )

   end subroutine destroy_pair_coupling
   ! ---------------------------------------------------------- 

   subroutine destroy_parent(parent)
     
     type(parent_data) :: parent
     INTEGER :: i

    

     if(associated(parent%pair)) then
        DO i = 1, SIZE(parent%pair)   
           !write(*,*) 'destroy parent%pair',i,'/',SIZE(parent%pair)
           CALL destroy_pair_coupling( parent%pair(i) )
        END DO
        !write(*,*) 'deallocate parent%pair'
        DEALLOCATE( parent%pair )
     end if


     if(associated(parent%ion)) then
        DO i = 1, SIZE(parent%ion)
           call destroy_ion_orbit(parent%ion(i))
        END DO
        DEALLOCATE( parent%ion )
     end if


   end subroutine destroy_parent

   ! ---------------------------------------------------------- 

   subroutine destroy_material(mat_data)

     TYPE(material_data) :: mat_data
     INTEGER :: i
     

     if(associated(mat_data%parent)) then
        DO i = 1, SIZE(mat_data%parent)
           !print*, '(debug) destroy parent',i
           call destroy_parent( mat_data%parent(i) ) 
        END DO
        DEALLOCATE( mat_data%parent )
     end if

     !................................................................
     ! We must be careful here because for simple and binary (1 parent)
     ! the mat_data%ion is obtained by type-copy.
     ! In this case (nr_parents == 1) the deallocation must be avoided
     !................................................................
     if (associated(mat_data%ion)) then
        if(sum(mat_data%nr_parents).ne.1) then 
           DO i = 1, SIZE(mat_data%ion)
              !print*, '(debug) destroy ion',i
              call destroy_ion_orbit( mat_data%ion(i) )
           END DO
        end if
        DEALLOCATE( mat_data%ion ) 
     end if
      
     !print*, '(debug) deallocate nr_parents'
     if(associated(mat_data%nr_parents)) then
        DEALLOCATE( mat_data%nr_parents)
     end if

     !print*, '(debug) deallocate bowing'
     if(associated(mat_data%bowing))  DEALLOCATE( mat_data%bowing )

   end subroutine destroy_material

   ! ---------------------------------------------------------- 
   

   ! ---------------------------------------------------------- 
   subroutine destroy_structure(str)

     TYPE(TStructure) :: str

     if(associated(str%atomtypes)) deallocate(str%atomtypes)
     if(associated(str%coord)) deallocate(str%coord)     
     if(associated(str%nn_list)) deallocate(str%nn_list)     
     if(associated(str%mat)) deallocate(str%mat)     
     if(associated(str%atomtype)) deallocate(str%atomtype)     
     if(associated(str%indexa)) deallocate(str%indexa)     
     if(associated(str%inv_indexa)) deallocate(str%inv_indexa)     

   end subroutine destroy_structure


    ! ---------------------------------------------------------- 

   subroutine nullify_structure(str)

     TYPE(TStructure) :: str

     nullify(str%atomtypes)
     nullify(str%coord)     
     nullify(str%nn_list)     
     nullify(str%mat)     
     nullify(str%atomtype)     
     nullify(str%indexa)     
     nullify(str%inv_indexa)     

   end subroutine nullify_structure

   subroutine bcast_struct(str)
     type(Tstructure) :: str
    
     
     !CHARACTER(MST)                       :: gen_filename
     !INTEGER                              :: n_atoms
     !CHARACTER(1)                         :: C_S     
     !INTEGER                              :: nr_atomtypes
     !CHARACTER(2), DIMENSION(:), POINTER  :: atomtypes
     !INTEGER                              :: n_danglings
     !INTEGER                              :: hydro_type
     !REAL ( dp ), DIMENSION(:,:), POINTER :: coord     
     !INTEGER, DIMENSION(:,:), POINTER     :: nn_list
     !INTEGER, DIMENSION(:), POINTER       :: mat              
     !INTEGER, DIMENSION(:), POINTER       :: atomtype
     !INTEGER, DIMENSION(:), POINTER       :: indexa
     !INTEGER, DIMENSION(:), POINTER       :: inv_indexa
     !REAL ( dp ), DIMENSION(3)            :: start
     !REAL ( dp ), DIMENSION(3,3)          :: prim

     ! Bcast fisrt size 
     ! then ALLOCATE
     ! then send data

   end subroutine bcast_struct

   subroutine bcast_matdata(matdata)
     type(material_data) :: matdata
    

     !INTEGER                                         :: nr_parents
     !TYPE (parent_data),     DIMENSION( : ), POINTER :: parent
     !TYPE (ion_orbit),       DIMENSION( : ), POINTER :: ion  ! Alloy averaged   
     !REAL ( dp ),            DIMENSION( : ), POINTER :: bowing
     !REAL ( dp ),            DIMENSION( 2 )          :: gap_corr
     !REAL ( dp ),            DIMENSION( 3 )          :: lat_par
     !REAL ( dp )                                     :: u_par
     !REAL ( dp )                                     :: int_str_par
     !INTEGER                                         :: n_pair
     !INTEGER                                         :: max_order
     !CHARACTER ( LEN = 20 )                          :: cry
     !CHARACTER ( LEN = 20 )                          :: type
     !CHARACTER ( LEN = 20 )                          :: name
     !LOGICAL                                         :: alloy
     !LOGICAL                                         :: alloy_random
     !LOGICAL                                         :: extra_cpl
   
   end subroutine bcast_matdata  

end module type_defs
