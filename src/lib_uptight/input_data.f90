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
!=============================================================================
!
!        Module "input_data" - (c) Jerome Gleize - August 2003
!
!=============================================================================
!
! contains all routines used to read the data needed for the ETB simulation.
!
! => read_parse
! => read_data
! => read_ETB_orbitals
! => read_ETB_couplings
! => read_lattice
!
!=============================================================================
!
! All the data are collected in a parse file, which contains :
!
! => the geometry input
! => the chemical input 
! => the simulation input
!
!=============================================================================

MODULE input_data
  
  !===========================================================================
  !
  !  parameters : latt_max, latt_dim, kvectorsdim_D, max_ind_D, prec
  !
  !  and all types and global variables
  !
  !===========================================================================

  USE precision
  USE globals, only : LST
  USE type_defs

  USE errors 
  USE input_output
  USE exceptions

  !===========================================================================
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_data

  !===========================================================================
  
CONTAINS

  !===========================================================================
  !
  ! Subroutine read_data 
  !
  !===========================================================================
  !
  ! Read all data needed to run the ETB simulation.
  !
  ! NOTE : all input and output are done via global variables
  !        declared in the parameter module "parameter.f90"
  !
  !===========================================================================
    
  SUBROUTINE read_data2(n_mat, mat_data, work_path, database_path)

    INTEGER :: n_mat
    TYPE(material_data), DIMENSION(:), POINTER :: mat_data
    CHARACTER(LST) :: work_path, database_path

    !=========================================================================
    ! Local variables
    !=========================================================================

    INTEGER :: i_mat, i_par
    INTEGER :: i_ion, n_ion
    INTEGER :: i_pair, n_pair, i_k, i_state, err, file_num, nr_parents
    LOGICAL :: error

    CHARACTER (LEN = 500) :: file_name
    CHARACTER (LEN = 500) :: read_format

    !=========================================================================
    ! Read from the data files
    !=========================================================================

    !WRITE( *, '( /, "***", //, "READ THE MATERIAL FILES", / )' )
    
    !=========================================================================
    ! Loop on the materials
    !=========================================================================
    nr_parents = sum( mat_data(i_mat)%nr_parents )

    DO i_mat = 1, n_mat
       
       !----------------------------------------------------------------------
       ! Loop on the parents for the current material
       !----------------------------------------------------------------------
       
       DO i_par = 1, nr_parents 
          
          !-------------------------------------------------------------------
          ! Open the current data file, with screen ouput and error handling
          !-------------------------------------------------------------------

          file_name = TRIM( work_path ) // &
               TRIM( mat_data( i_mat )%parent( i_par )%data_file )
          
          INQUIRE(FILE= file_name, EXIST= error)
          IF (error) THEN
            CALL open_file( file_name, file_num, "read", format_flag = .TRUE. )
          ELSE
            file_name = TRIM( database_path ) // &
               TRIM( mat_data( i_mat )%parent( i_par )%data_file )
            CALL open_file( file_name, file_num, "read", format_flag = .TRUE. )
          ENDIF

          !-------------------------------------------------------------------
          ! Read the type of parent material
          !-------------------------------------------------------------------
          
          read_format = ' '
          read_format = 'material type = '
          
          CALL label_in_file( file_num, read_format, error, 'no' )
          IF ( error ) CALL read_errors( 'missing', file_name, read_format )
          
          READ ( file_num, *, IOSTAT = err ) &
               mat_data( i_mat )%parent( i_par )%type
          IF ( err .NE. 0 ) CALL read_errors( 'not character string', &
               file_name, read_format )

          !-------------------------------------------------------------------
          ! Read the crystal structure of parent material
          !-------------------------------------------------------------------
          
          read_format = ' '
          read_format = 'basis type = '
       
          CALL label_in_file( file_num, read_format, error, 'no' )
          IF ( error ) CALL read_errors( 'missing', file_name, read_format )
          
          READ ( file_num, *, IOSTAT = err ) &
               mat_data( i_mat )%parent( i_par )%cry
          IF ( err .NE. 0 ) CALL read_errors( 'not character string', &
               file_name, read_format )
          
          !-------------------------------------------------------------------
          ! Parent type testing : get the number of ions in the current parent
          !-------------------------------------------------------------------
          n_ion = int_in_file( file_name, file_num, 'number of ions' ) 

          mat_data( i_mat )%parent( i_par )%n_ion = n_ion

          !-------------------------------------------------------------------
          ! Allocate each element of the parent's ion list
          !-------------------------------------------------------------------
          CALL allocate_parent_ion(mat_data(i_mat)%parent(i_par), n_ion)

          !-------------------------------------------------------------------
          ! Read the number of distinct pairs in the current parent
          !-------------------------------------------------------------------
          
          n_pair = int_in_file( file_name, file_num, 'number of pairs' )
          
          mat_data( i_mat )%parent( i_par )%n_pair = n_pair 
          
          IF ( n_pair .NE. mat_data( i_mat )%parent(1)%n_pair ) &
               call throw_init_exception(ERR_DB_PAIR)                
          
          !-------------------------------------------------------------------
          ! Allocate each element of pair
          !-------------------------------------------------------------------
          
          CALL allocate_parent_pair( mat_data(i_mat)%parent(i_par), n_pair)


          !-------------------------------------------------------------------
          ! Read the elements of the ion and pair data lists
          !-------------------------------------------------------------------

          !write(*,*) 'Reading ETB orbitals ...' 

          CALL read_ETB_orbitals(file_num,&
                                   mat_data( i_mat )%parent( i_par )%ion )
                
          !-------------------------------------------------------------------
          ! Here we add the valence band absolute values to offset         
          !-------------------------------------------------------------------
          do i_ion=1, n_ion 

              !...............................................................
              ! Add the offset to the ETB onsite energy
              !...............................................................
              mat_data(i_mat)%parent(i_par)%ion(i_ion)%offset = &
                        mat_data(i_mat)%parent(i_par)%ion(i_ion)%offset + &
                        mat_data(i_mat)%parent(i_par)%e_v 

              mat_data(i_mat)%parent(i_par)%ion(i_ion)%energy = &
                         mat_data(i_mat)%parent(i_par)%ion(i_ion)%energy +&
                         mat_data(i_mat)%parent(i_par)%ion(i_ion)%offset
          

          enddo                          

          !write(*,*) 'Reading ETB couplings ...'     

          !!CALL read_ETB_couplings( file_num, i_mat, i_par, mat_data )
          
          !IF ( TRIM( mat_data( i_mat )%parent( i_parent )%type ) .NE. &
          !     'extra_TB_couplings' ) THEN
          
          CLOSE( file_num )
          
          !-------------------------------------------------------------------
          ! End of parent material loop
          !-------------------------------------------------------------------
          
       END DO

       !-------------------------------------------------------------
       ! CHECKS CONSISTENCY IN PARENT LATTICE STRUCTURES 
       !-------------------------------------------------------------

       DO i_par = 1, nr_parents 

          IF (TRIM(mat_data( i_mat )%parent( i_par )%cry) .ne. &
               TRIM(mat_data( i_mat )%parent( 1 )%cry) ) THEN

             WRITE(*,*) "WARNING: incompatible parent lattices"

          END IF

       END DO

       mat_data(i_mat)%cry  = mat_data( i_mat )%parent( 1 )%cry

    END DO

  END SUBROUTINE read_data2              

  !===========================================================================
  !
  ! Subroutine read_data 
  !
  !===========================================================================
  !
  ! Read all data needed to run the ETB simulation.
  !
  ! NOTE : all input and output are done via global variables
  !        declared in the parameter module "parameter.f90"
  !
  !===========================================================================
    
  SUBROUTINE read_data(mat_data, work_path, database_path)

    TYPE(material_data) :: mat_data
    CHARACTER(LST) :: work_path, database_path

    !=========================================================================
    ! Local variables
    !=========================================================================

    INTEGER :: i_mat, i_par
    INTEGER :: i_ion, n_ion
    INTEGER :: i_pair, n_pair, i_k, i_state, nr_parents
    INTEGER :: err, file_num, n_path
    LOGICAL :: error

    CHARACTER (LEN = 500) :: file_name
    CHARACTER (LEN = 500) :: read_format
    CHARACTER (LEN = 10) :: scheme
    CHARACTER (LEN = 400), allocatable :: strarray(:)

    !=========================================================================
    ! Read from the data files
    !=========================================================================
 
    ! database_path can have the unix format path1:path2:...
    n_path = count(transfer(trim(database_path), 'a', len(trim(database_path))) == ':')
    n_path = n_path + 1
    allocate(strarray(n_path))

    i_mat = 1
    i_state = 1

    DO i_k = 1, len(trim(database_path))
      IF (database_path(i_k:i_k) == ':') THEN
        strarray(i_mat) = database_path(i_k-i_state+1:i_k-1)
        i_mat = i_mat + 1
        i_state = 1
      ELSE
        i_state = i_state + 1
      ENDIF
    ENDDO
    strarray(i_mat) = database_path(i_k-i_state+1:i_k)

    
       
    !----------------------------------------------------------------------
    ! Loop on the parents for the current material
    !----------------------------------------------------------------------
    nr_parents = sum(mat_data%nr_parents)
    
    DO i_par = 1, nr_parents 
       
       !-------------------------------------------------------------------
       ! Open the current data file, with screen ouput and error handling
       !-------------------------------------------------------------------


       DO i_k = 1, n_path
         
         file_name = TRIM( strarray(i_k) ) // &
            TRIM( mat_data%parent( i_par )%data_file )
         
         INQUIRE(FILE= file_name, EXIST= error)
         IF (error) THEN
           CALL open_file( file_name, file_num, "read", format_flag = .TRUE. )
           EXIT
         ENDIF

       ENDDO
         
       IF (.not. error) THEN
         CALL throw_init_exception(ERR_FILE_OPEN)
       ENDIF
       
       
       !-------------------------------------------------------------------
       ! Read the scheme of .etb file ('tan' or 'jancu')
       !-------------------------------------------------------------------
       
       read_format = ' '
       read_format = 'scheme = '
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       !If say nothing about scheme, 'jancu' is default
       IF ( error ) mat_data%parent( i_par )%scheme = 'jancu'
       
       IF (.NOT. error) THEN
         READ ( file_num, *, IOSTAT = err ) &
              mat_data%parent( i_par )%scheme
         IF ( err .NE. 0 ) CALL read_errors( 'not character string', &
                                            file_name, read_format )
       ENDIF
       scheme = TRIM(mat_data%parent( i_par )%scheme)
       
       IF (scheme .NE. 'jancu' .AND. scheme .NE. 'tan') THEN
         CALL throw_init_exception(ERR_HAM_UNSCHE)
       ENDIF
       
       !-------------------------------------------------------------------
       ! Read the type of parent material
       !-------------------------------------------------------------------
       
       read_format = ' '
       read_format = 'material type = '
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', file_name, read_format )
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_par )%type
       IF ( err .NE. 0 ) CALL read_errors( 'not character string', &
            file_name, read_format )

       !-------------------------------------------------------------------
       ! Read the crystal structure of parent material
       !-------------------------------------------------------------------
       
       read_format = ' '
       read_format = 'basis type = '
    
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', file_name, read_format )
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_par )%cry
       IF ( err .NE. 0 ) CALL read_errors( 'not character string', &
            file_name, read_format )

       !-------------------------------------------------------------------
       ! Read the reference axis of parent material
       !-------------------------------------------------------------------
       
       read_format = ' '
       read_format = 'reference axis = '

       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( .not. error ) THEN

         READ ( file_num, *, IOSTAT = err ) &
              mat_data%parent( i_par )%ref_axis
         IF ( err .NE. 0 ) CALL read_errors( 'not a 3-vector', &
              file_name, read_format )
       ELSE
         ! assign defaults, if reference axis is not present
         mat_data%parent( i_par )%ref_axis = (/ 0.d0, 0.d0, 1.d0 /)

         SELECT CASE ( mat_data%parent( i_par )%cry )

         CASE( 'wurtzite' )

            mat_data%parent( i_par )%ref_axis = (/ 1.d0, 1.d0, 1.d0 /)

         CASE( 'zinc-blende' )

            mat_data%parent( i_par )%ref_axis = (/ 0.d0, 0.d0, 1.d0 /)

         END SELECT
       END IF

       !-------------------------------------------------------------------
       ! Parent type testing : get the number of ions in the current parent
       !-------------------------------------------------------------------
       n_ion = int_in_file( file_name, file_num, 'number of ions' ) 

       mat_data%parent( i_par )%n_ion = n_ion

       !-------------------------------------------------------------------
       ! Allocate each element of the parent's ion list
       !-------------------------------------------------------------------
       CALL allocate_parent_ion(mat_data%parent(i_par), n_ion)

       !-------------------------------------------------------------------
       ! Read the number of distinct pairs in the current parent
       !-------------------------------------------------------------------
       
       n_pair = int_in_file( file_name, file_num, 'number of pairs' )
       
       mat_data%parent( i_par )%n_pair = n_pair 
       
       IF ( n_pair .NE. mat_data%parent(1)%n_pair ) &
            call throw_init_exception(ERR_DB_PAIR)                
       
       !-------------------------------------------------------------------
       ! Allocate each element of pair
       !-------------------------------------------------------------------
       
       CALL allocate_parent_pair( mat_data%parent(i_par), n_pair)


       !-------------------------------------------------------------------
       ! Read the elements of the ion and pair data lists
       !-------------------------------------------------------------------

       !write(*,*) 'Reading ETB orbitals ...' 

       CALL read_ETB_orbitals(file_num,&
                                mat_data%parent(i_par)%ion )
                                
       !-------------------------------------------------------------------
       ! Here we add the valence band absolute values to offset         
       !-------------------------------------------------------------------
       do i_ion=1, n_ion 

           !...............................................................
           ! Add the offset to the ETB onsite energy
           !...............................................................
           mat_data%parent(i_par)%ion(i_ion)%offset = &
                     mat_data%parent(i_par)%ion(i_ion)%offset + &
                     mat_data%parent(i_par)%e_v 
           ! Not add offset to onsite energies here,
           ! will be handled later in TB_ham.f90
           !mat_data%parent(i_par)%ion(i_ion)%energy = &
           !           mat_data%parent(i_par)%ion(i_ion)%energy +&
           !           mat_data%parent(i_par)%ion(i_ion)%offset
       

       enddo                          

       !write(*,*) 'Reading ETB couplings ...'     

       CALL read_ETB_couplings( file_num, mat_data, i_par, scheme )
       
       IF (scheme .EQ. 'tan') THEN 
         CALL read_ETB_ons_corr( file_num, mat_data, i_par )

         CALL read_ETB_intra( file_num, mat_data, i_par )
       ENDIF
       
       !IF ( TRIM( mat_data( i_mat )%parent( i_parent )%type ) .NE. &
       !     'extra_TB_couplings' ) THEN
       
       CLOSE( file_num )
       
       !-------------------------------------------------------------------
       ! End of parent material loop
       !-------------------------------------------------------------------
       
    END DO

    !-------------------------------------------------------------
    ! CHECKS CONSISTENCY IN PARENT LATTICE STRUCTURES 
    !-------------------------------------------------------------

    DO i_par = 1, nr_parents 

       IF (TRIM(mat_data%parent( i_par )%cry) .ne. &
            TRIM(mat_data%parent( 1 )%cry) ) THEN

          WRITE(*,*) "WARNING: incompatible parent lattices"

       END IF

    END DO

    IF (nr_parents .gt. 0) THEN
      mat_data%cry  = mat_data%parent( 1 )%cry
    END IF

    !-------------------------------------------------------------
    ! CHECKS CONSISTENCY IN THE CHOSEN SCHEME OF PARENT MATERIALS
    !-------------------------------------------------------------

    DO i_par = 1, nr_parents 

       IF (TRIM(mat_data%parent( i_par )%scheme) .ne. &
            TRIM(mat_data%parent( 1 )%scheme) ) THEN
          CALL throw_init_exception(ERR_SCHE_INCON)

       END IF

    END DO

    IF (nr_parents .gt. 0) THEN
      mat_data%scheme  = mat_data%parent( 1 )%scheme
    END IF


  END SUBROUTINE read_data              
  
  !===========================================================================
  !
  ! Subroutine read_ETB_orbitals( ion )
  !
  !===========================================================================
  !
  ! Read the onsite ETB parameters from a parent material data file :
  !
  ! Since the input is not redundant, this routine repeats the onsite ETB data
  ! according to the degeneracies of each atomic orbital.
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> ion : pointer to a ion data list
  !
  !===========================================================================
  
  SUBROUTINE read_ETB_orbitals( file_num, ion )

    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER :: file_num
    TYPE(ion_orbit), DIMENSION( : ), POINTER :: ion

    !=========================================================================
    ! Local variables
    !=========================================================================

    REAL      ( dp ), DIMENSION( : ), ALLOCATABLE :: energy
    REAL      ( dp ), DIMENSION( : ), ALLOCATABLE :: so_energy
    CHARACTER (STATELEN),       DIMENSION( : ), ALLOCATABLE :: orb

    CHARACTER ( LEN = 1 ),      DIMENSION( : ), ALLOCATABLE :: so_state
    INTEGER,                    DIMENSION( : ), ALLOCATABLE :: n_orb, n_so
    INTEGER,                    DIMENSION( : ), ALLOCATABLE :: deg
    
    REAL      ( dp ) :: offset
    
    INTEGER :: i_ion, n_ion
    INTEGER :: i_state, n_states
    INTEGER :: i_deg, i_energy, i_orb, i_so
    INTEGER :: err
    LOGICAL :: error
    CHARACTER (LEN = 500) :: read_format

    !=========================================================================
    ! Get the number of ions in the parent material
    !=========================================================================

    n_ion = SIZE(ion)

    !=========================================================================
    ! Allocate and read n_orb : ( n_ion ) dimension array
    ! --> number of distinct orbitals involved per ion 
    !=========================================================================
    
    ALLOCATE( n_orb( n_ion ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
         'read_ETB_orbitals', 'n_orb' )
    
    n_orb = 0
        
    read_format = ' '
    read_format = 'number of orbitals per ion ='

    CALL label_in_file( file_num, read_format, error, 'no' )
    IF ( error ) CALL read_errors( 'missing', 'database', read_format )
    
    READ( file_num, * ) n_orb

    !=========================================================================
    ! Allocate and read n_so : ( n_ion ) dimension array
    ! --> number of distinct spin-orbit couplings per ion 
    !=========================================================================
    
    ALLOCATE( n_so( n_ion ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
         'read_ETB_orbitals', 'n_so' )
    
    n_so = 0
    
    read_format = ' '
    read_format = 'number of spin-orbit parameters per ion ='
    
    CALL label_in_file( file_num, read_format, error, 'no' )
    IF ( error ) CALL read_errors( 'missing', 'database', read_format )
    
    READ( file_num, * ) n_so
        
    !=========================================================================
    ! Read the ions' names
    !=========================================================================
       
    read_format =  ' '             
    read_format = 'ions ='

    CALL label_in_file( file_num, read_format, error, 'no' )
    IF ( error ) CALL read_errors( 'missing', 'database', read_format )
    
    READ( file_num, * ) ion( : )%name

    !=========================================================================
    !
    ! Read the onsite and spin-orbit ETB parameters :
    !
    ! Loop on the ions in the parent material
    !
    !=========================================================================

    DO i_ion = 1, n_ion

       !----------------------------------------------------------------------
       ! Read the valence of each atom for neighbour list
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = TRIM( ion( i_ion )%name ) // ' valence ='

       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       READ( file_num, * ) ion( i_ion )%valence


       !----------------------------------------------------------------------
       ! Allocate orb : n_orb( i_ion ) dimension array
       ! --> labels of the atomic orbitals for the current ion ( i_ion )
       !----------------------------------------------------------------------
       
       ALLOCATE( orb( n_orb( i_ion ) ), STAT = err )
       IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
            'read_ETB_orbitals', 'orb' )

       orb = ' '
       
       !----------------------------------------------------------------------
       ! Allocate energy : n_orb( i_ion ) dimension array
       ! --> ETB onsite energies for the current ion ( i_ion )
       !----------------------------------------------------------------------

       ALLOCATE( energy( n_orb( i_ion ) ), STAT = err )
       IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
            'read_ETB_orbitals', 'energy' )
       
       energy = 0.0D0

       !----------------------------------------------------------------------
       ! Allocate deg : n_orb( i_ion ) dimension array
       ! --> degeneracies of the atomic orbitals for the current ion ( i_ion )
       !----------------------------------------------------------------------

       ALLOCATE( deg( n_orb( i_ion ) ), STAT = err )
       IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
            'read_ETB_orbitals', 'deg' )
          
       deg = 0

       !----------------------------------------------------------------------
       ! Read the bare onsite data
       !-------------------------------------------------------------------
       
       read_format = ' '
       read_format = TRIM( ion( i_ion )%name ) // ' onsite'

       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       READ( file_num, * ) ( orb( i_orb ), energy( i_orb ), deg( i_orb ), &
             i_orb = 1, n_orb( i_ion ) )

       !----------------------------------------------------------------------
       ! Read the band offset
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = TRIM( ion( i_ion )%name ) // ' offset ='

       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       READ( file_num, * ) ion( i_ion )%offset

       !----------------------------------------------------------------------
       ! Read the band offset strain factor for d orbitals
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = TRIM( ion( i_ion )%name ) // ' b_d parameter ='

       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       READ( file_num, * ) ion( i_ion )%b_d

       !----------------------------------------------------------------------
       ! Relativistic data
       !----------------------------------------------------------------------
       
       !-------------------------------------------------------------------
       ! Allocate so_state : n_so( i_ion ) dimension array
       ! --> label of the spin-orbit couplings for the current ion ( i_ion )
       !-------------------------------------------------------------------
       
       ALLOCATE( so_state( n_so( i_ion ) ), STAT = err )
       IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
            'read_ETB_orbitals', 'so_state' )
       
       so_state = ' '
       
       !-------------------------------------------------------------------
       ! Allocate so_energy : n_so( i_ion ) dimension array
       ! --> spin-orbit coupling energies for the current ion ( i_ion )
       !-------------------------------------------------------------------
       
       ALLOCATE( so_energy( n_so( i_ion ) ), STAT = err )
       IF ( err .NE. 0 ) CALL alloc_error( 'input_data', &
            'read_ETB_orbitals', 'so_energy' )
       
       so_energy = 0.0D0
       
       !-------------------------------------------------------------------
       ! Read the bare spin - orbit data
       !-------------------------------------------------------------------
       
       read_format =  ' '
       read_format = TRIM( ion( i_ion )%name ) // ' spin orbit' 
       
       CALL label_in_file( file_num, read_format, error )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       READ( file_num, * ) &
            ( so_state( i_so ), so_energy( i_so ),&
              i_so = 1, n_so( i_ion ) )
       
       !-------------------------------------------------------------------
       ! End of relativistic data
       !-------------------------------------------------------------------
       !======================================================================
       ! Calculate total number of atomic orbitals on the current ion
       !======================================================================
       
       n_states = SUM( deg )
       
       !======================================================================
       ! Allocate the arrays of the onsite data list
       !======================================================================
       !----------------------------------------------------------------------
       ! Allocate ion( i_ion )%energy : n_states dimension array
       ! --> ETB onsite energies for the current ion ( i_ion )
       !----------------------------------------------------------------------
       
       CALL allocate_ion(ion(i_ion), n_states, n_so(i_ion))
       
       ! In addition, allocate the ion(i_ion)%deg for degeneracy info
       ALLOCATE( ion(i_ion)%deg( SIZE(deg) ), STAT = err )
       IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 
       ion(i_ion)%deg = 0   
          
       !----------------------------------------------------------------------
       !
       ! Store the ETB onsite data in memory :
       !
       ! --> the data read from the file are not redundant but we must
       !     store them in memory in a redundant form :
       !     * the degeneracy is stored
       !
       !     * we use the degeneracy of the atomic orbitals to store
       !       each energy as many times as required
       !
       !     * an atomic state is associated with each energy
       !       to label the vectors of the ETB basis
       !       ( function "states_sym" )
       !
       !-------------------------------------------------------------------
       
       ion( i_ion )%deg(1:n_orb(i_ion)) = deg(1:n_orb(i_ion))

       ion( i_ion )%state = states_sym( orb , n_states )
       
       i_state = 1
       i_energy = 1
       
       !...................................................................
       ! Loop on the atomic orbitals associated with the current ion
       !...................................................................

       DO WHILE ( i_state .LE. n_states )

          !................................................................
          ! Loop on the number of orbitals sharing the same onsite energy
          !................................................................
          
          DO i_deg = 1, deg( i_energy )
             
             ion( i_ion )%energy( i_state ) = energy( i_energy )
             i_state = i_state + 1
             
          END DO
                    
          i_energy = i_energy + 1
          
       END DO
       
       !----------------------------------------------------------------------
       ! Deallocate the onsite work arrays
       !----------------------------------------------------------------------
       
       DEALLOCATE( orb, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'input_data', 'read_ETB_orbitals', 'orb' )

       DEALLOCATE( energy, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'input_data', 'read_ETB_orbitals', 'energy' )
       
       DEALLOCATE( deg, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'input_data', 'read_ETB_orbitals', 'deg' )

       !----------------------------------------------------------------------
       ! Store the relativistic data
       !----------------------------------------------------------------------
       
       ion( i_ion )%so_energy = so_energy
       ion( i_ion )%so_state = so_state
       
       DEALLOCATE( so_state, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'input_data', &
            'read_ETB_orbitals', 'so_state' )
       
       DEALLOCATE( so_energy, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'input_data', &
            'read_ETB_orbitals', 'so_energy' )
             
       !----------------------------------------------------------------------
       ! End of loop on the ions in the parent material
       !----------------------------------------------------------------------

    END DO

    
    !=========================================================================
    ! Deallocate work arrays
    !=========================================================================

    DEALLOCATE( n_orb, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'input_data', 'read_ETB_orbitals', 'n_orb' )
    
    DEALLOCATE( n_so, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'input_data', 'read_ETB_orbitals', 'n_so' )
    
    !=========================================================================
    
  END SUBROUTINE read_ETB_orbitals
  



  !===========================================================================
  !
  ! Subroutine read_ETB_ons_corr( file_num, mat_data, i_parent)
  !
  !===========================================================================
  !
  ! Read the onsite and SO corrections from a parent material data file
  !
  ! The corrections are pair-dependent and specific for each orbital type
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> mat_data : data of the material
  !   --> i_parent : index of a parent of material
  !
  !===========================================================================
  
  SUBROUTINE read_ETB_ons_corr( file_num, mat_data, i_parent)

    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: file_num, i_parent
    TYPE(material_data) :: mat_data 
        
    !=========================================================================
    ! Local variables
    !=========================================================================
    
    INTEGER :: n_pair, n_I, n_so, n_O
    INTEGER :: i_pair, i_I, i_so
    INTEGER, DIMENSION(:), ALLOCATABLE :: deg
    INTEGER :: err
    CHARACTER (LEN = 500) :: read_format
    LOGICAL :: error
    !=========================================================================
    ! Get the number of pairs in the current parent
    !=========================================================================

    n_pair = mat_data%parent( i_parent )%n_pair 
        
    !=========================================================================
    ! Loop on the number of pairs
    !=========================================================================
    
    DO i_pair = 1, n_pair

       !----------------------------------------------------------------------
       ! Onsite corrections
       !----------------------------------------------------------------------
       ! Select the relevant data list in the current material file
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = "onsite correction " // &
                     TRIM( int_to_char( i_pair ) ) // " ="
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       !----------------------------------------------------------------------
       ! Read the corresponding pair names, coupling elements number and order
       !----------------------------------------------------------------------
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_parent )%pair( i_pair )%name, &
            n_I, mat_data%parent( i_parent )%pair( i_pair )%order

       IF ( err .NE. 0 ) CALL read_errors( 'bad pairs', 'database', 'pairs' )
       
       IF ( mat_data%parent( i_parent )%pair( i_pair )%order .NE. &
            mat_data%parent( 1 )%pair( i_pair )%order ) &
            call throw_init_exception(ERR_DB_PAIR)
       
       !----------------------------------------------------------------------
       ! Allocate data arrays for onsite corrections
       !----------------------------------------------------------------------   
       CALL allocate_ons_corr(mat_data%parent(i_parent)%pair(i_pair), n_I)
       
       !----------------------------------------------------------------------
       ! Read the onsite corrections for the current pair ( i_pair )
       !----------------------------------------------------------------------
       
       ! look for e.g.  "Ga As coupling 1"

       read_format = ' '
       read_format = &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 1 ) ) // ' ' // &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 2 ) ) // &
            ' onsite correction ' // TRIM( int_to_char( mat_data&
            %parent( i_parent )%pair( i_pair )%order ) )
       
       CALL label_in_file( file_num, read_format, error )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       DO i_I = 1, n_I
          
          READ( file_num, * ) &
               mat_data%parent(i_parent)%pair(i_pair)%ons_corr(i_I), &
               mat_data%parent(i_parent)%pair(i_pair)%fac_I(i_I), &
               mat_data%parent(i_parent)%pair(i_pair)%l_I(i_I)

       END DO
       
       !----------------------------------------------------------------------
       ! SO corrections
       !----------------------------------------------------------------------
       ! Select the relevant data list in the current material file
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = "SO correction " // &
                     TRIM( int_to_char( i_pair ) ) // " ="
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       !----------------------------------------------------------------------
       ! Read the corresponding pair names, coupling elements number and order
       !----------------------------------------------------------------------
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_parent )%pair( i_pair )%name, &
            n_so, mat_data%parent( i_parent )%pair( i_pair )%order

       IF ( err .NE. 0 ) CALL read_errors( 'bad pairs', 'database', 'pairs' )
       
       IF ( mat_data%parent( i_parent )%pair( i_pair )%order .NE. &
            mat_data%parent( 1 )%pair( i_pair )%order ) &
            call throw_init_exception(ERR_DB_PAIR)
       
       !----------------------------------------------------------------------
       ! Allocate data arrays for the SO corrections
       !----------------------------------------------------------------------   
       CALL allocate_so_corr(mat_data%parent(i_parent)%pair(i_pair), n_so)
       
       !----------------------------------------------------------------------
       ! Read the SO correction for the current pair ( i_pair )
       !----------------------------------------------------------------------
       
       ! look for e.g.  "Ga As coupling 1"

       read_format = ' '
       read_format = &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 1 ) ) // ' ' // &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 2 ) ) // &
            ' SO correction ' // TRIM( int_to_char( mat_data&
            %parent( i_parent )%pair( i_pair )%order ) )
       
       CALL label_in_file( file_num, read_format, error )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
          
       DO i_so = 1, n_so
          
          READ( file_num, * ) &
               mat_data%parent(i_parent)%pair(i_pair)%so_corr(i_so), &
               mat_data%parent(i_parent)%pair(i_pair)%fac_so(i_so)

       END DO
       
       !----------------------------------------------------------------------
       ! Offset correction
       !----------------------------------------------------------------------
       ! Select the relevant data list in the current material file
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = "offset correction " // &
                     TRIM( int_to_char( i_pair ) ) // " ="
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       !----------------------------------------------------------------------
       ! Read the corresponding pair names, coupling elements number and order
       !----------------------------------------------------------------------
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_parent )%pair( i_pair )%name, &
            n_O, mat_data%parent( i_parent )%pair( i_pair )%order

       IF ( err .NE. 0 ) CALL read_errors( 'bad pairs', 'database', 'pairs' )
       
       IF ( mat_data%parent( i_parent )%pair( i_pair )%order .NE. &
            mat_data%parent( 1 )%pair( i_pair )%order ) &
            call throw_init_exception(ERR_DB_PAIR)
       
       !----------------------------------------------------------------------
       ! Read the offset correction for the current pair ( i_pair )
       !----------------------------------------------------------------------
       
       ! look for e.g.  "Ga As coupling 1"

       read_format = ' '
       read_format = &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 1 ) ) // ' ' // &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 2 ) ) // &
            ' offset correction ' // TRIM( int_to_char( mat_data&
            %parent( i_parent )%pair( i_pair )%order ) )
       
       CALL label_in_file( file_num, read_format, error )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       READ( file_num, * ) &
            mat_data%parent(i_parent)%pair(i_pair)%fac_O, &
            mat_data%parent(i_parent)%pair(i_pair)%l_O, &
            mat_data%parent(i_parent)%pair(i_pair)%delta_d
       
    ENDDO
    
  END SUBROUTINE read_ETB_ons_corr
  


  !===========================================================================
  !
  ! Subroutine read_ETB_intra( file_num, mat_data, i_parent)
  !
  !===========================================================================
  !
  ! Read the parameters for intracoupling from a parent material data file
  !
  ! The intracouplings are pair-dependent
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> mat_data : data of the material
  !   --> i_parent : index of a parent of material
  !
  !===========================================================================
   
  SUBROUTINE read_ETB_intra( file_num, mat_data, i_parent)

    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: file_num, i_parent
    TYPE(material_data) :: mat_data 
        
    !=========================================================================
    ! Local variables
    !=========================================================================
    
    INTEGER :: n_pair, n_intra
    INTEGER :: i_pair, i_intra
    INTEGER :: err
    CHARACTER (LEN = 500) :: read_format
    LOGICAL :: error
    !=========================================================================
    ! Get the number of pairs in the current parent
    !=========================================================================
    
    n_pair = mat_data%parent( i_parent )%n_pair 
        
    !=========================================================================
    ! Loop on the number of pairs
    !=========================================================================
    
    DO i_pair = 1, n_pair

       !----------------------------------------------------------------------
       ! Intracoupling parameters
       !----------------------------------------------------------------------
       ! Select the relevant data list in the current material file
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = "intra " // &
                     TRIM( int_to_char( i_pair ) ) // " ="
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       !----------------------------------------------------------------------
       ! Read the corresponding pair names, coupling elements number and order
       !----------------------------------------------------------------------
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_parent )%pair( i_pair )%name, &
            n_intra, mat_data%parent( i_parent )%pair( i_pair )%order

       IF ( err .NE. 0 ) CALL read_errors( 'bad pairs', 'database', 'pairs' )
       
       IF ( mat_data%parent( i_parent )%pair( i_pair )%order .NE. &
            mat_data%parent( 1 )%pair( i_pair )%order ) &
            call throw_init_exception(ERR_DB_PAIR)
 
       !----------------------------------------------------------------------
       ! Allocate data arrays for the intracoupling parameters
       !----------------------------------------------------------------------   
       CALL allocate_intra(mat_data%parent(i_parent)%pair(i_pair), n_intra)

       !----------------------------------------------------------------------
       ! Read the intracoupling parameters for the current pair ( i_pair )
       !----------------------------------------------------------------------
       
       ! look for e.g.  "Ga As coupling 1"

       read_format = ' '
       read_format = &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 1 ) ) // ' ' // &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 2 ) ) // &
            ' intra ' // TRIM( int_to_char( mat_data&
            %parent( i_parent )%pair( i_pair )%order ) )
       
       CALL label_in_file( file_num, read_format, error )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
          
       DO i_intra = 1, n_intra
          
          READ( file_num, * ) &
               mat_data%parent(i_parent)%pair(i_pair)%intra(i_intra), &
               mat_data%parent(i_parent)%pair(i_pair)%fac_C(i_intra)

       END DO

    ENDDO

  END SUBROUTINE read_ETB_intra
  


  !===========================================================================
  !
  ! Subroutine read_ETB_couplings( i_mat, i_parent )
  !
  !===========================================================================
  !  
  ! NOTE : data are read from a data file which has been opened in
  ! another routine : its file numerical handle must be provided.
  !
  !___________________________________________________________________________
  !
  ! * input arguments :
  !
  !   --> i_mat    : numerical label to the current material
  !   --> i_parent : numerical label to the current parent material
  !
  !
  !===========================================================================
  
  SUBROUTINE read_ETB_couplings(file_num, mat_data, i_parent, scheme)
    
    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: file_num, i_parent
    TYPE(material_data) :: mat_data 
    CHARACTER (LEN = 5) :: scheme
        
    !=========================================================================
    ! Local variables
    !=========================================================================
    
    INTEGER :: n_pair, n_cpl
    INTEGER :: i_pair, i_cpl
    INTEGER :: err
    CHARACTER (LEN = 500) :: read_format
    LOGICAL :: error
    !=========================================================================
    ! Get the number of pairs in the current parent
    !=========================================================================
    
    n_pair = mat_data%parent( i_parent )%n_pair 
        
    !=========================================================================
    ! Loop on the number of pairs
    !=========================================================================
    
    DO i_pair = 1, n_pair

       !----------------------------------------------------------------------
       ! Select the relevant data list in the current material file
       !----------------------------------------------------------------------
       
       read_format = ' '
       read_format = "pair " // TRIM( int_to_char( i_pair ) ) // " ="
       
       CALL label_in_file( file_num, read_format, error, 'no' )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
       
       !----------------------------------------------------------------------
       ! Read the corresponding pair names, coupling elements number and order
       !----------------------------------------------------------------------
       
       READ ( file_num, *, IOSTAT = err ) &
            mat_data%parent( i_parent )%pair( i_pair )%name, &
            n_cpl, mat_data%parent( i_parent )%pair( i_pair )%order

       IF ( err .NE. 0 ) CALL read_errors( 'bad pairs', 'database', 'pairs' )
       
       IF ( mat_data%parent( i_parent )%pair( i_pair )%order .NE. &
            mat_data%parent( 1 )%pair( i_pair )%order ) &
            call throw_init_exception(ERR_DB_PAIR)
       
       !----------------------------------------------------------------------
       ! check database consistency (same n_cpl far all parents)
       !----------------------------------------------------------------------   
       IF ( i_parent.GT.1 .AND. n_cpl .NE. &
            SIZE( mat_data%parent( 1 )%pair( i_pair )%cpl ) ) &
            call throw_init_exception(ERR_DB_PAIR)
 
       !----------------------------------------------------------------------
       ! Allocate data arrays for the ETB couplings
       !----------------------------------------------------------------------   
       CALL allocate_pair(mat_data%parent(i_parent)%pair(i_pair), n_cpl)

       !----------------------------------------------------------------------
       ! Read the ETB coupling data for the current pair ( i_pair )
       !----------------------------------------------------------------------
       
       ! look for e.g.  "Ga As coupling 1"

       read_format = ' '
       read_format = &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 1 ) ) // ' ' // &
            TRIM( mat_data%parent( i_parent )&
            %pair( i_pair )%name( 2 ) ) // &
            ' coupling ' // TRIM( int_to_char( mat_data&
            %parent( i_parent )%pair( i_pair )%order ) )
       
       CALL label_in_file( file_num, read_format, error )
       IF ( error ) CALL read_errors( 'missing', 'database', read_format )
          
       IF (scheme .EQ. 'tan') THEN

         DO i_cpl = 1, n_cpl
            
            READ( file_num, * ) &
                 mat_data%parent(i_parent)%pair(i_pair)%cpl(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%energy(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%scaling(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%fac_P(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%fac_S(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%fac_Q(i_cpl)
  
         END DO

       ELSEIF (scheme .EQ. 'jancu') THEN

         DO i_cpl = 1, n_cpl
            
            READ( file_num, * ) &
                 mat_data%parent(i_parent)%pair(i_pair)%cpl(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%energy(i_cpl), &
                 mat_data%parent(i_parent)%pair(i_pair)%scaling(i_cpl)
  
         END DO

       ENDIF
       
       !----------------------------------------------------------------------
       ! Read the reference distance for the current pair ( i_pair )
       !----------------------------------------------------------------------
 
       IF ( fuzzy ) THEN       

          read_format = ' '
          read_format = &
               TRIM( mat_data%parent( i_parent )&
               %pair( i_pair )%name( 1 ) ) // ' ' // &
               TRIM( mat_data%parent( i_parent )&
               %pair( i_pair )%name( 2 ) ) // &
               ' pair reference distance ' // &
               TRIM( int_to_char( mat_data%parent( i_parent )&
               %pair( i_pair )%order ) ) // ' ='
                   
          CALL label_in_file( file_num, read_format, error, 'no' )
          IF ( error ) CALL read_errors( 'missing', 'database', read_format )
          
          READ( file_num, * ) mat_data%parent( i_parent )&
               %pair( i_pair )%dist_ref

       END IF

       !----------------------------------------------------------------------
       ! End of loop on pair number
       !----------------------------------------------------------------------
       
    END DO
    
  END SUBROUTINE read_ETB_couplings

  !===========================================================================

  !===========================================================================
  !
  ! Subroutine "states_sym"
  !
  !===========================================================================

  FUNCTION states_sym( orb, n_states )
    
    !=========================================================================

    ! input arguments :

    CHARACTER ( LEN = * ), DIMENSION( : ), INTENT( IN ) :: orb

    INTEGER, INTENT( IN ) :: n_states

    !_________________________________________________________________________

    ! output result :
    
    CHARACTER ( LEN = 5 ), DIMENSION( n_states ) :: states_sym

    !_________________________________________________________________________

    ! local variables :

    INTEGER :: i_orb, i_state

    !=========================================================================

    i_state = 1

    DO i_orb = 1, SIZE( orb )

       !______________________________________________________________________

       IF ( ( orb( i_orb ) .EQ. "s" ) .OR. ( orb( i_orb ) .EQ. "s*" ) ) THEN

          states_sym( i_state ) = TRIM( orb( i_orb ) )

          i_state = i_state + 1  
   
          !-------------------------------------------------------------------
          
       ELSEIF ( ( orb( i_orb ) .EQ. "px" ) .OR. &
                ( orb( i_orb ) .EQ. "py" ) .OR. &
                ( orb( i_orb ) .EQ. "pz" )  )   THEN
          
          states_sym( i_state ) =  orb( i_orb ) 
          i_state = i_state + 1  
          
          ! -------------------------------------------------------------------

       ELSEIF ( orb( i_orb ) .EQ. "p" ) THEN

          states_sym( i_state ) = "px"
          states_sym( i_state + 1 ) = "py"
          states_sym( i_state + 2 ) = "pz"
          i_state = i_state + 3

          !-------------------------------------------------------------------

       ELSEIF ( orb( i_orb ) .EQ. "d" ) THEN

          states_sym( i_state ) = "dxy"
          states_sym( i_state + 1 ) = "dyz"
          states_sym( i_state + 2 ) = "dzx"
          states_sym( i_state + 3 ) = "dx2y2"
          states_sym( i_state + 4 ) = "dz2r2"
          i_state = i_state + 5
          
          !-------------------------------------------------------------------
       ELSEIF ( ( orb( i_orb ) .EQ. "dxy" ) .OR. &
                ( orb( i_orb ) .EQ. "dyz" ) .OR. &
                ( orb( i_orb ) .EQ. "dzx" ) .OR. &
                ( orb( i_orb ) .EQ. "dx2y2" ) .OR. &
                ( orb( i_orb ) .EQ. "dz2r2" ) ) THEN
          
          states_sym( i_state ) = TRIM(  orb( i_orb ) )
          i_state = i_state + 1  

       ELSEIF ( orb( i_orb ) .EQ. "d12" ) THEN
          
          states_sym( i_state ) = "dx2y2"
          states_sym( i_state + 1 ) = "dz2r2"
          i_state = i_state + 2
          
          !___________________________________________________________________
          
       ELSEIF ( orb( i_orb ) .EQ. "d15" ) THEN
          
          states_sym( i_state ) = "dxy"
          states_sym( i_state + 1 ) = "dyz"
          states_sym( i_state + 2 ) = "dzx"
          i_state = i_state + 3
          
          !-------------------------------------------------------------------
          
       END IF
       
       !______________________________________________________________________
       
    END DO
    
    !=========================================================================
    
  END FUNCTION states_sym

END MODULE input_data
