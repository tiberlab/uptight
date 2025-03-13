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
!        Module "alloys" - (c) Jerome Gleize - August 2002
!
!=============================================================================
!
! contains routines used to insert regions in the supercell :
!
!_____________________________________________________________________________
!
!=============================================================================

MODULE alloys

  !===========================================================================
  !
  !  parameters :
  !

  USE precision 
  USE type_defs 
  USE input_output
  USE errors
  USE exceptions
  
  !===========================================================================

  IMPLICIT NONE
  PRIVATE 
  !===========================================================================
  
  PUBLIC :: init_mat_ion
  PUBLIC :: random_alloy_common, random_alloy_uncommon

  PRIVATE :: ternary_states_common, ternary_states_uncommon
  PRIVATE :: alloy_states_nomixing


CONTAINS

  !===========================================================================
  ! SUBROUTINE init_mat_ion is used to initialize the mat_data%ion fields
  ! 
  !  This takes the on-site elements from the parent materials and copy
  !  all infos into the mat%ion field.
  !  VCA averages are done where appropriate
  !===========================================================================
  SUBROUTINE init_mat_ion(mat_data)

    TYPE(material_data) :: mat_data

    !----------------------------------------------------------------------
    ! Material type testing : get the number of ions in the current material
    !----------------------------------------------------------------------
    REAL(dp), DIMENSION(:), POINTER :: content
    LOGICAL :: alloy_random, average
    INTEGER :: i_parent, i_state, n_ion, i_ion, np, err, i_st

    np = SUM(mat_data%nr_parents)

    ALLOCATE( content(np) )

    content(1:np) = mat_data%parent(1:np)%content
    alloy_random = mat_data%alloy_random 
    average = mat_data%average_cat

    mat_data%ref_axis = mat_data%parent( 1 )%ref_axis

    SELECT CASE ( TRIM( mat_data%type ) )
       
       !-------------------------------------------------------------------
       ! Simple compound
       !-------------------------------------------------------------------
       
    CASE ( 'simple' )
       
       n_ion = 1
        
       CALL allocate_mat_ion( mat_data, n_ion )
       
       mat_data%ion = mat_data%parent( 1 )%ion

       mat_data%lat_par = mat_data%parent( 1 )%lat_par
       mat_data%u_par = mat_data%parent( 1 )%u_par
       mat_data%int_str_par = mat_data%parent( 1 )%int_str_par
       
       !-------------------------------------------------------------------
       ! Binary compound
       !-------------------------------------------------------------------
       
    CASE ( 'binary' )
       
       n_ion = 2
       
       CALL allocate_mat_ion( mat_data, n_ion )

       SELECT CASE ( np )

       CASE (1)
          !................................................................
          !
          ! Case 1 : the binary is a crystal with one parent material 
          !          --> GaAs, InP, ... like
          !................................................................

          mat_data%ion = mat_data%parent( 1 )%ion
          mat_data%lat_par = mat_data%parent( 1 )%lat_par
          mat_data%u_par = mat_data%parent( 1 )%u_par
          mat_data%int_str_par = mat_data%parent( 1 )%int_str_par

       CASE (2)   
          !................................................................
          !
          ! Case 2 : binary alloy (SiGe) made of two simple compounds
          !          --> SiGe like  -->  contains 2 parents (Si and Ge)
          !          --> needs an extra coupling Si-Ge 
          !          treated with an extra coupling like a fake ternary
          !................................................................
          i_ion = 0
          
          CALL ternary_states_uncommon( mat_data%ion, &
               mat_data%parent, i_ion, content, &
               mat_data%bowing(1), alloy_random, average )         

       END SELECT
       
    CASE ( 'ternary' )
       !-------------------------------------------------------------------
       ! Ternary alloy : AlGaAs or GaAsP
       !-------------------------------------------------------------------  
       !...................................................................
       !
       ! Case 1 : the ternay is made of two binaries
       !          ( with a common atom ) --> AlGaAs, AlGaN, ...
       !
       ! Common atoms are treated separately :
       !
       ! --> the common atom must be averaged in any case
       !
       ! --> the other ( "uncommon" ) atoms are averaged only for VCA 
       !
       !...................................................................
       n_ion = 3
       
       CALL allocate_mat_ion( mat_data, n_ion )

       i_ion = 1   ! Assign position 1 to common ion            

       CALL ternary_states_common( mat_data%ion, &
            mat_data%parent, i_ion, content, &
            mat_data%bowing(1), alloy_random )
       
       ! i_ion is now 1 so uncommon ion will go second.
       CALL ternary_states_uncommon( mat_data%ion, &
            mat_data%parent, i_ion, content, &
            mat_data%bowing(1), alloy_random, average )
       
       ! check reference axis (must be the same)
       if (any ((mat_data%parent(1)%ref_axis .ne. mat_data%parent(2)%ref_axis)) ) then
         write(*,*) "UPT ERROR: reference axes, if provided, must be equal in all parent materials"
         call throw_init_exception(ERR_GENERAL)
       end if
       !................................................................
       !
       ! Average the lattice ( correct for VCA / crude for random )
       !
       !...................................................................
       
       ! **************************************
       !  for any ternary ( both algaas or  sige)  !!!!!
       ! **************************************
       
       mat_data%lat_par = &
            content(1) * mat_data%parent(1)%lat_par &
            + content(2) * mat_data%parent(2)%lat_par
       
       mat_data%u_par = &
            content(1) * mat_data%parent(1)%u_par &
            +content(2) * mat_data%parent(2)%u_par
       
       mat_data%int_str_par = &
            content(1) * mat_data%parent(1)%int_str_par &
            + content(2) * mat_data%parent(2)%int_str_par
       
       !-------------------------------------------------------------------
       ! Quaternary alloy
       !-------------------------------------------------------------------
       
    CASE ( 'quaternary' )
       
       n_ion = 4
       
       CALL allocate_mat_ion( mat_data, n_ion )
       
       ! 2016-04-12: for simplicity we do not do any averageing here, relying
       !             on conistent TB paramters

       CALL alloy_states_nomixing( mat_data%ion, mat_data%parent )
       
       mat_data%lat_par = 0
       mat_data%u_par = 0
       mat_data%int_str_par = 0

       DO i_parent = 1, np
         mat_data%lat_par = mat_data%lat_par &
            + content(i_parent) * mat_data%parent(i_parent)%lat_par
       
         mat_data%u_par = mat_data%u_par &
            + content(i_parent) * mat_data%parent(i_parent)%u_par 
       
         mat_data%int_str_par = mat_data%int_str_par &
            + content(i_parent) * mat_data%parent(i_parent)%int_str_par 
       END DO
       
       !-------------------------------------------------------------------
       ! End of material type testing
       !-------------------------------------------------------------------
       
    END SELECT

    DEALLOCATE(content)
    
  END SUBROUTINE init_mat_ion
 
  !===========================================================================
  !
  ! Subroutine " ternary_states_common"
  !
  ! CONVENTION : the common atom is assigned to number 1 in the ion data list
  !
  !===========================================================================
  
  SUBROUTINE ternary_states_common( ion, parent, i_ion, content, bowing, &
                                                                alloy_random )

    !=========================================================================

    ! Input arguments :

    TYPE (parent_data), DIMENSION( : ), POINTER :: parent
    TYPE (ion_orbit),   DIMENSION( : ), POINTER :: ion

    REAL ( dp ), DIMENSION(:), POINTER :: content
    REAL ( dp ),  INTENT( IN ) :: bowing
    LOGICAL, INTENT(IN) :: alloy_random
    INTEGER, INTENT( IN ) :: i_ion  ! defines the position of common ion

    !_________________________________________________________________________

    ! local variables :

    CHARACTER(ATOMLEN), DIMENSION( 2 ) :: name

    INTEGER, DIMENSION( 2 ) :: n_so, n_orb

    INTEGER :: i_ion_parent, i_st, err

    err = 0
    
    IF ( SIZE(parent(1)%ion) .NE. SIZE(parent(2)%ion) ) err = 1

    IF ( err .ne. 0 ) call throw_init_exception(ERR_GENERAL)


    !=========================================================================
    ! Find the common ion
    !=========================================================================

    DO i_ion_parent = 1, SIZE( parent( 1 )%ion )

       name( 1 ) = parent( 1 )%ion( i_ion_parent )%name       
       name( 2 ) = parent( 2 )%ion( i_ion_parent )%name

       !IF ( relat ) THEN

       n_so( 1 ) = SIZE( parent( 1 )%ion( i_ion_parent )%so_state )
       n_so( 2 ) = SIZE( parent( 2 )%ion( i_ion_parent )%so_state )

       !END IF

       n_orb( 1 ) = SIZE( parent( 1 )%ion( i_ion_parent )%state )
       n_orb( 2 ) = SIZE( parent( 2 )%ion( i_ion_parent )%state )

       !----------------------------------------------------------------------
       ! Atoms test : this must be a common atom
       !----------------------------------------------------------------------

       IF ( name( 1 ) .EQ. name( 2 ) ) THEN

          !-------------------------------------------------------------------
          ! Average parameters
          !-------------------------------------------------------------------

          !IF ( relat ) THEN
          !   write (*,*)  'n_so( 1 ), n_so( 2 )', n_so( 1 ),n_so( 2 )
           
          IF ( ( n_so( 1 ) .NE. n_so( 2 ) ) &
               .OR. ( n_orb( 1 ) .NE. n_orb( 2 ) ) ) &
               call throw_init_exception(ERR_ALLOY_TBBLK)

          IF ( ( n_orb( 1 ) .NE. n_orb( 2 ) ) ) &
               call throw_init_exception(ERR_ALLOY_TBBLK)
          
          ion( i_ion )%name = name( 1 )


          CALL allocate_ion(ion(i_ion), n_orb(1), n_so(1))

          !-------------------------------------------------------------------
          ! Store parameters 
          !-------------------------------------------------------------------

          ion( i_ion )%ind_ref = parent( 1 )%ion( i_ion_parent )%ind_ref

          ion( i_ion )%ind_ref_inv = parent( 1 )%ion( i_ion_parent )%ind_ref_inv

          ion( i_ion )%valence = parent( 1 )%ion( i_ion_parent )%valence

          ion( i_ion )%state = parent( 1 )%ion( i_ion_parent )%state
   

          ion( i_ion )%so_state = parent( 1 )%ion( i_ion_parent )%so_state

          !-------------------------------------------------------------------
          ! VCA average of common ion properties 
          !-------------------------------------------------------------------

          ion( i_ion )%offset = &
               content(1) * parent(1)%ion( i_ion_parent )%offset &
               +content(2) * parent(2)%ion( i_ion_parent )%offset

          ion( i_ion )%energy = &
               content(1) * parent(1)%ion( i_ion_parent )%energy &
               + content(2) * parent(2)%ion( i_ion_parent )%energy

          ion( i_ion )%so_energy = &
               content(1) * parent(1)%ion( i_ion_parent )%so_energy &
               +content(2) * parent(2)%ion( i_ion_parent )%so_energy

          ion( i_ion )%b_d = &
               content(1) * parent(1)%ion( i_ion_parent )%b_d &
               +content(2) * parent(2)%ion( i_ion_parent )%b_d 
          
          !-------------------------------------------------------------------
          ! Add bowing (if any) 
          !-------------------------------------------------------------------

          ion( i_ion )%energy = ion( i_ion )%energy &
                  - bowing * &
                  SQRT( content(1) * content(2) ) * &
                  ( parent( 1 )%ion( i_ion_parent )%energy - &
                    parent( 2 )%ion( i_ion_parent )%energy )

          ion( i_ion )%so_energy = ion( i_ion )%so_energy &
                  - bowing * &
                  SQRT( content(1) * content(2) ) * &
                  ( parent( 1 )%ion( i_ion_parent )%so_energy - &
                    parent( 2 )%ion( i_ion_parent )%so_energy )

       END IF

       !----------------------------------------------------------------------
       ! End loop on parent ions
       !----------------------------------------------------------------------

    END DO

    !=========================================================================

  END SUBROUTINE ternary_states_common
  

  !===========================================================================
  !
  ! Subroutine "ternary_states_uncommon"
  !
  !===========================================================================

  SUBROUTINE ternary_states_uncommon( ion, parent, i_ion, content, bowing, &
                                                       alloy_random, average )

    !=========================================================================

    ! Input arguments :

    TYPE (ion_orbit),   DIMENSION( : ), POINTER :: ion
    TYPE (parent_data), DIMENSION( : ), POINTER :: parent

    REAL ( dp ),  DIMENSION(:), POINTER :: content
    REAL ( dp ),  INTENT( IN ) :: bowing
    LOGICAL, INTENT(IN) :: alloy_random
    LOGICAL, INTENT(IN) :: average 
    !_________________________________________________________________________

    ! Input - output arguments :

    INTEGER, INTENT( INOUT ) :: i_ion

    !_________________________________________________________________________

    ! local variables :

    CHARACTER(ATOMLEN), DIMENSION( 2 ) :: name

    INTEGER, DIMENSION( 2 ) :: n_so, n_orb

    INTEGER :: i_ion_parent, i_parent, i_state, err

    !=========================================================================
    ! Loop on parent ions
    !=========================================================================

    DO i_ion_parent = 1, SIZE( parent( 1 )%ion )

       name( 1 ) = parent( 1 )%ion( i_ion_parent )%name  
       name( 2 ) = parent( 2 )%ion( i_ion_parent )%name

       !IF ( relat ) THEN

       n_so( 1 ) = SIZE( parent( 1 )%ion( i_ion_parent )%so_state )
       n_so( 2 ) = SIZE( parent( 2 )%ion( i_ion_parent )%so_state )

       !END IF

       n_orb( 1 ) = SIZE( parent( 1 )%ion( i_ion_parent )%state )
       n_orb( 2 ) = SIZE( parent( 2 )%ion( i_ion_parent )%state )

       !----------------------------------------------------------------------
       ! Atoms test : this should not be a common atom
       !----------------------------------------------------------------------

       IF ( name( 1 ) .NE. name( 2 ) ) THEN

          !-------------------------------------------------------------------
          ! Loop on parent material
          !-------------------------------------------------------------------

          DO i_parent = 1, 2

             !----------------------------------------------------------------
             ! Increment ion counter (i_ion after call to common is = 1)
             !----------------------------------------------------------------

             i_ion = i_ion + 1  ! Assign positions to uncommon atoms

             ion( i_ion )%name = name( i_parent )

             !-------------------------------------------------------------
             ! Allocate ion( i_ion )%so_state :
             ! contains the spin-orbit orbital's symbols for ion "i_ion"
             !-------------------------------------------------------------
             CALL allocate_ion(ion(i_ion), n_orb(i_parent), n_so(i_parent)) 
            
             !----------------------------------------------------------------
             ion( i_ion )%ind_ref = &
                  parent( i_parent )%ion( i_ion_parent )%ind_ref
             
             ion( i_ion )%ind_ref_inv = &
                  parent( i_parent )%ion( i_ion_parent )%ind_ref_inv

             ion( i_ion )%state = &
                  parent( i_parent )%ion( i_ion_parent )%state

             ion( i_ion )%valence = &
                  parent( i_parent )%ion( i_ion_parent )%valence

             ion( i_ion )%so_state = &
                  parent( i_parent )%ion( i_ion_parent )%so_state

             !---------------------------------------------------------------
             ! Random alloy : store onsite elements
             !----------------------------------------------------------------
             IF ( alloy_random ) THEN
               
             
                if (average) then

                    !---------------------------------------------------------------
                    ! Take VCA average of onsites and offsets
                    ! offset contains alignment of Ev to 0.0
                    ! random ionic fluctuations come from 'pot_atom.dat'
                    !----------------------------------------------------------------
                    ion( i_ion )%offset = &
                       content(1) * parent(1)%ion( i_ion_parent )%offset &
                       + content(2) * parent(2)%ion( i_ion_parent )%offset
                  
                    ion( i_ion )%energy = &
                       content(1) * parent(1)%ion( i_ion_parent )%energy &
                       + content(2) * parent(2)%ion( i_ion_parent )%energy
                else
                    !---------------------------------------------------------------
                    ! Do not average offsets
                    ! offsets contains alignment of Ev to 0.0 and band offsets
                    !----------------------------------------------------------------

                    ion( i_ion )%offset = &
                       parent( i_parent )%ion( i_ion_parent )%offset
     
                    ion( i_ion )%energy = &
                       parent( i_parent )%ion( i_ion_parent )%energy
     
                end if      

                ion( i_ion )%so_energy = &
                     parent( i_parent )%ion( i_ion_parent )%so_energy
                
                ion( i_ion )%b_d = &
                     parent( i_parent )%ion( i_ion_parent )%b_d 
                
             ELSE
                
                !-------------------------------------------------------------
                ! VCA : average and add bowing
                !-------------------------------------------------------------
                ion( i_ion )%offset = &
                     content(1) * parent(1)%ion( i_ion_parent )%offset &
                     + content(2) * parent(2)%ion( i_ion_parent )%offset
                
                ion( i_ion )%energy = &
                     content(1) * parent(1)%ion( i_ion_parent )%energy &
                     + content(2) * parent(2)%ion( i_ion_parent )%energy
            
                 
                ion( i_ion )%energy = ion( i_ion )%energy &
                     - bowing * &
                     SQRT( content(1) * content(2) ) &
                     * ( parent( 1 )%ion( i_ion_parent )%energy - &
                     parent( 2 )%ion( i_ion_parent )%energy )

                ion( i_ion )%so_energy = &
                     content(1) * parent(1)%ion( i_ion_parent )%so_energy &
                     + content(2) * parent(2)%ion( i_ion_parent )%so_energy
                
                ion( i_ion )%so_energy = ion( i_ion )%so_energy &
                     - bowing * &
                     SQRT( content(1) * content(2) ) * &
                     ( parent( 1 )%ion( i_ion_parent )%so_energy - &
                       parent( 2 )%ion( i_ion_parent )%so_energy )
                
                ion( i_ion )%b_d = &
                     content(1) * parent(1)%ion( i_ion_parent )%b_d &
                     + &
                     content(2) * parent(2)%ion( i_ion_parent )%b_d
                
                 
                !-------------------------------------------------------------
                ! End of random test
                !-------------------------------------------------------------
                
             END IF

             !----------------------------------------------------------------
             ! End of non-common atoms loop
             !----------------------------------------------------------------

          END DO

          !-------------------------------------------------------------------
          ! End of atoms test
          !-------------------------------------------------------------------

       END IF

       !----------------------------------------------------------------------
       ! End loop on parent ions
       !----------------------------------------------------------------------

    END DO

    !=========================================================================

  END SUBROUTINE ternary_states_uncommon



  !===========================================================================
  !
  ! Subroutine "alloy_states_nomixing"
  !
  ! setup onsite parameters without doing any mixing, therefore not needing
  ! distinction between 'common' and 'uncommon' ions
  !
  !===========================================================================

  SUBROUTINE alloy_states_nomixing( ion, parent )

    !=========================================================================

    ! Input arguments :

    TYPE (ion_orbit),   DIMENSION( : ), POINTER :: ion
    TYPE (parent_data), DIMENSION( : ), POINTER :: parent

    !_________________________________________________________________________

    ! local variables :

    CHARACTER(ATOMLEN) :: name

    INTEGER :: n_so, n_orb

    INTEGER :: i_ion, i_ion_parent, i_parent, i, err

    LOGICAL :: exists


    ! counter for distinct ions
    i_ion = 0

    !=========================================================================
    ! Loop on parents
    !=========================================================================

    DO i_parent = 1, SIZE(parent)

      !print *, parent(i_parent)%name 

      !=======================================================================
      ! Loop on ions
      !=======================================================================

      DO i_ion_parent = 1, SIZE( parent( i_parent )%ion )

        name = parent( i_parent )%ion( i_ion_parent )%name  

        n_so = SIZE( parent( i_parent )%ion( i_ion_parent )%so_state )

        n_orb = SIZE( parent( i_parent )%ion( i_ion_parent )%state )

        !----------------------------------------------------------------------
        ! Check if it has already been added to the ions container
        !----------------------------------------------------------------------

        exists = .FALSE.
        DO i = 1, i_ion
          IF (name .EQ. ion(i)%name) THEN
            exists = .TRUE.
          END IF
        END DO

        IF (exists .EQV. .TRUE.) THEN
          CYCLE
        END IF
        !print *, "   ", name



        i_ion = i_ion + 1  ! Assign positions to uncommon atoms

        ion( i_ion )%name = name

        !-------------------------------------------------------------
        ! Allocate ion( i_ion )%so_state :
        ! contains the spin-orbit orbital's symbols for ion "i_ion"
        !-------------------------------------------------------------
        CALL allocate_ion(ion(i_ion), n_orb, n_so) 

        !---------------------------------------------------------------
        ! Random alloy : store onsite elements
        !----------------------------------------------------------------

        ion( i_ion )%ind_ref_inv = &
        parent( i_parent )%ion( i_ion_parent )%ind_ref_inv

        ion( i_ion )%ind_ref = &
        parent( i_parent )%ion( i_ion_parent )%ind_ref

        ion( i_ion )%state = &
        parent( i_parent )%ion( i_ion_parent )%state

        ion( i_ion )%valence = &
        parent( i_parent )%ion( i_ion_parent )%valence

        ion( i_ion )%so_state = &
        parent( i_parent )%ion( i_ion_parent )%so_state



        ion( i_ion )%offset = &
        parent( i_parent )%ion( i_ion_parent )%offset

        ion( i_ion )%energy = &
        parent( i_parent )%ion( i_ion_parent )%energy


        ion( i_ion )%so_energy = &
        parent( i_parent )%ion( i_ion_parent )%so_energy

        ion( i_ion )%b_d = &
        parent( i_parent )%ion( i_ion_parent )%b_d 


        !-------------------------------------------------------------------
        ! End loop on ions
        !-------------------------------------------------------------------

      END DO

      !----------------------------------------------------------------------
      ! End loop on parents 
      !----------------------------------------------------------------------

    END DO

    !=========================================================================

  END SUBROUTINE alloy_states_nomixing




  !===========================================================================

  SUBROUTINE random_alloy_common( atom_name, n_atom, conc, ref_name )
    
    !=========================================================================
    
    ! input arguments :

    INTEGER,                               INTENT( IN ) :: n_atom
    
    REAL ( dp ),                INTENT( INOUT ) :: conc
    
    CHARACTER(ATOMLEN), DIMENSION( : ), INTENT( IN ) :: ref_name
    
    !=========================================================================
    
    ! input - output arguments :
    
    CHARACTER(ATOMLEN), DIMENSION( n_atom ), INTENT( INOUT ) :: atom_name
    
    !=========================================================================
    
    ! local variables :
    
    REAL ( dp ), DIMENSION( 2 ) :: alloy_test
    REAL ( dp )                 :: rdm
            
    INTEGER, DIMENSION( : ), ALLOCATABLE :: index

    INTEGER, DIMENSION( n_atom ) :: atom_index
    LOGICAL, DIMENSION( n_atom ) :: atom_test

    INTEGER :: i_atom, n_index, i_random, n_random
    INTEGER :: j, random, random_test, random_limit, err
      
    !=========================================================================
    
    ! initialize the cation lattice :

    !=========================================================================

   IF ( equiv( conc, 1.0d0, 1.d-10, fuzzy ) ) THEN  
       
       !______________________________________________________________________
       
       ! if we have the pure material 2 : set all cation names to material 2's
       
       WHERE ( atom_name .EQ. ref_name( 2 ) ) atom_name = ref_name( 3 )
       
       !______________________________________________________________________
       
    ELSEIF ( equiv( conc, 0.0d0, 1.d0-10, fuzzy ) ) THEN  
       
       !______________________________________________________________________
       
       ! if we have the pure material 1 : set all cation names to material 1's
       
       WHERE ( atom_name .EQ. ref_name( 3 ) ) atom_name = ref_name( 2 )
       
       !______________________________________________________________________
       
    ELSE

       !______________________________________________________________________
       
       ! otherwhise : set all cation names to material 1's
       
       WHERE ( atom_name .EQ. ref_name( 3 ) ) atom_name = ref_name( 2 )
       
       !----------------------------------------------------------------------
       
       ! select "index" : cation indexes (to perform the random distribution)
       
       atom_test = .FALSE.

       atom_index = (/ ( i_atom, i_atom = 1, n_atom ) /)

       WHERE ( atom_name .EQ. ref_name( 2 ) ) atom_test = .TRUE.
       
       n_index = COUNT( atom_test )

       ALLOCATE( index( n_index ) )
       
       index = PACK( atom_index, atom_test )  !  numeri  dei soli atomi di Ga (cation)

       !______________________________________________________________________

    END IF
        
    !=========================================================================

    IF ( .NOT. equiv( conc, 1.0d0, 1.d0-10, fuzzy ) .AND. &
         .NOT. equiv( conc, 0.0d0, 1.d0-10, fuzzy ) ) THEN
       
       !______________________________________________________________________
       
       ! perform random cation distribution when the material is not pure :
       
       !______________________________________________________________________

       ! calculate "n_random" : the number of cations to replace :

       alloy_test( 1 ) = REAL( INT( conc * n_index ) )
       alloy_test( 2 ) = REAL( INT( conc * n_index ) + 1 )

       alloy_test = ABS( alloy_test / REAL( n_index ) - conc )
       
       IF ( alloy_test( 1 ) .LT. alloy_test( 2 ) ) THEN
          n_random = INT( conc * n_index )
       ELSE
          n_random = INT( conc * n_index ) + 1
       END IF


          !PRINT*, "n Random,n_index ", n_random, n_index  !  0 e  0 per  SiGe !!!??

       conc = DBLE( n_random ) / DBLE( n_index )




      ! PRINT*
      ! PRINT*, "Random alloy treatment"
      ! PRINT*
      ! PRINT*, "Concentration = ", conc
      ! PRINT*
      ! PRINT*, "Total number of cations = ", n_index
      ! PRINT*
      ! PRINT*, "Replaced cations = ", n_random
      ! PRINT*
       
       !______________________________________________________________________
       
       i_random = 1
       
       DO WHILE ( i_random .LE. n_random )
          
          !-------------------------------------------------------------------
          
          j = 0
          
          CALL RANDOM_NUMBER( rdm )
          random_limit = INT( rdm * 99 + 1 )          
          
          CALL RANDOM_NUMBER( rdm )
          random_test = INT( rdm * random_limit + 1 )
          
          DO WHILE ( j .LE. random_test )        
             j =  j + 1
             CALL RANDOM_NUMBER( rdm )
          END DO
          
          !-------------------------------------------------------------------
          
          random = INT( rdm * ( n_index - 1 ) + 1 )
          
          IF ( index( random ) .NE. 0 ) THEN
             
             i_random = i_random + 1
             atom_name( index( random ) ) = ref_name( 3 )
             index( random ) = 0
             
          END IF
          
          !-------------------------------------------------------------------
          
       END DO
       
       DEALLOCATE( index )
       
       !______________________________________________________________________
       
    END IF
    
    !=========================================================================
    
  END SUBROUTINE random_alloy_common
  
  !===========================================================================
 

! ------------------------------------------------------------------------------------ 
  !                            !!!! ATTENZIONE !!!!
! *************************************************************************************
! random_alloy_sige    per  random substitution  of Ge in Si supercell  ->  SiGe alloy
! *************************************************************************************


  !===========================================================================

  SUBROUTINE random_alloy_uncommon( atom_name, n_atom, conc, ref_name )
    
    !=========================================================================
    
    ! input arguments :

    INTEGER,                               INTENT( IN ) :: n_atom
    
    REAL ( dp ),                INTENT( INOUT ) :: conc
    
    CHARACTER(ATOMLEN), DIMENSION( : ), INTENT( IN ) :: ref_name
    
    !=========================================================================
    
    ! input - output arguments :
    
    CHARACTER(ATOMLEN), DIMENSION( n_atom ), INTENT( INOUT ) :: atom_name
    
    !=========================================================================
    
    ! local variables :
    
    REAL ( dp ), DIMENSION( 2 ) :: alloy_test
    REAL ( dp )                 :: rdm
            
    INTEGER, DIMENSION( : ), ALLOCATABLE :: index

    INTEGER, DIMENSION( n_atom ) :: atom_index
    LOGICAL, DIMENSION( n_atom ) :: atom_test

    INTEGER :: i_atom, n_index, i_random, n_random
    INTEGER :: j, random, random_test, random_limit, err
      
    !=========================================================================
    
    ! initialize the cation lattice :

    !=========================================================================

!   PRINT*, " ref_name(1 , 2,  3 )  " , ref_name( 1 ), "   ", ref_name( 2 ),  "    ",ref_name( 3 )

!!!  **** prima era : ref_name = AsGaAl = stringa contenente tutti gli  ioni del mat alloy  
!!  ora :  ref_name = SiGe ->   ref_name( 1 ) = Si , ref_name( 2 ) = Ge  !!!!!!! 
    
   IF ( equiv( conc, 1.0d0, 1.d0-10, fuzzy ) ) THEN  
       
       !______________________________________________________________________
  !  concentr. x(Ge) = 1  ->  tutti  atomi di mat 2 =  Ge 
       
       ! if we have the pure material 2 : set all cation names to material 2's
       
!!$       WHERE ( atom_name .EQ. ref_name( 2 ) ) atom_name = ref_name( 3 )

       WHERE ( atom_name .EQ. ref_name( 1 ) ) atom_name = ref_name( 2 )

       
       !______________________________________________________________________
       
    ELSEIF ( equiv( conc, 0.0d0, 1.d0-10, fuzzy ) ) THEN  

  !  concentr. x(Ge) = 0  ->  tutti  atomi di mat 2 =  Si 
 
       !______________________________________________________________________
       
       ! if we have the pure material 1 : set all cation names to material 1's
       
!!$       WHERE ( atom_name .EQ. ref_name( 3 ) ) atom_name = ref_name( 2 )

              WHERE ( atom_name .EQ. ref_name( 2 ) ) atom_name = ref_name( 1 )

       !______________________________________________________________________
       
    ELSE

       !______________________________________________________________________
       
       ! otherwhise : set all cation names to material 1's
!    x (Ge) /= 0 , 1  ->     tutti  atomi  =  Si
       
!!$       WHERE ( atom_name .EQ. ref_name( 3 ) ) atom_name = ref_name( 2 )

       WHERE ( atom_name .EQ. ref_name( 2 ) ) atom_name = ref_name( 1 )

       
       !----------------------------------------------------------------------
       
       ! select "index" : cation indexes (to perform the random distribution)
       
       atom_test = .FALSE.

       atom_index = (/ ( i_atom, i_atom = 1, n_atom ) /)

   !      PRINT*, "ref name 1  , 2 , 3 " ,ref_name( 1 ), ref_name( 2 ), ref_name( 3 )
!  PRINT*, "atom_name" ,atom_name

!!$       WHERE ( atom_name .EQ. ref_name( 2 ) ) atom_test = .TRUE.

! lista degli atomi = Si 
       WHERE ( atom_name .EQ. ref_name( 1 ) ) atom_test = .TRUE.


       
       n_index = COUNT( atom_test )

 !PRINT*, " n_index= num Si  atoms =  " , n_index
       
       ALLOCATE( index( n_index ) )
       
       index = PACK( atom_index, atom_test )  !  numeri  dei soli atomi di Si (cation)

 !PRINT*, " index  =  " , index  


       !______________________________________________________________________

    END IF
        
    !=========================================================================

    IF ( .NOT. equiv( conc, 1.0d0, 1.d0-10, fuzzy ) .AND. &
         .NOT. equiv( conc, 0.0d0, 1.d0-10, fuzzy ) ) THEN
       
       !______________________________________________________________________
       
       ! perform random cation distribution when the material is not pure :
       
       !______________________________________________________________________

       ! calculate "n_random" : the number of cations to replace :

!  PRINT*, "n_index " ,n_index

       alloy_test( 1 ) = REAL( INT( conc * n_index ) )
       alloy_test( 2 ) = REAL( INT( conc * n_index ) + 1 )

   ! PRINT*, "conc,alloy_test  " ,conc, alloy_test

       alloy_test = ABS( alloy_test / REAL( n_index ) - conc )
       
       IF ( alloy_test( 1 ) .LT. alloy_test( 2 ) ) THEN
          n_random = INT( conc * n_index )
       ELSE
          n_random = INT( conc * n_index ) + 1
       END IF


    !      PRINT*, "n Random,n_index ", n_random, n_index  !  0 e  0 per  SiGe !!!??

       conc = DBLE( n_random ) / DBLE( n_index )




     !  PRINT*
      ! PRINT*, "Random alloy treatment"
      ! PRINT*
      ! PRINT*, "Concentration = ", conc
      ! PRINT*
      ! PRINT*, "Total number of cations = ", n_index
      ! PRINT*
      ! PRINT*, "Replaced cations = ", n_random
      ! PRINT*
       
       !______________________________________________________________________
       
       i_random = 1
       
       DO WHILE ( i_random .LE. n_random )
          
          !-------------------------------------------------------------------
          
          j = 0
          
          CALL RANDOM_NUMBER( rdm )
          random_limit = INT( rdm * 99 + 1 )          
          
          CALL RANDOM_NUMBER( rdm )
          random_test = INT( rdm * random_limit + 1 )
          
          DO WHILE ( j .LE. random_test )        
             j =  j + 1
             CALL RANDOM_NUMBER( rdm )
          END DO
          
          !-------------------------------------------------------------------
          
          random = INT( rdm * ( n_index - 1 ) + 1 )
          
          IF ( index( random ) .NE. 0 ) THEN
             
             i_random = i_random + 1

!!$             atom_name( index( random ) ) = ref_name( 3 )

! random substitution  of Si  with  Ge
 
             atom_name( index( random ) ) = ref_name( 2 )


             index( random ) = 0
             
          END IF
          
          !-------------------------------------------------------------------
          
       END DO
       
       DEALLOCATE( index )
       
       !______________________________________________________________________
       
    END IF
    
    !=========================================================================
    
  END SUBROUTINE random_alloy_uncommon
  
  !===========================================================================


END MODULE alloys




