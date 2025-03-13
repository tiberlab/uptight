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
!                 Module "lists" - (c) Jerome Gleize - August 2003
!
!=============================================================================
!
! contains procedures related to pointer - directed lists :
!
!_____________________________________________________________________________
!
! => create_list
! => add_list
! => delete_list
! => write_list
!
!=============================================================================

MODULE lists
  
  !===========================================================================
  
  !USE TB_ham_parameter
  USE type_defs
  USE list_types
  USE input_output, only : int_to_char, int_format
  USE errors
  
  !===========================================================================
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: create_list, add_list, delete_list, write_list

  !===========================================================================
  !
  ! NOTE : the procedure are generic so that different types of data can be
  !        handled by a single procedure call.
  !
  !===========================================================================
  INTERFACE create_list
     MODULE PROCEDURE create_list_int
     MODULE PROCEDURE create_list_dp
     MODULE PROCEDURE create_list_names
  END INTERFACE

  !___________________________________________________________________________
  
  INTERFACE add_list
     MODULE PROCEDURE add_list_int
     MODULE PROCEDURE add_list_dp 
     MODULE PROCEDURE add_list_names
  END INTERFACE

  !___________________________________________________________________________
  
  INTERFACE delete_list
     MODULE PROCEDURE delete_list_int
     MODULE PROCEDURE delete_list_dp
     MODULE PROCEDURE delete_list_names
  END INTERFACE

  !___________________________________________________________________________

  INTERFACE write_list
     MODULE PROCEDURE write_list_int
     MODULE PROCEDURE write_list_dp
     MODULE PROCEDURE write_list_names
     MODULE PROCEDURE write_list_ion_orbit
     MODULE PROCEDURE write_list_pair_coupling
  END INTERFACE

  !===========================================================================
  
CONTAINS


  !===========================================================================
  !
  ! Subroutines "create_list_*" : create a given list.
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => new - pointer to a new data to be inserted in the list.
  ! => start - pointer to the first element in the list.
  ! => current - pointer to the current element in the list.
  !
  !===========================================================================
    
  SUBROUTINE create_list_int( new, start, current )
    
    !_________________________________________________________________________
    
    TYPE (ints), POINTER :: new, start, current

    !_________________________________________________________________________
    
    start => new
    current => new

    !_________________________________________________________________________
    
  END SUBROUTINE create_list_int

  !===========================================================================

  SUBROUTINE create_list_dp( new, start, current )

    !_________________________________________________________________________
    
    TYPE (dprec), POINTER :: new, start, current

    !_________________________________________________________________________
    
    start => new
    current => new

    !_________________________________________________________________________
    
  END SUBROUTINE create_list_dp

  !===========================================================================

  SUBROUTINE create_list_names( new, start, current )
    
    !_________________________________________________________________________
        
    TYPE (names), POINTER :: new, start, current    

    !_________________________________________________________________________
        
    start => new
    current => new

    !_________________________________________________________________________
        
  END SUBROUTINE create_list_names

  !===========================================================================



  !===========================================================================
  !
  ! Subroutines "add_list_*" : add a new element to a given list.
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => new - pointer to a new data to be inserted in the list.
  ! => current - pointer to the current element in the list.
  !
  !===========================================================================

  SUBROUTINE add_list_int( new, current )

    !_________________________________________________________________________
    
    TYPE (ints), POINTER :: new, current

    !_________________________________________________________________________
        
    current%next => new
    current => current%next

    !_________________________________________________________________________
    
  END SUBROUTINE add_list_int

  !===========================================================================
  
  SUBROUTINE add_list_dp( new, current )

    !_________________________________________________________________________
    
    TYPE (dprec), POINTER :: new, current

    !_________________________________________________________________________
        
    current%next => new
    current => current%next

    !_________________________________________________________________________
    
  END SUBROUTINE add_list_dp

  !===========================================================================

  SUBROUTINE add_list_names( new, current )
        
    !_________________________________________________________________________

    TYPE (names), POINTER :: new, current

    !_________________________________________________________________________
        
    current%next => new
    current => current%next

    !_________________________________________________________________________
    
  END SUBROUTINE add_list_names

  !===========================================================================
  


  !===========================================================================
  !
  ! Subroutines "delete_list_*" : deallocate a given list.
  ! 
  !===========================================================================

  SUBROUTINE delete_list_int( start )

    !_________________________________________________________________________
    
    TYPE (ints), POINTER :: start, current

    !_________________________________________________________________________
    
    DO WHILE ( ASSOCIATED( start%next ) )
       
       current => start%next
       NULLIFY( start%next )
       DEALLOCATE( start%nbr )
       DEALLOCATE( start )
       start => current
    
    END DO

    DEALLOCATE( start%nbr )
    DEALLOCATE( start )

    !_________________________________________________________________________
    
  END SUBROUTINE delete_list_int

  !===========================================================================


  SUBROUTINE delete_list_dp( start )

    !_________________________________________________________________________
    
    TYPE (dprec), POINTER :: start, current

    !_________________________________________________________________________
    
    DO WHILE ( ASSOCIATED( start%next ) )
       
       current => start%next
       NULLIFY( start%next )
       DEALLOCATE( start%nbr )
       DEALLOCATE( start )
       start => current
    
    END DO

    DEALLOCATE( start%nbr )
    DEALLOCATE( start )

    !_________________________________________________________________________
    
  END SUBROUTINE delete_list_dp

  !===========================================================================

  SUBROUTINE delete_list_names( start )

    !_________________________________________________________________________
    
    TYPE (names), POINTER :: start, current

    !_________________________________________________________________________
    
    DO WHILE ( ASSOCIATED( start%next ) )
       
       current => start%next
       NULLIFY( start%next )
       DEALLOCATE( start%name )
       DEALLOCATE( start )
       start => current
    
    END DO

    DEALLOCATE( start%name )
    DEALLOCATE( start )

    !_________________________________________________________________________
    
  END SUBROUTINE delete_list_names



  !===========================================================================
  !
  ! Subroutines "write_list_*" : output a given list to screen.
  ! 
  !===========================================================================

   SUBROUTINE write_list_int( start )

    !_________________________________________________________________________
    
    TYPE (ints), POINTER :: start, test

    !_________________________________________________________________________
    
    test => start

    DO WHILE ( ASSOCIATED( test ) )
    
       WRITE ( *, * ) test%nbr
       test => test%next
    
    END DO

    !_________________________________________________________________________
    
  END SUBROUTINE write_list_int

  !===========================================================================

   SUBROUTINE write_list_dp( start )

    !_________________________________________________________________________
    
    TYPE (dprec), POINTER :: start, test

    !_________________________________________________________________________
    
    test => start

    DO WHILE ( ASSOCIATED( test ) )
    
       WRITE ( *, * ) test%nbr
       test => test%next
    
    END DO

    !_________________________________________________________________________
    
  END SUBROUTINE write_list_dp

  !===========================================================================

  SUBROUTINE write_list_names( start )

    !_________________________________________________________________________
    
    TYPE (names), POINTER :: start, test

    !_________________________________________________________________________
    
    test => start

    DO WHILE ( ASSOCIATED( test ) )

       WRITE ( *, * ) test%name
       test => test%next
  
    END DO

    !_________________________________________________________________________
    
  END SUBROUTINE write_list_names


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  SUBROUTINE write_list_ion_orbit( ion_data )

    !=========================================================================
    
    ! input arguments :

    TYPE (ion_orbit), DIMENSION( : ), POINTER :: ion_data
    
    !_________________________________________________________________________

    ! local variables :

    INTEGER :: i_ion, i_state, i_so, n_state, n_so
    LOGICAL :: test

    CHARACTER (LEN = 350) :: write_format
    !=========================================================================
    
    write_format = ' '
    write_format = '( /, "***", //, "ONSITE TB PARAMETERS" )'

    WRITE ( *, write_format )

    !=========================================================================

    DO i_ion = 1, SIZE( ion_data )

       !______________________________________________________________________

       ! write current ion name :

       write_format = ' '
       write_format = '( /, "---" )'

       WRITE ( *, write_format )

       write_format = ' '
       write_format = &
            '( /, "ion ' // TRIM( int_to_char( i_ion ) ) // ' = ", a2, / )'
       
       WRITE ( *, write_format ) ion_data( i_ion )%name

       !______________________________________________________________________

       ! write current ion's states and onsite energies :

       test = ASSOCIATED( ion_data( i_ion )%state )
       test = test .AND. ASSOCIATED( ion_data( i_ion )%energy )

       IF ( test ) THEN

          !-------------------------------------------------------------------

          write_format = ' '
          write_format = '( /, " TB orbital states and energies ( eV ) :", / )'
       
          WRITE ( *, write_format )

          n_state = SIZE( ion_data( i_ion )%state )
       
          WRITE ( *, '( 2x, a6, 5x, f10.6 )' ) &
               ( ion_data( i_ion )%state( i_state ), &
               ion_data( i_ion )%energy( i_state ), i_state = 1, n_state )
       
          !-------------------------------------------------------------------

       END IF

       !______________________________________________________________________
       
       
       test = ASSOCIATED( ion_data( i_ion )%so_state )
       test = test .AND. ASSOCIATED( ion_data( i_ion )%so_energy )
       
       IF ( test ) THEN
          
          !................................................................
          
          n_so = SIZE( ion_data( i_ion )%so_state )
          
          write_format = ' '
          write_format = &
               '( /, " Spin-orbit orbital and energy ( eV ) :", / )'
          
          WRITE ( *, write_format )
          
          WRITE ( *, '( 2x, a1, " states", 5x, f10.6 )' ) ( &
               ion_data( i_ion )%so_state( i_so ), &
               ion_data( i_ion )%so_energy( i_so ), i_so = 1, n_so )
          
          !................................................................
          
       END IF

       
       !______________________________________________________________________
     
    END DO

    !=========================================================================
    
  END SUBROUTINE write_list_ion_orbit


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  SUBROUTINE write_list_pair_coupling( pair_data )

    !=========================================================================
    
    ! input arguments :

    TYPE (pair_coupling), DIMENSION( : ), POINTER :: pair_data

    !_________________________________________________________________________

    ! local variables :

    INTEGER :: i_pair, i_cpl, n_cpl
    LOGICAL :: test
    CHARACTER (LEN = 350) :: write_format
    !=========================================================================
    
    
    write_format = ' '
    write_format = '( /, "***", //, "TB COUPLING PARAMETERS" )'

    WRITE ( *, write_format )

    !=========================================================================

    DO i_pair = 1, SIZE( pair_data )
       
       !______________________________________________________________________

       ! write current pair name :

       write_format = ' '
       write_format = '( /, "---" )'

       WRITE ( *, write_format )

       write_format = ' '
       write_format = &
            '( /, "pair ' // TRIM( int_to_char( i_pair ) ) // &
            ' = ", a2, x, a2 / )'
       
       WRITE ( *, write_format ) pair_data( i_pair )%name


       ! write current pair order and reference distance :
       
       write_format = ' '
       write_format = &
            '( /,' // TRIM( int_format( pair_data( i_pair )%order ) ) // &
            '" nearest neighbours" )'

       WRITE ( *, write_format ) pair_data( i_pair )%order

       write_format = ' '
       write_format = &
            '( /, "reference distance = ", f10.6, " Angstrom" )'

       WRITE ( *, write_format ) pair_data( i_pair )%dist_ref

       !______________________________________________________________________

       test = ASSOCIATED( pair_data( i_pair )%cpl )
       test = test .AND. ASSOCIATED( pair_data( i_pair )%energy )
       test = test .AND. ASSOCIATED( pair_data( i_pair )%scaling )

       IF ( test ) THEN

          !-------------------------------------------------------------------

          write_format = ' '
          write_format = '( /, " TB couplings energies ( eV ) and ' // &
               'distance scaling exponent :", / )'
       
          WRITE ( *, write_format )
          
          n_cpl = SIZE( pair_data( i_pair )%cpl )
       
          WRITE ( *, '( 2x, a5, 5x, f10.6, 5x, f10.6 )' ) &
               ( pair_data( i_pair )%cpl( i_cpl ), &
               pair_data( i_pair )%energy( i_cpl ), &
               pair_data( i_pair )%scaling( i_cpl ), i_cpl = 1, n_cpl )
          
          !-------------------------------------------------------------------
          
       END IF

       !______________________________________________________________________

      

    END DO

    !_________________________________________________________________________
    
  END SUBROUTINE write_list_pair_coupling



END MODULE lists

!=============================================================================
