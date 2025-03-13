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
!               Module "errors" - (c) Jerome Gleize - August 2003
!
!=============================================================================
!
! contains routines related to errors and warnings :
!
! => read_errors
! => alloc_error
! => dealloc_error
! => lattice_errors
!
!=============================================================================

MODULE errors

  !===========================================================================
  USE exceptions
  USE input_output, only : int_format

  !===========================================================================

  IMPLICIT NONE
  PRIVATE
  !===========================================================================

  PUBLIC :: read_errors, alloc_error, dealloc_error, lattice_errors


CONTAINS

  !===========================================================================
  !
  ! Subroutine read_errors( error_type, file, list )
  !
  !===========================================================================
  !
  ! Displays an error message to the output ( screen ) 
  ! 
  ! The errors considered here are all syntax errors which may occur
  ! while reading data from an input file
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> error_type :  type of error which has occured.
  !   --> file       :  name of the input file involved.
  !
  !---------------------------------------------------------------------------
  !
  ! * optional input arguments :
  !
  !   --> list : label of the data which the program has failed to read.
  !  
  !===========================================================================
  
  SUBROUTINE read_errors( error_type, file, list )

    !=========================================================================
    ! Input arguments
    !=========================================================================

    CHARACTER ( LEN = * ), INTENT( IN )           :: file, error_type
    CHARACTER ( LEN = * ), INTENT( IN ), OPTIONAL :: list
    
    !=========================================================================
    ! Local variables
    !=========================================================================

    CHARACTER ( LEN = 500 ) :: write_format

    !=========================================================================
    ! Select the error type
    !=========================================================================

    SELECT CASE( error_type )

       !======================================================================
       
       CASE ( 'missing' )        ! Missing data label
       
       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list ''", a, "'' missing.", / )'
       
       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )
       
       !======================================================================
       
       CASE ( 'not integer' )    ! Failure to read an integer value
       
       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /, &
            &"An integer value is required for data ''", a, "'' .", / )'
       
       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================
      
       CASE ( 'not logical' )    ! Failure to read a logical flag       
       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /, &
            &"A logical flag is required for data ''", a, "'' .", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'not character string' )  ! Failure to read a character string

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "A character string is required for data ''", a, "''.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )
       
       !======================================================================
       
       CASE ( 'bad pairs' )      ! Wrong syntax for the pair data

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list ''", a, "'' should have the form :", /,&
            & "pair index <blank> names of pair atoms <blank>&
            & number of TB couplings.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad index' )      ! Index overflow in a numbered data list

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "An element index in data list ''", a, "'' is too large.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad compound type' )       ! Wrong compound type

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", //, &
            & "Data list should have the form :  ''", &
            & a, " type of compound'',", /, &
            & "where ''type of compound'' must be one of :", //, &
            & "''simple''     --> simple compound ( e.g. Si )", /, &
            & "''binary''     --> binary compound ( e.g. GaAs )", /, &
            & "''ternary''    --> ternary alloy ( e.g. SiGe or AlGaAs )", /, &
            & "''quaternary'' --> quaternary alloy ( e.g. AlGaAsSb )" )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad parent type' )         ! Wrong parent type

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", //, &
            & "Data list should have the form :  ''", &
            & a, " type of parent'',", /, &
            & "where ''type of parent'' must be one of :", //, &
            & "''simple''     --> simple compound ( e.g. Si )", /, &
            & "''binary''     --> binary compound ( e.g. GaAs )" )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad alloy' )

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "The number of parent materials in alloys cannot exceed 2.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad file' )

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list should have the form :  ''", &
            & a, " data file name''.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad point' )

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list should have the form :   ''", &
            & a, " coordinates of point''.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================
              
       CASE ( 'bad strains' )
       
       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list ''", a, "'' should have the form :", /,&
            & "real accuracy.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )

       !======================================================================

       CASE ( 'bad type' )

       write_format = ' '
       write_format = '( /, "DATA ERROR in file ''", a, "'' :", /,&
            & "wrong type in list ''", a, "''.", / )'
       
       WRITE ( *, TRIM( write_format ) ) TRIM( file ), TRIM( list )
       
       !======================================================================

    END SELECT

    !========================================================================

    call throw_init_exception(ERR_INPUT)

    !=========================================================================

  END SUBROUTINE read_errors



  !===========================================================================
  !
  ! Subroutine alloc_error( module_name, routine_name, name )
  !
  !===========================================================================
  !
  ! Display an error message if an allocation of memory went wrong
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> module_name  : name of the module where the allocation occurs
  !   --> routine_name : name of the routine where the allocation occurs
  !   --> name         : name of the variable whose allocation failed
  !
  !===========================================================================

  SUBROUTINE alloc_error( module_name, routine_name, name )

    !=========================================================================
    ! input arguments
    !=========================================================================

    CHARACTER ( LEN = * ), INTENT( IN ) :: module_name, routine_name, name

    !=========================================================================
    ! Local variables
    !=========================================================================

    CHARACTER ( LEN = 500 ) :: write_format

    !=========================================================================
    
    write_format = ' '
    write_format = '( /, "ALLOCATION ERROR in module ''", a,&
         & "'' - routine ''", a, "'' : ", /,&
         & "allocation of ''", a, "'' failed.", / )'
       
    WRITE ( *, TRIM( write_format ) ) &
         TRIM( module_name ), TRIM( routine_name ), TRIM( name )
 
    call throw_init_exception(ERR_ALLOC_ERR)

    !=========================================================================

  END SUBROUTINE alloc_error


  !===========================================================================
  !
  ! Subroutine dealloc_error( module_name, routine_name, name )
  !
  !===========================================================================
  !
  ! Display an error message if a deallocation of memory went wrong
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> module_name  : name of the module where the deallocation occurs
  !   --> routine_name : name of the routine where the deallocation occurs
  !   --> name         : name of the variable whose deallocation failed
  !
  !===========================================================================
  
  SUBROUTINE dealloc_error( module_name, routine_name, name )

    !=========================================================================
    ! input arguments
    !=========================================================================

    CHARACTER ( LEN = * ), INTENT( IN ) :: module_name, routine_name, name

    !=========================================================================
    ! Local variables
    !=========================================================================

    CHARACTER ( LEN = 500 ) :: write_format

    !=========================================================================
    
    write_format = ' '
    write_format = '( /, "DEALLOCATION ERROR in module ''", a,&
         & "'' - routine ''", a, "'' : ", /,&
         & "deallocation of ''", a, "'' failed.", / )'
       
    WRITE ( *, TRIM( write_format ) ) &
         TRIM( module_name ), TRIM( routine_name ), TRIM( name )
 
    call throw_init_exception(ERR_ALLOC_ERR)    

    !=========================================================================

  END SUBROUTINE dealloc_error



  !===========================================================================
  !
  ! Subroutine lattice_errors( error_type, param, value )
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> error_type   :  type of error encountered
  !
  !---------------------------------------------------------------------------
  !
  ! * optional input arguments
  !
  !   --> param   :  name of the parameter to change
  !   --> value   :  present value of the parameter
  !
  !===========================================================================

  SUBROUTINE lattice_errors( error_type, param, value )

    !=========================================================================
    ! Input arguments
    !=========================================================================

    CHARACTER ( LEN = * ),           INTENT( IN ) :: error_type
    CHARACTER ( LEN = * ), OPTIONAL, INTENT( IN ) :: param
    INTEGER,               OPTIONAL, INTENT( IN ) :: value

    !=========================================================================
    ! Local variables
    !=========================================================================

    CHARACTER ( LEN = 500 ) :: write_format

    !=========================================================================

    SELECT CASE ( error_type )

       !======================================================================

       CASE ( 'exceed limit' )

       write_format = ' '
       write_format = &
            '( /, "PARAMETER ERROR - file ''parameter.f90'' :", /,&
            & "Parameter ''", a, "'' may be too small for the number of&
            & nearest neighbours and/or stars requested.", /,&
            & "Increase its value beyond the value ", ' &
            & // TRIM( int_format( value ) ) // ', "." )'

       WRITE ( *, TRIM( write_format ) ) TRIM( param ), value

       !======================================================================

       CASE ( 'zero distance' )

       write_format = ' '
       write_format = '( /, "ERROR : at least one distance in a star should&
            & be non zero. Check data." )'

       !======================================================================

    END SELECT

    call throw_init_exception(ERR_GENERAL)

    !=========================================================================

  END SUBROUTINE lattice_errors

  !===========================================================================

END MODULE errors
