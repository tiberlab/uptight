!=============================================================================
!
!            Module "input_output" - (c) Jerome Gleize - August 2003
!
!=============================================================================
!
! contains routines related to file handling :
!
! => open_file
! => label_in_file
! => int_in_file
! => flag_in_file
! => int_format
! => int_to_char
!
!=============================================================================

MODULE input_output

  !===========================================================================
    
  IMPLICIT NONE

  !===========================================================================
  
CONTAINS
  
  !===========================================================================
  !
  ! Subroutine open_file( file_name, file_num, operation, &
  !     format_flag, output_flag, replace_flag )
  !
  !===========================================================================
  !
  ! Opens a file with error handling and screen output.
  !
  ! Since many files are opened by the "UPTIGHT" software, this routine
  ! makes sure that a distinct numerical handle is associated with each.
  !
  ! The standard options of the FORTRAN 90 "open" function are included. 
  !
  !===========================================================================
  !
  ! * input arguments :
  !
  !   --> file_name : name of the file to open
  !   --> operation : operation to perform on the file : "read" or "write"
  !   --> format_flag : TRUE => formatted file / FALSE => unformatted file
  !
  !---------------------------------------------------------------------------
  !
  ! * optional input arguments :
  !
  !   --> output_flag  : TRUE => screen output / FALSE => no screen output
  !   --> replace_flag : for writing operations only
  !                      TRUE  => if the file already exists, replace it 
  !                      FALSE => if the file already exists, write at the end
  !
  !---------------------------------------------------------------------------
  !
  ! * output arguments :
  !
  !   --> file_num : numerical handle to the file
  !
  !       it is determined by the program to avoid horrible errors
  !
  !===========================================================================
  !
  ! NOTE : This module does not call for other modules and can be used
  !        as it is in for other applications. 
  !
  !===========================================================================
  
  SUBROUTINE open_file( file_name, file_num, operation, &
       format_flag, output_flag, replace_flag )
    
    !=========================================================================
    ! Input arguments
    !=========================================================================
    
    CHARACTER ( LEN = * ), INTENT( IN ) :: file_name, operation
    
    LOGICAL, INTENT( IN ) :: format_flag
    
    LOGICAL, INTENT( IN ), OPTIONAL :: replace_flag, output_flag
    
    !=========================================================================
    ! Output arguments
    !=========================================================================
    
    INTEGER, INTENT( OUT ) :: file_num
    
    !=========================================================================
    ! Local variables
    !=========================================================================

    CHARACTER ( LEN = 500 ) :: write_format

    INTEGER :: err

    LOGICAL :: output, replace, open_test, exist_test
    
    !=========================================================================
    ! Initialization
    !=========================================================================
  
    output = .TRUE.
    IF ( PRESENT( output_flag ) ) output = output_flag

    replace = .TRUE.
    IF ( PRESENT( replace_flag ) ) replace = replace_flag
    
    !=========================================================================
    ! Screen output
    !=========================================================================

    IF ( output ) &
         WRITE ( *, '( /, "Open file : ", a, " ..." )' ) TRIM( file_name )
           
    !=========================================================================
    ! Check if file is not already opened
    !=========================================================================
    open_test = .TRUE.
    
    INQUIRE( FILE = file_name, OPENED = open_test )
   
    IF ( open_test ) &
         WRITE ( *, '( /, "ERROR : file ''", a, "'' is already opened." )' ) &
         TRIM( file_name )
    !=========================================================================
    ! Find a numerical handle to associate to the file
    !=========================================================================
    
    file_num = 0
    exist_test = .TRUE.
    
    DO WHILE ( exist_test )
       
       file_num = file_num + 1
       INQUIRE( UNIT = file_num, OPENED = exist_test )
   
    END DO
    
    !=========================================================================
    ! Open the file
    !=========================================================================

    SELECT CASE ( TRIM( operation ) )

       !======================================================================

       CASE ( 'read' )       ! Read
       
       !----------------------------------------------------------------------
       ! From a formatted file
       !----------------------------------------------------------------------
       IF ( format_flag ) THEN
          
          OPEN( UNIT = file_num, FILE = file_name, &
               FORM = 'FORMATTED', STATUS = 'OLD', IOSTAT = err )
          
          !-------------------------------------------------------------------
          ! From an unformatted file
          !-------------------------------------------------------------------
          
       ELSE 
          OPEN( UNIT = file_num, FILE = file_name, &
               FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = err )
       END IF

       !======================================================================
             
       CASE ( 'write' )       ! Write

       !----------------------------------------------------------------------
       ! From a formatted file
       !----------------------------------------------------------------------
       
       IF ( format_flag ) THEN
          
          IF ( replace ) THEN
             
             OPEN( UNIT = file_num, FILE = file_name, STATUS = 'REPLACE', &
                  FORM = 'FORMATTED', IOSTAT = err )
          
          ELSE
             
             OPEN( UNIT = file_num, FILE = file_name, STATUS = 'UNKNOWN', &
                  POSITION = "APPEND", FORM = 'FORMATTED', IOSTAT = err )
             
          END IF
          
       ELSE
          
          !-------------------------------------------------------------------
          ! From an unformatted file
          !-------------------------------------------------------------------
          
          IF ( replace ) THEN
             
             OPEN( UNIT = file_num, FILE = file_name, STATUS = 'REPLACE', &
                  FORM = 'UNFORMATTED', IOSTAT = err )
             
          ELSE
             
             OPEN( UNIT = file_num, FILE = file_name, STATUS = 'UNKNOWN', &
                  POSITION = "APPEND", FORM = 'UNFORMATTED', IOSTAT = err )
             
          END IF
          
       END IF
    
       !======================================================================
   
    END SELECT
    
    !=========================================================================
    ! Error message
    !=========================================================================

    IF ( err .NE. 0 ) THEN

       write_format = ' '
       write_format = '( /, "Opening of file ''", a, "'' failed.", /,  &
            &"Check file name / path." )'
       WRITE ( *, TRIM( write_format ) ) TRIM( file_name )
       STOP

    END IF

    !=========================================================================

  END SUBROUTINE open_file

  
  !===========================================================================
  !
  ! Subroutine label_in_file( file_num, label, error, advance )
  !
  !===========================================================================
  !
  ! Find a label in a data file.
  !
  ! Used while reading data, to go at the relevant location in the file.
  !  
  ! Two syntaxes are possible : 
  !
  ! => a group of data follows the label after one line feed.
  ! => the data immediately follow the label on the same line.
  ! 
  ! The optional argument "advance" is used to choose between the two.
  !
  !===========================================================================
  !
  ! * input arguments :
  !
  !   --> file_num : numerical handle to the file.
  !   --> label    : label to find.
  !
  !---------------------------------------------------------------------------
  !
  ! * optional input arguments :
  !
  !   --> advance  : syntax control flag.
  !
  !---------------------------------------------------------------------------
  !
  ! * output arguments :
  !
  !   -->> error : error flag.
  !
  !===========================================================================

  SUBROUTINE label_in_file( file_num, label, error, advance )
    
    !=========================================================================
    ! input arguments 
    !=========================================================================

    CHARACTER ( LEN = * ), INTENT( IN )           :: label
    CHARACTER ( LEN = * ), INTENT( IN ), OPTIONAL :: advance
    
    INTEGER, INTENT( IN ) :: file_num
    
    !=========================================================================
    ! output arguments
    !=========================================================================
    
    LOGICAL, INTENT( OUT ) :: error

    !=========================================================================
    ! local variables
    !=========================================================================

    CHARACTER ( LEN = LEN( label ) ) :: list
    CHARACTER ( LEN = 1 )            :: char

    INTEGER :: err
    LOGICAL :: adv, blank_char

    !=========================================================================
    ! Initialization
    !=========================================================================

    REWIND( file_num )
    error = .FALSE.

    !=========================================================================
    ! Select syntax
    !=========================================================================

    adv = .TRUE.
    IF ( PRESENT( advance ) ) THEN
       IF ( TRIM( advance ) .EQ. 'no' ) adv = .FALSE.
    END IF
    
    !=========================================================================
    ! Start of label search
    !=========================================================================

    DO

       !======================================================================
       ! Search initialization
       !======================================================================

       list = ' '
       char = ' '
       blank_char = .FALSE.

       !======================================================================
       ! Read a string of the same size as the label, on the current line
       !======================================================================

       DO WHILE( LEN_TRIM( list ) .LT. LEN_TRIM( label ) )

          !-------------------------------------------------------------------
          !
          ! Read characters on the line until the next non-blank character.
          !
          ! If the end of the line is reached -> jump to label 10 ( cycle )
          ! If the end of the file is reached -> jump to label 20 ( break )
          !
          !-------------------------------------------------------------------

          DO WHILE( char .EQ. ' ' )
             
             READ( file_num, '( a )', ADVANCE = 'NO', EOR = 10, END = 20 ) char
             
             IF ( char .EQ. ' ' ) blank_char = .TRUE.
             
          END DO

          !-------------------------------------------------------------------
          ! Trim extra blank characters to avoid syntax errors :
          ! two non-blank characters are separated by a single blank
          !-------------------------------------------------------------------
          
          IF ( LEN_TRIM( list ) .LT. LEN_TRIM( label ) ) THEN
             
             IF ( blank_char ) THEN
                
                list = TRIM( list ) // ' ' // char
                
             ELSE
                
                list = TRIM( list ) // char
                
             END IF
             
             char = ' '
             
             blank_char = .FALSE.
             
          END IF
          
          !===================================================================
          ! End of string read
          !===================================================================

       END DO

       !======================================================================
       !       
       ! Compare the string to the label and stop if the label is found
       !
       ! Note : force a line feed if the data is on the following line
       !
       !======================================================================
      
       IF ( TRIM( list ) .EQ. TRIM( label ) ) THEN

          IF ( adv ) READ( file_num, *, END = 20 ) 
          EXIT
          
       END IF

       !======================================================================
       ! Go to the next line
       !======================================================================

       READ( file_num, *, END = 20 )

       !======================================================================
       ! End of label search
       !======================================================================

10  END DO

    !=========================================================================
    ! If the label is not found --> activate the error flag
    !=========================================================================

20  IF ( TRIM( list ) .NE. TRIM( label ) ) error = .TRUE.

    !=========================================================================

  END SUBROUTINE label_in_file


  !===========================================================================
  !
  ! Function int_in_file( file_name, file_num, label )
  !
  !===========================================================================
  !
  ! Read an integer value in a file with the syntax :
  !
  !              label = integer value,
  !
  ! where label is a character string associated to the data.
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> file_name  : name of the data file.
  !   --> file_num   : numerical handle to the data file.
  !   --> label      : label of the data to find.
  !
  !---------------------------------------------------------------------------
  !
  ! * output arguments
  !
  !   --> int_in_file : integer value read from the file.
  !
  !===========================================================================
  !
  ! NOTE : this function calls the subroutine "label_in_file" 
  !        and is actually a shortcut of the latter to read 
  !        a single integer value.
  !
  !===========================================================================
  
  FUNCTION int_in_file( file_name, file_num, label )

    !=========================================================================
    ! input arguments
    !=========================================================================

    CHARACTER ( LEN = * ), INTENT( IN ) :: file_name, label
    INTEGER,               INTENT( IN ) :: file_num
    
    !=========================================================================
    ! output result
    !=========================================================================
    
    INTEGER :: int_in_file

    !=========================================================================
    ! local variables
    !=========================================================================

    CHARACTER ( LEN = LEN_TRIM( label ) + 2 ) :: list

    CHARACTER ( LEN = 500 ) :: write_format

    INTEGER :: err

    LOGICAL :: error

    !=========================================================================
    ! Find the label in the file
    !=========================================================================

    list = label // ' ='

    CALL label_in_file( file_num, list, error, advance = 'no' )

    !=========================================================================
    ! If the label is not found --> display an error message
    !=========================================================================

    IF ( error ) THEN

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list ''", a, "'' missing.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file_name ), TRIM( label )

    END IF

    !=========================================================================
    ! Read the integer data
    !=========================================================================

    READ( file_num, *, IOSTAT = err  ) int_in_file

    !=========================================================================
    ! If the value is not an integer --> display an error message
    !=========================================================================
    
    IF ( err .NE. 0 ) THEN

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /, &
            &"An integer value is required for data ''", a, "'' .", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file_name ), TRIM( label )

       STOP

    END IF
         
    !=========================================================================

  END FUNCTION int_in_file


  !===========================================================================
  !
  ! Function flag_in_file( file_name, file_num, label )
  !
  !===========================================================================
  !
  ! Read a logical flag in a file with the syntax
  !
  !                   label = logical flag,
  !
  ! where label is a character string associated to the data.
  !
  !===========================================================================
  !
  ! * input arguments
  !
  !   --> file_name : name of the data file.
  !   --> file_num  : numerical handle to the data file.
  !   --> label     : label of the data to find.
  !
  !---------------------------------------------------------------------------
  !
  ! * output arguments
  !
  !   --> flag_in_file - logical vlaue to read from the file.
  !
  !===========================================================================
  !
  ! NOTE : this function calls the subroutine "label_in_file" 
  !        and is actually a shortcut of the latter to read 
  !        a single logical value.
  !
  !===========================================================================
  
  FUNCTION flag_in_file( file_name, file_num, label )

    !=========================================================================
    ! Input arguments
    !=========================================================================

    CHARACTER ( LEN = * ), INTENT( IN ) :: file_name, label
    INTEGER,               INTENT( IN ) :: file_num
    
    !=========================================================================
    ! Output result
    !=========================================================================
    
    LOGICAL :: flag_in_file

    !=========================================================================
    ! Local variables
    !=========================================================================

    CHARACTER ( LEN = LEN_TRIM( label ) + 2 ) :: list

    CHARACTER ( LEN = 500 ) :: write_format

    INTEGER :: err

    LOGICAL :: error

    !=========================================================================
    ! Find the label in the file 
    !=========================================================================

    list = label // ' ='

    CALL label_in_file( file_num, list, error, advance = 'no' )

    !=========================================================================
    ! If the label is not found --> display an error message
    !=========================================================================

    IF ( error ) THEN

       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /,&
            & "Data list ''", a, "'' missing.", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file_name ), TRIM( label )
       
    END IF
    
    !=========================================================================
    ! Read the logical data
    !=========================================================================

    READ( file_num, '( l10 )', IOSTAT = err  ) flag_in_file

    !=========================================================================
    ! If the value is not logical --> display an error message
    !=========================================================================
    
    IF ( err .NE. 0 ) THEN
       
       write_format = ' '
       write_format = '( /, "SYNTAX ERROR in file ''", a, "'' :", /, &
            &"A logical flag value is required for data ''", a, "'' .", / )'

       WRITE ( *, TRIM( write_format ) ) TRIM( file_name ), TRIM( label )

       STOP

    END IF
         
    !=========================================================================

  END FUNCTION flag_in_file


  !===========================================================================
  !
  ! Function int_format( i )
  !
  !===========================================================================
  !
  ! Selects a correct output format string for integer values 
  ! ... just for the sake of displaying a good looking output on the screen
  !
  !===========================================================================
  !
  ! * input arguments :
  !
  !   --> i : integer value
  !
  !---------------------------------------------------------------------------
  !
  ! * output value :
  !
  !   --> int_format : format strin
  !
  !===========================================================================
  
  FUNCTION int_format( i )
    
    !=========================================================================
    ! input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: i
    
    !=========================================================================
    ! output result
    !=========================================================================
    
    CHARACTER ( LEN = 2 ) :: int_format
    
    !=========================================================================
    ! Adjust the format string to the size of the integer value
    !=========================================================================

    SELECT CASE ( i )

       !======================================================================

       CASE ( 0 : 9 )
       int_format = "i1"
       
       CASE ( 10 : 99 )
       int_format = "i2"

       CASE ( 100 : 999 )
       int_format = "i3"
 
       CASE ( 1000 : 9999 )
       int_format = "i4"

       CASE ( 10000 : 99999 )
       int_format = "i5"

       CASE ( 100000 : 999999 )
       int_format = "i6"

       CASE ( 1000000 : 9999999 )
       int_format = "i7"

       CASE DEFAULT
       int_format = "i8"       

       !======================================================================

    END SELECT

    !=========================================================================

  END FUNCTION int_format

  !===========================================================================
  !
  ! Subroutine int_to_char( i )
  !
  !===========================================================================
  !
  ! Convert an integer to a character string for screen output purposes.
  !
  !===========================================================================
  !
  ! * input arguments :
  !
  !   --> i : integer value
  !
  !---------------------------------------------------------------------------
  !
  ! * output value :
  !
  !   --> int_to_char : character representation of the integer value
  ! 
  !===========================================================================

  FUNCTION int_to_char( i )

    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: i

    !=========================================================================
    ! Output result
    !=========================================================================

    CHARACTER ( LEN = 10 ) :: int_to_char

    !=========================================================================
    ! Local variables
    !=========================================================================

    INTEGER :: factor, i_digit, i_counter

    !=========================================================================
    ! Initialization
    !=========================================================================

    int_to_char = ' '
    i_digit  = 1
    i_counter = 1
    factor = 1

    !=========================================================================
    ! Compute the number of digits in the integer value
    !=========================================================================

    DO WHILE ( factor .LE. i )
       
       factor = factor * 10
       
    END DO

    i_digit = i

    !=========================================================================
    !
    ! Convert each digit to a character :
    !
    !
    !             0  0  0  2  5  0  0  0  0  0    -->  integer to convert
    !
    !            10  9  8  7  6  5  4  3  2  1    -->  digits rank
    !                      |
    !                      |_ i_digit = 2 500 000
    !                         factor  = 1 000 000
    !
    ! -> the factor is seet to the current digit  --> 1 000 000
    !                    
    ! -> the integer value "i_digit / factor" 
    !    gives the leading integer digit      -->  2 500 000 / 1 000 000 = 2  
    !
    ! -> it is converted in an ASCII character 
    !    by using the FORTRAN 90 function CHAR
    !
    ! -> i_digit is updated to the following digit --> 2 500 000 - 2 000 000
    !
    !=========================================================================

    DO WHILE ( factor .GT. 1 )
       
       factor = factor / 10
       int_to_char( 1 : i_counter ) = &
            TRIM( int_to_char ) // CHAR( 48 + i_digit / factor )
       i_digit = i_digit - ( i_digit / factor ) * factor
       i_counter = i_counter + 1

    END DO

    !=========================================================================

  END FUNCTION int_to_char

  !===========================================================================

END MODULE input_output
