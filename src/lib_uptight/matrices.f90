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
!        Module "matrices" - (c) Jerome Gleize - 2002
!
!=============================================================================

MODULE matrices

  !===========================================================================

  !USE parameter

  !===========================================================================
  USE exceptions
  IMPLICIT NONE
  PRIVATE
  !===========================================================================

  PUBLIC :: determinant, transp, inverse

CONTAINS

  !===========================================================================

  RECURSIVE FUNCTION determinant( mat ) RESULT( det )

    !=========================================================================
    
    ! input arguments :
    
    DOUBLE PRECISION, DIMENSION( :, : ), INTENT( IN ) :: mat
    
    !_________________________________________________________________________
    
    ! output result :
    
    DOUBLE PRECISION :: det

    !_________________________________________________________________________

    ! local variables :
    
    DOUBLE PRECISION, &
         DIMENSION( SIZE( mat, 1 ) - 1, SIZE( mat, 2 ) - 1 ) :: mat_i

    DOUBLE PRECISION :: parity

    INTEGER :: n_row, n_col
    INTEGER :: i_row, j_row, k_row

    !=========================================================================

    det = 0.0d0

    n_row = SIZE( mat, 1 )
    n_col = SIZE( mat, 2 )

    !=========================================================================

    IF ( n_row .NE. n_col ) call throw_solve_exception(ERR_MAT_NOTSQRE)
    
    !=========================================================================

    SELECT CASE ( n_row )

       !______________________________________________________________________

       CASE ( 1 )
       
       det = mat( 1, 1 )

       !______________________________________________________________________

       CASE ( 2 )

       det = mat( 1, 1 ) * mat( 2, 2 ) - mat( 2, 1 ) * mat( 1, 2 )

       !______________________________________________________________________

    CASE DEFAULT
       
       !______________________________________________________________________

       parity = - 1.0d0

       !______________________________________________________________________

       DO i_row = 1, n_row

          !-------------------------------------------------------------------

          parity = parity * ( -1.0d0 )

          k_row = 0

          !-------------------------------------------------------------------

          DO j_row = 1, n_row

             !................................................................

             IF ( j_row .EQ. i_row ) CYCLE

             !................................................................

             k_row = k_row + 1
             
             mat_i( k_row, : ) = mat( j_row, 2 : n_row )

             !................................................................

          END DO

          !-------------------------------------------------------------------

          det = det + parity * mat( i_row, 1 ) * determinant( mat_i )
          
          !-------------------------------------------------------------------
          
       END DO

       !______________________________________________________________________
       
    END SELECT
    
    !=========================================================================
    
  END FUNCTION determinant


  !===========================================================================


  FUNCTION transp( mat )

    !=========================================================================

    ! input arguments :
    
    DOUBLE PRECISION, DIMENSION( :, : ), INTENT( IN ) :: mat

    !_________________________________________________________________________

    ! output result :

    DOUBLE PRECISION, DIMENSION( SIZE( mat, 2 ), SIZE( mat, 1 ) ) :: transp

    !_________________________________________________________________________

    ! local variables

    INTEGER :: i_row

    !=========================================================================

    DO i_row = 1, SIZE( mat, 2 )

       transp( i_row, : ) = mat( :, i_row )

    END DO

    !=========================================================================

  END FUNCTION transp

  
  !===========================================================================


  FUNCTION inverse( mat )

    !=========================================================================

    ! input arguments :

    DOUBLE PRECISION, DIMENSION( :, : ), INTENT( IN ) :: mat

    !_________________________________________________________________________

    ! output result :

    DOUBLE PRECISION, DIMENSION( SIZE( mat, 1 ), SIZE( mat, 2 ) ) :: inverse

    !_________________________________________________________________________

    ! local variables :

    DOUBLE PRECISION, DIMENSION( SIZE( mat, 1 ), SIZE( mat, 2 ) ) :: mat_i
    
    DOUBLE PRECISION :: det

    INTEGER :: i_row, i_col, n_row, n_col

    !=========================================================================

    n_row = SIZE( mat, 1 )
    n_col = SIZE( mat, 2 )

    !_________________________________________________________________________

    IF ( n_row .NE. n_col ) call throw_solve_exception(ERR_MAT_NOTSQRE)

    !_________________________________________________________________________

    det = determinant( mat )

    IF ( ABS( det ) .LE. emach ) call throw_solve_exception(ERR_MAT_NOINV) 
    
    !=========================================================================

    SELECT CASE ( n_row )

       !______________________________________________________________________

       CASE ( 1 )

       inverse = 1 / mat( 1, 1 )
       
       !______________________________________________________________________

    CASE DEFAULT

       !______________________________________________________________________
       
       DO i_row = 1, n_row

          !.................................................................

          DO i_col = 1, n_col

             mat_i = mat
             mat_i( i_row, : ) = 0.0d0
             mat_i( :, i_col ) = 0.0d0
             mat_i( i_row, i_col ) = 1.0d0
             inverse( i_row, i_col ) = determinant( mat_i )
             
          END DO
          
          !.................................................................
          
       END DO
       
       inverse = transp( inverse / det )
       
       !-------------------------------------------------------------------
       
    END SELECT

    !=========================================================================

  END FUNCTION inverse

  !===========================================================================

END MODULE matrices

!=============================================================================
