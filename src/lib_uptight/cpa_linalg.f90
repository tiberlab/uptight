! cpa_linalg.f90
!
! Module wrapper cho LAPACK ZGETRF/ZGETRI dung de invert complex dense matrix.
! Duoc dung boi cpa_solver.f90 de tinh Green's function:
!     G(k,z) = [ z - Sigma(z) - H_K(k) ]^{-1}
!
MODULE cpa_linalg

  USE precision, ONLY : dp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: invert_complex_matrix

CONTAINS

  !===========================================================================
  ! Subroutine invert_complex_matrix
  !
  ! Invert ma tran vuong phuc A (n x n), in-place.
  !
  ! INPUT/OUTPUT:
  !   A(n,n)  : ma tran can invert (input). Bi ghi de boi A^{-1} (output).
  !   n       : kich thuoc ma tran
  !
  ! OUTPUT:
  !   ierr    : = 0 neu thanh cong
  !             > 0 neu ma tran singular (U(ierr,ierr) = 0 trong LU)
  !             < 0 neu tham so thu |ierr| bi loi (ZGETRF/ZGETRI)
  !===========================================================================

  SUBROUTINE invert_complex_matrix( A, n, ierr )

    INTEGER,                     INTENT( IN )    :: n
    COMPLEX( dp ), DIMENSION(n,n), INTENT( INOUT ) :: A
    INTEGER,                     INTENT( OUT )   :: ierr

    !=========================================================================
    ! Local variables
    !=========================================================================

    INTEGER, DIMENSION(n) :: ipiv

    COMPLEX( dp ), DIMENSION(:), ALLOCATABLE :: work
    COMPLEX( dp ) :: work_query(1)

    INTEGER :: lwork, info

    !=========================================================================
    ! LU decomposition : A = P * L * U
    !=========================================================================

    CALL ZGETRF( n, n, A, n, ipiv, info )

    IF ( info .NE. 0 ) THEN
       ierr = info
       RETURN
    END IF

    !=========================================================================
    ! Workspace query cho ZGETRI (LWORK = -1)
    !=========================================================================

    CALL ZGETRI( n, A, n, ipiv, work_query, -1, info )

    IF ( info .NE. 0 ) THEN
       ierr = info
       RETURN
    END IF

    lwork = INT( REAL( work_query(1), dp ) )
    lwork = MAX( lwork, 1 )

    ALLOCATE( work( lwork ) )

    !=========================================================================
    ! Tinh A^{-1} tu LU decomposition
    !=========================================================================

    CALL ZGETRI( n, A, n, ipiv, work, lwork, info )

    DEALLOCATE( work )

    ierr = info

  END SUBROUTINE invert_complex_matrix

END MODULE cpa_linalg
