!=============================================================================
!
!  1. subroutine setup_david_precond_jac     (sigma)
!  2. subroutine setup_david_precond_jac_real(sigma)
!
!  3. subroutine dealloc_david_precond_jac
!  4. subroutine dealloc_david_precond_jac_real
!
!  4. subroutine transform_sparse_to_ilu1D     (kindC)
!
!  5. subroutine transform_sparse_to_ilu2D     (kindC)
!  6. subroutine transform_sparse_to_ilu2D_real(kindC)
!
!  7. subroutine transform_sparse_to_ilu3D     (kindC)
!  8. subroutine transform_sparse_to_ilu3D_real(kindC)
!
!=============================================================================

MODULE mod_precond_ilu_davidson

  !===========================================================================
  !
  ! Note: Assumption: Hermitian matrix stored in module sparse_matrix.
  !                   Values also include multiplication vector diag_kpV.
  !
  !___________________________________________________________________________
  !
  ! index_rowV( 1 : n + 1 )
  ! value_rowM( 1 : 2, 1 : num_elements ),
  ! num_off_elements = size( ija ) - N - 1
  ! inv_diag_iluV( 1 .. n )  : D^-1
  !
  ! index_rowV:   1  |  2  |  3  |   .........         |  n  |  n + 1
  !
  !             row1   row2
  !
  !              |      |
  !              | /----/
  !              v v 
  !
  ! value_rowM:  position  ( 1 : num_elements ) marked by index_rowV
  !              contains in component1 : column
  !                       in component2 : position of the value in sa
  !
  !===========================================================================
  
  IMPLICIT NONE
  
  !===========================================================================
  
  COMPLEX ( KIND(0.0D0) ), DIMENSION( : ),      ALLOCATABLE :: &
       davidson_precond_jacV

  REAL    ( KIND(0.0D0) ), DIMENSION( : ),      ALLOCATABLE :: &
       davidson_precond_jacV_real, inv_diag_iluV

  INTEGER,                   DIMENSION( : ),    ALLOCATABLE :: index_rowV
  INTEGER,                   DIMENSION( :, : ), ALLOCATABLE :: value_rowM

  COMPLEX ( KIND(0.0D0) ) :: E_sigma

  !===========================================================================
 
CONTAINS

  !===========================================================================
 
  SUBROUTINE setup_david_precond_jac( sigma )
 
    !=========================================================================
    !
    ! Jacobi preconditioning for Jacobi-Davidson method.
    ! Complex version.
    !
    !=========================================================================

    USE sparse_matrix, ONLY : sa, ija

    !=========================================================================

    COMPLEX( KIND(0.0D0) ), INTENT( IN ) :: sigma

    !_________________________________________________________________________

    INTEGER :: istatus, n

    !=========================================================================

    n = ija( 1 ) - 2

    ALLOCATE( davidson_precond_jacV( n ), STAT = istatus )
    
    davidson_precond_jacV = 1.0D0 / ( sa( 1 : n ) - sigma )
  
    !=========================================================================

  END SUBROUTINE setup_david_precond_jac
 
  !===========================================================================

  SUBROUTINE setup_david_precond_jac_real( sigma )
 
    !=========================================================================
    !
    ! Jacobi preconditioning for Jacobi-Davidson method.
    ! Real version.
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    REAL ( KIND(0.0D0) ), INTENT( IN ) :: sigma

    !_________________________________________________________________________

    INTEGER :: istatus, n

    !=========================================================================

    n = ija( 1 ) - 2

    ALLOCATE( davidson_precond_jacV_real( n ), STAT = istatus )

    davidson_precond_jacV_real = 1.0D0 / ( sa( 1 : n ) - sigma )
 
    !=========================================================================

  END SUBROUTINE setup_david_precond_jac_real

  !===========================================================================

  SUBROUTINE dealloc_david_precond_jac

    !=========================================================================
    !
    ! Complex version.
    !
    !=========================================================================

    DEALLOCATE( davidson_precond_jacV )

    !=========================================================================
 
  END SUBROUTINE dealloc_david_precond_jac
 
  !===========================================================================

  SUBROUTINE dealloc_david_precond_jac_real

    !=========================================================================
    !
    ! Real version.
    !
    !=========================================================================
 
    DEALLOCATE( davidson_precond_jacV_real )

    !=========================================================================

  END SUBROUTINE dealloc_david_precond_jac_real

  !===========================================================================

  SUBROUTINE transform_sparse_to_ilu1D( kindC )

    !=========================================================================
    !
    ! NOTE: The matrix to be diagonalized is D^(1/2) H D^(1/2)
    !       In the setup, the matrix H is stored as sparse matrix.
    !       In davidson routines, the routine amul directly calls
    !       the sparse matrix-vector product. So the diagonal 
    !       matrices D^(1/2) , stored in diag_kpV, have to be 
    !       taken into account. 
    !
    !=========================================================================

    USE normalization_factors_sg1D, ONLY : diag_kpV
    USE sparse_matrix,              ONLY : sa, ija

    !=========================================================================

    INTEGER :: i, k, n
  
    CHARACTER( LEN = 2 ) :: kindC

    !=========================================================================
 
    SELECT CASE ( kindC )

       !______________________________________________________________________

       CASE ( 'kp' )

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n
          
          sa( i ) = sa( i ) * CMPLX( diag_kpV( i )**2 )
          
          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * CMPLX( diag_kpV( i ) ) * &
                  CMPLX( diag_kpV( ija( k ) ) )
          END DO
          
       END DO
       
       !______________________________________________________________________

    CASE default
 
       !______________________________________________________________________

       STOP 'error transform_sparse_to_ilu1D: no davidson with &
            & sg+magnetic in1D'

       !______________________________________________________________________

    END SELECT

    !=========================================================================
    
  END SUBROUTINE transform_sparse_to_ilu1D
  
  !===========================================================================
  
  SUBROUTINE transform_sparse_to_ilu2D( kindC )
 
    !=========================================================================
    !
    ! NOTE: The matrix to be diagonalized is D^(1/2) H D^(1/2).
    !       In the setup, the matrix H is stored as a sparse matrix.
    !       In Davidson routines, the routine "amul" directly calls
    !       the sparse matrix-vector product. So the diagonal 
    !       matrices D^(1/2) have to be taken into account. 
    !
    ! kindC = 'kp'  --> If it is a k.p matrix: D^(1/2) is stored in
    !                   diag_kpV.
    !       = 'sg'  --> If it is a single band Schroedinger matrix
    !                   with magnetic field: D^(1/2) is stored 
    !                   in diag_sgV.
    !
    ! Complex version.
    !
    !=========================================================================

    USE normalization_factors_sg2D, ONLY : diag_kpV, diag_sgV
    USE sparse_matrix             , ONLY : sa, ija

    !=========================================================================

    INTEGER :: i, k, n

    CHARACTER( LEN = 2 ) :: kindC
   
    !=========================================================================
 
    SELECT CASE ( kindC )
       
       !______________________________________________________________________

       CASE ( 'kp' )

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n

          sa( i ) = sa( i ) * CMPLX( diag_kpV( i )**2 )

          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * CMPLX( diag_kpV( i ) ) * &
                  CMPLX( diag_kpV( ija( k ) ) )
          END DO
         
       END DO

       !______________________________________________________________________

       CASE ( 'sg' )

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n

          sa( i ) = sa( i ) * CMPLX( diag_sgV( i )**2 )
         
          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * CMPLX( diag_sgV( i ) ) * &
                  CMPLX( diag_sgV( ija( k ) ) )
          END DO

       END DO

       !______________________________________________________________________

    CASE DEFAULT
 
       !______________________________________________________________________

       STOP 'Error transform_sparse_to_ilu2D: kindC not allowed.'

       !______________________________________________________________________

    END SELECT

    !=========================================================================
 
  END SUBROUTINE transform_sparse_to_ilu2D
  
  !===========================================================================
 
  SUBROUTINE transform_sparse_to_ilu2D_real( kindC )
 
    !=========================================================================
    !
    ! NOTE: The matrix to be diagonalized is D^(1/2) H D^(1/2).
    !       In the setup, the matrix H is stored as a sparse matrix.
    !       In davidson routines, the routine "amul" directly calls
    !       the sparse matrix-vector product. So the diagonal 
    !       matrices D^(1/2) have to be taken into account. 
    !
    ! kindC = 'kp'  --> Not useful in real version.
    !       = 'sg'  --> If it is a single band Schroedinger matrix
    !                   with magnetic field: D^(1/2) is stored 
    !                   in diag_sgV.
    !
    ! Real version.
    !
    !=========================================================================

    USE normalization_factors_sg2D, ONLY : diag_sgV
    USE sparse_matrix_real        , ONLY : sa, ija

    !=========================================================================

    INTEGER :: i, k, n

    CHARACTER( LEN = 2 ) :: kindC

    !=========================================================================
    
    SELECT CASE ( kindC )
    
       !______________________________________________________________________

       CASE ( 'kp' )
  
       !______________________________________________________________________

       STOP 'Error transform_sparse_to_ilu2D_real: &
            & Real version not useful for complex k.p matrix.'

       !______________________________________________________________________

       CASE ( 'sg' )

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n

          sa( i ) = sa( i ) * ( diag_sgV( i )**2 )

          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * diag_sgV( i ) * diag_sgV( ija( k ) )
          END DO

       END DO

       !______________________________________________________________________

    CASE DEFAULT
 
       !______________________________________________________________________

       STOP 'Error transform_sparse_to_ilu2D_real: kindC not allowed.'

       !______________________________________________________________________
       
    END SELECT

    !=========================================================================
 
  END SUBROUTINE transform_sparse_to_ilu2D_real
 
  !===========================================================================
 
  SUBROUTINE transform_sparse_to_ilu3D( kindC )
 
    !=========================================================================
    !
    ! NOTE: The matrix to be diagonalized is D^(1/2) H D^(1/2).
    !       In the setup, the matrix H is stored as a sparse matrix.
    !       In Davidson routines, the routine amul directly calls
    !       the sparse matrix-vector product. So the diagonal 
    !       matrices D^(1/2) have to be taken into account. 
    !
    ! kindC = 'kp'  --> If it is a k.p matrix: D^(1/2) is stored in
    !                   diag_kpV.
    !       = 'sg'  --> If it is a single band Schroedinger matrix
    !                   with magnetic field: D^(1/2) is stored 
    !                   in diag_sgV.
    !
    ! Complex version.
    !
    !=========================================================================
 
    USE normalization_factors_sg3D, ONLY : diag_kpV, diag_sgV
    USE sparse_matrix             , ONLY : sa, ija

    !=========================================================================

    INTEGER :: i, k, n

    CHARACTER( LEN = 2 ) :: kindC
   
    !=========================================================================
   
    SELECT CASE ( kindC )

       !______________________________________________________________________

       CASE ('kp')

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n

          sa( i ) = sa( i ) * CMPLX( diag_kpV( i )**2 )
          
          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * CMPLX( diag_kpV( i ) ) * &
                  CMPLX( diag_kpV( ija( k ) ) )
          END DO

       END DO

       !______________________________________________________________________

       CASE ( 'sg' )

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n

          sa( i ) = sa( i ) * CMPLX( diag_sgV( i )**2 )

          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * CMPLX( diag_sgV( i ) ) * &
                  CMPLX( diag_sgV( ija( k ) ) )
          END DO

       END DO

       !______________________________________________________________________

    CASE DEFAULT
 
       !______________________________________________________________________

       STOP 'Error transform_sparse_to_ilu3D: kindC not allowed.'

       !______________________________________________________________________

    END SELECT

    !=========================================================================
 
  END SUBROUTINE transform_sparse_to_ilu3D
 
  !===========================================================================
 
  SUBROUTINE transform_sparse_to_ilu3D_real( kindC )
 
    !=========================================================================
    !
    ! NOTE: The matrix to be diagonalized is D^(1/2) H D^(1/2).
    !       In the setup, the matrix H is stored as a sparse matrix.
    !       In davidson routines, the routine amul directly calls
    !       the sparse matrix-vector product. So the diagonal 
    !       matrices D^(1/2) have to be taken into account. 
    !
    ! kindC = 'kp'  --> If it is a k.p matrix: D^(1/2) is stored in
    !                   diag_kpV.
    !       = 'sg'  --> If it is a single band Schroedinger matrix
    !                   with magnetic field: D^(1/2) is stored 
    !                   in diag_sgV.
    !
    ! Real version.
    ! 
    !=========================================================================
    
    USE normalization_factors_sg3D, ONLY : diag_kpV, diag_sgV
    USE sparse_matrix_real        , ONLY : sa, ija
   
    !=========================================================================

    INTEGER               :: i, k, n
    CHARACTER ( LEN = 2 ) :: kindC

    !=========================================================================
 
    SELECT CASE ( kindC )

       !______________________________________________________________________

       CASE ( 'kp' )

       !______________________________________________________________________

       
       STOP 'Error transform_sparse_to_ilu3D_real: &
            &Real version not useful for complex k.p matrix.'

       !______________________________________________________________________

       CASE ( 'sg' )

       !______________________________________________________________________

       n = ija( 1 ) - 2

       DO i = 1, n

          sa( i ) = sa( i ) * ( diag_sgV(i)**2 )
         
          DO k = ija( i ), ija( i + 1 ) - 1
             sa( k ) = sa( k ) * diag_sgV( i ) * diag_sgV( ija(k) )
          END DO
         
       END DO
      
       !______________________________________________________________________

    CASE DEFAULT
    
       !______________________________________________________________________
  
       STOP 'Error transform_sparse_to_ilu3D: kindC not allowed.'
      
       !______________________________________________________________________

    END SELECT
   
    !=========================================================================
   
  END SUBROUTINE transform_sparse_to_ilu3D_real

  !===========================================================================
  
END MODULE mod_precond_ilu_davidson

!=============================================================================
