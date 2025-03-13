PROGRAM Driver

  !===========================================================================

  USE sparse_matrix      , ONLY : sa, ija 
  USE sparse_matrix_rout , ONLY : setup_sparse_matrix, &
       check_sparse_matrix, add_element_sparse_restr

  USE mod_precond_ilu_davidson
  USE mod_eigensolver_davidson   , ONLY : eigensolver_davidson

  !===========================================================================

  IMPLICIT NONE

  !===========================================================================

  COMPLEX( KIND( 0.0d0 ) ), DIMENSION( :, : ), POINTER :: eigenfunctions_pointM
  REAL   ( KIND( 0.0d0 ) ), DIMENSION( : )   , POINTER :: eigenval_pointV
  
  COMPLEX( KIND( 0.0d0 ) ) :: sigma, value
  
  REAL   ( KIND( 0.0d0 ) ) :: tolerance

  INTEGER( 4 ) :: N
  
  INTEGER :: num_ev, max_iterations, max_offdiag
  INTEGER :: istatus, i, j, row, column

  !===========================================================================

  !------------------------------------
  ! Parameters
  !------------------------------------

  N = 100                      ! matrix dimension
  max_offdiag = N * N / 2 + N  ! max_ofdiagonal
  num_ev = 5                   ! number of eigenvalues
  sigma  = ( 0.0, 0.0 )        ! value near which the eigenvalues are sought
  tolerance = 1e-5             ! MIN(SchroedingerkpResidual3D,min_precision_kp)
  max_iterations = 10000       ! SchroedingerkpIterations3D


  !------------------------------------
  ! Allocate sparse matrix.
  !------------------------------------

  CALL setup_sparse_matrix( N, max_offdiag )

  !------------------------------------
  ! Add elements
  !------------------------------------

  ! CALL add_element_sparse_restr( row, column, value )

  ! TEST

  DO i = 1, N

     DO j = i, N

        value = dcmplx( i ) + ( 0.D0, 1.D0 ) * dcmplx( j - i )
        CALL add_element_sparse_restr( i, j, value )

     END DO

  END DO


  !-----------------------------------
  ! Allocate eigenvalues eigen vectors
  !------------------------------------

  ALLOCATE( eigenval_pointV( num_ev ), &
       eigenfunctions_pointM( N, num_ev ), STAT = istatus )
  IF ( istatus .NE. 0 ) STOP &
       'Error calculate_kp_eigenvalues: Allocation failed.'

  !CALL transform_sparse_to_ilu3D('kp')

  CALL setup_david_precond_jac( sigma ) 

  CALL eigensolver_davidson( N, num_ev, eigenval_pointV, &
       eigenfunctions_pointM, sigma, N, max_iterations, tolerance )

  !===========================================================================

END PROGRAM Driver



