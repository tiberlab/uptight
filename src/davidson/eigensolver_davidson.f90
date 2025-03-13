!=============================================================================

MODULE mod_eigensolver_davidson

  !===========================================================================
  
  IMPLICIT NONE

  !===========================================================================

CONTAINS

  !===========================================================================

  SUBROUTINE eigensolver_davidson( ndim, num_ev, eigenvalueV, eigenvectorM, &
       sigma, num_ev_converged, max_iterations, eps )

    !=========================================================================
    !
    ! Sparse hermitian matrix in module sparse_matrix_rout must be defined.
    !
    ! ndim             ... dimension of matrix
    ! num_ev           ... number of wanted eigenvalues
    ! sigma            ... the value near which the eigenvalues are sought
    ! num_ev_converged ... number of converged eigenvalues
    ! max_iterations   ... max. number of iterations
    ! eps              ... precision
    !
    ! eigenvalueV ( 1..num_ev )
    ! eigenvectorM( 1..ndim, 1..num_ev )
    !
    !=========================================================================
    
    ! input arguments :

    COMPLEX( KIND( 0.0d0 ) ), INTENT( IN )  :: sigma
    REAL   ( KIND( 0.0d0 ) ), INTENT( IN )  :: eps
    INTEGER,                  INTENT( IN )  :: ndim, num_ev, max_iterations
    
    !_________________________________________________________________________
    
    ! output arguments :
    
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( :, : ), POINTER :: eigenvectorM
    REAL   ( KIND( 0.0d0 ) ), DIMENSION( : ),    POINTER :: eigenvalueV   
 
    INTEGER, INTENT( OUT ) :: num_ev_converged
    
    !_________________________________________________________________________

    ! local variables :

    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( :, : ), ALLOCATABLE :: eivec, zwork 
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ),    ALLOCATABLE :: alpha, beta


    INTEGER :: kmax, jmax, jmin, method, m, l, mxmv, maxstep, n, &
         order, testspace, lwork, istatus, k

    COMPLEX( KIND( 0.0d0 ) ) :: target
    REAL   ( KIND( 0.0d0 ) ) :: lock

    INTEGER, DIMENSION( num_ev ) :: orderV
    
    LOGICAL :: wanted
    
    !=========================================================================
    
    n         = ndim
    target    = sigma
    kmax      = num_ev

    !jmax      = MAX( 3 * kmax, 20 ) ! changed by S. Birner
    !jmin      = 2 * kmax            ! following recommendations in manual

    jmax     = MAX( 2 * kmax, 10 ) ! This must be used for huge matrices
    jmin      = 1.5 * kmax         ! ( e.g. dots in 3D )

    method    = 2  ! 1 = GMRES
                   ! 2 = BICGSTAB ( Linear equation solver, Ax = b )

    maxstep   = max_iterations ! 1000 could be a reasonable value

    m         = 30
    l         = 2
    mxmv      = 100
    
    lock      = 1.0d-9
    testspace = 2
    wanted    = .TRUE.
    order     = 0
    lwork     = 10 + 6*l + 5*jmax + 3*kmax

    PRINT *, "jmax, Dimension of matrix, No. of eigenvalues, lwork, maxstep"
    PRINT *, jmax, n, kmax, lwork, maxstep

    ALLOCATE( alpha( jmax ), beta( jmax ), STAT = istatus )
    IF ( istatus .NE. 0 ) STOP &
         'Error eigensolver_davidson: Allocation failed for alpha or beta.'

    ALLOCATE( eivec( n, kmax ), STAT = istatus )
    IF ( istatus .NE. 0 ) STOP &
         'Error eigensolver_davidson: Allocation failed for eivec.'

    ALLOCATE( zwork( n, lwork ), STAT = istatus )
    IF ( istatus .NE. 0 ) THEN
       PRINT *, "dim   =", n
       PRINT *, "lwork =", lwork
       STOP 'Error eigensolver_davidson: Allocation failed for zwork. &
            &Suggestion: Reduce jmax and jmin.'
    END IF

    alpha = CMPLX( 0.0d0 )
    beta  = CMPLX( 0.0d0 )
    eivec = eigenvectorM

    CALL JDQZ( alpha, beta, eivec, wanted, n, target, eps, &
         kmax, jmax, jmin, method, m, l, mxmv, maxstep, &
         lock, order, testspace, zwork, lwork )
    
    PRINT *, 'num_ev_converged', kmax
    num_ev_converged = kmax
 
    eigenvalueV( 1 : num_ev ) = &
         DBLE( alpha( 1 : num_ev ) / beta( 1 : num_ev ) )

    !------------------------------------------
    
    CALL eigensort_complex( ndim, num_ev, eigenvalueV, orderV )
   
    WRITE ( *, * ) "sorted eigenvalues:"
    WRITE ( *, '( 5( f9.4, x ) )' ) eigenvalueV
    
    DO k = 1, num_ev
       eigenvectorM( 1 : ndim, k ) = eivec( 1 : ndim, orderV( k ) )
    END DO

    DEALLOCATE( alpha, beta, eivec, zwork )
    
    !=========================================================================
    
  CONTAINS
    
    !=========================================================================

    SUBROUTINE eigensort_complex( ndim, max_eigval, eigenvalue_cV, orderV )

      !=======================================================================
      !
      ! Input:
      ! ------
      !
      ! ndim         ...   size of matrix
      ! max_eigval   ...   number of eigenvalues
      !
      ! Input/Output:
      ! -------------
      !
      ! eigenvalue_cV ...  complex, DIMENSION(max_eigval)
      !                    Contains complex eigenvalues 
      !                    |imaginary part| < resid ,
      !                    otherwise error occurs.
      !
      ! eigenvectorM ...   complex (1:ndim,max_eigval)
      !                    corresponding eigenvalues
      !
      !_______________________________________________________________________
      !
      ! Sorts real eigenvectors in descending order.
      !
      !=======================================================================

      ! input arguments :

      INTEGER, INTENT( IN ) :: ndim, max_eigval

      !_______________________________________________________________________

      ! output arguments :

      INTEGER, DIMENSION( max_eigval ), INTENT( INOUT ) :: orderV

      REAL( KIND( 0.0d0 ) ), DIMENSION( : ), POINTER :: eigenvalue_cV
      
      !_______________________________________________________________________
  
      ! local variables :

      REAL( KIND( 0.0d0 ) ) :: eigenvalue_help
      
      INTEGER                 :: i, j, k, zw
      INTEGER, DIMENSION( 1 ) :: jV
      
      !=======================================================================

      orderV = (/ ( i, i = 1, max_eigval ) /)

      DO i = 1, max_eigval - 1

         jV = MAXLOC( DBLE( eigenvalue_cV( i : max_eigval ) ) ) + i - 1
         j = jV( 1 )

         IF ( j .NE. 1 ) THEN

            ! swap i <-> j
        
            eigenvalue_help    = eigenvalue_cV( i )
            eigenvalue_cV( i ) = eigenvalue_cV( j )
            eigenvalue_cV( j ) = eigenvalue_help
           
            zw          = orderV( i )
            orderV( i ) = orderV( j )
            orderV( j ) = zw

         END IF
      
      END DO

      !=======================================================================

    END SUBROUTINE eigensort_complex

    !=========================================================================

  END SUBROUTINE eigensolver_davidson


  !===========================================================================
  

  SUBROUTINE eigensolver_davidson_real(ndim,num_ev,eigenvalueV,eigenvectorM, &
       sigma,num_ev_converged,max_iterations,eps)

    !----------------------------------------------------------------------
    !
    ! Sparse hermitian matrix in MODULE sparse_matrix_rout must be defined.
    !
    ! ndim             ... dimension of matrix
    ! num_ev           ... number of wanted eigenvalues
    ! sigma            ... the value near which the eigenvalues are sought
    ! num_ev_converged ... number of converged eigenvalues
    ! max_iterations   ... max. number of iterations
    ! eps              ... precision
    !
    ! eigenvalueV (1..num_ev)
    ! eigenvectorM(1..ndim,1..num_ev)
    !
    !---------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER                            ,INTENT(in)  :: ndim,num_ev
    REAL   (KIND(0.0d0)),DIMENSION(:)  ,POINTER     :: eigenvalueV
    REAL   (KIND(0.0d0)),DIMENSION(:,:),POINTER     :: eigenvectorM
    REAL   (KIND(0.0d0))               ,INTENT(in)  :: sigma
    INTEGER                            ,INTENT(out) :: num_ev_converged
    INTEGER                            ,INTENT(in)  :: max_iterations
    REAL   (KIND(0.0d0))               ,INTENT(in)  :: eps

    !--------------------------------------------------------------------
    
    INTEGER :: kmax,jmax,jmin,method,m,l,mxmv,maxstep,n, &
         order,testspace,lwork,istatus,k
    COMPLEX(KIND(0.0d0)),DIMENSION(:)  ,ALLOCATABLE :: alpha,beta
    COMPLEX(KIND(0.0d0)),DIMENSION(:,:),ALLOCATABLE :: eivec 
    LOGICAL                                         :: wanted
    COMPLEX(KIND(0.0d0))                            :: target
    REAL   (KIND(0.0d0))                            :: lock
    COMPLEX(KIND(0.0d0)),DIMENSION(:,:),ALLOCATABLE :: zwork
    INTEGER             ,DIMENSION(num_ev)          :: orderV
    
    !---------------------------------------------------------------------
    
    n         = ndim
    target    = CMPLX(sigma)
    kmax      = num_ev
    jmax      = MAX(3*kmax,20) ! changed by S. Birner following recommendations in manual
    !jmax      = MAX(2*kmax,10)
    jmin      = 2*kmax         ! changed by S. Birner following recommendations in manual
    !jmin      = 1.5*kmax 
    method    = 2  ! 1=GMRES, 2=BICGSTAB  (Linear equation solver, Ax=b)
    m         = 30
    l         = 2
    mxmv      = 100
    maxstep   = max_iterations ! 1000 could be a reasonable value
    lock      = 1d-9
    testspace = 2
    wanted    = .TRUE.
    order     = 0
    lwork     = 10 + 6*l + 5*jmax + 3*kmax
    
    PRINT *,"jmax, Dimension of matrix, No. of eigenvalues, lwork"
    PRINT *, jmax,n,kmax,lwork
    
    ALLOCATE(alpha(jmax),beta(jmax),STAT=istatus)
    IF (istatus/=0) STOP &
         'Error eigensolver_davidson: Allocation failed for alpha or beta.'
    
    ALLOCATE(eivec(n,kmax),STAT=istatus)
    IF (istatus/=0) STOP &
         'Error eigensolver_davidson: Allocation failed for eivec.'
    
    ALLOCATE(zwork(n,lwork),STAT=istatus)
    IF (istatus/=0) STOP &
         'Error eigensolver_davidson: Allocation failed zwork.'
    
    
    alpha = CMPLX(0.0d0)
    beta  = CMPLX(0.0d0)
    eivec = CMPLX(eigenvectorM)
    
    ! Although this subroutine is called "real" it actually uses the complex 
    ! Jacobi-Davidson method.
    CALL JDQZ_real(alpha,beta,eivec,wanted,n,target,eps,kmax,jmax,jmin,method,&
         m,l,mxmv,maxstep,lock,order,testspace,zwork,lwork)
    
    PRINT *,'num_ev_converged',kmax
    num_ev_converged = kmax
    
    eigenvalueV(1:num_ev) = DBLE(alpha(1:num_ev)/beta(1:num_ev))
    
    
    !------------------------------------------
    
    CALL eigensort_real(ndim,num_ev,eigenvalueV,orderV)
    PRINT *,"sorted eigenvalues:",eigenvalueV
    
    DO k=1,num_ev
       eigenvectorM(1:ndim,k) = REAL( eivec(1:ndim,orderV(k)) )
    END DO
    
    DEALLOCATE(alpha,beta,eivec,zwork)
    
    !------------------------------------------
    
  CONTAINS
    
    !----------------------------------------------------------------------
    SUBROUTINE eigensort_real(ndim,max_eigval,eigenvalue_cV,orderV)
      !---------------------------------------------------------------------
      !
      ! Input:
      ! ------
      !
      ! ndim         ...   size of matrix
      ! max_eigval   ...   number of eigenvalues
      !
      ! Input/Output:
      ! -------------
      !
      ! eigenvalue_cV ...  complex, DIMENSION(max_eigval)
      !                    Contains complex eigenvalues 
      !                    |imaginary part| < resid ,
      !                    otherwise error occurs.
      !
      ! eigenvectorM ...   complex (1:ndim,max_eigval)
      !                    corresponding eigenvalues
      !
      !------------------------------------------------
      !
      ! Sorts real eigenvectors in descending order.
      !
      !--------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER                                ,INTENT(in)     :: ndim,max_eigval
      REAL(KIND(0.0d0)),DIMENSION(:),POINTER                 :: eigenvalue_cV
      INTEGER          ,DIMENSION(max_eigval),INTENT(in out) :: orderV
      !---------------------------------------------
      INTEGER                                                :: i,j,k,zw
      INTEGER          ,DIMENSION(1)                         :: jV
      REAL(KIND(0.0d0))                                      :: eigenvalue_help
      
      !------------------------------------------------------------------

      orderV = (/ (i, i=1,max_eigval) /)

      DO i=1,max_eigval-1

         jV=MAXLOC(DBLE(eigenvalue_cV(i:max_eigval)))+i-1
         j=jV(1)

         IF (j/=1) THEN
            ! swap i <-> j
            eigenvalue_help  = eigenvalue_cV(i)
            eigenvalue_cV(i) = eigenvalue_cV(j)
            eigenvalue_cV(j) = eigenvalue_help
            zw        = orderV(i)
            orderV(i) = orderV(j)
            orderV(j) = zw
         END IF
      
      END DO
      
      !---------------------------------------------------
    END SUBROUTINE eigensort_real
    !---------------------------------------------------
    
    !---------------------------------------------------------------------
  END SUBROUTINE eigensolver_davidson_real
  !---------------------------------------------------------------------

  
  !---------------------------------------------------------------------
END MODULE mod_eigensolver_davidson
!---------------------------------------------------------------------
