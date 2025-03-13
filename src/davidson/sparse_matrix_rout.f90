!=============================================================================
!
! NOTE: complex, hermitian version
!
! ( Hermitian matrix: 
!
!    Take the complex conjugate of A, then take the tranpose of it;
!    if this is equal to A, then it is called hermitian.
!    The diagonal elements have to be real.
!    The real parts of the offdiagonal elements must be symmetric,
!    its imaginary parts must be symmetric but of opposite sign.      )
!
! 1. subroutine setup_sparse_matrix
!
! 2. subroutine add_element_sparse
!
! 3. subroutine check_sparse_matrix
!
! 4. subroutine add_element_sparse_restr
!
! 5. subroutine matvec_sparse
!
! 6. subroutine matTvec_sparse
!
! 7. subroutine deallocate_sparse_matrix
!
! 8. subroutine num_rec_v_sp_m
!    numerical_recipes_version_of_sparse_matrix
!    won't be needed, just for testing purposes
!
! 9. subroutine num_rec_v_sp_m_upper
!
!    numerical_recipes_version_of_sparse_matrix_for_upper_part
!    won't be needed, just for testing purposes
!
! birner: this routine is slightly different from the Numerical Recipes version
!
!          ija( n + 1 ) ist not the value as indicated in the book
!
!  I WILL HAVE TO CHANGE IT!!!!!!!!!!!!!!!!
!
!=============================================================================

MODULE sparse_matrix_rout
  
  !===========================================================================

  IMPLICIT NONE
  
  !===========================================================================
  
CONTAINS
  
  !===========================================================================
  
  SUBROUTINE setup_sparse_matrix( N, max_offdiag )
    
    !=========================================================================
    !
    ! Arrays in module sparse_matrix must be deallocated.
    ! Allocates them with temporary dimension N + 1 + max_offdiag.
    !
    !=========================================================================
    
    USE sparse_matrix, ONLY : sa, ija
    
    !=========================================================================

    IMPLICIT NONE
    
    !=========================================================================
    
    INTEGER, INTENT( IN ) :: N, max_offdiag
    
    INTEGER :: istatus
    
    !=========================================================================
    
    IF ( ( N .LT. 1 ) .OR. ( max_offdiag .LT. 0 ) ) STOP &
         'Error setup_sparse_matrix: Wrong input.'
    
    !=========================================================================
    
    PRINT *, "Allocate matrix: dim =", N, "max_offdiag =", max_offdiag
    
    ALLOCATE( ija( N + 1 + max_offdiag ), STAT = istatus )
    IF ( istatus .NE. 0 ) STOP &
         'Error setup_sparse_matrix: Allocation for ija failed.'
    
    ALLOCATE( sa( N + 1 + max_offdiag ), STAT = istatus )
    IF ( istatus .NE. 0 ) STOP &
         'Error setup_sparse_matrix: Allocation for sa failed.'
    
    !=========================================================================
    
    sa  = CMPLX( 0.0d0 )
    
    ija = 0
    ija( 1 : N ) = N + 2
    ija( N + 1 ) = N + 2 ! set to initial value
    
    ! ija( N + 1 ) = size( ija ) + 1 would be correct IF and only if 
    ! the number of max_offdiag is identical to the number of off-diagonal
    ! elements.
    
    ! ija( N + 1 ) = size( ija ) + 1  ! IS THIS REALLY CORRECT??? S.Birner
    ! Why this? It probably is (see above)
    ! BUT it is different to the other routines
    ! AND More Importantly it doesn't work!!!
    
    !=========================================================================
    
  END SUBROUTINE setup_sparse_matrix
  
  
  !===========================================================================
  !
  !
  !
  !===========================================================================
  
  SUBROUTINE add_element_sparse( row, column, value )
    
    !=========================================================================
    ! adds matrix element to sparse matrix
    ! only for diagonal and upper part of matrix
    !=========================================================================
    
    USE sparse_matrix, ONLY : sa, ija
    
    !=========================================================================
    
    IMPLICIT NONE
    
    !=========================================================================
    
    INTEGER,                INTENT( IN ) :: row, column
    COMPLEX( KIND(0.0D0) ), INTENT( IN ) :: value
    
    INTEGER :: N, dim, targ, num_el, k, n_unshift

    LOGICAL, DIMENSION( : ), ALLOCATABLE :: unshift

    !=========================================================================
    
    IF ( ( row .LT. 1 ) .OR. ( row .GT. ija( 1 ) - 2 ) .OR. &
         ( column .LT. 1 ) .OR. ( column .GT. ija( 1 ) - 2 ) ) &
         STOP 'Error add_element_sparse (rout): Indices out of range.'

    !_________________________________________________________________________
    
    IF ( row .GT. column ) STOP &
         'Error add_element_sparse (rout): Only upper part of matrix.'

    !_________________________________________________________________________

    IF ( row .EQ. column ) THEN
    
       !------ diagonal element ----------
       
       ! checks if matrix is free of complex diagonal elements
       ! which is a requirement of Hermititan matrices

       IF ( AIMAG( value ) .EQ. 0 ) THEN
                 
          sa( row ) = value
          
       ELSE
          
          STOP 'Error add_element_sparse (rout): &
               & Matrix is not Hermitian because diagonal element is complex.'
          
       END IF

       !----------------------------------------------------------------------
       
    ELSE
       
       !------ offdiagonal element -------
       
       !-- shift all elements in greater rows to the right --
       
       N = ija( 1 ) - 2
       
       dim = SIZE( ija )
       targ = ija( row + 1 )
       num_el = ija( row + 1 ) - ija( row )
       
       IF ( targ .GT. SIZE( ija ) ) STOP &
            'Error add_element_sparse (rout): Matrix inconsistent (-1).'
       
       IF ( num_el .LT. 0 ) STOP &
            'Error add_element_sparse (rout): Matrix inconsistent (-2).'
       
       IF ( ija( dim ) .NE. 0 ) STOP &
            'Error add_element_sparse (rout): Matrix inconsistent (-3).'
       
       IF ( targ .LT. dim ) THEN
          
          ija( row + 1 : N ) = ija( row + 1 : N ) + 1
          ija( targ + 1 : dim ) = ija( targ : dim - 1 )
          sa( targ + 1 : dim ) = sa( targ : dim - 1 )
          
       END IF
       
       !---- write new element -----------------------------

       IF ( num_el .EQ. 0 ) THEN
          
          targ = ija( row )
          ija( targ ) = column
          sa( targ ) = value
          
       ELSE
          
          targ = ija( row ) 
          
          IF ( targ .GT. SIZE( ija ) ) STOP &
               'Error add_element_sparse (rout): Matrix inconsistent (-4).'
          
          IF ( column .GT. ija( targ + num_el - 1 ) ) THEN
             
             k = num_el
             ija( targ + k ) = column
             sa( targ + k ) = value 

          ELSE

             IF ( ANY( column .EQ. ija( targ : targ + num_el - 1 ) ) ) &
                  STOP 'Error add_element_sparse (rout):  &
                  &Overwrite existing matrix element'

             ALLOCATE( unshift( num_el ) )
             unshift = .TRUE.

             print*, unshift

             WHERE( column .LT. ija( targ : targ + num_el - 1 ) ) &
                  unshift = .FALSE.

             n_unshift = COUNT( unshift )

             ija( targ + n_unshift + 1 : targ + num_el ) = &
                  ija( targ + n_unshift : targ + num_el - 1 )
             
             sa( targ + n_unshift + 1 : targ + num_el ) = &
                  sa( targ + n_unshift : targ + num_el - 1 )
             
             ija( targ + n_unshift ) = column
             sa( targ + n_unshift ) = value
             
             DEALLOCATE( unshift )
             
          END IF
                          
       END IF
       
       ! This is to set the value of ija( N + 1 ).
       ! Location N+1 of ija is one greater than the index in sa of the last
       ! off-diagonal element of the last row.
       ! (It can be read to determine the number of non-zero elements
       ! in the matrix, or the logical length of the arrays sa and ija.)
       
       ! ( Start value was ija( N + 1 ) = N + 2 )
       ! Initial Number of off-diagonal elements before CALL was
       ! ija( N + 1 ) - ija( 1 )
       ! Now we need one more
       
       ija( N + 1 ) = ija( N + 1 ) + 1

    END IF
    
    !=========================================================================
    
  END SUBROUTINE add_element_sparse
  
  
  !===========================================================================
  !
  !
  !===========================================================================

  SUBROUTINE check_sparse_matrix
    
    !=========================================================================
    !
    ! checks if sparse matrix storage scheme is correct, 
    ! only useful if max_offdiag = num_offdiag
    !
    !=========================================================================

    USE sparse_matrix , ONLY: sa, ija 

    !=========================================================================
    
    IMPLICIT NONE
    
    !=========================================================================
    
    INTEGER :: N, k
    
    !=========================================================================
    
    N = ija( 1 ) - 2  ! N = dimension of N x N matrix
    
    ! checks if zero elements were recorded on the right part
    ! of the sparse matrix storage scheme
    ! ( if max_offdiag = num_off_diag this test is useful,
    ! otherwise it always fails )        
    ! This check leads to an error if max_offdiag is an upper limit
    
    DO k = N + 2, SIZE( ija )
       
       IF ( ija( k ) .EQ. 0 ) STOP &
            'Error check_sparse_matrix: Matrix inconsistent (1).'
       
    END DO
    
    ! This check doesn't work if max_offdiag is an upper limit
    ! rather than the number of actual offdiagonal elements 
 
    DO k = 1, N
       
       IF ( ija( k ) .GT. SIZE( ija ) + 1 ) STOP &
            'Error check_sparse_matrix: Matrix inconsistent (2).'
       
       ! StefanB Why not ija( k ) > size( ija ) ????
       ! I think that is more suitable.
       
    END DO
    
    ! This check doesn't work if max_offdiag is an upper limit
    ! rather than the number of actual offdiagonal elements
    ! This value must be equal
    
    IF ( ija( N + 1 ) .NE. SIZE( ija ) + 1 ) &
         STOP 'Error check_sparse_matrix: Matrix inconsistent (3).'
    
    !=========================================================================
    
  END SUBROUTINE check_sparse_matrix
  
  
  !===========================================================================
  
  SUBROUTINE add_element_sparse_restr( row, column, value )
    
    !=========================================================================
    !
    ! Adds matrix element to sparse matrix.
    ! Only for diagonal and upper part of matrix.
    !
    !=========================================================================
    
    USE sparse_matrix, ONLY : sa, ija
    
    !=========================================================================
    
    IMPLICIT NONE
    
    !=========================================================================
    
    ! input arguments :
    
    COMPLEX ( KIND( 0.0d0 ) ), INTENT( IN ) :: value
    INTEGER,                   INTENT( IN ) :: row, column
    
    !_________________________________________________________________________
    
    ! local variables :
    
    INTEGER :: N, dim,targ, num_el, k
    LOGICAL :: writeL
    
    !=========================================================================

    IF ( ( row .LT. 1 ) .OR. ( row .GT. ija( 1 ) - 2 ) .OR. &
         ( column .LT. 1 ) .OR. ( column .GT. ija( 1 ) - 2 ) ) &
         STOP 'error add_element_sparse_restr (rout): indices out of range'
    
    IF ( row .GT. column ) STOP &
         'error add_element_sparse_restr (rout): only upper part of matrix'
    
    !=========================================================================
    
    IF ( row .EQ. column ) THEN
       
       !------ diagonal element ----------
       
       IF ( AIMAG( value ) .EQ. 0 ) THEN
          
          ! checks if matrix is free of complex diagonal elements
          ! which is a requirement of Hermititan matrices
          
          sa( row ) = value
          
       ELSE
          
          STOP 'Error add_element_sparse_restr (rout): &
               & Matrix is not Hermitian because diagonal element is complex.'
          
       END IF
       
    ELSE
       
       !------ offdiagonal element -------
       
       !-- shift all elements in greater rows to the right --
       
       N = ija( 1 ) - 2
       
       dim = SIZE( ija )
       targ = ija( row + 1 )
       num_el = ija( row + 1 ) - ija( row )
       
       IF ( targ .GT. SIZE( ija ) ) STOP &
            'error add_element_sparse_restr (rout): matrix inconsistent -1'
       
       IF ( num_el .LT. 0 ) STOP &
            'error add_element_sparse_restr (rout): matrix inconsistent -2'
       
       IF ( ija( dim ) .NE. 0 ) STOP &
            'error add_element_sparse_restr (rout): matrix inconsistent -3'
       
       ! -- note: in this special routine NO shift occurs --
      
       IF ( row .LT. N ) THEN
          
          ija( row + 1 ) = ija( row + 1 ) + 1
          
          !ija( targ + 1 ) = ija( targ )
          
       ELSE
          
          ija( N + 1 ) = ija( N + 1 ) + 1
          
       END IF
       
       IF ( ija( N + 1 ) .GT. SIZE( ija ) + 1 ) STOP '--d---'
       
       !---- write new element -----------------------------
       
       IF ( num_el .EQ. 0 ) THEN
          
          targ = ija( row )
          ija( targ ) = column
          sa( targ ) = value
          
       ELSE
          
          targ = ija( row ) 
          writeL = .FALSE.
          
          IF ( targ .GT. SIZE( ija ) ) STOP &
               'error add_element_sparse_restr (rout):matrix inconsistent-4'
          
          DO k = 0, num_el - 1
             
             IF ( column .EQ. ija( targ + k ) ) &
                  STOP 'error add_element_sparse_restr (rout): &
                  & overwrite existing matrix element!'
             
             IF ( column .LT. ija( targ + k ) ) THEN
                
                ija( targ + k + 1 : targ + num_el ) = &
                     ija( targ + k : targ + num_el - 1 )
                
                sa( targ + k + 1 : targ + num_el ) = &
                     sa( targ + k : targ + num_el - 1 )
                
                ija( targ + k ) = column
                sa( targ + k ) = value
                writeL = .TRUE.       
                
                EXIT

             END IF
             
          END DO
          
          IF ( .NOT.( writeL ) ) THEN
             
             k = num_el
             ija( targ + k ) = column
             sa( targ + k ) = value 
             
          END IF
          
       END IF
       
    END IF
    
    IF ( row .LT. ija( 1 ) - 2 ) ija( row + 2 ) = ija( row + 1 )
    
    !=========================================================================
    
  END SUBROUTINE add_element_sparse_restr
  
  
  !===========================================================================
  !
  !
  !
  !===========================================================================
  
  SUBROUTINE matvec_sparse( n, xV, bV )
    
    !=========================================================================
    !
    ! calculates   bV = Matrix * xV
    !
    !=========================================================================

    USE sparse_matrix , ONLY: sa, ija
    
    !========================================================================= 
    
    IMPLICIT NONE
    
    !=========================================================================
    
    INTEGER, INTENT( IN )  :: n
    
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ), INTENT( IN )  :: xV
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ), INTENT( OUT ) :: bV
    
    INTEGER :: i, k
    
    !=========================================================================
    
    IF ( ija( 1 ) - 2 .NE. n ) STOP &
         'error matvec_sparse: vector has wrong DIMENSION'

    bV = CMPLX( 0.0d0 )
    
    DO i = 1, n

       bV( i ) = bV( i ) + sa( i ) * xV( i )

       DO k = ija( i ), ija( i + 1 ) - 1
          
          bV( i ) = bV( i ) + sa( k ) * xV( ija( k ) )
          bV( ija( k ) ) = bV( ija( k ) ) + CONJG( sa( k ) ) * xV( i )
         
       END DO

    END DO
    
    !=========================================================================

  END SUBROUTINE matvec_sparse

  
  !===========================================================================
  !
  !
  !
  !
  !===========================================================================

  SUBROUTINE matTvec_sparse( n, xV, bV )

    !=========================================================================
    ! calculates   bV = Matrix ^ T* xV
    !=========================================================================

    USE sparse_matrix, ONLY : sa, ija

    !=========================================================================
    
    IMPLICIT NONE

    !=========================================================================
 
    INTEGER, INTENT( IN ) :: n

    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ), INTENT( IN )  :: xV
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ), INTENT( OUT ) :: bV

    INTEGER :: i, j, k

    !=========================================================================

    IF ( ija( 1 ) - 2 .NE. n ) STOP &
         'error matTvec_sparse: vector has wrong DIMENSION'

    bV = CMPLX( 0.0d0 )

    DO i = 1, n

       bV( i ) = bV( i ) + sa( i ) * xV( i ) !start with diagonal terms
       
    END DO
 
    DO i = 1, n
      
       DO k = ija( i ), ija( i + 1 ) - 1 !loop over off-diagonal terms
      
          j = ija( k ) 
          bV( j ) = bV( j ) + sa( k ) * xV( i )
          bV( i ) = bV( i ) + CONJG( sa( k ) ) * xV( j )

       END DO
       
    END DO
    
    !=========================================================================
    
  END SUBROUTINE matTvec_sparse

  
  !===========================================================================
  !
  !
  !===========================================================================

  SUBROUTINE num_rec_v_sp_m( a, n, np, thresh, nmax )

    !=========================================================================
    !
    ! written by Stefan Birner to test the sparse matrix storage system
    ! converts a square matrix a( 1 : n, 1 : n ) with physical dimension np
    ! into row-indexed sparse storage mode.
    ! only elements of a with magnitude >= thresh are retained.
    ! Output is in two linear arrays with physical dimension nmax
    ! ( an input parameter ) :
    ! sa( 1: ) contains array values, indexed by ija( 1: ).
    ! The logical size of sa and ija on output are both
    ! ija( ija( 1 ) - 1 ) - 1 ( see Numerial Recipes text )
    !
    !=========================================================================

    USE sparse_matrix, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    ! Numerical Recipes version

    INTEGER, INTENT( IN ) :: n, nmax, np

    REAL( KIND( 0.0d0 ) ), INTENT( IN )  :: thresh
 
    !INTEGER :: ija(nmax)

    COMPLEX( KIND( 0.0d0 ) ) :: a( np, np )

    !COMPLEX(KIND(0.0d0)) :: sa(nmax)

    INTEGER :: i, j, k

    !=========================================================================

    DO j = 1, n                   !store diagonal elements
       sa( j ) = a( j, j ) 
    END DO

    ija( 1 ) = n + 2       !index to 1st row off-diagonal element, if any
    k = n + 1
    DO i = 1, n                   !loop over rows
       DO j = 1, n                !loop over columns
          
          IF ( ABS( a( i, j ) ) .GE. thresh ) THEN

             IF ( i .NE. j ) THEN 

                ! store off-diagonal elements and their columns
               
                k = k + 1
      
                IF ( k .GT. nmax ) STOP &
                     'nmax too small in birner_version_of_sparse_matrix'
               
                sa( k ) = a( i, j )
                ija( k ) = j         
               
             END IF
         
          END IF
      
       END DO
      
       ija( i + 1 ) = k + 1   !As each row is completed, store index to next.
   
    END DO

    RETURN

    !=========================================================================
 
  END SUBROUTINE num_rec_v_sp_m

  !===========================================================================
  !
  !
  !===========================================================================
 
  SUBROUTINE num_rec_v_sp_m_upper( a, n, np, thresh, nmax )

    !=========================================================================
    !
    ! written by Stefan Birner to test the sparse matrix storage system
    ! converts a square matrix a( 1 : n, 1 : n ) with physical dimension np
    ! into row-indexed sparse storage mode.
    ! Only elements of a with magnitude >= thresh are retained.
    ! Output is in two linear arrays with physical dimension nmax
    ! ( an input parameter ):
    ! sa( 1: ) contains array values, indexed by ija( 1: ).
    ! The logical size of sa and ija on output are both
    ! ija( ija( 1 ) - 1 ) - 1 ( see Numerial Recipes text )
    !
    !=========================================================================
   
    USE sparse_matrix, ONLY : sa, ija
   
    !=========================================================================

    IMPLICIT NONE
   
    !=========================================================================

    ! Numerical Recipes version
   
    INTEGER, INTENT( IN ) :: n, nmax, np

    REAL( KIND( 0.0d0 ) ), INTENT( IN ) :: thresh

    !INTEGER :: ija(nmax)

    COMPLEX( KIND( 0.0d0 ) ) :: a( np, np )

    !COMPLEX(KIND(0.0d0)) :: sa(nmax)

    INTEGER :: i, j, k
   
    !=========================================================================

    DO j = 1, n                   !store diagonal elements
       sa( j ) = a( j, j ) 
    END DO

    ija( 1 ) = n + 2       !index to 1st row off-diagonal element, if any
    k = n + 1

    DO i = 1, n                   !loop over rows
     
       !ONLY upper part of the matrix version
       !do j=1,n                !loop over columns

       DO j = i, n                !loop over columns
        
          IF ( ABS( a( i, j ) ) .GE. thresh ) THEN

             IF ( i .NE. j ) THEN

                ! store off-diagonal elements and their columns
               
                k = k + 1
               
                IF ( k .GT. nmax ) STOP &
                     'nmax too small in birner_version_of_sparse_matrix'
               
                sa( k ) = a( i, j )
                ija( k ) = j         
           
             END IF

          END IF

       END DO

       ija( i + 1 ) = k + 1    !As each row is completed, store index to next.

    END DO

    RETURN

    !=========================================================================

  END SUBROUTINE num_rec_v_sp_m_upper

  
  !===========================================================================
  !
  !
  !===========================================================================

  SUBROUTINE deallocate_sparse_matrix

    !=========================================================================

    USE sparse_matrix, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    DEALLOCATE( sa, ija )
   
    !=========================================================================

  END SUBROUTINE deallocate_sparse_matrix

  !===========================================================================

END MODULE sparse_matrix_rout


!============================================================================
!
!
!
!============================================================================

MODULE SPARSKIT_COMPLEX

  USE sparse_matrix, ONLY : nrow, sparse_fmt, sa, ija    

  IMPLICIT NONE
  
  !===========================================================================
  
CONTAINS

  SUBROUTINE msr2msr_full()

    ! OPERATES ON sparse_matrix

    ! INTERNAL ARRAYS:

    INTEGER :: nzmax, i, ierr
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ), ALLOCATABLE, TARGET ::  sa2
    INTEGER, DIMENSION( : ), ALLOCATABLE, TARGET        :: ija2
    
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( : ), ALLOCATABLE :: sa3, sa4
    INTEGER, DIMENSION( : ), ALLOCATABLE :: iwk, ia2, ija3, ia3, ia4, ija4

    ! ALLOCATE NEW ARRAYS 

    if(sparse_fmt.EQ.'F' .OR. sparse_fmt.EQ.'f') return

    nzmax = 2*SIZE(sa) - nrow - 1

    WRITE(*,*) SIZE(sa), nzmax

    ALLOCATE( sa2( nzmax ), sa3( nzmax ) )
    ALLOCATE( ija2( nzmax ), ija3( nzmax ) )
    ALLOCATE( ia2(nrow+1), ia3(nrow+1) )

    ALLOCATE( iwk(nrow+1) )

    IF ( nrow .NE. ija( 1 ) - 2 ) STOP "Error msr2csr_full: Wrong dimension."

    WRITE(*,*) 'convert msr to csr'
    ! CREATE an intermediate csr matrix (actually hsr)
    CALL zmsrcsr( nrow, sa, ija,    sa2, ija2, ia2,     sa3, ia3 )
 
    DEALLOCATE(sa, ija)

    ALLOCATE( sa4( nzmax ), ija4( nzmax), ia4(nrow+1) )

    ! MAKE hermitian hsr into csr in place
    !SELECT CASE(sparse_fmt)
    !CASE('U','u')
    !
    !   WRITE(*,*) 'connot yet convert uhsr to csr'
    !   ! buggy routine ...
    !   !CALL zuhsrcsr( nrow, sa2, ija2, ia2, nzmax, sa2, ija2, ia2, iwk, ierr )
    !   deallocate(sa2,ija2,wk, iwk, ia2)
    !   return!
    !
    !CASE('L','l')
    !
    !   WRITE(*,*) 'convert lhsr to csr'
    !   CALL zlhsrcsr( nrow, sa2, ija2, ia2, nzmax, sa2, ija2, ia2, iwk, ierr )
    !
    !END SELECT
    !    zcsrcsc2(n,n,job,ipos,a,ja,ia,ao,jao,iao)

    CALL zcsrcsc2(nrow, nrow, 1, 1, sa2, ija2, ia2, sa3, ija3, ia3)
    do i=1,nzmax
       sa3(i)=conjg(sa3(i))
    enddo

    CALL zaplb(nrow, nrow, 1, sa2, ija2, ia2, sa3, ija3, ia3, &
                                   sa4, ija4, ia4, nzmax, iwk, ierr)

    if (ierr.ne.0) STOP 'ERROR in msr2msr_full'

    !CONVERT back to msr

    WRITE(*,*) 'convert csr to msr'
    CALL zcsrmsr( nrow, sa4, ija4, ia4,    sa2, ija2,   sa3, iwk )

    ! MAKE the sparse matrix to point to new sparse matrix
    
    sa => sa2
    ija => ija2

    WRITE(*,*) 'sa, ija redirected'

    sparse_fmt='F'

    DEALLOCATE( iwk, ia2, ija3, ija4, ia3, ia4, sa3, sa4 )


  END SUBROUTINE msr2msr_full
  !===========================================================================
  
  SUBROUTINE print_sparse_matrix_complex( NDIM, filename )
    
    INTEGER, INTENT( IN ) :: NDIM
    
    CHARACTER ( LEN = * ), INTENT( IN ) :: filename

    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( NDIM )            :: tmpcV, bcV
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( SIZE( ija ) - 1 ) ::  a_csr, ao_csr
    COMPLEX( KIND( 0.0d0 ) ), DIMENSION( 2 * SIZE( ija ) - 1 - NDIM ) :: ACSR
    
    INTEGER, DIMENSION( SIZE( ija ) - 1 ) :: ja_csr, jao_csr
    INTEGER, DIMENSION( NDIM + 1 )        :: iwk_csr, ia_csr, iao_csr, IACSR
    INTEGER, DIMENSION( 2 * SIZE( ija ) - 1 - NDIM ) :: JACSR
    
    INTEGER :: open_error
    
    ! For printing matrices
    
    CHARACTER( LEN = 72 ) :: title  = "TITLE"
    CHARACTER( LEN = 8 )  :: key    = "KEY"
    CHARACTER( LEN = 3 )  :: type1  = "TYP"
    CHARACTER( LEN = 2 )  :: guesol = "ab"
    
    ! GZ are proper options, ab is nothing
    
    IF ( NDIM .NE. ija( 1 ) - 2 ) STOP &
         "Error print_sparse_matrix_complex: Wrong dimension."
    
    PRINT *, sa
    
    !  print *,NDIM
    !  print *,size(tmpcV),size(bcV)
    !  print *,size(a_csr),size(ao_csr),size(ja_csr),size(jao_csr)
    !  print *,size(acsr),size(jacsr)
    
    ! Create CSR storage format from MSR storage format
    
    CALL zmsrcsr( NDIM, a_csr, ja_csr, ia_csr, tmpcV, iwk_csr )
    
    ! okay, now we only have the upper triangular part.

    ! We also need the lower complex conjugated part
    ! ( here the diagonal is zero )
    
    CALL csrcsc_complex( NDIM, a_csr, ja_csr, ia_csr, &
         ao_csr, jao_csr, iao_csr )

    ! Now we add the lower part to the upper part
    
    CALL aplb1_complex( NDIM, a_csr, ja_csr, ia_csr, &
         ao_csr, jao_csr, iao_csr, ACSR, JACSR, IACSR, 2 * SIZE( ija ) - 1 )
    
    ! We now print a H/B file of the matrix
    
    OPEN ( UNIT = 8, FILE = filename, STATUS = "REPLACE", &
         ACTION = "WRITE", IOSTAT = open_error )
    
    IF ( open_error .NE. 0 ) STOP "error during file writing"
   
    CALL prtmt( NDIM, NDIM, ACSR, JACSR, IACSR, &
         bcV, guesol, title, key, type1, 4, 2, 8 )
    
    !=========================================================================
    
  END SUBROUTINE print_sparse_matrix_complex
  
  !===========================================================================
  
END MODULE SPARSKIT_COMPLEX


