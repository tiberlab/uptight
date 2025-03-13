!=============================================================================
!
! note: real, symmetric version
!
! 1. SUBROUTINE setup_sparse_matrix_real
! 2. SUBROUTINE add_element_sparse_real
! 3. SUBROUTINE check_sparse_matrix_real
! 4. SUBROUTINE matvec_sparse_real
! 5. SUBROUTINE matvec_sparse_real_pointer
! 6. SUBROUTINE deallocate_sparse_matrix_real
! 7. SUBROUTINE add_element_sparse_real_restr
! 8. SUBROUTINE write2file_sparse_matrix_real
! 9. SUBROUTINE num_rec_v_sp_m_upper_real
!    numerical_recipes_version_of_sparse_matrix_for_upper_part_real
!    for testing purposes
!
!=============================================================================

MODULE sparse_matrix_rout_real

  !===========================================================================

  IMPLICIT NONE

  !===========================================================================
  
CONTAINS

  !===========================================================================
  
  SUBROUTINE setup_sparse_matrix_real( N, max_offdiag )
  
    !=========================================================================
    !
    ! arrays in MODULE sparse_matrix must be deallocated
    ! allocates them with temporary dimension N+1+max_offdiag
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija
    
    !=========================================================================

    IMPLICIT NONE

    !=========================================================================
    
    INTEGER, INTENT( IN ) :: N, max_offdiag
    INTEGER               :: istatus

    !=========================================================================
    
    IF ( ( N < 1 ) .OR. ( max_offdiag < 0 ) ) STOP &
         'error setup_sparse_matrix: wrong input'
    
    ALLOCATE( sa( N+1+max_offdiag ), ija( N+1+max_offdiag ), STAT = istatus )
    IF ( istatus .NE. 0 ) STOP &
         'Error setup_sparse_matrix_real: Allocation failed.'
    
    sa  = 0.0D0
    ija = 0
    
    ija( 1 : N+1 ) = N + 2
    
    !=========================================================================

  END SUBROUTINE setup_sparse_matrix_real
  

  !===========================================================================


  SUBROUTINE add_element_sparse_real( row, column, value )

    !=========================================================================
    !
    ! adds matrix element to sparse matrix
    ! only for diagonal and upper part of matrix
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    REAL( KIND( 0.0d0 ) ), INTENT( IN ) :: value
    INTEGER,               INTENT( IN ) :: row, column

    INTEGER :: N, dim, targ, num_el, k
    LOGICAL :: writeL

    !=========================================================================

    IF ( ( row .LT. 1 ) .OR. ( row .GT. ija(1) - 2 ) .OR. &
         ( column .LT. 1 ) .OR. ( column .GT. ija(1) - 2 ) ) &
         STOP 'error add_element_sparse: indices out of range'
    
    IF ( row .GT. column ) STOP &
         'Error add_element_sparse_real: Only upper part of matrix.'
    
    
    !------ diagonal element ----------
    
    IF ( row .EQ. column ) THEN
       
       sa(row) = value
       
    ELSE
       
       !------ offdiagonal element -------
       
       !-- shift all elements in greater rows to the right --
       
       N = ija( 1 ) - 2
       dim = SIZE( ija )
       targ = ija( row + 1 )
       num_el = ija( row + 1 ) - ija( row )

       IF ( targ .GT. SIZE( ija ) ) STOP &
            'Error add_element_sparse_real: Matrix inconsistent (-1).'
       IF ( num_el .LT. 0 ) STOP &
            'Error add_element_sparse_real: Matrix inconsistent (-2).'
       IF ( ija(dim) .NE. 0 ) STOP &
            'Error add_element_sparse_real: Matrix inconsistent (-3).'

       IF ( targ .LT. dim ) THEN

          ija( row + 1 : N )    = ija( row + 1 : N ) + 1
          ija( targ + 1 : dim ) = ija( targ : dim - 1 )
          sa( targ + 1 : dim )  = sa( targ : dim - 1 )

       END IF

       !---- write new element -----------------------------

       IF ( num_el .EQ. 0 ) THEN

          targ = ija( row )
          ija( targ ) = column
          sa( targ ) = value

       ELSE

          targ = ija( row )
          writeL = .FALSE.
          IF ( targ .GT. SIZE( ija ) ) STOP &
               'Error add_element_sparse_real:matrix inconsistent-4'

          DO k = 0, num_el - 1

             IF ( column .EQ. ija( targ + k ) ) STOP &
                  'Error add_element_sparse_real: &
                  & Overwrite existing matrix element!'

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

          IF ( .NOT. writeL ) THEN

             k = num_el
             ija( targ + k ) = column
             sa( targ + k ) = value 

          END IF

       END IF

       ! This is to set the value of ija(N+1).
       ! Location N+1 of ija is one greater than the index in sa of the last
       ! off-diagonal element of the last row.
       ! (It can be read to determine the number of non-zero elements 
       ! in the matrix, or the logical length of the arrays sa and ija.)

       ! (Start value was ija(N+1)=N+2)
       ! Initial Number of off-diagonal elements 
       !before CALL was ija(N+1)-ija(1)
       ! Now we need one more
       
       ija( N + 1 ) = ija( N + 1 ) + 1

    END IF

    !=========================================================================

  END SUBROUTINE add_element_sparse_real
 

  !===========================================================================
  

  SUBROUTINE check_sparse_matrix_real

    !=========================================================================
    !
    ! Checks if sparse matrix storage scheme is correct.
    ! Only useful if max_offdiag=num_offdiag.
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija
    
    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    INTEGER :: N, k

    !=========================================================================

    N = ija( 1 ) - 2 ! N = dimension of N x N matrix

    ! Checks if zero elements were recorded on the right part
    ! of the sparse matrix storage scheme.
    ! (if max_offdiag=num_off_diag this test is useful,
    ! otherwise it always fails)
    ! This check leads to an error if max_offdiag is an upper limit.

    DO k = N + 2, SIZE( ija )

       IF ( ija( k ) .EQ. 0 ) STOP &
            'Error check_sparse_matrix_real: Matrix inconsistent (1).'

    END DO

    ! This check doesn't work if max_offdiag is an upper limit
    ! rather than the number of actual offdiagonal elements 
    
    DO k = 1, N

       IF ( ija(k) .GT. SIZE(ija) + 1 ) STOP &
            'Error check_sparse_matrix_real: Matrix inconsistent (2).'
       
       ! StefanB Why not ija(k)>size(ija) ???? I think that is more suitable.

    END DO


    ! This check doesn't work if max_offdiag is an upper limit
    ! rather than the number of actual offdiagonal elements
    ! This value must be equal
    
    PRINT *, ija( N + 1 ), SIZE(ija) + 1
    IF ( ija( N + 1 ) .NE. SIZE(ija) + 1 ) STOP &
         'Error check_sparse_matrix_real: Matrix inconsistent (3).'

    !=========================================================================

  END SUBROUTINE check_sparse_matrix_real
  
  !===========================================================================

  SUBROUTINE matvec_sparse_real( n, xV, bV )

    !=========================================================================
    !
    ! calculates   bV = Matrix*xV
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    INTEGER                           , INTENT( IN )  :: n
    REAL ( KIND(0.0d0) ), DIMENSION(:), INTENT( IN )  :: xV
    REAL ( KIND(0.0d0) ), DIMENSION(:), INTENT( OUT ) :: bV

    INTEGER :: i, k

    !=========================================================================

    IF ( ija(1) - 2 .NE. n ) STOP &
         'Error matvec_sparse_real: Vector has wrong dimension.'

    bV = 0.0D0

    DO i = 1, n

       bV(i) = bV(i) + sa(i) * xV(i)

       DO k = ija(i), ija(i+1) - 1

          bV(i)      = bV(i)      + sa(k) * xV(ija(k))
          bV(ija(k)) = bV(ija(k)) + sa(k) * xV(i)

       END DO

    END DO

    !=========================================================================

  END SUBROUTINE matvec_sparse_real

  !===========================================================================
  
  SUBROUTINE matvec_sparse_real_pointer( n, xV, bV )
    
    !=========================================================================
    !
    ! calculates   bV = Matrix*xV
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    INTEGER                           , INTENT( IN ) :: n
    REAL ( KIND(0.0d0) ), DIMENSION(:), POINTER      :: xV
    REAL ( KIND(0.0d0) ), DIMENSION(:), POINTER      :: bV

    INTEGER :: i, k

    !=========================================================================

    IF ( ija(1) - 2 .NE. n ) STOP &
         'Error matvec_sparse_real_pointer: Vector has wrong dimension.'

    bV = 0.0D0

    DO i = 1, n

       bV(i) = bV(i) + sa(i) * xV(i)

       DO k = ija(i), ija(i+1) - 1

          bV(i) = bV(i) + sa(k) * xV(ija(k))
          bV(ija(k)) = bV(ija(k)) + sa(k) * xV(i)

       END DO

    END DO

    !=========================================================================

  END SUBROUTINE matvec_sparse_real_pointer
  
  !===========================================================================

  SUBROUTINE deallocate_sparse_matrix_real
  
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    DEALLOCATE( sa, ija )

    !=========================================================================

  END SUBROUTINE deallocate_sparse_matrix_real

  !===========================================================================

  SUBROUTINE num_rec_v_sp_m_upper_real( a, n, np, thresh, nmax )
  
    !=========================================================================
    !
    ! written by Stefan Birner to test the sparse matrix storage system
    ! converts a square matrix a(1:n,1:n) with physical dimension np
    ! into row-indexed sparse storage mode.
    ! Only elements of a with magnitude >=thresh are retained.
    ! Output is in two linear arrays with physical dimension nmax
    ! (an input parameter):
    ! sa(1:) contains array values, indexed by ija(1:).
    !The logical size of sa and ija on output are both ija(ija(1)-1)-1
    ! (see Numerial Recipes text)
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    !Numerical Recipes version
    
    INTEGER,INTENT( IN ) :: n, nmax, np
    
    !INTEGER :: ija(nmax)
    
    REAL( KIND(0.0d0) ) :: a( np, np )
    
    !complex(kind(0.0d0)) :: sa(nmax)
    
    REAL ( KIND(0.0d0) ), INTENT(IN ) :: thresh
    
    INTEGER :: i, j, k

    !=========================================================================

    DO j = 1, n                   !store diagonal elements
       sa(j) = a(j,j) 
    END DO

    ija(1) = n + 2    !index to 1st row off-diagonal element, if any
    k = n + 1
    
    DO i = 1, n                   !loop over rows
       
       !only upper part of the matrix version
       !do j=1,n                !loop over columns
       
       DO j = i, n                !loop over columns
          
          IF ( ABS( a(i,j) ) .GE. thresh ) THEN

             IF ( i .NE. j ) THEN
                ! store off-diagonal elements and their columns
                k = k + 1
                IF (k .GT. nmax ) STOP &
                     'nmax too small in birner_version_of_sparse_matrix'
                sa(k) = a(i,j)
                ija(k) = j         
             END IF

          END IF

       END DO

       ija(i+1) = k + 1       
       !As each row is completed, store index to next.

    END DO
    RETURN

    !=========================================================================

  END SUBROUTINE num_rec_v_sp_m_upper_real

  !===========================================================================

  SUBROUTINE add_element_sparse_real_restr(row,column,value)

    !=========================================================================
    !
    ! Adds matrix element to sparse matrix.
    ! Only for diagonal and upper part of matrix.
    ! The same as add_element_sparse_real with the only
    ! restriction that if elements of row n are written
    ! all rows>n are not defined, i.e. these values are not shifted
    ! to the right.
    !
    ! Note: Here the user must ensure that:
    !       - ALL rows are looped trough from 1 to N.
    !       - Always check with add_element_sparse_real.
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    INTEGER            , INTENT( IN ) :: row, column
    REAL( KIND(0.0d0) ), INTENT( IN ) :: value

    INTEGER :: N, dim, targ, num_el, k
    LOGICAL :: writeL

    !=========================================================================

    IF ( ( row .LT. 1 ) .OR. ( row .GT. ija(1)-2 ) .OR. &
         ( column .LT. 1 ) .OR. ( column .GT. ija(1)-2 ) ) &
         STOP 'Error add_element_sparse_real_restr: Indices out of range.'

    IF ( row .GT. column ) STOP &
         'Error add_element_sparse_real_restr: Only upper part of matrix.'


    !------ diagonal element ----------

    IF ( row .EQ. column ) THEN

       sa(row) = value

    ELSE

       !------ offdiagonal element -------
       
       !-- shift all elements in greater rows to the right --

       N = ija(1) - 2
       dim = SIZE(ija)
       targ = ija( row + 1 )
       num_el = ija( row + 1 ) - ija( row )

       IF ( targ .GT. SIZE(ija) ) STOP &
            'Error add_element_sparse_real_restr: Matrix inconsistent (-1).'
       IF ( num_el .LT. 0 )  STOP &
            'Error add_element_sparse_real_restr: Matrix inconsistent (-2).'
       IF ( ija(dim) .NE. 0 ) STOP &
            'Error add_element_sparse_real_restr: Matrix inconsistent (-3).'

       ! -- note: in this special routine NO shift occurs --

       IF ( row .LT. N ) THEN

          ija( row + 1 ) = ija( row + 1 ) + 1
          !   ija(targ+1)=ija(targ)

       ELSE

          ija( N + 1 ) = ija( N + 1 ) + 1

       END IF

       IF ( ija( N + 1 ) .GT. SIZE(ija) + 1 ) STOP '--d---'

       !---- write new element ------------------------------

       if ( num_el .EQ. 0 ) THEN

          targ = ija(row)
          ija(targ) = column
          sa(targ) = value

       ELSE 

          targ = ija(row)
          writeL = .FALSE.
          IF ( targ .GT. SIZE(ija) ) STOP &
               'Error add_element_sparse_real_restr: Matrix inconsistent (-4).'

          DO k = 0, num_el - 1

             IF ( column .EQ. ija( targ + k ) ) STOP &
                  'Error add_element_sparse_real_restr: &
                  &Overwrite existing matrix element!'

             IF ( column .LT. ija( targ + k ) ) THEN

                ija(targ+k+1:targ+num_el) = ija(targ+k:targ+num_el-1)
                sa (targ+k+1:targ+num_el) = sa (targ+k:targ+num_el-1)

                ija(targ+k) = column
                sa (targ+k) = value
                writeL = .TRUE.       

                EXIT

             END IF

          END DO

          IF ( .NOT. writeL ) THEN  ! if writeL == .false.

             k = num_el
             ija(targ+k) = column
             sa(targ+k) = value 

          END IF

       END if

    END IF

    IF ( row .LT. ija(1)-2 ) ija(row+2) = ija(row+1)

    !=========================================================================

  END SUBROUTINE add_element_sparse_real_restr
  

  !===========================================================================

  
  SUBROUTINE write2file_sparse_matrix_real(dim,WhichPart)
   
    !=========================================================================
    !
    ! Write sparse matrix stored in sa,ija format into a file.
    ! The format is like this:
    !
    !   row    column    value
    !   ...     ...       ...
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    INTEGER     , INTENT( IN ) :: dim
    CHARACTER(6), INTENT( IN ) :: WhichPart
    
    INTEGER :: i, j, k, row
    CHARACTER(29) :: filename
    INTEGER ijb
    REAL sb

    !=========================================================================

    IF ( ija(1)-2 .NE. dim ) STOP &
         'Error write2file_sparse_matrix_real: &
         &Matrix has inconsistent dimension.'

    SELECT CASE ( WhichPart )
       
       !-----------------------------------------------------------------------
       
       CASE ('DUpper') ! Diagonal and upper part
       filename='sparse_matrix_real_DUpper.out'

       OPEN  ( 8, file = filename )

       DO i = 1, dim
          
          WRITE ( 8, * ) i, i, sa(i)      ! write diagonal element

          DO k = ija(i), ija(i+1) - 1     ! loop over offdiagonal elements
             WRITE ( 8, * ) i, ija(k), sa(k)
          END DO

       END DO

       CLOSE (8)


       !-----------------------------------------------------------------------
       
       CASE ('DLower') ! Diagonal and lower part
       filename = 'sparse_matrix_real_DLower.out'
       
       !CALL Transpose_sparse_matrix_real(sb,ijb)
       
       OPEN ( 8, file = filename )

       DO i = dim, 1, -1

          DO k = ija(i), ija(i+1) - 1     ! loop over offdiagonal elements

             j = ija(k)
             WRITE ( 8, * ) i, ija(k), sa(k)
 
          END DO

          WRITE ( 8, * ) i, i, sa(i)           ! write diagonal element

       END DO

       CLOSE (8)

       !-----------------------------------------------------------------------

       CASE ('LoDiUp') ! Lower, diagonal and upper part
       filename = 'sparse_matrix_real_LoDiUp.out'

       OPEN ( 8, file = filename )
       DO i = 1, dim
          WRITE ( 8, * ) i, i, sa(i) ! write diagonal element
       END DO

       CLOSE (8)

       !-----------------------------------------------------------------------
       
       CASE ('Diagon') ! Diagonal only

       filename = 'sparse_matrix_real_Diagon.out'

       OPEN ( 8, file = filename )

       DO i = 1, dim
          WRITE ( 8, * ) i, i, sa(i) ! write diagonal element
       END DO

       CLOSE (8)

       !-----------------------------------------------------------------------
   
    CASE DEFAULT
    
       PRINT *,  WhichPart
       STOP 'Error write2file_sparse_matrix_real: WhichPart ill-defined.'
       
       !-----------------------------------------------------------------------
   
    END SELECT

    !=========================================================================
 
  END SUBROUTINE write2file_sparse_matrix_real
  
  !===========================================================================
  
  SUBROUTINE Transpose_sparse_matrix_real( sb, ijb )
  
    !=========================================================================
    !
    ! Construct the transpose of a sparse square matrix,
    ! from row-index sparse storage arrays sa
    ! and ija into arrays sb and ijb.
    !
    !=========================================================================

    USE sparse_matrix_real, ONLY : sa, ija

    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    INTEGER ijb(*)
    REAL    sb(*)

    !USES iindexx Version of indexx with all REAL variables changed to INTEGER.

    INTEGER j, jl, jm, jp, ju, k, m, n2, noff, inc, iv
    REAL    v

    !=========================================================================

    n2 = ija(1)   ! Linear size of matrix plus 2.

    DO j = 1, n2-2 ! Diagonal elements.
       sb(j) = sa(j)
    END DO

    !(ija(n2-1)-n2    ,ija(n2),ijb(n2))

    CALL iindexx( ija(n2-1)-ija(1), ija(n2:n2), ijb(n2) )

    ! Index all off-diagonal elements by their columns.

    jp = 0

    DO k = ija(1), ija(n2-1)-1 ! Loop over output off-diagonal elements.
       m = ijb(k) + n2 - 1 ! Use index table to store by (former) columns.
       sb(k) = sa(m)
       DO j = jp + 1, ija(m) ! Fill in the index to any omitted rows.
          ijb(j) = k
       END DO
       jp = ija(m)
       ! Use bisection to .nd which row element m is in and put that
       
       jl = 1      ! into ijb(k). 
       ju = n2-1
5      IF ( ju-jl .GT. 1 ) THEN
          jm = ( ju+jl ) / 2
          IF ( ija(jm) .GT. m ) THEN
             ju = jm
          ELSE
             jl = jm
          END IF
          GOTO 5
       END IF
       ijb(k) = jl
    END DO

    DO j = jp+1, n2-1
       ijb(j) = ija(n2-1)
    END DO ! Make a final pass to sort each row by Shell sort algorithm.
    
    DO j = 1, n2-2
       jl = ijb(j+1)-ijb(j)
       noff = ijb(j)-1
       inc = 1
1      inc = 3*inc+1
       IF ( inc .LE. jl ) GOTO 1
2      CONTINUE
       inc = inc/3
       DO k = noff+inc+1, noff+jl
          iv = ijb(k)
          v = sb(k)
          m = k
3         IF ( ijb(m-inc) .GT. iv ) THEN
             ijb(m) = ijb(m-inc)
             sb(m) = sb(m-inc)
             m = m-inc
             IF ( m-noff .LE. inc ) GOTO 4
             GOTO 3
          END IF
4         ijb(m) = iv
          sb(m) = v
       END DO
       IF ( inc .GT. 1 ) GOTO 2
    END DO
    RETURN

  CONTAINS

    !--------------------------------------------------------------------------

    SUBROUTINE iindexx(n,arr,indx)
     
      !------------------------------------------------------------------------
      
      INTEGER n,indx(n),M,NSTACK
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      
      !Indexes an array arr(1:n), i.e., outputs the array indx(1:n)
      ! such that arr(indx(j))
      !is in ascending order for j = 1, 2, . . . ,N.
      ! The input quantities n and arr are not changed.
      
      INTEGER i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK)
      INTEGER a

      DO j = 1, n
         indx(j) = j
      END DO

      jstack = 0
      l = 1
      ir = n
1     IF ( ir-l .LT. M ) THEN
         
         DO j = l+1, ir
            indxt = indx(j)
            a = arr(indxt)
            DO i = j-1, l, -1
               IF ( arr( indx(i) ) .LE. a ) GOTO 2
               indx(i+1) = indx(i)
            END DO
            i = l-1
2           indx(i+1) = indxt
         END DO

         IF ( jstack .EQ. 0 ) RETURN
         ir = istack(jstack)
         l = istack(jstack-1)
         jstack = jstack-2

      ELSE

         k = (l+ir)/2
         itemp = indx(k)
         indx(k) = indx(l+1)
         indx(l+1) = itemp
         IF ( arr( indx(l) ) .GT. arr( indx(ir) ) ) THEN
            itemp = indx(l)
            indx(l) = indx(ir)
            indx(ir) = itemp
         END IF

         IF ( arr( indx(l+1) ) .GT. arr( indx(ir) ) ) THEN
            itemp = indx(l+1)
            indx(l+1) = indx(ir)
            indx(ir) = itemp
         END IF
         
         IF ( arr( indx(l) ) .GT. arr( indx(l+1) ) ) THEN
            itemp = indx(l)
            indx(l) = indx(l+1)
            indx(l+1) = itemp
         END IF

         i = l+1
         j = ir
         indxt = indx(l+1)
         a = arr(indxt)

3        CONTINUE

         i = i+1
         IF ( arr( indx(i) ) .LT. a ) GOTO 3

4        CONTINUE

         j = j-1
         IF ( arr( indx(j) ) .GT. a ) GOTO 4
        
         IF ( j .LT. i ) GOTO 5
         
         itemp = indx(i)
         indx(i) = indx(j)
         indx(j) = itemp
         GOTO 3

5        indx(l+1) = indx(j)
         indx(j) = indxt
         jstack = jstack+2
         IF ( jstack .GT. NSTACK ) STOP 'NSTACK too small in indexx'
      
         IF ( ir-i+1 .GE. j-l ) THEN
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j-1
         ELSE
            istack(jstack) = j-1
            istack(jstack-1) = l
            l = i
         END IF
      
      END IF
      
      GOTO 1
      
      !------------------------------------------------------------------------
      
    END SUBROUTINE iindexx
  
    !=========================================================================

  END SUBROUTINE Transpose_sparse_matrix_real

  !===========================================================================

END MODULE sparse_matrix_rout_real

