c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                     UNARY SUBROUTINES MODULE                         c
c----------------------------------------------------------------------c
c     contents:                                                        c
c----------------------------------------------------------------------c
c     
c     --> submat    : extracts a submatrix from a sparse matrix. 
c     --> filter    : filters matrix elements according to their magnitude.
c     --> filterm   : same as above, but for the MSR format.
c     --> csort     : sorts the elements in increasing order of columns.
c     --> clncsr    : clean up the CSR format matrix
c     --> transp    : in-place transposition routine
c     ( see also csrcsc in formats )
c     --> copmat    : copy of a CSR matrix into another CSR matrix
c     --> msrcop    : copies a MSR matrix into a MSR matrix
c     --> getelm    : returns a( i, j ) for any CSR-stored ( i, j )
c     --> getdia    : extracts a specified diagonal from a matrix.
c     --> getl      : extracts lower triangular part.
c     --> getu      : extracts upper triangular part.
c     --> levels    : gets the level scheduling structure
c     for lower triangular matrices.
c     --> amask     : extracts  C = A mask M.
c     --> rperm     : permutes the rows of a matrix ( B = P A )
c     --> cperm     : permutes the columns of a matrix ( B = A Q )
c     --> dperm     : permutes both the rows and columns of a matrix
c     ( B = P A Q )
c     --> dperm1    : general extraction routine ( extracts arbitrary rows )
c     --> dmperm    : symmetric permutation of row and column ( B = PAP' )
c     in MSR format
c     --> dvperm    : permutes a real vector ( in-place )
c     --> ivperm    : permutes an integer vector ( in-place )
c     --> retmx     : returns the max absolute value in each row of the matrix
c     --> diapos    : returns the positions of diagonal elements in A.
c     --> extbdg    : extracts the main diagonal blocks of a matrix.
c     --> getbwd    : returns the bandwidth information on a matrix.
c     --> blkfnd    : finds the block-size of a matrix.
c     --> blkchk    : checks if a given integer is the block size of A.
c     --> infdia    : obtains information on the diagonals of A.
c     --> amubdg    : gets number of nonzeros in each row of A * B.
c     --> aplbdg    : gets number of nonzeros in each row of A + B.
c     --> rnrms     : computes the norms of the rows of A.
c     --> cnrms     : computes the norms of the columns of A.
c     --> roscal    : scales the rows of a matrix by their norms.
c     --> coscal    : scales the columns of a matrix by their norms.
c     --> addblk    : Adds a matrix B into a block of A.
c     --> get1up    : Collects the first elements of each row
c     of the upper triangular portion of a matrix.
c     --> xtrows    : extracts given rows from a matrix in CSR format.
c     --> csrkvstr  :  Finds block row partitioning of a CSR matrix.
c     --> csrkvstc  :  Finds block column partitioning of a CSR matrix.
c     --> kvstmerge : Merges block partitionings, for conformal row/col pattern
c
c     CSR = Compresser Sparse Row
c     MSR = Modified Sparse Row
c
c----------------------------------------------------------------------c
      
      SUBROUTINE submat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)

c-----------------------------------------------------------------------

      INTEGER n, job, i1, i2, j1, j2, nr, nc
      INTEGER ia( * ), ja( * ), jao( * ), iao( * )
      DOUBLE PRECISION a( * ), ao( * )

c-----------------------------------------------------------------------
c
c     extracts the submatrix A(i1:i2,j1:j2) and puts the result in 
c     matrix ao,iao,jao
c
c---- In place: ao,jao,iao may be the same as a,ja,ia.
c
c-----------------------------------------------------------------------
c     on input
c-----------------------------------------------------------------------
c
c n	= row dimension of the matrix 
c i1,i2 = two integers with i2 .GE. i1 indicating the range of rows to be
c          extracted. 
c j1,j2 = two integers with j2 .GE. j1 indicating the range of columns 
c         to be extracted.
c         * There is no checking whether the input values for i1, i2, j1,
c           j2 are between 1 and n. 
c a,
c ja,
c ia    = matrix in compressed sparse row format. 
c
c job	= job indicator: if job .ne. 1 then the real values in a are not
c         extracted, only the column indices (i.e. data structure) are.
c         otherwise values as well as column indices are extracted...
c
c-----------------------------------------------------------------------
c     on output
c-----------------------------------------------------------------------
c
c nr	= number of rows of submatrix 
c nc	= number of columns of submatrix 
c	  * IF either of nr or nc is nonpositive the code will quit.
c
c ao,
c jao,iao = extracted matrix in general sparse format with jao containing
c	the column indices,and iao being the pointer to the beginning 
c	of the row,in arrays a,ja.
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
     
      nr = i2 - i1 + 1
      nc = j2 - j1 + 1
     
      IF ( ( nr .LE. 0 ) .OR. ( nc .LE. 0 ) ) RETURN
     
      klen = 0
     
c     simple procedure. proceeds row-wise...
     
      DO i = 1, nr

         ii = i1 + i - 1
         k1 = ia( ii )
         k2 = ia( ii+1 ) - 1
         iao( i ) = klen + 1

         DO k = k1, k2

            j = ja( k )
            
            IF ( ( j .GE. j1 ) .AND. ( j .LE. j2 ) ) THEN
               klen = klen + 1
               IF ( job .EQ. 1 ) ao( klen ) = a( k )
               jao( klen ) = j - j1 + 1
            END IF
            
         END DO

      END DO

      iao( nr + 1 ) = klen + 1
      
      RETURN

c------------end-of submat---------------------------------------------- 

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE filter(n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)

c-----------------------------------------------------------------------

      DOUBLE PRECISION a(*), b(*), drptol
      INTEGER ja(*), jb(*), ia(*), ib(*), n, job, len, ierr

c-----------------------------------------------------------------------
c
c     This module removes any elements whose absolute value
c     is small from an input matrix A and puts the resulting
c     matrix in B.  The input parameter job selects a definition
c     of small.
c
c-----------------------------------------------------------------------
c on entry:
c-----------------------------------------------------------------------
c
c  n	 = integer. row dimension of matrix
c  job   = integer. used to determine strategy chosen by caller to
c         drop elements from matrix A. 
c          job = 1  
c              Elements whose absolute value is less than the
c              drop tolerance are removed.
c          job = 2
c              Elements whose absolute value is less than the 
c              product of the drop tolerance and the Euclidean
c              norm of the row are removed. 
c          job = 3
c              Elements whose absolute value is less that the
c              product of the drop tolerance and the largest
c              element in the row are removed.
c 
c drptol = real. drop tolerance used for dropping strategy.
c a	
c ja
c ia     = input matrix in compressed sparse format
c len	 = integer. the amount of space available in arrays b and jb.
c
c-----------------------------------------------------------------------
c on return:
c-----------------------------------------------------------------------
c
c b	
c jb
c ib    = resulting matrix in compressed sparse format. 
c 
c ierr	= integer. containing error message.
c         ierr .EQ. 0 indicates normal RETURN
c         ierr .GT. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr THEN contains the row number where filter stopped.
c
c-----------------------------------------------------------------------
c note:
c------ This module is in place. (b,jb,ib can ne the same as 
c       a, ja, ia in which case the result will be overwritten).
c
c----------------------------------------------------------------------c
c           contributed by David Day,  Sep 19, 1989.                   c
c----------------------------------------------------------------------c

c     local variables

      DOUBLE PRECISION norm, loctol
      INTEGER index, row, k, k1, k2

c-----------------------------------------------------------------------

      Index = 1
      
      DO row = 1, n

         k1 = ia( row )
         k2 = ia( row + 1 ) - 1
         ib( row ) = index
         GOTO ( 100, 200, 300 ) job

 100     norm = 1.0D0
         GOTO 400
         
 200     norm = 0.0D0
         DO k = k1, k2
            norm = norm + a( k ) * a( k )
         END DO
         norm = SQRT( norm )
         GOTO 400

 300     norm = 0.0D0
         DO k = k1, k2
            IF ( ABS( a( k ) )  .GT. norm ) THEN
               norm = ABS( a( k ) )
            END IF
         END DO
         
 400     loctol = drptol * norm

         DO k = k1, k2
            
            IF ( ABS( a( k ) ) .GT. loctol ) THEN
               
               IF ( index .GT. len ) THEN
                  ierr = row
                  RETURN
               END IF
               
               b( index ) =  a( k )
               jb( index ) = ja( k )
               index = index + 1
               
            END IF

         END DO

      END DO
      
      ib( n + 1 ) = index

      RETURN

c--------------------end-of-filter -------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE filterm (n,job,drop,a,ja,b,jb,len,ierr)

c-----------------------------------------------------------------------

      DOUBLE PRECISION a(*),b(*),drop
      INTEGER ja(*),jb(*),n,job,len,ierr

c-----------------------------------------------------------------------
c
c     This subroutine removes any elements whose absolute value
c     is small from an input matrix A. Same as filter but
c     uses the MSR format.
c
c-----------------------------------------------------------------------
c on entry:
c-----------------------------------------------------------------------
c
c  n	 = integer. row dimension of matrix
c  job   = integer. used to determine strategy chosen by caller to
c         drop elements from matrix A. 
c          job = 1  
c              Elements whose absolute value is less than the
c              drop tolerance are removed.
c          job = 2
c              Elements whose absolute value is less than the 
c              product of the drop tolerance and the Euclidean
c              norm of the row are removed. 
c          job = 3
c              Elements whose absolute value is less that the
c              product of the drop tolerance and the largest
c              element in the row are removed.
c 
c drop = real. drop tolerance used for dropping strategy.
c a	
c ja     = input matrix in Modifief Sparse Row format
c len	 = integer. the amount of space in arrays b and jb.
c
c-----------------------------------------------------------------------
c on return:
c-----------------------------------------------------------------------
c
c b, jb = resulting matrix in Modifief Sparse Row format
c 
c ierr	= integer. containing error message.
c         ierr .EQ. 0 indicates normal RETURN
c         ierr .GT. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr THEN contains the row number where filter stopped.
c
c----------------------------------------------------------------------c
c note:
c------ This module is in place. (b,jb can ne the same as 
c       a, ja in which case the result will be overwritten).
c
c----------------------------------------------------------------------c
c           contributed by David Day,  Sep 19, 1989.                   c
c----------------------------------------------------------------------c

c local variables

      DOUBLE PRECISION norm, loctol
      INTEGER index, row, k, k1, k2

c----------------------------------------------------------------------c

      index = n + 2

      DO row = 1, n

         k1 = ja( row )
         k2 = ja( row + 1 ) - 1
         jb( row ) = index
         GOTO ( 100, 200, 300 ) job

 100     norm = 1.0D0
         GOTO 400

 200     norm = a( row )**2
         DO k = k1, k2
            norm = norm + a( k ) * a( k )
         END DO
         norm = SQRT( norm )
         GOTO 400

 300     norm = ABS( a( row ) )
         DO k = k1, k2
            norm = MAX( ABS( a( k ) ), norm )
         END DO

 400     loctol = drop * norm

         DO k = k1, k2
            IF ( ABS( a( k ) ) .GT. loctol ) THEN
               IF ( index .GT. len ) THEN
                  ierr = row
                  RETURN
               END IF
               b( index ) =  a( k )
               jb( index ) = ja( k )
               index = index + 1
            END IF
         END DO

      END DO
      
      jb( n + 1 ) = index
      
      RETURN

c--------------------end-of-filterm-------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE csort( n, a, ja, ia, iwork, values ) 

c-----------------------------------------------------------------------

      LOGICAL values
      INTEGER n, ja(*), ia( n + 1 ), iwork(*) 
      DOUBLE PRECISION a(*) 

c-----------------------------------------------------------------------
c
c This routine sorts the elements of  a matrix (stored in Compressed
c Sparse Row Format) in increasing order of their column indices within 
c each row. It uses a form of bucket sort with a cost of O(nnz) where
c nnz = number of nonzero elements. 
c requires an INTEGER work array of length 2*nnz.  
c
c-----------------------------------------------------------------------
c on entry:
c-----------------------------------------------------------------------
c
c n     = the row dimension of the matrix
c a     = the matrix A in compressed sparse row format.
c ja    = the array of column indices of the elements in array a.
c ia    = the array of pointers to the rows. 
c iwork = INTEGER work array of length max ( n+1, 2*nnz ) 
c         where nnz = (ia(n+1)-ia(1))  ) .
c values= logical indicating whether or not the REAL values a(*) must 
c         also be permuted. IF (.not. values) THEN the array a is not
c         touched by csort and can be a dummy array. 
c
c----------------------------------------------------------------------- 
c on return:
c-----------------------------------------------------------------------
c
c the matrix stored in the structure a, ja, ia is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork(1:nnz) contains the permutation used  to rearrange the elements.
c
c----------------------------------------------------------------------- 
c Y. Saad - Feb. 1, 1991.
c-----------------------------------------------------------------------

c     local variables

      INTEGER i, k, j, ifirst, nnz, next  

c     count the number of elements in each column
      
      DO i = 1, n + 1
         iwork( i ) = 0
      END DO

      DO i = 1, n
         DO k = ia( i ), ia( i + 1 ) - 1
            j = ja( k ) + 1
            iwork( j ) = iwork( j ) + 1
         END DO
      END DO

c     compute pointers from lengths.

      iwork( 1 ) = 1
      DO i = 1, n
         iwork( i + 1 ) = iwork( i ) + iwork( i + 1 )
      END DO
 
c     get the positions of the nonzero elements in order of columns.

      ifirst = ia( 1 )
      nnz = ia( n + 1 ) - ifirst
      DO i = 1, n
         DO k = ia( i ), ia( i + 1 ) - 1
            j = ja( k )
            next = iwork( j )
            iwork( nnz + next ) = k
            iwork( j ) = next + 1
         END DO
      END DO

c convert to coordinate format

      DO i = 1, n
         DO k = ia( i ), ia( i + 1 ) - 1
            iwork( k ) = i
         END DO
      END DO

c loop to find permutation: for each element find the correct 
c position in (sorted) arrays a, ja. Record this in iwork. 

      DO k = 1, nnz
        
         ko = iwork( nnz + k )
         irow = iwork( ko )
         next = ia( irow )

c     the current element should go in next position in row. iwork
c     records this position. 
         
         iwork( ko ) = next
         ia( irow )  = next + 1

      END DO

c     perform an in-place permutation of the  arrays.

      CALL ivperm( nnz, ja( ifirst ), iwork )
      IF ( values ) CALL dvperm( nnz, a( ifirst ), iwork )
      
c     reshift the pointers of the original matrix back.
      
      DO i = n, 1, -1
         ia( i + 1 ) = ia( i )
      END DO

      ia( 1 ) = ifirst

      RETURN

c---------------end-of-csort--------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE clncsr( job, value2, nrow, a, ja, ia, indu, iwk )

c-----------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER job, nrow, value2

c     .. Array Arguments ..

      INTEGER ia(nrow+1), indu(nrow), iwk(nrow+1), ja(*)
      DOUBLE PRECISION  a(*)


c-----------------------------------------------------------------------
c
c     This routine performs two tasks to clean up a CSR matrix
c     -- remove duplicate/zero entries,
c     -- perform a partial ordering, new order lower triangular part,
c        main diagonal, upper triangular part.
c
c     On entry:
c
c     job   = options
c         0 -- nothing is done
c         1 -- eliminate duplicate entries, zero entries.
c         2 -- eliminate duplicate entries and perform partial ordering.
c         3 -- eliminate duplicate entries, sort the entries in the
c              increasing order of clumn indices.
c
c     value2  -- 0 the matrix is pattern only (a is not touched)
c                1 matrix has values too.
c     nrow    -- row dimension of the matrix
c     a,ja,ia -- input matrix in CSR format
c
c     On return:
c     a,ja,ia -- cleaned matrix
c     indu    -- pointers to the beginning of the upper triangular
c                portion IF job > 1
c
c     Work space:
c     iwk     -- INTEGER work space of size nrow+1
c
c-----------------------------------------------------------------------

c     .. Local Scalars ..

      INTEGER i, j, k, ko, ipos, kfirst, klast
      DOUBLE PRECISION  tmp

c-----------------------------------------------------------------------

      IF ( job .LE. 0 ) RETURN

c-----------------------------------------------------------------------
c
c     .. eliminate duplicate entries --
c     array INDU is used as marker for existing indices, it is also the
c     location of the entry.
c     IWK is used to stored the old IA array.
c     matrix is copied to squeeze out the space taken by the duplicated
c     entries.
c
c-----------------------------------------------------------------------
 
      DO i = 1, nrow
         indu( i ) = 0
         iwk( i ) = ia( i )
      END DO

      iwk( nrow + 1 ) = ia( nrow + 1 )
      k = 1
      
      DO i = 1, nrow
      
         ia( i ) = k
         ipos = iwk( i )
         klast = iwk( i + 1 )
         
 100     IF ( ipos .LT. klast ) THEN

            j = ja( ipos )
            
            IF ( indu( j ) .EQ. 0 ) THEN
               
c     .. new entry ..
               
               IF ( value2 .NE. 0 ) THEN
                  
                  IF ( a( ipos ) .NE. 0.0D0 ) THEN
                     indu( j ) = k
                     ja( k ) = ja( ipos )
                     a( k ) = a( ipos )
                     k = k + 1
                  END IF
                  
               ELSE
                  
                  indu( j ) = k
                  ja( k ) = ja( ipos )
                  k = k + 1
                  
               END IF
               
            ELSE IF ( value2 .NE. 0 ) THEN
               
c     .. duplicate entry ..
               
               a( indu( j ) ) = a( indu( j ) ) + a( ipos )
               
            END IF
            
            ipos = ipos + 1
            
            GOTO 100

         END IF

c     .. remove marks before working on the next row ..

         DO ipos = ia( i ), k - 1
            indu( ja( ipos ) ) = 0
         END DO

      END DO

      ia( nrow + 1 ) = k
      IF ( job .LE. 1 ) RETURN

c     .. partial ordering ..
c     split the matrix into strict upper/lower triangular
c     parts, INDU points to the the beginning of the upper part.

      DO i = 1, nrow
      
         klast = ia( i + 1 ) - 1
         kfirst = ia( i )

 130     IF ( klast .GT. kfirst ) THEN
            
            IF ( ( ja( klast ) .LT. i ) .AND.
     &           ( ja( kfirst ) .GE. i ) ) THEN

c     .. swap klast with kfirst ..

               j = ja( klast )
               ja( klast ) = ja( kfirst )
               ja( kfirst ) = j
         
               IF ( value2 .NE. 0 ) THEN
                  tmp = a( klast )
                  a( klast ) = a( kfirst )
                  a( kfirst ) = tmp
               END IF

            END IF

            IF ( ja( klast ) .GE. i ) klast = klast - 1
            IF ( ja( kfirst ) .LT. i ) kfirst = kfirst + 1
            GOTO 130
            
         END IF
         
         IF ( ja( klast ) .LT. i ) THEN
            indu( i ) = klast + 1
         ELSE
            indu( i ) = klast
         END IF

      END DO

      IF ( job .LE. 2 ) RETURN

c     .. order the entries according to column indices
c     burble-sort is used

      DO i = 1, nrow

         DO ipos = ia( i ), indu( i ) - 1

            DO j = indu( i ) - 1, ipos + 1, -1

               k = j - 1

               IF ( ja( k ) .GT. ja( j ) ) THEN

                  ko = ja( k )
                  ja( k ) = ja( j )
                  ja( j ) = ko
                  IF ( value2 .NE. 0 ) THEN
                     tmp = a( k )
                     a( k ) = a( j )
                     a( j ) = tmp
                  END IF

               END IF

            END DO

         END DO

         DO ipos = indu( i ), ia( i + 1 ) - 1
           
            DO j = ia( i + 1 ) - 1, ipos + 1, -1

               k = j - 1

               IF ( ja( k ) .GT. ja( j ) ) THEN

                  ko = ja( k )
                  ja( k ) = ja( j )
                  ja( j ) = ko

                  IF ( value2 .NE. 0 ) THEN
                     tmp = a( k )
                     a( k ) = a( j )
                     a( j ) = tmp
                  END IF

               END IF

            END DO

         END DO

      END DO

      RETURN

c---- end of clncsr ----------------------------------------------------

      END

c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 

      SUBROUTINE copmat( nrow, a, ja, ia, ao, jao, iao, ipos, job )

c----------------------------------------------------------------------- 

      DOUBLE PRECISION a(*), ao(*)
      INTEGER nrow, ia(*), ja(*), jao(*), iao(*), ipos, job
      
c----------------------------------------------------------------------
c
c copies the matrix a, ja, ia, into the matrix ao, jao, iao. 
c
c----------------------------------------------------------------------
c on entry:
c----------------------------------------------------------------------- 
c
c nrow	= row dimension of the matrix 
c a,
c ja,
c ia    = input matrix in compressed sparse row format. 
c ipos  = integer. indicates the position in the array ao, jao
c         where the first element should be copied. Thus 
c         iao(1) = ipos on RETURN. 
c job   = job indicator. IF (job .NE. 1) the values are not copies 
c         (i.e., pattern only is copied in the form of arrays ja, ia).
c
c-----------------------------------------------------------------------
c on return:
c-----------------------------------------------------------------------
c
c ao,
c jao,
c iao   = output matrix containing the same data as a, ja, ia.
c
c-----------------------------------------------------------------------
c           Y. Saad, March 1990. 
c-----------------------------------------------------------------------

c local variables

      INTEGER kst, i, k

      kst    = ipos - ia( 1 )
      DO i = 1, nrow + 1
         iao( i ) = ia( i ) + kst
      END DO

      DO k = ia( 1 ), ia( nrow + 1 ) - 1
         jao( kst + k ) = ja( k )
      END DO

      IF ( job .NE. 1 ) RETURN
      DO k = ia( 1 ), ia( nrow + 1 ) - 1
         ao( kst + k ) = a( k )
      END DO

      RETURN

c--------end-of-copmat -------------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE msrcop ( nrow, a, ja, ao, jao, job )

c-----------------------------------------------------------------------

      DOUBLE PRECISION a(*), ao(*)
      INTEGER nrow, ja(*), jao(*), job

c----------------------------------------------------------------------
c
c copies the MSR matrix a, ja, into the MSR matrix ao, jao 
c
c----------------------------------------------------------------------
c on entry:
c-----------------------------------------------------------------------
c
c nrow	= row dimension of the matrix 
c a,ja  = input matrix in Modified compressed sparse row format. 
c job   = job indicator. Values are not copied IF job .NE. 1 
c
c-----------------------------------------------------------------------       
c on return:
c-----------------------------------------------------------------------
c
c ao, jao   = output matrix containing the same data as a, ja.
c
c-----------------------------------------------------------------------
c           Y. Saad, 
c-----------------------------------------------------------------------

c local variables

      INTEGER i, k

      DO i = 1, nrow + 1
         jao( i ) = ja( i )
      END DO

      DO  k = ja( 1 ), ja( nrow + 1 ) - 1
         jao( k ) = ja( k )
      END DO

      IF ( job .NE. 1 ) RETURN
      
      DO k = ja( 1 ), ja( nrow + 1 ) - 1
         ao( k ) = a( k )
      END DO

      DO k = 1, nrow
         ao( k ) = a( k )
      END DO

      RETURN

c--------end-of-msrcop -------------------------------------------------

      END

c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 

      DOUBLE PRECISION FUNCTION getelm( i, j, a, ja, ia, iadd, sorted )

c-----------------------------------------------------------------------
c
c     purpose:
c     -------- 
c     this function returns the element a(i,j) of a matrix a, 
c     for any pair (i,j).  the matrix is assumed to be stored 
c     in compressed sparse row (csr) format. getelm performs a
c     binary search in the case where it is known that the elements 
c     are sorted so that the column indices are in increasing order. 
c     also returns (in iadd) the address of the element a(i,j) in 
c     arrays a and ja when the search is successsful (zero IF not).
c----- 
c     first contributed by noel nachtigal (mit). 
c     recoded jan. 20, 1991, by y. saad [in particular
c     added handling of the non-sorted case + the iadd output] 
c
c-----------------------------------------------------------------------
c
c     parameters:
c     ----------- 
c
c----------------------------------------------------------------------- 
c on entry: 
c----------------------------------------------------------------------- 
c
c     i      = the row index of the element sought (input).
c     j      = the column index of the element sought (input).
c     a      = the matrix a in compressed sparse row format (input).
c     ja     = the array of column indices (input).
c     ia     = the array of pointers to the rows' data (input).
c     sorted = logical indicating whether the matrix is knonw to 
c              have its column indices sorted in increasing order 
c              (sorted=.true.) or not (sorted=.false.).
c              (input). 
c
c----------------------------------------------------------------------- 
c on return:
c----------------------------------------------------------------------- 
c
c     getelm = value of a(i,j). 
c     iadd   = address of element a(i,j) in arrays a, ja IF found,
c              zero IF not found. (output) 
c
c     note: the inputs i and j are not checked for validity. 
c
c-----------------------------------------------------------------------
c     noel m. nachtigal october 28, 1990 -- youcef saad jan 20, 1991.
c----------------------------------------------------------------------- 
 
      INTEGER i, ia(*), iadd, j, ja(*)
      DOUBLE PRECISION a(*)
      LOGICAL sorted

c     local variables.

      INTEGER ibeg, iend, imid, k

c     initialization 

      iadd = 0
      getelm = 0.0D0
      ibeg = ia( i )
      iend = ia( i + 1 ) - 1

c     case where matrix is not necessarily sorted
     
      IF ( .NOT. sorted ) THEN

c     scan the row - exit as soon as a(i,j) is found

         DO k = ibeg, iend
            IF ( ja( k ) .EQ. j ) THEN
               iadd = k
               GOTO 20
            END IF
         END DO
    
c     end unsorted case. begin sorted case
     
      ELSE

c     begin binary search.   compute the middle index.
   
 10      imid = ( ibeg + iEND ) / 2
         
c     test if  found
    
         IF ( ja( imid ) .EQ. j ) THEN
            iadd = imid
            GOTO 20
         END IF

         IF ( ibeg .GE. iend ) GOTO 20
     
c     else update the interval bounds. 
     
         IF ( ja( imid ) .GT. j ) THEN
            iend = imid - 1
         ELSE
            ibeg = imid + 1
         END IF
         GOTO 10
    
c     end both cases
     
      END IF
     
 20   IF ( iadd .NE. 0 ) getelm = a( iadd )

      RETURN

c--------end-of-getelm--------------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)

c-----------------------------------------------------------------------

      DOUBLE PRECISION diag(*), a(*)
      INTEGER nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)

c-----------------------------------------------------------------------
c
c this subroutine extracts a given diagonal from a matrix stored in csr 
c format. the output matrix may be transformed with the diagonal removed
c from it IF desired (as indicated by job.) 
c
c----------------------------------------------------------------------- 
c
c our definition of a diagonal of matrix is a vector of length nrow
c (always) which contains the elements in rows 1 to nrow of
c the matrix that are contained in the diagonal offset by ioff
c with respect to the main diagonal. IF the diagonal element
c falls outside the matrix THEN it is defined as a zero entry.
c thus the proper definition of diag(*) with offset ioff is 
c
c     diag(i) = a(i,ioff+i) i=1,2,...,nrow
c     with elements falling outside the matrix being defined as zero.
c 
c------------------------------------------------------------------------ 
c on entry:
c------------------------------------------------------------------------ 
c
c nrow	= integer. the row dimension of the matrix a.
c ncol	= integer. the column dimension of the matrix a.
c job   = integer. job indicator.  IF job = 0 THEN
c         the matrix a, ja, ia, is not altered on RETURN.
c         IF job.NE.0  THEN getdia will remove the entries
c         collected in diag from the original matrix.
c         this is done in place.
c
c a,ja,
c    ia = matrix stored in compressed sparse row a,ja,ia,format
c ioff  = INTEGER,containing the offset of the wanted diagonal
c	  the diagonal extracted is the one corresponding to the
c	  entries a(i,j) with j-i = ioff.
c	  thus ioff = 0 means the main diagonal
c
c------------------------------------------------------------------------ 
c on return:
c------------------------------------------------------------------------ 
c
c len   = number of nonzero elements found in diag.
c         (len .LE. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
c
c diag  = DOUBLE PRECISION array of length nrow containing the wanted diagonal.
c	  diag contains the diagonal (a(i,j),j-i = ioff ) as defined 
c         above. 
c
c idiag = INTEGER array of  length len, containing the poisitions 
c         in the original arrays a and ja of the diagonal elements
c         collected in diag. a zero entry in idiag(i) means that 
c         there was no entry found in row i belonging to the diagonal.
c         
c a, ja,
c    ia = IF job .NE. 0 the matrix is unchanged. otherwise the nonzero
c         diagonal entries collected in diag are removed from the 
c         matrix and therefore the arrays a, ja, ia will change.
c	  (the matrix a, ja, ia will contain len fewer elements) 
c 
c----------------------------------------------------------------------c
c     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c 
c----------------------------------------------------------------------c

c     local variables

      INTEGER istart, max, iend, i, kold, k, kdiag, ko
     
      istart = max( 0, -ioff )
      iend = min( nrow, ncol - ioff )
      len = 0

      DO i = 1, nrow
         idiag( i ) = 0
         diag( i ) = 0.0D0
      END DO

c     extract  diagonal elements

      DO 6 i = istart + 1, iend
         DO k = ia( i ), ia( i + 1 ) - 1
            IF ( ja( k ) - i .EQ. ioff ) THEN
               diag( i )= a( k )
               idiag( i ) = k
               len = len + 1
               GOTO 6
            END IF
         END DO
 6    CONTINUE
      
      IF ( ( job .EQ. 0 ) .OR. ( len .EQ.0 ) ) RETURN

c     remove diagonal elements and rewind structure
 
      ko = 0
      DO i = 1, nrow
         kold = ko
         kdiag = idiag( i )
         DO k = ia( i ), ia( i + 1 ) - 1
            IF ( k .NE. kdiag ) THEN
               ko = ko + 1
               a( ko ) = a( k )
               ja( ko ) = ja( k )
            END IF
         END DO
         ia( i ) = kold + 1
      END DO

c     redefine ia(nrow+1)

      ia( nrow + 1 ) = ko + 1

      RETURN

c------------end-of-getdia----------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE transp( nrow, ncol, a, ja, ia, iwk, ierr )

c-----------------------------------------------------------------------

      INTEGER nrow, ncol, ia(*), ja(*), iwk(*), ierr
      DOUBLE PRECISION a(*)

c------------------------------------------------------------------------
c In-place transposition routine.
c------------------------------------------------------------------------
c
c this SUBROUTINE transposes a matrix stored in compressed sparse row 
c format. the transposition is done in place in that the arrays a,ja,ia
c of the transpose are overwritten onto the original arrays.
c
c------------------------------------------------------------------------
c on entry:
c-----------------------------------------------------------------------
c
c nrow	= integer. The row dimension of A.
c ncol	= integer. The column dimension of A.
c a	= REAL array of size nnz (number of nonzero elements in A).
c         containing the nonzero elements 
c ja	= INTEGER array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= INTEGER of size n+1, where n = max(nrow,ncol). On entry
c         ia(k) contains the position in a,ja of  the beginning of 
c         the k-th row.
c
c iwk	= INTEGER work array of same length as ja.
c
c-----------------------------------------------------------------------
c on return:
c-----------------------------------------------------------------------
c
c ncol	= actual row dimension of the transpose of the input matrix.
c         Note that this may be .LE. the input value for ncol, in
c         case some of the last columns of the input matrix are zero
c         columns. In the case where the actual number of rows found
c         in transp(A) exceeds the input value of ncol, transp will
c         RETURN without completing the transposition. see ierr.
c a,
c ja,
c ia	= contains the transposed matrix in compressed sparse
c         row format. The row dimension of a, ja, ia is now ncol.
c
c ierr	= integer. error message. IF the number of rows for the
c         transposed matrix exceeds the input value of ncol,
c         THEN ierr is  set to that number and transp quits.
c         Otherwise ierr is set to 0 (normal RETURN).
c
c-----------------------------------------------------------------------
c Note: 
c----- 1) IF you DO not need the transposition to be done in place
c         it is preferrable to use the conversion routine csrcsc 
c         (see conversion routines in formats).
c      2) the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use csrcsc
c         IF you want them sorted.
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c  modified Oct. 11, 1989.                                             c
c----------------------------------------------------------------------c

c local variables

      DOUBLE PRECISION t, t1

      ierr = 0
      nnz = ia( nrow + 1 ) - 1

c     determine column dimension

      jcol = 0

      DO k = 1, nnz
         jcol = max( jcol, ja( k ) )
      END DO

      IF ( jcol .GT. ncol ) THEN
         ierr = jcol
         RETURN
      END IF

c     convert to coordinate format. use iwk for row indices.

      ncol = jcol

      DO i = 1, nrow
         DO k = ia( i ), ia( i + 1 ) - 1
            iwk( k ) = i
         END DO
      END DO

c     find pointer array for transpose. 

      DO i = 1, ncol + 1
         ia( i ) = 0
      END DO

      DO k = 1, nnz
         i = ja( k )
         ia( i + 1 ) = ia( i + 1 ) + 1
      END DO

      ia( 1 ) = 1

c------------------------------------------------------------------------

      DO i = 1, ncol
         ia( i + 1 ) = ia( i ) + ia( i + 1 )
      END DO

c     loop for a cycle in chasing process. 
     
      init = 1
      k = 0

 5    t = a( init )
      i = ja( init )
      j = iwk( init )
      iwk( init ) = -1
      
c------------------------------------------------------------------------

 6    k = k + 1

c     current row number is i.  determine  where to go.
 
      l = ia( i )

c     save the chased element.
 
      t1 = a( l )
      inext = ja( l )

c     then occupy its location.

      a( l )  = t
      ja( l ) = j

c     update pointer information for next element to be put in row i.
 
      ia( i ) = l + 1

c     determine  next element to be chased

      IF ( iwk( l ) .LT. 0 ) GOTO 65
      t = t1
      i = inext
      j = iwk( l )
      iwk( l ) = -1
      IF ( k .LT. nnz ) GOTO 6
      GOTO 70

 65   init = init + 1
      IF ( init .GT. nnz ) GOTO 70
      IF ( iwk( init ) .LT. 0 ) GOTO 65

c     restart chasing --	

      GOTO 5
 
 70   CONTINUE
      
      DO i = ncol, 1, -1
         ia( i + 1 ) = ia( i )
      END DO

      ia( 1 ) = 1

      RETURN

c------------------end-of-transp ----------------------------------------

      END 

c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 

      SUBROUTINE getl( n, a, ja, ia, ao, jao, iao )

c------------------------------------------------------------------------
 
      INTEGER n, ia(*), ja(*), iao(*), jao(*)
      DOUBLE PRECISION a(*), ao(*)

c------------------------------------------------------------------------
c
c this subroutine extracts the lower triangular part of a matrix 
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia IF desired.
c
c------------------------------------------------------------------------ 
c on input:
c------------------------------------------------------------------------ 
c
c n     = dimension of the matrix a.
c a, ja, 
c    ia = matrix stored in compressed sparse row format.
c On return:
c ao, jao, 
c    iao = lower triangular matrix (lower part of a) 
c	stored in a, ja, ia, format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 ) 
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getl will overwrite the result on a, ja, ia.
c
c------------------------------------------------------------------------

c local variables

      DOUBLE PRECISION t
      INTEGER ko, kold, kdiag, k, i

c inititialize ko (pointer for output matrix)

      ko = 0
      
      DO i = 1, n

         kold = ko
         kdiag = 0
         
         DO 71 k = ia( i ), ia( i + 1 ) -1
            IF ( ja( k )  .GT. i ) GOTO 71
            ko = ko + 1
            ao( ko ) = a( k )
            jao( ko ) = ja( k )
            IF (ja( k )  .EQ. i ) kdiag = ko
 71      CONTINUE

         IF ( ( kdiag .EQ. 0 ) .OR. ( kdiag .EQ. ko ) ) GOTO 72

c     exchange
         
         t = ao( kdiag )
         ao( kdiag ) = ao( ko )
         ao( ko ) = t
         
         k = jao( kdiag )
         jao( kdiag ) = jao( ko )
         jao( ko ) = k
 72      iao( i ) = kold + 1
         
      END DO

c     redefine iao(n+1)

      iao( n + 1 ) = ko + 1

      RETURN

c----------end-of-getl ------------------------------------------------- 

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE getu( n, a, ja, ia, ao, jao, iao )

c-----------------------------------------------------------------------

      INTEGER n, ia(*), ja(*), iao(*), jao(*)
      DOUBLE PRECISION a(*), ao(*)

c------------------------------------------------------------------------
c
c this subroutine extracts the upper triangular part of a matrix 
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia IF desired.
c
c-----------------------------------------------------------------------
c on input:
c-----------------------------------------------------------------------
c
c n     = dimension of the matrix a.
c a, ja, 
c    ia = matrix stored in a, ja, ia, format
c On return:
c ao, jao, 
c    iao = upper triangular matrix (upper part of a) 
c	stored in compressed sparse row format
c
c-----------------------------------------------------------------------
c
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 ) 
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getu will overwrite the result on a, ja, ia.
c
c------------------------------------------------------------------------

c local variables

      DOUBLE PRECISION t
      INTEGER ko, k, i, kdiag, kfirst

      ko = 0

      DO i = 1, n
      
         kfirst = ko + 1
         kdiag = 0
         
         DO 71 k = ia( i ), ia( i + 1 ) - 1
         
            IF ( ja( k )  .LT. i ) GOTO 71
            ko = ko + 1
            ao( ko ) = a( k )
            jao( ko ) = ja( k )
            IF ( ja( k ) .EQ. i ) kdiag = ko

 71      CONTINUE

         IF ( ( kdiag .EQ. 0 ) .OR. ( kdiag .EQ. kfirst ) ) GOTO 72

c     exchange

         t = ao( kdiag )
         ao( kdiag ) = ao( kfirst )
         ao( kfirst ) = t
     
         k = jao( kdiag )
         jao( kdiag ) = jao( kfirst )
         jao( kfirst ) = k
 72      iao( i ) = kfirst

      END DO

c     redefine iao(n+1)

      iao( n + 1 ) = ko + 1

      RETURN

c----------end-of-getu ------------------------------------------------- 

      END

c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 

      SUBROUTINE levels ( n, jal, ial, nlev, lev, ilev, levnum )

c----------------------------------------------------------------------- 

      INTEGER jal(*),ial(*), levnum(*), ilev(*), lev(*)

c-----------------------------------------------------------------------
c
c levels gets the level structure of a lower triangular matrix 
c for level scheduling in the parallel solution of triangular systems
c strict lower matrices (e.g. unit) as well matrices with their main 
c diagonal are accepted. 
c
c-----------------------------------------------------------------------
c on entry:
c----------------------------------------------------------------------- 
c
c n        = integer. The row dimension of the matrix
c jal, ial = 
c 
c----------------------------------------------------------------------- 
c on return:
c----------------------------------------------------------------------- 
c
c nlev     = integer. number of levels found
c lev      = INTEGER array of length n containing the level
c            scheduling permutation.
c ilev     = INTEGER array. pointer to beginning of levels in lev.
c            the numbers lev(i) to lev(i+1)-1 contain the row numbers
c            that belong to level number i, in the level scheduling
c            ordering. The equations of the same level can be solved
c            in parallel, once those of all the previous levels have
c            been solved.
c----------------------------------------------------------------------- 
c work arrays:
c-----------------------------------------------------------------------
c 
c levnum   = INTEGER array of length n (containing the level numbers
c            of each unknown on RETURN)
c
c-----------------------------------------------------------------------
 
      DO i = 1, n
         levnum( i ) = 0
      END DO

c     compute level of each node --

      nlev = 0
      
      DO i = 1, n
         levi = 0
         DO j = ial( i ), ial( i + 1 ) - 1
            levi = MAX( levi, levnum( jal( j ) ) )
         END DO
         levi = levi + 1
         levnum( i ) = levi
         nlev = MAX( nlev, levi )
      END DO

c-------------set data structure  --------------------------------------
     
      DO j = 1, nlev + 1
         ilev( j ) = 0
      END DO

c------count  number   of elements in each level ----------------------- 

      DO j = 1, n
         i = levnum( j ) + 1
         ilev( i ) = ilev( i ) + 1
      END DO

c---- set up pointer for  each  level ---------------------------------- 

      ilev( 1 ) = 1
      DO j = 1, nlev
         ilev( j + 1 ) = ilev( j ) + ilev( j + 1 )
      END DO

c-----determine elements of each level -------------------------------- 

      DO j = 1, n
         i = levnum( j )
         lev( ilev( i ) ) = j
         ilev( i ) = ilev( i ) + 1
      END DO

c     reset pointers backwards

      DO j = nlev, 1, -1
         ilev( j + 1 ) = ilev( j )
      END DO

      ilev( 1 ) = 1

      RETURN

c----------end-of-levels------------------------------------------------ 

      END

c-----------------------------------------------------------------------
c---------------------------------------------------------------------

      SUBROUTINE amask( nrow, ncol, a, ja, ia, jmask, imask,
     &     c,jc, ic, iw, nzmax, ierr )

c---------------------------------------------------------------------

      DOUBLE PRECISION a(*), c(*)
      INTEGER ia( nrow + 1 ), ja(*), jc(*)
      INTEGER ic( nrow + 1 ), jmask(*), imask( nrow + 1 )
      LOGICAL iw( ncol )

c-----------------------------------------------------------------------
c
c This subroutine builds a sparse matrix from an input matrix by 
c extracting only elements in positions defined by the mask jmask, imask
c
c-----------------------------------------------------------------------
c On entry:
c-----------------------------------------------------------------------
c
c nrow  = integer. row dimension of input matrix 
c ncol	= integer. Column dimension of input matrix.
c
c a,
c ja,
c ia	= matrix in Compressed Sparse Row format
c
c jmask,
c imask = matrix defining mask (pattern only) stored in compressed
c         sparse row format.
c
c nzmax = length of arrays c and jc. see ierr.
c
c----------------------------------------------------------------------- 
c On return:
c-----------------------------------------------------------------------
c
c a, ja, ia and jmask, imask are unchanged.
c
c c
c jc, 
c ic	= the output matrix in Compressed Sparse Row format.
c 
c ierr  = integer. serving as error message.c
c         ierr = 1  means normal RETURN
c         ierr .GT. 1 means that amask stopped when processing
c         row number ierr, because there was not enough space in
c         c, jc according to the value of nzmax.
c
c-----------------------------------------------------------------------
c work arrays:
c-----------------------------------------------------------------------
c
c iw	= logical work array of length ncol.
c
c-----------------------------------------------------------------------
c note: 
c------ the  algorithm is in place: c, jc, ic can be the same as 
c a, ja, ia in which cas the code will overwrite the matrix c
c on a, ja, ia
c
c-----------------------------------------------------------------------

      ierr = 0
      len = 0
      DO j = 1, ncol
         iw( j ) = .FALSE.
      END DO

c     unpack the mask for row ii in iw

      DO ii = 1, nrow

c     save pointer in order to be able to do things in place

         DO k = imask( ii ), imask( ii + 1 ) - 1
            iw( jmask( k ) ) = .TRUE.
         END DO

c     add umasked elemnts of row ii

         k1 = ia( ii )
         k2 = ia( ii + 1 ) - 1
         ic( ii ) = len + 1
         DO k = k1, k2
            j = ja( k )
            IF ( iw( j ) ) THEN
               len = len + 1
               IF ( len .GT. nzmax ) THEN
                  ierr = ii
                  RETURN
               END IF
               jc( len ) = j
               c( len ) = a( k )
            END IF
         END DO

         DO k = imask( ii ), imask( ii + 1 ) - 1
            iw( jmask( k ) ) = .FALSE.
         END DO

      END DO
      
      ic( nrow + 1 ) = len + 1

      RETURN

c-----end-of-amask -----------------------------------------------------

      END

c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 

      SUBROUTINE rperm( nrow, a, ja, ia, ao, jao, iao, perm, job )

c----------------------------------------------------------------------- 

      INTEGER nrow, ja(*), ia( nrow + 1 )
      INTEGER jao(*), iao( nrow + 1 ), perm( nrow ), job
      DOUBLE PRECISION a(*), ao(*) 

c-----------------------------------------------------------------------
c
c this subroutine permutes the rows of a matrix in CSR format. 
c rperm  computes B = P A  where P is a permutation matrix.  
c the permutation P is defined through the array perm: for each j, 
c perm(j) represents the destination row number of row number j. 
c Youcef Saad -- recoded Jan 28, 1991.
c
c-----------------------------------------------------------------------
c on entry:
c----------------------------------------------------------------------- 
c
c n 	= dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm 	= INTEGER array of length nrow containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix. 
c         ---> a(i,j) in the original matrix becomes a(perm(i),j) 
c         in the output  matrix.
c
c job	= INTEGER indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of REAL values ao and
c                       the array iao).
c 		job .NE. 1 :  ignore REAL values.
c                     (in which case arrays a and ao are not needed nor
c                      used).
c
c----------------------------------------------------------------------- 
c on return: 
c----------------------------------------------------------------------- 
c 
c ao, jao, iao = input matrix in a, ja, ia format
c
c----------------------------------------------------------------------- 
c
c note : 
c        IF (job.NE.1)  THEN the arrays a and ao are not used.
c
c----------------------------------------------------------------------c
c           Y. Saad, May  2, 1990                                      c
c----------------------------------------------------------------------c
 
      LOGICAL values

      values = ( job .EQ. 1 )
     
c     determine pointers for output matix 
     
      DO j = 1, nrow
         i = perm( j )
         iao( i + 1 ) = ia( j + 1 ) - ia( j )
      END DO

c get pointers from lengths

      iao(1) = 1
      DO j = 1, nrow
         iao( j + 1 ) = iao( j + 1 ) + iao( j )
      END DO

c copying 

      DO ii = 1, nrow

c old row = ii  -- new row = iperm(ii) -- ko = new pointer
        
         ko = iao( perm( ii ) )
         DO k = ia( ii ), ia( ii + 1 ) - 1
            jao( ko ) = ja( k )
            IF ( values ) ao( ko ) = a( k )
            ko = ko + 1
         END DO

      END DO

      RETURN

c---------end-of-rperm ------------------------------------------------- 

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE cperm( nrow, a, ja, ia, ao, jao, iao, perm, job )

c-----------------------------------------------------------------------

      INTEGER nrow, ja(*), ia(nrow+1), jao(*), iao(nrow+1), perm(*), job
      DOUBLE PRECISION a(*), ao(*) 

c-----------------------------------------------------------------------
c
c this subroutine permutes the columns of a matrix a, ja, ia.
c the result is written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on RETURN 
c      a(i,j) becomes a(i,perm(j)) in new matrix 
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
c 
c-----------------------------------------------------------------------
c on entry:
c-----------------------------------------------------------------------
c
c nrow 	= row dimension of the matrix
c
c a, ja, ia = input matrix in csr format. 
c
c perm	= INTEGER array of length ncol (number of columns of A
c         containing the permutation array  the columns: 
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c job	= INTEGER indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of REAL values ao and
c                       the array iao).
c 		job .NE. 1 :  ignore REAL values ao and ignore iao.
c
c-----------------------------------------------------------------------
c on return: 
c-----------------------------------------------------------------------
c
c ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
c
c-----------------------------------------------------------------------
c Notes:
c------- 
c 1. IF job=1 THEN ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same. 
c 3. IF the matrix is initially sorted (by increasing column number) 
c    THEN ao,jao,iao  may not be on RETURN. 
c 
c----------------------------------------------------------------------c

c local parameters:

      INTEGER k, i, nnz

      nnz = ia( nrow + 1 ) - 1

      DO k = 1, nnz
         jao( k ) = perm( ja( k ) )
      END DO

c     done with ja array. RETURN IF no need to touch values.

      IF ( job .NE. 1 ) RETURN

c     else get new pointers -- and copy values too.
      
      DO i = 1, nrow + 1
         iao( i ) = ia( i )
      END DO

      DO k = 1, nnz
         ao( k ) = a( k )
      END DO

      RETURN

c---------end-of-cperm-------------------------------------------------- 

      END

c----------------------------------------------------------------------- 

      SUBROUTINE dperm( nrow, a, ja, ia, ao, jao, iao, perm,
     &     qperm, job )

c----------------------------------------------------------------------- 

      INTEGER nrow, ja(*), ia(nrow+1), jao(*), iao(nrow+1), perm(nrow),
     &     qperm(*), job
      DOUBLE PRECISION a(*), ao(*) 

c-----------------------------------------------------------------------
c
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices. 
c P maps row i into row perm(i) and Q maps column j into column qperm(j): 
c      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
c In the particular case where Q is the transpose of P (symmetric 
c permutation of A) THEN qperm is not needed. 
c note that qperm should be of length ncol (number of columns) but this
c is not checked. 
c
c-----------------------------------------------------------------------
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c-----------------------------------------------------------------------
c
c n 	= dimension of the matrix
c a, ja, 
c    ia = input matrix in a, ja, ia format
c perm 	= INTEGER array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .LE. 2) 
c
c qperm	= same thing for the columns. This should be provided only
c         IF job=3 or job=4, i.e., only in the case of a nonsymmetric
c	  permutation of rows and columns. Otherwise qperm is a dummy
c
c job	= INTEGER indicating the work to be done:
c * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c 		job = 2 permute matrix ignoring REAL values.
c * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q 
c 		job = 3	permute a, ja, ia into ao, jao, iao 
c 		job = 4 permute matrix ignoring REAL values.
c
c-----------------------------------------------------------------------
c on return: 
c-----------------------------------------------------------------------
c
c ao, jao, iao = input matrix in a, ja, ia format
c
c in case job .EQ. 2 or job .EQ. 4, a and ao are never referred to 
c and can be dummy arguments. 
c
c-----------------------------------------------------------------------
c Notes:
c------- 
c  1) algorithm is in place 
c  2) column indices may not be sorted on RETURN even  though they may be 
c     on entry.
c----------------------------------------------------------------------c

c local variables 

      INTEGER locjob, mod

c     locjob indicates whether or not REAL values must be copied. 
     
      locjob = mod( job, 2 )

c     permute rows first 
 
      CALL rperm( nrow, a, ja, ia, ao, jao, iao, perm, locjob )

c     then permute columns

      locjob = 0

      IF ( job .LE. 2 ) THEN
         CALL cperm( nrow,ao,jao,iao,ao,jao,iao,perm,locjob )
      ELSE
         CALL cperm( nrow,ao,jao,iao,ao,jao,iao,qperm,locjob )
      END IF 
     
      RETURN

c-------end-of-dperm----------------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE dperm1( i1, i2, a, ja, ia, b, jb, ib, 
     &     perm, ipos, job )

c-----------------------------------------------------------------------

      INTEGER i1, i2, job, ja(*), ia(*), jb(*), ib(*), perm(*)
      DOUBLE PRECISION a(*), b(*)

c----------------------------------------------------------------------- 
c     general submatrix extraction routine.
c----------------------------------------------------------------------- 
c
c     extracts rows perm(i1), perm(i1+1), ..., perm(i2) (in this order) 
c     from a matrix (doing nothing in the column indices.) The resulting 
c     submatrix is constructed in b, jb, ib. A pointer ipos to the
c     beginning of arrays b,jb,is also allowed (i.e., nonzero elements
c     are accumulated starting in position ipos of b, jb). 
c
c-----------------------------------------------------------------------
c Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991 / modified for PSPARSLIB 
c Sept. 1997.. 
c-----------------------------------------------------------------------
c on entry: 
c-----------------------------------------------------------------------
c 
c n 	= dimension of the matrix
c a,ja,
c   ia  = input matrix in CSR format
c perm 	= INTEGER array of length n containing the indices of the rows
c         to be extracted. 
c
c job   = job indicator. IF (job .NE.1) values are not copied (i.e.,
c         only pattern is copied).
c
c-----------------------------------------------------------------------
c on return: 
c-----------------------------------------------------------------------
c
c b,ja,
c ib   = matrix in csr format. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1) 
c     contain the value and column indices respectively of the nnz
c     nonzero elements of the permuted matrix. thus ib(1)=ipos.
c
c-----------------------------------------------------------------------
c Notes:
c------- 
c  algorithm is NOT in place 
c-----------------------------------------------------------------------
 
c local variables

      INTEGER ko, irow, k
      LOGICAL values

      values = ( job .EQ. 1 )
      ko = ipos
      ib( 1 ) = ko

      DO i = i1, i2
         irow = perm( i )
         DO k = ia( irow ), ia( irow + 1 ) - 1
            IF ( values ) b( ko ) = a( k )
            jb( ko ) = ja( k )
            ko = ko + 1
         END DO
         ib( i - i1 + 2 ) = ko
      END DO

      RETURN

c--------end-of-dperm1--------------------------------------------------

      END

c-----------------------------------------------------------------------

      SUBROUTINE dperm2( i1, i2, a, ja, ia, b, jb, ib,
     &     cperm, rperm, istart, ipos, job )

c-----------------------------------------------------------------------

      INTEGER i1,i2,job,istart,ja(*),ia(*),jb(*),ib(*),cperm(*),rperm(*)
      DOUBLE PRECISION a(*), b(*)

c----------------------------------------------------------------------- 
c     general submatrix permutation/ extraction routine.
c----------------------------------------------------------------------- 
c
c     extracts rows rperm(i1), rperm(i1+1), ..., rperm(i2) and does an
c     associated column permutation (using array cperm). The resulting
c     submatrix is constructed in b, jb, ib. For added flexibility, the
c     extracted elements are put in sequence starting from row 'istart' 
c     of B. In addition a pointer ipos to the beginning of arrays b,jb,
c     is also allowed (i.e., nonzero elements are accumulated starting in
c     position ipos of b, jb). In most applications istart and ipos are 
c     equal to one. However, the generality adds substantial flexiblity.
c     EXPLE: (1) to permute msr to msr (excluding diagonals) 
c     call dperm2 (1,n,a,ja,ja,b,jb,jb,rperm,rperm,1,n+2) 
c            (2) To extract rows 1 to 10: define rperm and cperm to be
c     identity permutations (rperm(i)=i, i=1,n) and THEN
c            call dperm2 (1,10,a,ja,ia,b,jb,ib,rperm,rperm,1,1) 
c            (3) to achieve a symmetric permutation as defined by perm: 
c            call dperm2 (1,10,a,ja,ia,b,jb,ib,perm,perm,1,1) 
c            (4) to get a symmetric permutation of A and appEND the
c            resulting data structure to A's data structure (useful!) 
c            call dperm2 (1,10,a,ja,ia,a,ja,ia(n+1),perm,perm,1,ia(n+1))
c-----------------------------------------------------------------------
c Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c-----------------------------------------------------------------------
c
c n 	= dimension of the matrix
c i1,i2 = extract rows rperm(i1) to rperm(i2) of A, with i1<i2.
c
c a,ja,
c   ia  = input matrix in CSR format
c cperm = INTEGER array of length n containing the permutation arrays
c	  for the columns: cperm(i) is the destination of column j, 
c         i.e., any column index ja(k) is transformed into cperm(ja(k)) 
c
c rperm	=  permutation array for the rows. rperm(i) = origin (in A) of
c          row i in B. This is the reverse permutation relative to the
c          ones used in routines cperm, dperm,.... 
c          rows rperm(i1), rperm(i1)+1, ... rperm(i2) are 
c          extracted from A and stacked into B, starting in row istart
c          of B. 
c istart= starting row for B where extracted matrix is to be added.
c         this is also only a pointer of the be beginning address for
c         ib , on RETURN. 
c ipos  = beginning position in arrays b and jb where to start copying 
c         elements. Thus, ib(istart) = ipos. 
c
c job   = job indicator. IF (job .NE.1) values are not copied (i.e.,
c         only pattern is copied).
c
c-----------------------------------------------------------------------
c on return:
c-----------------------------------------------------------------------
c
c b,ja,
c ib   = matrix in csr format. positions 1,2,...,istart-1 of ib 
c     are not touched. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1) 
c     contain the value and column indices respectively of the nnz
c     nonzero elements of the permuted matrix. thus ib(istart)=ipos.
c
c-----------------------------------------------------------------------
c Notes:
c------- 
c  1) algorithm is NOT in place 
c  2) column indices may not be sorted on RETURN even  though they 
c     may be on entry.
c-----------------------------------------------------------------------
c 
c local variables

      INTEGER ko, irow, k 
      LOGICAL values

c-----------------------------------------------------------------------

      values = ( job .EQ. 1 )
      ko = ipos
      ib(istart) = ko
      DO i = i1, i2
         irow = rperm( i )
         DO k = ia( irow ), ia( irow + 1 ) - 1
            IF ( values ) b( ko ) = a( k )
            jb( ko ) = cperm( ja( k ) )
            ko = ko + 1
         END DO
         ib( istart + i - i1 + 1 ) = ko
      END DO

      RETURN

c--------end-of-dperm2--------------------------------------------------

      END

c-----------------------------------------------------------------------

      SUBROUTINE dmperm( nrow, a, ja, ao, jao, perm, job )

c-----------------------------------------------------------------------

      INTEGER nrow, ja(*), jao(*), perm(nrow), job
      DOUBLE PRECISION a(*), ao(*) 

c-----------------------------------------------------------------------
c
c This routine performs a symmetric permutation of the rows and 
c columns of a matrix stored in MSR format. i.e., it computes 
c B = P A transp(P), where P, is  a permutation matrix. 
c P maps row i into row perm(i) and column j into column perm(j): 
c      a(i,j)    becomes   a(perm(i),perm(j)) in new matrix
c (i.e.  ao(perm(i),perm(j)) = a(i,j) ) 
c calls dperm. 
c
c-----------------------------------------------------------------------
c Y. Saad, Nov 15, 1991. 
c-----------------------------------------------------------------------
c on entry: 
c-----------------------------------------------------------------------
c
c n 	= dimension of the matrix
c a, ja = input matrix in MSR format. 
c perm 	= INTEGER array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .LE. 2) 
c
c job	= INTEGER indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c 		job = 2 permute matrix ignoring REAL values.
c
c----------------------------------------------------------------------
c on return: 
c-----------------------------------------------------------------------
c
c ao, jao = output matrix in MSR. 
c
c in case job .EQ. 2 a and ao are never referred to and can be dummy 
c arguments. 
c
c-----------------------------------------------------------------------
c Notes:
c------- 
c  1) algorithm is NOT in place 
c  2) column indices may not be sorted on RETURN even  though they may be 
c     on entry.
c
c----------------------------------------------------------------------c

c     local variables
     
      INTEGER n1, n2

      n1 = nrow + 1
      n2 = n1 + 1
     
      CALL dperm( nrow,a,ja,ja,ao(n2),jao(n2),jao,perm,perm,job )
     
      jao( 1 ) = n2
      DO j = 1, nrow
         ao( perm( j ) ) = a( j )
         jao( j + 1 ) = jao( j + 1 ) + n1
      END DO

      RETURN

c--------end-of-dmperm--------------------------------------------------

      END

c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 

      SUBROUTINE dvperm( n, x, perm ) 

c-----------------------------------------------------------------------

      INTEGER n, perm( n ) 
      DOUBLE PRECISION x( n )

c-----------------------------------------------------------------------
c
c     this subroutine performs an in-place permutation of a real vector x 
c     according to the permutation array perm(*), i.e., on return, 
c     the vector x satisfies,
c
c	x( perm( j ) ) :== x( j ), j = 1, 2, .., n
c
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> n 	= length of vector x.
c     --> perm 	= integer array of length n containing the permutation  array.
c     --> x	= input vector
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c     
c     --> x = vector x permuted according to x( perm(*) ) :=  x(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c

c     local variables 

      DOUBLE PRECISION tmp, tmp1

c-----------------------------------------------------------------------

      init = 1
      tmp = x( init )
      ii = perm( init )
      perm( init  )= -perm( init )
      k = 0
     
c     loop
 
 6    k = k + 1

c     save the chased element --

      tmp1 = x( ii ) 
      x( ii ) = tmp
      next = perm( ii ) 
      IF ( next .LT. 0 ) GOTO 65
     
c     test for end

      IF ( k .GT. n ) GOTO 101
      tmp = tmp1
      perm( ii ) = - perm( ii )
      ii = next 

c     end loop 

      GOTO 6

c     reinitilaize cycle --

 65   init = init + 1
      IF ( init .GT. n ) GOTO 101
      IF ( perm( init ) .LT. 0 ) GOTO 65
      tmp = x(init )
      ii = perm( init )
      perm( init ) = -perm( init )
      GOTO 6
     
 101  CONTINUE
     
      DO j = 1, n
         perm( j ) = -perm( j )
      END DO
     
      RETURN

c-------------------end-of-dvperm--------------------------------------- 

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE ivperm( n, ix, perm ) 

c-----------------------------------------------------------------------

      INTEGER n, perm( n ), ix( n )

c-----------------------------------------------------------------------
c
c     this subroutine performs an in-place permutation of an integer vector 
c     ix according to the permutation array perm(*), i.e., on return, 
c     the vector x satisfies,
c
c     ix( perm( j ) ) :== ix( j ), j = 1, 2, .., n
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> n = length of vector x.
c
c     --> perm = integer array of length n containing the permutation  array.
c
c     --> ix = input vector
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c 
c     --> ix = vector x permuted according to ix( perm(*) ) :=  ix(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c

c     local variables
      
      INTEGER tmp, tmp1

c-----------------------------------------------------------------------

      init = 1
      tmp = ix( init )
      ii = perm( init )
      perm( init ) = -perm( init )
      k = 0
     
c     loop
 
 6    k = k + 1
      
c     save the chased element --
 
      tmp1 = ix( ii ) 
      ix( ii ) = tmp
      next = perm( ii ) 
      IF ( next .LT. 0 ) GOTO 65
     
c     test for end

      IF ( k .GT. n ) GOTO 101
      tmp = tmp1
      perm( ii ) = - perm( ii )
      ii = next 

c     end loop 

      GOTO 6

c     reinitilaize cycle --

 65   init = init + 1
      IF ( init .GT. n ) GOTO 101
      IF ( perm( init ) .LT. 0)  GOTO 65
      tmp = ix( init )
      ii = perm( init )
      perm( init ) = -perm( init )
      GOTO 6
     
 101  CONTINUE

      DO j = 1, n
         perm( j ) = -perm( j )
      END DO
     
      RETURN

c-------------------end-of-ivperm--------------------------------------- 

      END

c-----------------------------------------------------------------------  
c-----------------------------------------------------------------------  

      SUBROUTINE retmx( n, a, ja, ia, dd )

c-----------------------------------------------------------------------

      DOUBLE PRECISION a(*), dd(*)
      INTEGER n, ia(*), ja(*)

c-----------------------------------------------------------------------
c
c     Returns in dd(*) the max absolute value of elements in row *.
c     Used for scaling purposes. Superseded by rnrms  .
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> n = dimension of A
c
c     --> a,ja,ia = matrix stored in compressed sparse row format
c
c-----------------------------------------------------------------------
c     on output :
c-----------------------------------------------------------------------
c
c     --> dd = DOUBLE PRECISION array of length n. 
c     
c     entry dd(i) contains the element of row i that has the largest
c     absolute value.
c     Moreover the sign of dd is modified such that it is the
c     same as that of the diagonal element in row i.
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c

c     local variables
 
      INTEGER k2, i, k1, k
      DOUBLE PRECISION t, t1, t2

c-----------------------------------------------------------------------

c     initialize 

      k2 = 1

      DO i = 1, n
         
         k1 = k2
         k2 = ia( i + 1 ) - 1
         t = 0.0D0

         DO k = k1, k2

            t1 = ABS( a( k ) )
            IF ( t1 .GT. t ) t = t1

            IF ( ja( k ) .EQ. i ) THEN 
            
               IF ( a( k ) .GE. 0.0D0 ) THEN 
                  t2 = a( k ) 
               ELSE 
                  t2 = - a( k )
               END IF

            END IF

         END DO
      
         dd( i ) =  t2 * t

c     we do not invert diag

      END DO

      RETURN

c---------end of retmx -------------------------------------------------

      END

c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------

      SUBROUTINE diapos( n, ja, ia, idiag ) 

c-----------------------------------------------------------------------

      INTEGER ia( n + 1 ), ja(*), idiag( n ) 

c-----------------------------------------------------------------------
c
c     this subroutine returns the positions of the diagonal elements of a
c     sparse matrix a, ja, ia, in the array idiag.
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c     
c     --> n = integer - row dimension of the matrix a.
c
c     --> a, ja, ia = matrix stored in CSR format; a array skipped.
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c
c     --> idiag  = integer array of length n. 
c
c     The i-th entry of idiag points to the diagonal element a( i, i )
c     in the arrays a, ja.
c     ( i.e., a( idiag( i ) ) = element A( i, i ) of matrix A )
c     
c     If no diagonal element is found the entry is set to 0.
c
c----------------------------------------------------------------------c
c     Y. Saad, March, 1990
c----------------------------------------------------------------------c
      
      DO i = 1, n 
         idiag( i ) = 0
      END DO
     
c     sweep through data structure. 
     
      DO i = 1, n
         DO k = ia( i ), ia( i + 1 ) -1
            IF ( ja( k ) .EQ. i ) idiag( i ) = k
         END DO
      END DO

      RETURN

c----------- -end-of-diapos---------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE dscaldg( n, a, ja, ia, diag, job )

c-----------------------------------------------------------------------
      
      DOUBLE PRECISION a(*), diag(*), t
      INTEGER ia(*), ja(*)

c-----------------------------------------------------------------------
c
c     scales rows by diag where diag is either given ( job = 0 )
c     or to be computed :
c
c     job = 1 : scale row i by  +/- max | a( i, j ) | and put inverse of 
c     scaling factor in diag( i ), where +/- is the sign of a( i, i ).
c
c     job = 2 : scale by 2-norm of each row..
c
c     If diag( i ) = 0, then diag( i ) is replaced by one ( no scaling ).
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
 
      GOTO ( 12, 11, 10 ) job + 1

 10   DO j = 1, n

         k1 = ia( j )
         k2 = ia( j + 1 ) - 1
         t = 0.0D0
         
         DO k = k1, k2
            t = t + a( k ) * a( k )
         END DO
         
         diag( j ) = SQRT( t )
         
      END DO

      GOTO 12
            
 11   CONTINUE
            
      CALL retmx( n, a, ja, ia, diag )
            
 12   DO j = 1, n

         IF ( diag( j ) .NE. 0.0D0 ) THEN 
            diag( j ) = 1.0D0 / diag( j )
         ELSE 
            diag( j ) = 1.0D0
         END IF
         
      END DO
      
      DO i = 1, n

         t = diag( i )

         DO k = ia( i ), ia( i + 1 ) -1
            a( k ) = a( k ) * t
         END DO

      END DO

      RETURN
 
c--------end of dscaldg -----------------------------------------------

      END

c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 

      SUBROUTINE extbdg( n, a, ja, ia, bdiag, nblk, ao, jao, iao )

c-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION ( a-h, o-z )
      DOUBLE PRECISION bdiag(*), a(*), ao(*)
      INTEGER ia(*), ja(*), jao(*), iao(*) 

c-----------------------------------------------------------------------
c     
c     this subroutine extracts the main diagonal blocks of a 
c     matrix stored in compressed sparse row format and puts the result
c     into the array bdiag and the remainder in ao, jao, iao.
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> n = integer. The row dimension of the matrix a.
c
c     --> a, ja, ia = matrix stored in csr format
c     
c     --> nblk = dimension of each diagonal block.
c
c     The diagonal blocks are stored in CSR format
c     i.e., we store in succession the i nonzeros of the i-th row
c     after those of row number i - 1 ..
c
c-----------------------------------------------------------------------
c     on return:
c----------------------------------------------------------------------- 
c
c     --> bdiag = DOUBLE PRECISION array of size ( n x nblk ) containing the diagonal
c     blocks of A on return
c
c     --> ao, jao, iao = remainder of the matrix stored in CSR format.
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c

      m = 1 + ( n - 1 ) / nblk

c     this version is sequential -- there is a more parallel version
c     that goes through the structure twice ....

      ltr =  ( ( nblk - 1 ) * nblk ) / 2 
      l = m * ltr

      DO i = 1, l
         bdiag( i ) = 0.0D0
      END DO

      ko = 0
      kb = 1
      iao( 1 ) = 1

      DO jj = 1, m

         j1 = ( jj - 1 ) * nblk + 1
         j2 =  min0( n, j1 + nblk - 1 )

         DO j = j1, j2

            DO i = ia( j ), ia( j + 1 ) -1
         
               k = ja( i )

               IF ( k .LT. j1 ) THEN

                  ko = ko + 1
                  ao( ko ) = a( i )
                  jao( ko ) = k

               ELSE IF (k .LT. j) THEN
                  
c     kb = (jj-1)*ltr+((j-j1)*(j-j1-1))/2+k-j1+1
c     bdiag(kb) = a(i)

                  bdiag( kb + k- j1 ) = a( i )
                  
               END IF
               
            END DO

            kb = kb + j - j1
            iao( j + 1 ) = ko + 1
            
         END DO
         
      END DO
      
      RETURN
      
c---------end-of-extbdg------------------------------------------------- 

      END

c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 

      SUBROUTINE getbwd( n, a, ja, ia, ml, mu )

c-----------------------------------------------------------------------
c
c     gets the bandwidth of lower part and upper part of A.
c     does not assume that A is sorted.
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> n = integer - the row dimension of the matrix
c     --> a, ja, ia = matrix in CSR format.
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c
c     --> ml = integer. The bandwidth of the strict lower part of A
c     --> mu = integer. The bandwidth of the strict upper part of A 
c
c-----------------------------------------------------------------------
c
c     Notes:
c
c     ml and mu are allowed to be negative on return.
c     This may be useful since it will tell us whether a band is confined 
c     in the strict  upper/lower triangular part. 
c     Indeed, the definitions of ml and mu are
c
c     ml = max ( ( i - j ) s.t. a( i, j ) .NE. 0  )
c     mu = max ( ( j - i ) s.t. a( i, j ) .NE. 0  )
c
c----------------------------------------------------------------------c
c     Y. Saad, Sep. 21 1989                                            c
c----------------------------------------------------------------------c
 
      DOUBLE PRECISION a(*) 
      INTEGER ja(*), ia( n + 1 ), ml, mu, ldist, i, k 

c----------------------------------------------------------------------- 

      ml = - n
      mu = - n
      
      DO i = 1, n
         DO k = ia( i ), ia( i + 1 ) - 1 
            ldist = i - ja( k )
            ml = max( ml, ldist )
            mu = max( mu, -ldist )
         END DO
      END DO

      RETURN

c---------------end-of-getbwd ------------------------------------------ 

      END

c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------

      SUBROUTINE blkfnd( nrow, ja, ia, nblk )

c-----------------------------------------------------------------------
c
c     This routine attemptps to determine whether or not the input
c     matrix has a block structure and finds the blocks size if it does.
c     A block matrix is one which is comprised of small square dense blocks.
c     If there are zero elements within the square blocks and the original
c     data structure takes these zeros into account, then blkchk may fail 
c     to find the correct block size. 
c
c----------------------------------------------------------------------- 
c     on entry
c-----------------------------------------------------------------------
c
c     --> nrow = integer equal to the row dimension of the matrix.  
c
c     --> ja   = integer array containing the column indices
c     of the nonzero entries of the matrix stored by row.
c
c     --> ia   = integer array of length nrow + 1 containing the pointers 
c     beginning of each row in array ja.
c		
c     --> nblk = integer containing the assumed value of nblk if job = 0
c
c-----------------------------------------------------------------------
c     on return
c-----------------------------------------------------------------------
c
c     --> nblk  = integer containing the value found for nblk when job = 1.
c     If imsg .NE. 0 this value is meaningless however.
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c

      INTEGER ia( nrow + 1 ), ja(*) 

c-----------------------------------------------------------------------
c
c     First part of code will find candidate block sizes.
c     Criterion used here is a simple one: scan rows and  determine groups 
c     of rows that have the same length and such that the first column 
c     number and the last column number are identical.
c 
c-----------------------------------------------------------------------
      
      minlen = ia( 2 ) - ia( 1 )
      irow   = 1

      DO i = 2, nrow
         len = ia( i + 1 ) - ia( i )
         IF ( len .LT. minlen ) THEN
            minlen = len 
            irow = i
         END IF
      END DO

c---- candidates are all dividers of minlen

      nblk = 1
      IF ( minlen .LE. 1 ) RETURN

      DO 99 iblk = minlen, 1, -1

         IF ( MOD( minlen, iblk ) .NE. 0 ) GOTO 99
         
         len = ia( 2 ) - ia( 1 )
         len0 = len
         jfirst = ja( 1 ) 
         jlast = ja( ia( 2 ) - 1 )
         
         DO jrow = irow + 1, irow + nblk - 1
         
            i1 = ia( jrow )
            i2 = ia( jrow + 1 ) - 1
            len = i2 + 1 - i1
            jf = ja( i1 )
            jl = ja( i2 ) 
            IF ( ( len .NE. len0 ) .OR. ( jf .NE. jfirst ) .OR. 
     &           ( jl .NE. jlast ) ) GOTO 99
            
         END DO
         
c     check for this candidate ----
         
         CALL blkchk( nrow, ja, ia, iblk, imsg )
         IF ( imsg .EQ. 0 ) THEN
            
c     block size found
            
            nblk = iblk
            RETURN
            
         END IF 
         
 99   CONTINUE

c--------end-of-blkfnd ------------------------------------------------- 
      
      END
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE blkchk( nrow, ja, ia, nblk, imsg )

c----------------------------------------------------------------------- 
c     
c     This routine checks whether the input matrix is a block
c     matrix with block size of nblk. A block matrix is one which is
c     comprised of small square dense blocks. If there are zero
c     elements within the square blocks and the data structure
c     takes them into account then blkchk may fail to find the
c     correct block size.
c     
c-----------------------------------------------------------------------
c     on entry
c-----------------------------------------------------------------------
c     
c     --> nrow = integer equal to the row dimension of the matrix.
c     
c     --> ja   = integer array containing the column indices
c     of the nonzero entries of the matrix stored by row.
c     
c     --> ia   = integer array of length nrow + 1 containing the pointers 
c     beginning of each row in array ja.
c     
c     --> nblk = integer containing the value of nblk to be checked. 
c     
c-----------------------------------------------------------------------
c     on return
c-----------------------------------------------------------------------
c     
c     --> imsg = integer containing a message with the following meaning :
c
c     imsg = 0 means that the output value of nblk is a correct block size.
c     imsg < 0 means nblk not correct block size :
c
c     imsg = -1 : nblk does not divide nrow 
c     imsg = -2 : a starting element in a row is at wrong position
c     ( j .NE. mult * nblk + 1 )
c     imsg = -3 : nblk does divide a row length
c     imsg = -4 : an element is isolated outside a block or
c     two rows in same group have different lengths
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
   
      INTEGER ia( nrow + 1 ), ja(*) 

c----------------------------------------------------------------------
c
c     First part of code will find candidate block sizes.
c     This is not guaranteed to work, so a check is done at the end.
c     The criterion used here is a simple one:
c     scan rows and determine groups of rows that have the same length
c     and such that the first column number and the last column number 
c     are identical. 
c
c----------------------------------------------------------------------
      
      imsg = 0
      IF ( nblk .LE. 1 ) RETURN

      nr = nrow / nblk
      IF ( nr * nblk .NE. nrow ) GOTO 101

c--   main loop --------------------------------------------------------- 

      irow = 1

      DO ii = 1, nr

c     i1 = starting position for group of nblk rows in original matrix

         i1 = ia( irow )
         j2 = i1

c     lena = length of each row in that group  in the original matrix

         lena = ia( irow + 1 ) - i1

c     len = length of each block-row in that group in the output matrix

         len = lena / nblk
         IF ( len * nblk .NE. lena ) GOTO 103

c     for each row
     
         DO i = 1, nblk
         
            irow = irow + 1
            IF ( ia( irow ) - ia( irow - 1 ) .NE. lena ) GOTO 104
     
c     for each block

            DO k = 0, len - 1

               jstart = ja( i1 + nblk * k ) - 1
               IF ( ( jstart / nblk ) * nblk .NE. jstart ) GOTO 102

c     for each column
     
               DO j = 1, nblk
                  
                  IF ( jstart + j .NE. ja( j2 ) )  GOTO 104
                  j2 = j2 + 1
               
               END DO

            END DO

         END DO

      END DO

c     went through all loops successfully:

      RETURN
 
 101  imsg = -1
      RETURN
 
 102  imsg = -2
      RETURN
 
 103  imsg = -3
      RETURN
 
 104  imsg = -4
      RETURN

c----------------end of chkblk ----------------------------------------- 

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE infdia( n, ja, ia, ind, idiag ) 

c-----------------------------------------------------------------------

      INTEGER ia(*), ind(*), ja(*)

c-----------------------------------------------------------------------
c
c     obtains information on the diagonals of A. 
c
c-----------------------------------------------------------------------
c
c     this subroutine finds the lengths of each of the 2 * n - 1
c     diagonals of A.
c     it also outputs the number of nonzero diagonals found. 
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> n = dimension of the matrix a.
c
c     --> a (not needed here ), ja, ia = matrix stored in CSR format
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c
c     --> idiag = integer - number of nonzero diagonals found. 
c 
c     --> ind   = integer array of length at least 2 * n - 1.
c
c     The k-th entry in ind contains the number of nonzero elements
c     in the diagonal number k, the numbering being from the lowermost
c     diagonal ( bottom-left ).
c     In other words ind( k ) = length of diagonal whose offset wrt
c     the main diagonal is = - n + k.
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c

      n2 = n + n - 1

      DO i = 1, n2
         ind( i ) = 0
      END DO

      DO i = 1, n
         DO k = ia( i ), ia( i + 1 ) - 1
            j = ja( k )
            ind( n + j - i ) = ind( n + j - i ) + 1
         END DO
      END DO
      
c     count the nonzero ones.

      idiag = 0 
      
      DO k = 1, n2
         IF ( ind( k ) .NE. 0 ) idiag = idiag + 1
      END DO
      
      RETURN

c------end-of-infdia ---------------------------------------------------
      
      END

c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------

      SUBROUTINE amubdg( nrow, ncol, ncolb, ja, ia, jb, ib,
     &     ndegr, nnz, iw ) 

c-----------------------------------------------------------------------

      INTEGER ja(*), jb(*), ia( nrow + 1 ), ib( ncol + 1 )
      INTEGER ndegr( nrow ), iw( ncolb ) 

c-----------------------------------------------------------------------
c
c     gets the number of nonzero elements in each row of A*B and the total 
c     number of nonzero elements in A*B.
c     
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c     
c     --> nrow = integer - row dimension of matrix A
c     --> ncol = integer - column dimension of matrix A
c     ( = row dimension of matrix B ).
c
c     --> ncolb = integer - the colum dimension of the matrix B.
c
c     --> ja, ia = row structure of input matrix A :
c
c     ja = column indices of the nonzero elements of A stored by rows.
c     ia = pointer to beginning of each row  in ja.
c 
c     --> jb, ib = row structure of input matrix B :
c
c     jb = column indices of the nonzero elements of A stored by rows.
c     ib = pointer to beginning of each row  in jb.
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c
c     --> ndegr = integer array of length nrow containing the degrees
c     ( i.e., the number of nonzeros in each row ) of the matrix A * B 
c				
c     --> nnz = total number of nonzero elements found in A * B
c
c-----------------------------------------------------------------------
c     work arrays:
c-----------------------------------------------------------------------
c
c     --> iw = integer work array of length ncolb. 
c
c-----------------------------------------------------------------------

      DO k = 1, ncolb 
         iw( k ) = 0
      END DO
      
      DO k = 1, nrow
         ndegr( k ) = 0
      END DO

c     method used :
c     
c     Transp( A ) * A = sum [over i = 1, nrow]  a( i )^T a(i)
c     where a( i ) = i-th row of A.
c     We must be careful not to add the elements already accounted for.

      DO ii = 1, nrow 

c     for each row of A
     
         ldg = 0 

c     end-of-linked list
     
         last = -1 
         DO j = ia( ii ), ia( ii + 1 ) - 1

c     row number to be added:
     
            jr = ja( j )
            
            DO k = ib( jr ), ib( jr + 1 ) - 1
               
               jc = jb( k ) 

               IF ( iw( jc ) .EQ. 0 ) THEN 

c     add one element to the linked list 
     
                  ldg = ldg + 1
                  iw( jc ) = last 
                  last = jc

               END IF

            END DO

         END DO

         ndegr( ii ) = ldg
     
c     reset iw to zero

         DO k = 1, ldg 
            j = iw( last ) 
            iw( last ) = 0
            last = j

         END DO

      END DO
     
      nnz = 0

      DO ii = 1, nrow 
         nnz = nnz + ndegr( ii )
      END DO
     
      RETURN

c---------------end-of-amubdg ------------------------------------------

      END

c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------

      SUBROUTINE aplbdg( nrow, ncol, ja, ia, jb, ib, ndegr, nnz, iw ) 

c-----------------------------------------------------------------------

      INTEGER ja(*), jb(*), ia( nrow + 1 ), ib( nrow + 1 )
      INTEGER iw( ncol ), ndegr( nrow ) 

c-----------------------------------------------------------------------
c
c     gets the number of nonzero elements in each row of A + B
c     and the total number of nonzero elements in A + B. 
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> nrow = integer - The row dimension of A and B
c     --> ncol = integer - The column dimension of A and B.
c
c     --> a, ja, ia = Matrix A in CSR format.
c 
c     --> b, jb, ib = Matrix B in CSR format.
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c
c     --> ndegr = integer array of length nrow containing the degrees
c     ( i.e., the number of nonzeros in each row ) of the matrix A + B.
c				
c     --> nnz = total number of nonzero elements found in A + B
c
c-----------------------------------------------------------------------
c     work arrays:
c-----------------------------------------------------------------------
c
c     --> iw = integer work array of length equal to ncol. 
c
c-----------------------------------------------------------------------

      DO k = 1, ncol 
         iw( k ) = 0
      END DO
      
      DO k = 1, nrow
         ndegr( k ) = 0
      END DO

      DO ii = 1, nrow 
         
         ldg = 0 
         
c     end-of-linked list
         
         last = -1 

c     row of A
         
         DO j = ia( ii ), ia( ii + 1 ) - 1 
            
            jr = ja( j ) 
            
c     add element to the linked list 
            
            ldg = ldg + 1
            iw( jr ) = last 
            last = jr
            
         END DO
         
c     row of B
         
         DO j = ib( ii ), ib( ii + 1 ) - 1
            
            jc = jb( j )
            
            IF ( iw( jc ) .EQ. 0 ) THEN 
               
c     add one element to the linked list 
               
               ldg = ldg + 1
               iw( jc ) = last 
               last = jc
               
            END IF
            
         END DO
         
c     done with row ii. 

         ndegr( ii ) = ldg

c     reset iw to zero
     
         DO k = 1, ldg 
            j = iw( last ) 
            iw( last ) = 0
            last = j
         END DO

      END DO
      
      nnz = 0

      DO ii = 1, nrow 
         nnz = nnz + ndegr( ii )
      END DO

      RETURN

c----------------end-of-aplbdg -----------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE rnrms( nrow, nrm, a, ja, ia, diag ) 

c-----------------------------------------------------------------------
      
      DOUBLE PRECISION a(*), diag( nrow ), scal 
      INTEGER ja(*), ia( nrow + 1 ) 

c-----------------------------------------------------------------------
c
c gets the norms of each row of A. ( choice of three norms )
c
c-----------------------------------------------------------------------
c
c     on entry:
c
c-----------------------------------------------------------------------
c
c     --> nrow = integer - The row dimension of A
c     
c     --> nrm  = integer - norm indicator.
c
c     nrm = 1  means 1-norm
c     nrm = 2  means the 2-nrm
c     nrm = 0  means max norm
c
c     --> a, ja, ia = Matrix A in CSR format.
c
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c     
c     --> diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------------

      DO ii = 1, nrow
     
         scal = 0.0D0
         k1 = ia( ii )
         k2 = ia( ii + 1 ) - 1

         IF ( nrm .EQ. 0 ) THEN
            
            DO k = k1, k2
               scal = MAX( scal, ABS( a( k ) ) ) 
            END DO
         
         ELSEIF ( nrm .EQ. 1 ) THEN
            
            DO k = k1, k2
               scal = scal + ABS( a( k ) ) 
            END DO

         ELSE
         
            DO k = k1, k2
               scal = scal + a( k )**2
            END DO

         END IF 

         IF ( nrm .EQ. 2 ) scal = SQRT( scal )
 
         diag( ii ) = scal

      END DO

      RETURN

c-------------end-of-rnrms----------------------------------------------

      END 

c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------

      SUBROUTINE cnrms( nrow, nrm, a, ja, ia, diag ) 

c-----------------------------------------------------------------------

      DOUBLE PRECISION a(*), diag( nrow ) 
      INTEGER ja(*), ia( nrow + 1 ) 

c-----------------------------------------------------------------------
c
c     gets the norms of each column of A. ( choice of three norms )
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> nrow	= integer - The row dimension of A
c
c     --> nrm   = integer - norm indicator.
c
c     nrm = 1  means 1-norm
c     nrm = 2  means the 2-nrm
c     nrm = 0  means max norm
c
c     --> a, ja, ia = Matrix A in CSR format.
c 
c-----------------------------------------------------------------------
c     on return:
c-----------------------------------------------------------------------
c
c
c     --> diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------------

      DO k = 1, nrow 
         diag( k ) = 0.0D0
      END DO

      DO ii = 1, nrow
         k1 = ia( ii )
         k2 = ia( ii + 1 ) - 1
         DO k = k1, k2
            j = ja( k )
 
c     update the norm of each column
      
            IF ( nrm .EQ. 0 ) THEN
               diag( j ) = MAX( diag( j ), ABS( a( k ) ) ) 
            ELSEIF ( nrm .EQ. 1 ) THEN
               diag( j ) = diag( j ) + ABS( a( k ) ) 
            ELSE
               diag( j ) = diag( j ) + a( k )**2
            END IF 

         END DO

      END DO

      IF ( nrm .NE. 2 ) RETURN

      DO k = 1, nrow
         diag( k ) = SQRT( diag( k ) )
      END DO

      RETURN

c------------end-of-cnrms-----------------------------------------------

      END 

c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 

      SUBROUTINE roscal( nrow, job, nrm, a, ja, ia, diag, b,
     &     jb, ib, ierr ) 

c----------------------------------------------------------------------- 

      DOUBLE PRECISION a(*), b(*), diag( nrow ) 
      INTEGER nrow, job, nrm, ja(*), jb(*)
      INTEGER ia( nrow + 1 ), ib( nrow + 1 ), ierr 

c-----------------------------------------------------------------------
c
c     scales the rows of A such that their norms are one on return
c     3 choices of norms: 1-norm, 2-norm, max-norm.
c
c-----------------------------------------------------------------------
c     on entry:
c----------------------------------------------------------------------- 
c
c     --> nrow = integer - The row dimension of A
c
c     --> job = integer - job indicator.
c
c     Job = 0 means get array b only
c     job = 1 means get b, and the integer arrays ib, jb.
c
c     --> nrm = integer - norm indicator.
c
c     nrm = 1  means 1-norm
c     nrm = 2  means the 2-nrm
c     nrm = 0  means max norm
c
c     --> a, ja, ia = Matrix A in CSR format.
c 
c----------------------------------------------------------------------- 
c     on return:
c----------------------------------------------------------------------- 
c
c     --> diag = diagonal matrix stored as a vector containing the matrix
c     by which the rows have been scaled, i.e., on return B = Diag*A.
c     
c     --> b, jb, ib = resulting matrix B in CSR format.
c	    
c     --> ierr  = error message.
c
c     ierr = 0     : Normal return.
c     ierr = i > 0 : Row number i is a zero row.
c
c----------------------------------------------------------------------- 
c
c     Notes:
c
c     1)        The column dimension of A is not needed. 
c     2)        algorithm in place ( B can take the place of A ).
c
c----------------------------------------------------------------------- 

      CALL rnrms( nrow, nrm, a, ja, ia, diag )
      
      ierr = 0
      
      DO j = 1, nrow
         IF ( diag( j ) .EQ. 0.0D0 ) THEN
            ierr = j 
            RETURN
         ELSE
            diag( j ) = 1.0D0 / diag( j )
         END IF
      END DO
c     birner
c     call diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      RETURN

c-------end-of-roscal---------------------------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE coscal( nrow, job, nrm, a, ja, ia, diag,
     &     b, jb, ib, ierr ) 

c-----------------------------------------------------------------------
      
      DOUBLE PRECISION a(*), b(*), diag( nrow ) 
      INTEGER nrow, job, ja(*), jb(*)
      INTEGER ia( nrow + 1 ), ib( nrow + 1 ),ierr 

c-----------------------------------------------------------------------
c
c     scales the columns of A such that their norms are one on return.
c     result matrix written on b, or overwritten on A.
c     3 choices of norms: 1-norm, 2-norm, max-norm. in place.
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> nrow	= integer - The row dimension of A
c
c     --> job   = integer - job indicator.
c
c     job = 0 means get array b only
c     job = 1 means get b, and the integer arrays ib, jb.
c
c     --> nrm   = integer - norm indicator.
c
c     nrm = 1 means 1-norm
c     nrm = 2 means the 2-nrm
c     nrm = 0 means max norm
c
c     --> a, ja, ia = Matrix A in CSR format.
c 
c-----------------------------------------------------------------------
c on return:
c-----------------------------------------------------------------------
c
c     --> diag = diagonal matrix stored as a vector containing the matrix
c     by which the columns have been scaled, i.e., on return B = A * Diag
c
c     --> b, jb, ib = resulting matrix B in CSR format.
c
c     --> ierr = error message.
c
c     ierr = 0     : Normal return.
c     ierr = i > 0 : Column number i is a zero row.
c
c-----------------------------------------------------------------------
c
c     Notes:
c
c     1)     The column dimension of A is not needed. 
c     2)     algorithm in place ( B can take the place of A ).
c
c-----------------------------------------------------------------------
 
      CALL cnrms( nrow, nrm, a, ja, ia, diag )
      
      ierr = 0
      
      DO j = 1, nrow
         IF ( diag( j ) .EQ. 0.0D0 ) THEN
            ierr = j 
            RETURN
         ELSE
            diag( j ) = 1.0D0 / diag( j )
         END IF
      END DO
c     birner
c     call amudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      RETURN

c--------end-of-coscal-------------------------------------------------- 

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE addblk( nrowa, ncola, a, ja, ia, ipos, jpos, job,
     &     nrowb, ncolb, b, jb, ib, nrowc, ncolc, c, jc, ic, nzmx,
     &     ierr )

c-----------------------------------------------------------------------

c     implicit none

      INTEGER nrowa, nrowb, nrowc, ncola, ncolb, ncolc, ipos, jpos
      INTEGER nzmx, ierr, job
      INTEGER ja( 1 : * ), ia( 1 : * ), jb( 1 : * ), ib( 1 : * )
      INTEGER jc( 1 : * ), ic( 1 : * )
      DOUBLE PRECISION a( 1 : * ), b( 1 : * ), c( 1 : * )

c-----------------------------------------------------------------------
c
c     This subroutine adds a matrix B into a submatrix of A whose 
c     ( 1, 1 ) element is located in the starting position ( ipos, jpos ). 
c     The resulting matrix is allowed to be larger than A ( and B ), 
c     and the resulting dimensions nrowc, ncolc will be redefined 
c     accordingly upon return.  
c     The input matrices are assumed to be sorted, i.e. in each row
c     the column indices appear in ascending order in the CSR format.
c
c-----------------------------------------------------------------------
c     on entry:
c-----------------------------------------------------------------------
c
c     --> nrowa      = number of rows in A.
c     --> bcola      = number of columns in A.
c     --> a, ja, ia  = Matrix A in CSR format with entries sorted
c     --> nrowb      = number of rows in B.
c     --> ncolb      = number of columns in B.
c     --> b, jb, ib  = Matrix B in CSR format with entries sorted
c     --> nzmax      = integer - length of the arrays c and jc.
c
c     addblk will stop if the number of nonzero elements in the matrix C
c     exceeds nzmax. See ierr.
c
c----------------------------------------------------------------------- 
c on return:
c-----------------------------------------------------------------------
c
c     --> nrowc      = number of rows in C.
c     --> ncolc      = number of columns in C.
c     --> c, jc, ic  = resulting matrix C in CSR format with entries
c                      sorted ascendly in each row. 
c	    
c     --> ierr	     = integer - serving as error message. 
c
c     ierr = 0 means normal return,
c     ierr > 0 means that addblk stopped while computing the i-th row of C
c     with i = ierr, because the number of elements in C exceeds nzmax.
c
c-----------------------------------------------------------------------
c
c     Notes: 
c
c     this will not work if any of the two input matrices is not sorted
c
c-----------------------------------------------------------------------

c     local variables

      LOGICAL values
      INTEGER i, j1, j2, ka, kb, kc, kamax, kbmax
      
c-----------------------------------------------------------------------

      values = ( job .NE. 0 ) 
      ierr = 0
      nrowc = MAX( nrowa, nrowb + ipos - 1 )
      ncolc = MAX( ncola, ncolb + jpos - 1 )
      kc = 1
      kbmax = 0
      ic( 1 ) = kc

      DO i = 1, nrowc

         IF ( i .LE. nrowa ) THEN
            ka = ia( i )
            kamax = ia( i + 1 ) - 1
         ELSE
            ka = ia( nrowa + 1 )
         END IF
         
         IF ( ( i .GE. ipos ) .AND. 
     &        ( ( i - ipos ) .LE. nrowb ) ) THEN
            kb = ib( i - ipos + 1 )
            kbmax = ib( i - ipos + 2 ) - 1 
         ELSE
            kb = ib( nrowb + 1 )
         END IF

c     a do-while type loop -- goes through all the elements in a row.

 20      CONTINUE 
         
         IF ( ka .LE. kamax ) THEN
            j1 = ja( ka )
         ELSE
            j1 = ncolc + 1
         END IF
         
         IF ( kb .LE. kbmax ) THEN 
            j2 = jb( kb ) + jpos - 1
         ELSE 
            j2 = ncolc + 1
         END IF
         
c     if there are more elements to be added.
         
         IF ( ( ka .LE. kamax .OR. kb .LE. kbmax ) .AND.
     &        ( j1 .LE. ncolc .OR. j2 .LE. ncolc ) ) THEN
            
c     three cases
            
            IF ( j1 .EQ. j2 ) THEN
 
               IF ( values ) c( kc ) = a( ka ) + b( kb )
               jc( kc ) = j1
               ka = ka + 1
               kb = kb + 1
               kc = kc + 1

            ELSE IF ( j1 .LT. j2 ) THEN

               jc( kc ) = j1
               IF ( values ) c( kc ) = a( ka )
               ka = ka + 1
               kc = kc + 1

            ELSE IF ( j1 .GT. j2 ) THEN

               jc( kc ) = j2
               IF ( values ) c( kc ) = b( kb )
               kb = kb + 1
               kc = kc + 1

            END IF

            IF ( kc .GT. nzmx ) GOTO 999
            GOTO 20

         END IF

         ic( i + 1 ) = kc

      END DO

      RETURN

 999  ierr = i 

      RETURN

c------------------------end-of-addblk-----------------------------------

      END

c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 
      
      SUBROUTINE get1up( n, ja, ia, ju )

c-----------------------------------------------------------------------

      INTEGER  n, ja(*), ia(*), ju(*)

c----------------------------------------------------------------------
c
c     obtains the first element of each row of the upper triangular part
c     of a matrix. Assumes that the matrix is already sorted.
c
c-----------------------------------------------------------------------
c     parameters
c-----------------------------------------------------------------------
c     input
c-----------------------------------------------------------------------
c
c     --> ja = integer array containing the column indices of aij
c     --> ia = pointer array.
c
c     ia( j ) contains the position of the beginning of row j in ja
c 
c-----------------------------------------------------------------------
c     output 
c-----------------------------------------------------------------------
c
c     --> ju = integer array of length n.
c
c     ju( i ) is the address in ja of the first element of the uper
c     triangular part of of A ( including the diagonal. Thus if row i
c     does have a nonzero diagonal element then ju( i ) will point to it ).
c     This is a more general version of diapos.
c     
c-----------------------------------------------------------------------

c     local variables
      
      INTEGER i, k 

c-----------------------------------------------------------------------
     
      DO 5 i = 1, n

         ju( i ) = 0
         k = ia( i ) 
         
 1       CONTINUE
         
         IF ( ja( k ) .GE. i ) THEN
         
            ju( i ) = k
            GOTO 5
            
         ELSEIF ( k .LT. ia( i + 1 ) - 1 ) THEN
 
c     go try next element in row 
            
            k = k + 1
            GOTO 1
            
         END IF 
         
 5    CONTINUE

      RETURN

c------------------------end-of-get1up---------------------------------

      END 

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      SUBROUTINE xtrows( i1, i2, a, ja, ia, ao, jao, iao, iperm, job )

c----------------------------------------------------------------------

      INTEGER i1, i2, ja(*), ia(*), jao(*), iao(*), iperm(*), job
      DOUBLE PRECISION a(*), ao(*) 

c-----------------------------------------------------------------------
c
c     This subroutine extracts given rows from a matrix in CSR format. 
c     Specifically, rows number iperm( i1 ), iperm( i1 + 1 ), ..., iperm( i2 )
c     are extracted and put in the output matrix ao, jao, iao, in CSR
c     format.  Not in place. 
c
c-----------------------------------------------------------------------
c
c     Youcef Saad -- coded Feb 15, 1992. 
c     
c-----------------------------------------------------------------------
c     on entry:
c----------------------------------------------------------------------
c
c     --> i1,i2 = two integers indicating the rows to be extracted.
c
c     xtrows will extract rows iperm( i1 ), iperm( i1 + 1 ), ..., iperm( i2 )
c     from original matrix and stack them in output matrix 
c     ao, jao, iao in csr format
c
c     --> a, ja, ia = input matrix in csr format
c
c     --> iperm = integer array of length nrow containing
c     the reverse permutation array for the rows.
c
c     row number iperm( j ) in permuted matrix PA used to be row number j
c     in unpermuted matrix : a( i, j ) in the permuted matrix was
c     a( iperm( i ), j ) in the inout matrix.
c
c     --> job = integer indicating the work to be done:
c
c     job .ne. 1 : get structure only of output matrix,,
c     i.e., ignore real values. ( in which case arrays a 
c     and ao are not used nor accessed ).
c     
c     job = 1 : get complete data structure of output matrix. 
c     ( i.e., including arrays ao and iao ).
c
c----------------------------------------------------------------------
c     on return: 
c----------------------------------------------------------------------
c 
c     ao, jao, iao = input matrix in a, ja, ia format
c
c----------------------------------------------------------------------
c
c     note :
c 
c     If job .ne. 1 then the arrays a and ao are not used.
c
c----------------------------------------------------------------------c
c
c     Y. Saad, revised May  2, 1990                              
c
c----------------------------------------------------------------------c

c     Local variables

      LOGICAL values

c----------------------------------------------------------------------c

      values = ( job .EQ. 1 )

c     copying 

      ko = 1
      iao( 1 ) = ko
      
      DO j = i1, i2

c     ii = iperm( j ) is the index of old row to be copied.

         ii = iperm( j )

         DO k = ia( ii ), ia( ii + 1 ) - 1
            jao( ko ) = ja( k )
            IF ( values ) ao( ko ) = a( k )
            ko = ko + 1
         END DO
         
         iao( j - i1 + 2 ) = ko
         
      END DO

      RETURN

c------------------------end-of-xtrows----------------------------------

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE csrkvstr( n, ia, ja, nr, kvstr )

c-----------------------------------------------------------------------

      INTEGER n, ia( n + 1 ), ja( * ), nr, kvstr( * )

c-----------------------------------------------------------------------
c
c     Finds block row partitioning of matrix in CSR format.
c
c-----------------------------------------------------------------------
c     On entry:
c-----------------------------------------------------------------------
c
c     n       = number of matrix scalar rows
c     ia,ja   = input matrix sparsity structure in CSR format
c
c-----------------------------------------------------------------------
c     On return:
c-----------------------------------------------------------------------
c
c     nr      = number of block rows
c     kvstr   = first row number for each block row
c
c-----------------------------------------------------------------------
c
c     Notes:
c
c     Assumes that the matrix is sorted by columns.
c     This routine does not need any workspace.
c
c-----------------------------------------------------------------------

c     local variables
      
      INTEGER i, j, jdiff

c-----------------------------------------------------------------------

      nr = 1
      kvstr( 1 ) = 1

      DO i = 2, n

         jdiff = ia( i + 1 ) - ia( i )

         IF ( jdiff .EQ. ia( i ) - ia( i - 1 ) ) THEN

            DO j = ia( i ), ia( i + 1 ) - 1
             
               IF ( ja( j ) .NE. ja( j - jdiff ) ) THEN
               
                  nr = nr + 1
                  kvstr( nr ) = i
                  GOTO 299
                  
               END IF
               
            END DO
            
 299        CONTINUE
            
         ELSE
            
            nr = nr + 1
            kvstr( nr ) = i
            
         END IF

      END DO

      kvstr( nr + 1 ) = n + 1
      
      RETURN
      
c------------------------end-of-csrkvstr--------------------------------
      
      END
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE csrkvstc( n, ia, ja, nc, kvstc, iwk )

c-----------------------------------------------------------------------

      INTEGER n, ia( n + 1 ), ja( * ), nc, kvstc( * ), iwk( * )

c-----------------------------------------------------------------------
c
c     Finds block column partitioning of matrix in CSR format.
c
c-----------------------------------------------------------------------
c     On entry:
c-----------------------------------------------------------------------
c
c     n       = number of matrix scalar rows
c     ia,ja   = input matrix sparsity structure in CSR format
c
c-----------------------------------------------------------------------
c     On return:
c-----------------------------------------------------------------------
c
c     nc      = number of block columns
c     kvstc   = first column number for each block column
c
c-----------------------------------------------------------------------
c     Work space:
c-----------------------------------------------------------------------
c
c     iwk(*) of size equal to the number of scalar columns plus one.
c        Assumed initialized to 0, and left initialized on return.
c
c-----------------------------------------------------------------------
c
c     Notes:
c
c     Assumes that the matrix is sorted by columns.
c
c-----------------------------------------------------------------------

c     local variables

      INTEGER i, j, k, ncol

c-----------------------------------------------------------------------

c     use ncol to find maximum scalar column number

      ncol = 0

c     mark the beginning position of the blocks in iwk

      DO i = 1, n
      
         IF ( ia( i ) .LT. ia( i + 1 ) ) THEN
      
            j = ja( ia( i ) )
            iwk( j ) = 1
         
            DO k = ia( i ) + 1, ia( i + 1 ) - 1
               j = ja( k )
               IF ( ja( k - 1 ) .NE. j - 1 ) THEN
                  iwk( j ) = 1
                  iwk( ja( k - 1 ) + 1 ) = 1
               END IF
            END DO
            
            iwk( j + 1 ) = 1
            ncol = max0( ncol, j )
         
         END IF
      
      END DO

      nc = 1
      kvstc( 1 ) = 1
      DO i = 2, ncol + 1
         IF ( iwk( i ) .NE. 0 ) THEN
            nc = nc + 1
            kvstc( nc ) = i
            iwk( i ) = 0
         END IF
      END DO
      nc = nc - 1

      RETURN

c------------------------end-of-csrkvstc--------------------------------
      
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE kvstmerge( nr, kvstr, nc, kvstc, n, kvst )

c-----------------------------------------------------------------------

      INTEGER nr, kvstr( nr + 1 ), nc, kvstc( nc + 1 ), n, kvst( * )

c-----------------------------------------------------------------------
c
c     Merges block partitionings, for conformal row/col pattern.
c
c-----------------------------------------------------------------------
c     On entry:
c-----------------------------------------------------------------------
c     
c     nr,nc   = matrix block row and block column dimension
c     kvstr   = first row number for each block row
c     kvstc   = first column number for each block column
c     
c-----------------------------------------------------------------------
c     On return:
c-----------------------------------------------------------------------
c     
c     n       = conformal row/col matrix block dimension
c     kvst    = conformal row/col block partitioning
c     
c-----------------------------------------------------------------------
c     
c     Notes:
c     
c     If matrix is not square, this routine returns without warning.
c     
c-----------------------------------------------------------------------
      
c-----local variables
      
      INTEGER i, j
      
c-----------------------------------------------------------------------
      
      IF ( kvstr( nr + 1 ) .NE. kvstc( nc + 1 ) ) RETURN

      i = 1
      j = 1
      n = 1

  200 IF ( i .GT. nr + 1 ) THEN
         kvst( n ) = kvstc( j )
         j = j + 1
      ELSEIF ( j .GT. nc + 1 ) THEN
         kvst( n ) = kvstr( i )
         i = i + 1
      ELSEIF ( kvstc( j ) .EQ. kvstr( i ) ) THEN
         kvst( n ) = kvstc( j )
         j = j + 1
         i = i + 1
      ELSEIF ( kvstc( j ) .LT. kvstr( i ) ) THEN
         kvst( n ) = kvstc( j )
         j = j + 1
      ELSE
         kvst( n ) = kvstr( i )
         i = i + 1
      END IF

      n = n + 1

      IF ( ( i .LE. nr + 1 ) .OR. ( j .LE. nc + 1 ) ) GOTO 200

      n = n - 2

      RETURN

c------------------------end-of-kvstmerge-------------------------------

      END

c-----------------------------------------------------------------------
