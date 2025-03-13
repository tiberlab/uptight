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
!        Module "sparse_numrec" - (c) Martin Persson - April 2005
!
!=============================================================================
!
! Contains routines for the handling of sparse matrixes, the matrix is stored
! in the row-indexed sparse storage mode as described in Numerical Recipes F77 or C.
! The standard numerical recipes functions are found in the numrec/ catalog
! The Fortran 90 choice for sparse matrixes is NOT used as it is less suited for
! the 'shift_ham' routine and is less memory conservative.
!_____________________________________________________________________________
!
! sprs_matrix        loads matrix from file into the sparse format
! sprs_element       finds the element i,j in the sparse matrix
! sprs_shift     adds a constant value to all diagonal elements of the sparse matrix
! sprs_ax             multiplies sparse matrix (M,Mij) with a vector v to its right
! sprs_screen        outputs matrix to screen row by row
!
!_____________________________________________________________________________


MODULE sparse_numrec

!===========================================================================
  USE mpi_globals
  USE precision
  USE errors
  USE exceptions
  USE sparse_matrix
  USE rcm_module
!===========================================================================

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sprs_matrix, sprs_matrix_mem, sprs_element, sprs_shift
  PUBLIC :: compute_ovr
  PUBLIC :: sprs_ax
  PUBLIC :: sprs_screen, set_sprs_element, test_eigvect, test_eigvect_csr
  PUBLIC :: sprs_scale, premult, mult
  PUBLIC :: reorder

  INTERFACE sprs_shift
     module procedure sprs_shift_ex 
     module procedure sprs_shift_csr
     module procedure sprs_shift_csr_real     
  END INTERFACE

  INTERFACE sprs_ax
     module procedure sprs_ax_ex 
     module procedure sprs_ax_csr
     module procedure sprs_ax_csr_real
  END INTERFACE

  INTERFACE test_eigvect
     module procedure test_eigvect_ex 
     module procedure test_eigvect_csr
     module procedure test_eigvect_csr_real
  END INTERFACE

  INTERFACE sprs_element
     module procedure sprs_element_ex 
     module procedure sprs_element_csr
  END INTERFACE
  
  INTERFACE set_sprs_element
     module procedure set_sprs_element_ex 
     module procedure set_sprs_element_csr
  END INTERFACE

  INTERFACE sprs_matrix
     module procedure sprs_matrix_ex 
     module procedure sprs_matrix_csr
  END INTERFACE
  
  INTERFACE premult
     module procedure zpremultcsr
  END INTERFACE
  
  INTERFACE mult
     module procedure zmultcsr
  END INTERFACE

CONTAINS
!----------------------------------------------------------------------------
! Allocates memory for the sparse martix (M,Mij)
! Reads matrix from file into the sparse matrix format
!
! N size of matrix
! M vector containing matrix element values
! Mij vector containing index information
! file_value file number for matrix elements
! file_index file number for matrix indexes
! n_max number of non-zero elements
!----------------------------------------------------------------------------
  SUBROUTINE sprs_matrix_ex(M, Mij, file_value,file_index,N,i_max,v_max,block_size)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: file_value, file_index, N, v_max, i_max, block_size
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:), POINTER :: Mij
    COMPLEX (dp), DIMENSION(:), POINTER :: M
!-------------------------------------------------
    INTEGER :: row, col, i_read, v_read, index, index_r, err, &
                                             master_i, temp_index_i, temp_index_v
    INTEGER, DIMENSION(block_size) :: temp_i
    COMPLEX (dp), DIMENSION(block_size) :: temp_v
!Allocalte memory for sparce matrix M and Mij
    LOGICAL :: flag

!write (*,*) 'reading sparse matrix'

    IF(ASSOCIATED(M)) THEN
       DEALLOCATE( M, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'Mi' )
    END IF

    ALLOCATE( M(i_max + 2 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'M' )
    
    IF(ASSOCIATED(Mij)) THEN
       DEALLOCATE( Mij, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'Mij' )
    END IF

    ALLOCATE( Mij(i_max + 2 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'Mij' )


    i_read   = 0
    v_read   = 0
    index    = N + 1   ! Actual row index. Initialize index as NROWs + 1
    index_r  = 0       ! row index counter
    master_i = 1
    temp_index_v = block_size + 1  ! initialization that force initial reading
    temp_index_i = block_size + 1  ! same story


! Read row by row from file
    DO WHILE ( master_i .LE. v_max)
       
! CHECKS IF IT IS NECESSARY TO READ ANOTHER BUFFER LINE FROM INDEX FILE
       IF (temp_index_i .GT. block_size) THEN
! CHECKS IF THE NUMBER OF REMAINING INDECES IS GREATER THAN BUFFER SIZE
          IF(i_max-i_read*block_size .GE. block_size) THEN
             READ ( file_index ) temp_i(1 : block_size) !Index
          ELSE 
             READ ( file_index ) temp_i(1 : i_max-i_read*block_size) !Index
          END IF
! UPDATES NUBER OF LINES READ FROM INDEX FILE
          i_read = i_read + 1
          temp_index_i = 1
       END IF

! CHECKS IF IT IS NECESSARY TO READ ANOTHER BUFFER LINE FROM VALUE FILE
       IF (temp_index_v .GT. block_size) THEN
! CHECKS IF THE NUMBER OF REMAINING VALUES IS GREATER THAN BUFFER SIZE
          IF(v_max-v_read*block_size .GE. block_size) THEN
             READ ( file_value ) temp_v(1 : block_size) !Values
          ELSE 
             READ ( file_value ) temp_v(1 : v_max-v_read*block_size) !Values
          END IF
! UPDATES THE NUMBER OF LINES READ FROM FILE
          v_read = v_read + 1
          temp_index_v = 1
       END IF

       IF(temp_index_i .GT. block_size) write (*,*) 'ERROR', temp_index_i,i_max-i_read*block_size
 
! CHECKS IF THE INDEX IS NEGATIVE. THIS IS THEN A ROW INDEX
       IF ( temp_i(temp_index_i) .LT. 0 ) THEN
 
          row = -temp_i(temp_index_i)
! UPDATES ROW index (marks the beginning of each row)
          index_r = index + 1
! INCREMENT index TO READ NEXT col INDEX
          temp_index_i = temp_index_i + 1
! INITIALIZE THE DIAGONAL ELEMENT FOR A POSSIBLE 0 VALUE
          M(row) = 0.0D0
          Mij(row) = index_r   

! CHECKS IF IT IS NECESSARY TO READ ANOTHER BUFFER LINE FROM INDEX FILE
          IF (temp_index_i .GT. block_size) THEN
! CHECKS IF THE NUMBER OF REMAINING INDECES IS GREATER THAN BUFFER SIZE
             IF(i_max-i_read*block_size .GE. block_size) THEN
                READ ( file_index ) temp_i(1 : block_size) !Index
             ELSE 
                READ ( file_index ) temp_i(1 : i_max-i_read*block_size) !Index
             END IF
! UPDATES NUBER OF LINES READ FROM INDEX FILE
             i_read = i_read + 1
             temp_index_i = 1
          END IF

! READS THE col INDEX AND INCREMENT INDEX TO READ NEXT col INDEX
          col = temp_i(temp_index_i)
          if(col.lt.0) call throw_solve_exception(ERR_HAM_ZEROLN) 
          temp_index_i = temp_index_i + 1

       ELSE ! IF IT IS NOT THE BEGINNING OF A NEW ROW, CONTINUE READING col INDECES

          col = temp_i(temp_index_i)
          temp_index_i = temp_index_i + 1

       END IF
!WRITE (*,*) row, col, master_i, v_max
 
! PUT ELEMENTS ON THE MODIFIED CSR SPARSE MATRIX
! THE DIAGONAL ELEMENTS ARE THE FIRST N ELEMENTS OF M, IN row ORDER.
! CHECKS IF THE VALUE IS DIAGONAL
       IF ( row .EQ. col ) THEN
          M(row) = temp_v(temp_index_v)
          Mij(row) = index_r
! increase value counter index ready for next value
          temp_index_v = temp_index_v + 1
! THE DIAGONAL ELEMENT IS FOUND -> FLAG IS SET TO .true.
!flag = .true.
       ELSE  ! OFF DIAGONAL ELEMENTS:

          index = index + 1
          M(index) = temp_v(temp_index_v)
          temp_index_v = temp_index_v + 1        
          Mij(index) = col
          
       END IF
       
       master_i = master_i + 1
    END DO
! WRITE (*,*) row, col, master_i, v_max

! Mark last used index
    Mij(row + 1) = index + 1

  END SUBROUTINE sprs_matrix_ex

!----------------------------------------------------------------------------
! Allocates memory for the sparse martix (M,Mij)
! Reads matrix from file into the sparse matrix format
!
! N size of matrix
! M vector containing matrix element values
! Mi vector containing rowpointers information
! Mj vector containing colind infos
! file_value file number for matrix elements
! file_index file number for matrix indexes
! n_max number of non-zero elements
!----------------------------------------------------------------------------
  SUBROUTINE sprs_matrix_csr(M, Mj, Mi, file_value,file_index,N,i_max,v_max,block_size)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: file_value, file_index, N, v_max, i_max, block_size
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:), POINTER :: Mj       ! COLIND
    INTEGER, DIMENSION(:), POINTER :: Mi       ! ROWPNT
    COMPLEX (dp), DIMENSION(:), POINTER :: M
!-------------------------------------------------
    INTEGER :: row, col, i_read, v_read, index, index_r, err, &
                                             master_i, temp_index_i, temp_index_v
    INTEGER, DIMENSION(block_size) :: temp_i
    COMPLEX (dp), DIMENSION(block_size) :: temp_v
!Allocalte memory for sparce matrix M and Mij
    LOGICAL :: flag

!write (*,*) 'reading sparse matrix'

    IF(ASSOCIATED(M)) THEN
       DEALLOCATE( M, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'M' )
    END IF

    ALLOCATE( M(i_max + 2 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'M' )
    
    IF(ASSOCIATED(Mi)) THEN
       DEALLOCATE( Mi, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'Mi' )
    END IF

    ALLOCATE( Mi(N + 1 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'Mi' )

    IF(ASSOCIATED(Mj)) THEN
       DEALLOCATE( Mj, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'Mj' )
    END IF

    ALLOCATE( Mj(i_max + 2 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'Mj' )


    i_read   = 0
    v_read   = 0
    index    = 0
    index_r  = 0       ! row index counter
    master_i = 1
    temp_index_v = block_size + 1  ! initialization that force initial reading
    temp_index_i = block_size + 1  ! same story


! Read row by row from file
    DO WHILE ( master_i .LE. v_max)
       
! CHECKS IF IT IS NECESSARY TO READ ANOTHER BUFFER LINE FROM INDEX FILE
       IF (temp_index_i .GT. block_size) THEN
! CHECKS IF THE NUMBER OF REMAINING INDECES IS GREATER THAN BUFFER SIZE
          IF(i_max-i_read*block_size .GE. block_size) THEN
             READ ( file_index ) temp_i(1 : block_size) !Index
          ELSE 
             READ ( file_index ) temp_i(1 : i_max-i_read*block_size) !Index
          END IF
! UPDATES NUBER OF LINES READ FROM INDEX FILE
          i_read = i_read + 1
          temp_index_i = 1
       END IF

! CHECKS IF IT IS NECESSARY TO READ ANOTHER BUFFER LINE FROM VALUE FILE
       IF (temp_index_v .GT. block_size) THEN
! CHECKS IF THE NUMBER OF REMAINING VALUES IS GREATER THAN BUFFER SIZE
          IF(v_max-v_read*block_size .GE. block_size) THEN
             READ ( file_value ) temp_v(1 : block_size) !Values
          ELSE 
             READ ( file_value ) temp_v(1 : v_max-v_read*block_size) !Values
          END IF
! UPDATES THE NUMBER OF LINES READ FROM FILE
          v_read = v_read + 1
          temp_index_v = 1
       END IF

       IF(temp_index_i .GT. block_size) write (*,*) 'ERROR', temp_index_i,i_max-i_read*block_size
 
! CHECKS IF THE INDEX IS NEGATIVE. THIS IS THEN A ROW INDEX
       IF ( temp_i(temp_index_i) .LT. 0 ) THEN
 
          row = -temp_i(temp_index_i)
! UPDATES ROW index (marks the beginning of each row)
          index_r = index + 1
! INCREMENT index TO READ NEXT col INDEX
          temp_index_i = temp_index_i + 1
! INITIALIZE THE DIAGONAL ELEMENT FOR A POSSIBLE 0 VALUE
!M(index_r) = 0.0D0
          Mi(row) = index_r   

! CHECKS IF IT IS NECESSARY TO READ ANOTHER BUFFER LINE FROM INDEX FILE
          IF (temp_index_i .GT. block_size) THEN
! CHECKS IF THE NUMBER OF REMAINING INDECES IS GREATER THAN BUFFER SIZE
             IF(i_max-i_read*block_size .GE. block_size) THEN
                READ ( file_index ) temp_i(1 : block_size) !Index
             ELSE 
                READ ( file_index ) temp_i(1 : i_max-i_read*block_size) !Index
             END IF
! UPDATES NUBER OF LINES READ FROM INDEX FILE
             i_read = i_read + 1
             temp_index_i = 1
          END IF

! READS THE col INDEX AND INCREMENT INDEX TO READ NEXT col INDEX
          col = temp_i(temp_index_i)
          if(col.lt.0) call throw_solve_exception(ERR_HAM_ZEROLN) 
          temp_index_i = temp_index_i + 1

       ELSE ! IF IT IS NOT THE BEGINNING OF A NEW ROW, CONTINUE READING col INDECES

          col = temp_i(temp_index_i)
          temp_index_i = temp_index_i + 1

       END IF
!WRITE (*,*) row, col, master_i, v_max
 
! PUT ELEMENTS ON THE MODIFIED CSR SPARSE MATRIX
! THE DIAGONAL ELEMENTS ARE THE FIRST N ELEMENTS OF M, IN row ORDER.
! CHECKS IF THE VALUE IS DIAGONAL
!IF ( row .EQ. col ) THEN
!   M(row) = temp_v(temp_index_v)
!   Mij(row) = index_r
!   ! increase value counter index ready for next value
!   temp_index_v = temp_index_v + 1
!   ! THE DIAGONAL ELEMENT IS FOUND -> FLAG IS SET TO .true.
!   !flag = .true.
!ELSE  ! OFF DIAGONAL ELEMENTS:

          index = index + 1
          M(index) = temp_v(temp_index_v)
          temp_index_v = temp_index_v + 1        
          Mj(index) = col
          
!END IF
       
       master_i = master_i + 1
    END DO
! WRITE (*,*) row, col, master_i, v_max

! Mark last used index
    Mi(row + 1) = index + 1

  END SUBROUTINE sprs_matrix_csr

!*********************************************************************************
!----------------------------------------------------------------------------
! Allocates memory for the sparse martix (M,Mij)
! Reads matrix from memory into the sparse matrix format
!
! N size of matrix
! M vector containing matrix element values
! Mij vector containing index information
! file_value file number for matrix elements
! file_index file number for matrix indexes
! n_max number of non-zero elements
!----------------------------------------------------------------------------
  SUBROUTINE sprs_matrix_mem(N, M, Mij, values, indeces, i_max, v_max)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: N, v_max, i_max
    COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: values
    INTEGER, DIMENSION(:), ALLOCATABLE :: indeces
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:), POINTER :: Mij
    COMPLEX (dp), DIMENSION(:), POINTER :: M
!-------------------------------------------------
    INTEGER :: row, col, i_read, v_read, index, index_r, err, master_i, index_v

!Allocalte memory for sparce matrix M and Mij
    LOGICAL :: flag

!write (*,*) 'reading sparse matrix'

    IF(ASSOCIATED(M)) THEN
       DEALLOCATE( M, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'M' )
    END IF

    ALLOCATE( M(i_max + 2 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'M' )
    
    IF(ASSOCIATED(Mij)) THEN
       DEALLOCATE( Mij, STAT = err )
       IF ( err .NE. 0 ) CALL dealloc_error( 'UPTIGHT', 'main', 'Mij' )
    END IF

    ALLOCATE( Mij(i_max + 2 ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparse_numrec', 'sparse_matrix', 'Mij' )

    i_read   = 0
    v_read   = 0
    index    = N + 1   ! Actual row index. Initialize index as NROWs + 1
    index_r  = 0       ! row index counter
    master_i = 1
    index_v  = 1       ! initialization

! Read row by row from file
    DO WHILE ( master_i .LE. v_max)
       
! CHECKS IF THE INDEX IS NEGATIVE. THIS IS THEN A ROW INDEX
       IF ( indeces(master_i) .LT. 0 ) THEN
 
          row = -indeces(master_i)
! UPDATES ROW index (marks the beginning of each row)
          index_r = index + 1
! INCREMENT index TO READ NEXT col INDEX
          master_i = master_i + 1
! FLAG CONTROLS THE POSSIBILITY THAT THE DIAGONAL ELEMENT IS 0
! FLAG is .false. UNTIL AN INDEX col == row IS FOUND
          flag = .false.

          col = indeces(master_i)
          if(col.lt.0) call throw_solve_exception(ERR_HAM_ZEROLN)
          master_i = master_i + 1          

       ELSE ! IF IT IS NOT THE BEGINNING OF A NEW ROW, CONTINUE READING col INDECES

          col = indeces(master_i)
          master_i = master_i + 1 

       END IF
 
! PUT ELEMENTS ON THE MODIFIED CSR SPARSE MATRIX
! THE DIAGONAL ELEMENTS ARE THE FIRST N ELEMENTS OF M, IN row ORDER.
! CHECKS IF THE VALUE IS DIAGONAL
       IF ( row .EQ. col ) THEN

          M(row) = values(index_v)
          Mij(row) = index_r
! increase value counter index ready for next value
          index_v = index_v + 1
! THE DIAGONAL ELEMENT IS FOUND -> FLAG IS SET TO .true.
          flag = .true.

       ELSE 
! INITIALIZE THE DIAGONAL ELEMENT FOR A POSSIBLE 0 VALUE
          IF ( .NOT.flag ) THEN 
             M(row) = 0.0D0
             Mij(row) = index_r   !
             flag = .true.        ! initialization only at first round
          END IF

! OFF DIAGONAL ELEMENTS:
          index = index + 1
          M(index) = values(index_v)
          index_v = index_v + 1          
          Mij(index) = col
          
       END IF
       
       master_i = master_i + 1
    END DO
! WRITE (*,*) row, col, master_i, v_max

! Mark last used index
    Mij(row + 1) = index + 1

  END SUBROUTINE sprs_matrix_mem

!*********************************************************************************
!*********************************************************************************

! Returns value of element M(row, col) from sparse matrix
!
  COMPLEX (dp) FUNCTION sprs_element_ex(M, Mij,row,col)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: row, col
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:) :: Mij
    COMPLEX (dp), DIMENSION(:) :: M
!-------------------------------------------------
    INTEGER :: index, ibeg, iend, imid

    IF (row .EQ. col) THEN
       sprs_element_ex = M(row)
    ELSE

       index = 0
       ibeg = Mij( row )
       iend = Mij( row + 1 ) - 1

       DO WHILE (iend.GE.ibeg)

          imid = (ibeg + iend) / 2

          IF (Mij(imid).eq.col) THEN
             index = imid
             exit
          END IF

          IF ( Mij(imid) .GT. col) then
             iend = imid - 1
          ELSE
             ibeg = imid + 1
          END IF

       END DO


!DO WHILE ( index .LT. iend )
!   IF ( Mij(index) .EQ. col ) EXIT
!   index=index+1
!END DO

       IF ( index .EQ. 0 ) THEN
          sprs_element_ex = 0.d0
       ELSE
          sprs_element_ex = M(index)
       END IF

    END IF

  END FUNCTION sprs_element_ex

!*********************************************************************************
!*********************************************************************************

! Returns value of element M(row, col) from sparse matrix
!
  COMPLEX (dp) FUNCTION sprs_element_csr(M, colind, rowpnt, row, col)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: row, col
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt    
    COMPLEX (dp), DIMENSION(:) :: M
!-------------------------------------------------
    INTEGER :: index, ibeg, iend, imid


    index = 0
    ibeg = rowpnt( row )
    iend = rowpnt( row + 1 ) - 1

    DO WHILE (iend.GE.ibeg)

       imid = (ibeg + iend) / 2
       
       IF ( colind(imid) .eq. col ) THEN
          index = imid
          exit
       END IF
       
       IF ( colind(imid) .GT. col) then
          iend = imid - 1
       ELSE
          ibeg = imid + 1
       END IF
       
    END DO

    IF ( index .EQ. 0 ) THEN
       sprs_element_csr = 0.d0
    ELSE
       sprs_element_csr = M(index)
    END IF
    
  END FUNCTION sprs_element_csr

!*********************************************************************************
!*********************************************************************************

  SUBROUTINE set_sprs_element_ex(M, Mij, row, col, val)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: row, col
    COMPLEX (dp) :: val
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:) :: Mij
    COMPLEX (dp), DIMENSION(:) :: M
!-------------------------------------------------
    INTEGER :: index, ibeg, iend, imid

    IF (row .EQ. col) THEN
       M(row) = val
    ELSE

       index = 0
       ibeg = Mij(row)
       iend = Mij(row+1) - 1
       
       DO WHILE (iend.GE.ibeg)

          imid = (ibeg + iend) / 2

          IF ( Mij(imid) .eq. col) THEN
             index = imid
             exit
          END IF

          IF ( Mij(imid) .GT. col) then
             iend = imid - 1
          ELSE
             ibeg = imid + 1
          END IF

       END DO

       IF (index .EQ. 0) THEN
          write(*,*) 'Element',row,col,'not found'
       ELSE
          M(index) = val
       END IF

    END IF

  END SUBROUTINE set_sprs_element_ex
!*********************************************************************************
!*********************************************************************************

  SUBROUTINE set_sprs_element_csr(M, colind, rowpnt, row, col, val)

!--------------------------------------------------
! IN data
    INTEGER, INTENT( IN ) :: row, col
    COMPLEX (dp) :: val
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt
    COMPLEX (dp), DIMENSION(:) :: M
!-------------------------------------------------
    INTEGER :: index, ibeg, iend, imid

    index = 0
    ibeg = rowpnt(row)
    iend = rowpnt(row+1) - 1
    
    DO WHILE (iend.GE.ibeg)
       
       imid = (ibeg + iend) / 2
       
       IF (colind(imid) .eq. col) THEN
          index = imid
          exit
       END IF
       
       IF (colind(imid) .GT. col) then
          iend = imid - 1
       ELSE
          ibeg = imid + 1
       END IF
       
    END DO

    IF (index .EQ. 0) THEN
       write(*,*) 'Element',row,col,'not found'
    ELSE
       M(index) = val
    END IF

  END SUBROUTINE set_sprs_element_csr
  !*************************************************************************

  !**************************************************************************
  SUBROUTINE sprs_shift_ex(M, Mij, shift)
    !--------------------------------------------------
    ! IN data
    REAL (dp), INTENT( IN ) :: shift
    !-------------------------------------------------
    ! OUT data
    INTEGER, DIMENSION(:) :: Mij
    COMPLEX (dp), DIMENSION(:) :: M
    !-------------------------------------------------

    M(1:Mij(1)-2) = M(1:Mij(1)-2) + CMPLX(shift, 0.0d0)

  END SUBROUTINE sprs_shift_ex

  !------------------------------------------------------------------------  
  SUBROUTINE sprs_shift_csr(M,colind,rowpnt,shift)
    !--------------------------------------------------
    ! IN data
    REAL (dp), INTENT( IN ) :: shift
    !-------------------------------------------------
    ! OUT data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt    
    COMPLEX (dp), DIMENSION(:) :: M
    !-------------------------------------------------

    integer :: i, k, iend, nrow, str, stp, row_index

    nrow = size(rowpnt)-1

    DO i = shift_init , shift_end
       str=rowpnt(i)
       stp=rowpnt(i+1)-1

       DO k = str, stp
          if (colind(k) .eq. i) then
             M(k) = M(k) + cmplx(shift,0.0D0) 
             exit
          endif
       end do
    end do
    
  END SUBROUTINE sprs_shift_csr

  !------------------------------------------------------------------------
  SUBROUTINE sprs_shift_csr_real(M,colind,rowpnt,shift)
    !--------------------------------------------------
    ! IN data
    REAL (dp), INTENT( IN ) :: shift
    !-------------------------------------------------
    ! OUT data
    REAL(dp), DIMENSION(:) :: M
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt    
    !-------------------------------------------------

    integer :: i, k, iend, nrow, str, stp, row_index

    nrow = size(rowpnt)-1

    DO i = shift_init , shift_end
       str=rowpnt(i)
       stp=rowpnt(i+1)-1

       DO k = str, stp
          if (colind(k) .eq. i) then
             M(k) = M(k) + cmplx(shift,0.0D0) 
             exit
          endif
       end do
    end do
    
  END SUBROUTINE sprs_shift_csr_real

  !------------------------------------------------------------------------
  SUBROUTINE compute_ovr(M,colind,rowpnt, col_ind_low, col_ind_hig)

    INTEGER :: col_ind_low
    INTEGER :: col_ind_hig
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt    
    COMPLEX (dp), DIMENSION(:) :: M
    !-------------------------------------------------

    integer :: i, k, iend, nrow, str, stp, row_index

    nrow = size(rowpnt)-1
    
    col_ind_low = colind(1)
    col_ind_hig = colind(1)
    
    DO i = shift_init , shift_end
       str=rowpnt(i)
       stp=rowpnt(i+1)-1
       DO k = str, stp
          if (colind(k) .LE. col_ind_low) col_ind_low = colind(k)
          if (colind(k) .GE. col_ind_hig) then 
             col_ind_hig = colind(k) 
             row_index = i
          endif
       end do
       
    end do
     
  END SUBROUTINE COMPUTE_OVR

!*********************************************************************************
  SUBROUTINE sprs_scale(M,colind,rowpnt,shift)

!--------------------------------------------------
! IN data
    REAL (dp), INTENT( IN ) :: shift
!-------------------------------------------------
! OUT data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt    
    COMPLEX (dp), DIMENSION(:) :: M
!-------------------------------------------------

    integer :: i, k, iend

    iend = size(rowpnt)-1

    do i = 1, iend
       do  k = rowpnt(i),rowpnt(i+1) -1
          if (colind(k) .eq. i) then
             M(k) = M(k) * cmplx(shift,0.d0) 
             exit
          endif
       end do
    end do
    
  END SUBROUTINE sprs_scale



!---------------------------------------------------------------------------------
! Returns multiplication of matrix (M,Mij) with a vector x to its right where
!---------------------------------------------------------------------------------
  SUBROUTINE sprs_ax_ex(M, Mij, sp_fmt, x, ax)
!--------------------------------------------------
! IN data
    INTEGER, DIMENSION(:) :: Mij
    COMPLEX ( dp ), DIMENSION(:) :: M, x
    CHARACTER(1), INTENT(IN) :: sp_fmt
!--------------------------------------------------
! OUTPUT data
!COMPLEX (dp ), DIMENSION(SIZE(x)) :: sprs_ax
    COMPLEX (dp ), DIMENSION(:) :: ax
!--------------------------------------------------
    INTEGER i,k,str,stp,p, nrow
    COMPLEX(dp) :: tmp1, tmp2


    IF ( Mij(1).NE. SIZE(x)+2) THEN 
       write (*,*) 'mismatched vector and matrix in sprsax',Mij(1)-2,SIZE(x)
       call throw_solve_exception(ERR_GENERAL)
    END IF

    nrow = Mij(1)-2

    SELECT CASE( TRIM(sp_fmt) )

       CASE( 'U', 'L', 'u', 'l' )
!------------------------------------------------!
! UPPER OR LOWER MULTIPLICATIONS
!------------------------------------------------
!Multiplication of the matrix diagonal
          DO i = 1,nrow
             ax(i) = M(i)*x(i)
          END DO

          DO i = 1,nrow
! Row-vector multiplication
             str=Mij(i)
             stp=Mij(i+1)-1

             DO k = str,stp
                
                p = Mij(k)
                tmp2 = M(k)
                ax(i) = ax(i)     + tmp2 * x(p)
                ax(p) = ax(p) + CONJG(tmp2)*x(i)
                
             END DO
          END DO
          
       CASE ( 'F', 'f' )
!------------------------------------------------!
! FULL MATRIX MULTIPLICATIONS
!------------------------------------------------
!Multiplication of the matrix diagonal

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,tmp1,str,stp)
! Matrix*vector: i is the row index over which we parallelize
!$OMP DO
          DO i = 1,nrow
!Multiplication of the matrix diagonal
             ax(i) = M(i)*x(i)
!Row-vector multiplication
             tmp1=(0.d0,0.d0)
             str=Mij(i)
             stp=Mij(i+1)-1
             DO k = str, stp
 
                tmp1     = tmp1    + M(k) * x(Mij(k))

             END DO
             ax(i) = ax(i) + tmp1
          END DO
!$OMP END DO

!$OMP END PARALLEL
          
    END SELECT

       
  END SUBROUTINE sprs_ax_ex
!*********************************************************************************

!---------------------------------------------------------------------------------
! Returns multiplication of matrix (M,Mij) with a vector x to its right where
!---------------------------------------------------------------------------------
  SUBROUTINE sprs_ax_csr(M, colind, rowpnt, sp_fmt, x, ax)
    !--------------------------------------------------
    ! IN data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt
    COMPLEX ( dp ), DIMENSION(:) :: M
    COMPLEX ( dp ), DIMENSION(:) ::  x
    CHARACTER(1), INTENT(IN) :: sp_fmt
    !--------------------------------------------------
    ! OUTPUT data
    COMPLEX (dp ), DIMENSION(:) :: ax
    !--------------------------------------------------
    INTEGER i,k,w,str,stp,p, nrow, nwarps, nlines
    COMPLEX(dp) :: tmp1 !, tmp2, tmp3
    INTEGER :: id
    INTEGER, EXTERNAL :: omp_get_thread_num, omp_get_num_threads

    nrow = SIZE(rowpnt) - 1

    SELECT CASE( TRIM(sp_fmt) )

       CASE( 'U', 'L', 'u', 'l' )
        !------------------------------------------------!
        ! UPPER OR LOWER MULTIPLICATIONS
        !------------------------------------------------
          ax = (0.d0, 0.d0)

          DO i = shift_init , shift_end
             ! Row-vector multiplication
             str=rowpnt(i)
             stp=rowpnt(i+1)-1

             DO k = str,stp
                
                p = colind(k)
                tmp1 = M(k)
                ax(i) = ax(i)     + tmp1 * x(p)
                if (i.ne.p) ax(p) = ax(p) + CONJG(tmp1)*x(i)
                
             END DO
          END DO
          
       CASE ( 'F', 'f' )
          
          !------------------------------------------------!
          ! FULL MATRIX MULTIPLICATIONS
          !------------------------------------------------
          ! Matrix*vector: i is the row index over which we parallelize
          
          !$OMP PARALLEL  SHARED(M,rowpnt,colind,x,ax,shift_init, shift_end) PRIVATE(i,k,tmp1,str,stp)
          !$OMP DO
          DO i = shift_init, shift_end
             str=rowpnt(i)
             stp=rowpnt(i+1)-1
             tmp1 = (0.0d0, 0.0d0)
             DO k = str, stp
                tmp1 =  tmp1  +  M(k) * x(colind(k))
             END DO
             
             ax(i) = tmp1
             
          END DO
          !$OMP END DO
          !$OMP END PARALLEL
          
       END SELECT


     END SUBROUTINE sprs_ax_csr

     !---------------------------------------------------------------------------------
     ! Returns multiplication of matrix (M,Mij) with a vector x to its right where
     !---------------------------------------------------------------------------------
     SUBROUTINE sprs_ax_csr_real(M, colind, rowpnt, sp_fmt, x, ax)
       !--------------------------------------------------
       ! IN data
       INTEGER, DIMENSION(:) :: colind
       INTEGER, DIMENSION(:) :: rowpnt
       REAL ( dp ), DIMENSION(:) :: M
       REAL ( dp ), DIMENSION(:) ::  x
       CHARACTER(1), INTENT(IN) :: sp_fmt
       !--------------------------------------------------
       ! OUTPUT data
       REAL (dp ), DIMENSION(:) :: ax
       !--------------------------------------------------
       INTEGER i,k,str,stp,p, nrow
       REAL(dp) :: tmp1, tmp2
       
       
       !IF ( SIZE(rowpnt)-1 .NE. SIZE(x) ) THEN
       !   write (*,*) 'mismatched vector and matrix in sprsax',&
       !        SIZE(rowpnt)-1,SIZE(x)
       !   call throw_solve_exception(ERR_GENERAL)
       !END IF
       
       nrow = SIZE(rowpnt) - 1
       
       SELECT CASE( TRIM(sp_fmt) )
          
       CASE( 'U', 'L', 'u', 'l' )
          !------------------------------------------------!
          ! UPPER OR LOWER MULTIPLICATIONS
          !------------------------------------------------
          ax = (0.d0, 0.d0)
          
          DO i = 1,nrow
             ! Row-vector multiplication
             str=rowpnt(i)
             stp=rowpnt(i+1)-1
             
             DO k = str,stp
                
                p = colind(k)
                tmp2 = M(k)
                ax(i) = ax(i)     + tmp2 * x(p)
                if (i.ne.p) ax(p) = ax(p) + (tmp2)*x(i)
                
             END DO
          END DO
          
       CASE ( 'F', 'f' )
          !------------------------------------------------!
          ! FULL MATRIX MULTIPLICATIONS
          !------------------------------------------------
          !Multiplication of the matrix diagonal
          
          !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,tmp1,str,stp)
          
          !$OMP DO
          DO i = 1,nrow
             ax(i) = (0.d0, 0.d0)
          END DO
          !$OMP END DO NOWAIT
          ! Matrix*vector: i is the row index over which we parallelize
          !$OMP DO
          DO i = 1,nrow
             ! Row-vector multiplication
             tmp1=(0.d0,0.d0)
             str=rowpnt(i)
             stp=rowpnt(i+1)-1
             DO k = str, stp
                tmp1     = tmp1    + M(k) * x(colind(k))
             END DO
             ax(i) = ax(i) + tmp1
          END DO
          !$OMP END DO
          
          !$OMP END PARALLEL
          
       END SELECT
       
       
     END SUBROUTINE sprs_ax_csr_real

     !---------------------------------------------------------------------------------
     ! Outputs sparse matrix row by row to the screen first element is diagonal element
     !---------------------------------------------------------------------------------
     SUBROUTINE sprs_screen(M, Mj, Mi)
       !--------------------------------------------------
       ! IN data
       INTEGER, DIMENSION(:) :: Mi
       INTEGER, DIMENSION(:) :: Mj
       COMPLEX ( dp ), DIMENSION(:) :: M
       !-------------------------------------------------
       INTEGER i
       DO WHILE (i .LT. SIZE(Mi)-1 )
          WRITE (*,*) ,'[',i,Mj(Mi(i):Mi(i+1)-1),']'
          WRITE (*,*) ,' ',M(Mi(i):Mi(i+1)-1)
          WRITE (*,*)
       END DO
     END SUBROUTINE sprs_screen
     
     
     !***************************************************************************
     !***************************************************************************
  
!---------------------------------------------------------------------------
! Multiplies a vector x with the Hamiltonian (M,Mij) and returns the norm of
! the new vetor divided by the norm of the old vector
!---------------------------------------------------------------------------
  REAL ( dp ) FUNCTION test_eigvect_ex(M, Mij, sp_fmt, x)
    
!--------------------------------------------------
! IN data
    INTEGER, DIMENSION(:) :: Mij
    COMPLEX ( dp ), DIMENSION(:) :: M, x
    CHARACTER(1) :: sp_fmt

! LOCAL
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE :: y
!COMPLEX ( dp ), DIMENSION(:), POINTER :: yp
    INTEGER :: err

! avoids using the stack
    ALLOCATE( y(size(x)), STAT=err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparce_numrec', 'test_eigenvect', 'y' )

!yp => y

    call sprs_ax_ex(M, Mij, sp_fmt, x, y)
    test_eigvect_ex = DOT_PRODUCT(x,y)/DOT_PRODUCT(x,x)

    DEALLOCATE(y)

  END FUNCTION test_eigvect_ex

  !---------------------------------------------------------------------------
  ! Multiplies a vector x with the Hamiltonian (M,Mij) and returns the norm of
  ! the new vetor divided by the norm of the old vector
  !---------------------------------------------------------------------------
  REAL ( dp ) FUNCTION test_eigvect_csr(M, colind, rowpnt, sp_fmt, x)
    
    !--------------------------------------------------
    ! IN data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt
    COMPLEX ( dp ), DIMENSION(:) :: M
    COMPLEX ( dp ), DIMENSION(:) :: x
    CHARACTER(1) :: sp_fmt
    
    ! LOCAL
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE :: y
    !    COMPLEX ( dp ), DIMENSION(:), POINTER :: yp
    INTEGER :: err, size_local_Mi, j
    real (dp) :: temp1, temp2, temp1ar, temp2ar
    
    ! avoids using the stack
    ALLOCATE( y(size(x)), STAT=err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparce_numrec', 'test_eigenvect', 'y' )
    
    !yp => y
    
    call sprs_ax(M, colind, rowpnt, sp_fmt, x, y)
    temp1= DOT_PRODUCT( x(shift_init:shift_end), y(shift_init:shift_end) )
    temp2= DOT_PRODUCT( x(shift_init:shift_end), x(shift_init:shift_end) )

#ifdef UPT_MPI
    call mpi_allreduce(temp1, temp1ar, 1 , MPI_double_precision, MPI_SUM, upt_comm,ierr)
    call mpi_allreduce(temp2, temp2ar, 1 , MPI_double_precision, MPI_SUM, upt_comm,ierr)
    test_eigvect_csr = temp1ar/temp2ar
#else
    test_eigvect_csr = temp1/temp2
#endif
    
    
    DEALLOCATE(y)
    
  END FUNCTION test_eigvect_csr


!---------------------------------------------------------------------------
! Multiplies a vector x with the Hamiltonian (M,Mij) and returns the norm of
! the new vetor divided by the norm of the old vector
!---------------------------------------------------------------------------
  REAL ( dp ) FUNCTION test_eigvect_csr_real(M, colind, rowpnt, sp_fmt, x)
    
!--------------------------------------------------
! IN data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt
    REAL ( dp ), DIMENSION(:) :: M
    REAL ( dp ), DIMENSION(:) :: x
    CHARACTER(1) :: sp_fmt

! LOCAL
    REAL ( dp ), DIMENSION(:), ALLOCATABLE :: y
!    COMPLEX ( dp ), DIMENSION(:), POINTER :: yp
    INTEGER :: err

! avoids using the stack
    ALLOCATE( y(size(x)), STAT=err )
    IF ( err .NE. 0 ) CALL alloc_error( 'sparce_numrec', 'test_eigenvect', 'y' )

!yp => y

    call sprs_ax_csr_real(M, colind, rowpnt, sp_fmt, x, y)
    test_eigvect_csr_real = DOT_PRODUCT(x,y)/DOT_PRODUCT(x,x)

    DEALLOCATE(y)

  END FUNCTION test_eigvect_csr_real


!*****************************************************************
!  ZMULT CSR COMPUTE SPACE NEEDED
!Input:
!A_csr: primo fattore in formato CSR
!B_csr: secondo fattore in formato  CSR
!B_ncol: numero di colonne in A_csr e B_csr
!
!
!
!****************************************************************
  FUNCTION zpremultcsr(A_csr,B_csr) RESULT(nnz)
    type(CSR) :: A_csr,B_csr,C_csr
    integer, DIMENSION(:), ALLOCATABLE :: iw 
    integer, DIMENSION(:), ALLOCATABLE :: colind 
    integer, DIMENSION(:), ALLOCATABLE :: rowpnt 
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval 
    integer :: ierr,B_ncol,nnz,k

    IF (A_csr%ncol.NE.B_csr%nrow) THEN
       WRITE(*,*) 'ERROR (zmult_csr): matrices don''t match';
       WRITE(*,*) 'A%ncol=',A_csr%ncol,'B%nrow=',B_csr%nrow
       call throw_solve_exception(ERR_GENERAL)       
    ENDIF

    IF ((A_csr%nnz.EQ.0).OR.(B_csr%nnz.EQ.0)) THEN
       nnz = 0
    ELSE

     DO k = 2, 6  

       B_ncol=B_csr%ncol
!Alloca le parti di C_csr di interesse
       allocate(nzval(1))
       nnz = k * ( A_csr%Mi(A_csr%nrow+1)-1 + B_csr%Mi(A_csr%nrow+1)-1)
       allocate(colind(nnz), stat=ierr)
       allocate(rowpnt(A_csr%nrow+1), stat=ierr)
       allocate(iw(B_ncol),stat=ierr)

         if (ierr.ne.0) then
            write(*,*) 'ERROR (zmultcsr) ALLOCATION ERROR'
            call throw_solve_exception(ERR_ALLOC_ERR)
         end if  
       rowpnt=0

       call zamub(A_csr%nrow,B_ncol,0,A_csr%M,A_csr%Mj,A_csr%Mi,&
            B_csr%M,B_csr%Mj,B_csr%Mi, nzval,colind,rowpnt,&
            nnz,iw,ierr)

       nnz=rowpnt(A_csr%nrow+1)-1
         deallocate(nzval,rowpnt,colind,iw)

       if (ierr.eq.0) exit 
         if (ierr.ne.0 .and. k.eq.6) then
            write(*,*) 'ERROR in premult: maximum size exceeded'
            call throw_solve_exception(ERR_GENERAL)     
         endif

        END DO

    END IF

  END FUNCTION zpremultcsr



!*****************************************************************
!  ZMULT CSR  MULTIPLY 2 CSR MATRICES : OUTPUT ALREADY ALLOCATED
!Input:
!A_csr: primo fattore in formato CSR
!B_csr: secondo fattore in formato  CSR
!B_ncol: numero di colonne in A_csr e B_csr
!C_csr: risultato in formato CSR (l'allocazione esatta viene eseguita
!nella subroutine
!
!****************************************************************
  SUBROUTINE zmultcsr(A_csr,B_csr,C_csr)
    type(CSR) :: A_csr,B_csr,C_csr
    integer, DIMENSION(:), ALLOCATABLE :: iw 
    integer, DIMENSION(:), ALLOCATABLE :: colind 
    integer, DIMENSION(:), ALLOCATABLE :: rowpnt 
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval 
    integer :: ierr,B_ncol,nnz

    IF (A_csr%ncol.NE.B_csr%nrow) THEN
       WRITE(*,*) 'ERROR (zmult_csr): matrices don''t match';
       WRITE(*,*) 'A%ncol=',A_csr%ncol,'B%nrow=',B_csr%nrow
       call throw_solve_exception(ERR_GENERAL)       
    ENDIF

    allocate(iw(B_ncol),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'ERROR (zmultcsr) ALLOCATION ERROR'
         call throw_solve_exception(ERR_ALLOC_ERR)
    end if       
    rowpnt=0

! Prodotto C_csr=A_csr*B_csr
    call zamub(A_csr%nrow,B_ncol,1,A_csr%M,A_csr%Mj,A_csr%Mi,&
            B_csr%M,B_csr%Mj,B_csr%Mi, C_csr%M, C_csr%Mj, C_csr%Mi,&
            C_csr%nnz, iw, ierr)

    if (ierr.ne.0) then
       write(*,*) 'Error in amub subroutine: exceeding C%nnz dimension'
       call throw_solve_exception(ERR_GENERAL)
    end if
        
    deallocate(iw)  


  END SUBROUTINE zmultcsr


!-----------------------------------------------------------------------
! performs the matrix by matrix product C = A B
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A = row dimension of C
! ncol  = integer. The column dimension of B = column dimension of C
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!!
! b,
! jb,
! ib    =  Matrix B in compressed sparse row format.
!
! nzmax = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr  = integer. serving as error message.
!         ierr = 0 means normal return,
!!        ierr .gt. 0 means that amub stopped while computing the
!!        i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw    = integer work array of length equal to the number of
!         columns in A.
! Note:
!-------
!   The row dimension of B is not needed. However there is no checking
!   on the condition that ncol(A) = nrow(B).
!
!-----------------------------------------------------------------------
  SUBROUTINE zamub (nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr) 
    complex(dp), dimension(:) :: a, b, c 
    integer, dimension(:) :: ja,jb,jc,ib,ic
    integer :: ia(nrow+1),iw(ncol)
    integer :: nrow,ncol,job,nzmax,ierr

    integer :: ii,jj,k,ka,kb,len,jcol,jpos    
    complex(dp) :: scal 
    logical :: values

    values = (job .ne. 0) 
    len = 0
    ic(1) = 1 
    ierr = 0
! initialize array iw.
    do jj=1, ncol
       iw(jj) = 0
    enddo
 
    do ii=1, nrow 
       do ka=ia(ii), ia(ii+1)-1 
          if (values) scal = a(ka)
          jj   = ja(ka)
          do kb=ib(jj),ib(jj+1)-1
             jcol = jb(kb)
             jpos = iw(jcol)
             if (jpos .eq. 0) then
                len = len+1
                if (len .gt. nzmax) then
                   ierr = ii
                   return
                endif
                jc(len) = jcol
                iw(jcol)= len
                if (values) c(len)  = scal*b(kb)
             else
                if (values) c(jpos) = c(jpos) + scal*b(kb)
             endif
          end do   
       end do
       do k=ic(ii), len
           iw(jc(k)) = 0
       end do
       ic(ii+1) = len+1
    enddo

  END SUBROUTINE zamub

!------------------------------------------------------------------------
! In-place transposition routine.
!------------------------------------------------------------------------
! this subroutine transposes a matrix stored in compressed sparse row
! format. the transposition is done in place in that the arrays a,ja,ia
! of the transpose are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!---------
! nrow  = integer. The row dimension of A.
! ncol  = integer. The column dimension of A.
! a     = real array of size nnz (number of nonzero elements in A).
!         containing the nonzero elements
! ja    = integer array of length nnz containing the column positions
!         of the corresponding elements in a.
! ia    = integer of size n+1, where n = max(nrow,ncol). On entry
!         ia(k) contains the position in a,ja of  the beginning of
!         the k-th row.
!
! iwk   = integer work array of same length as ja.
!
! on return:
! ----------
!
! ncol  = actual row dimension of the transpose of the input matrix.
!          Note that this may be .le. the input value for ncol, in
!          case some of the last columns of the input matrix are zero
!          columns. In the case where the actual number of rows found
!          in transp(A) exceeds the input value of ncol, transp will
!         return without completing the transposition. see ierr.
! a,
! ja,
! ia    = contains the transposed matrix in compressed sparse
!          row format. The row dimension of a, ja, ia is now ncol.
!
! ierr  = integer. error message. If the number of rows for the
!         transposed matrix exceeds the input value of ncol,
!          then ierr is  set to that number and transp quits.
!         Otherwise ierr is set to 0 (normal return).
!
!c Note:
!----- 1) If you do not need the transposition to be done in place
!         it is preferrable to use the conversion routine csrcsc
!         (see conversion routines in formats).
!      2) the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use csrcsc
!         if you want them sorted.
!----------------------------------------------------------------------c
  SUBROUTINE ztransp (nrow,ncol,a,ja,ia,iwk,ierr)
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr
      complex(dp) a(*) 

      complex(dp) t, t1
      integer nnz, jcol, k, i, init, j, l, inext
      ierr = 0
      nnz = ia(nrow+1)-1

! determine column dimension
      jcol = 0
      do k=1, nnz
         jcol = max(jcol,ja(k))
      enddo
      if (jcol .gt. ncol) then
         ierr = jcol
         return
      endif
!
! convert to coordinate format. use iwk for row indices.
!
      ncol = jcol
     
      do i=1,nrow
         do k=ia(i),ia(i+1)-1
            iwk(k) = i
         end do 
      end do
! find pointer array for transpose.
      do i=1,ncol+1
         ia(i) = 0
      end do
      do k=1,nnz
         i = ja(k)
         ia(i+1) = ia(i+1)+1
      enddo 
      ia(1) = 1 
!--------------------------
      do i=1,ncol
         ia(i+1) = ia(i) + ia(i+1)
      enddo

! loop for a cycle in chasing process.
      init = 1
      k = 0
      do while (init .le. nnz) 
        t = a(init)
        i = ja(init)
        j = iwk(init)
        iwk(init) = -1
!-----------------------------------------------------------
        do while (k .lt. nnz)
          k = k+1
! current row number is i.  determine  where to go.
          l = ia(i)
! save the chased element.
          t1 = a(l)
          inext = ja(l)
! then occupy its location.
          a(l)  = t
          ja(l) = j
! update pointer information for next element to be put in row i.
          ia(i) = l+1
! determine  next element to be chased
          if (iwk(l) .lt. 0) exit 
          t = t1
          i = inext
          j = iwk(l)
          iwk(l) = -1
       enddo
        
       do while (iwk(init).lt.0) 
          init = init+1
          if (init.gt.nnz) exit
       enddo

      enddo      

      do i=ncol,1,-1 
         ia(i+1) = ia(i)
      enddo
      ia(1) = 1

  END SUBROUTINE ztransp

!======================================================================
! REORDERS A MATRIX USING RCM SORTING
!
!======================================================================
  subroutine reorder(mat, perm)
    type(CSR) :: mat
    
    type(CSR) :: P, Tmp

    integer, dimension(:) :: perm
    !integer, dimension(:), allocatable :: perm
    integer :: i, nrow, ierr
    integer, dimension(:), allocatable :: iw

    nrow=mat%nrow

    !allocate(perm(nrow))

    call genrcm(nrow, mat%nnz, mat%Mi, mat%Mj, perm)

    call create_matrix(P,nrow,nrow,nrow)

    do i=1,nrow
       P%M(i)=1
       P%Mj(i)=perm(i)
       P%Mi(i)=i
    enddo
    P%Mi(nrow+1)=nrow+1

    
    call create_matrix(Tmp,nrow,nrow,mat%nnz)

    allocate(iw(nrow+1))

    call zamub(nrow,nrow,1,P%M,P%Mj,P%Mi,mat%M,mat%Mj,mat%Mi,&
               Tmp%M,Tmp%Mj,Tmp%Mi,mat%nnz,iw,ierr) 

    call ztransp(nrow,nrow,P%M,P%Mj,P%Mi,iw,ierr)

    call zamub(nrow,nrow,1,Tmp%M,Tmp%Mj,Tmp%Mi,&
               P%M,P%Mj,P%Mi,mat%M,mat%Mj,mat%Mi,mat%nnz,iw,ierr) 

    call destroy_matrix(P)
    call destroy_matrix(Tmp)
    !deallocate(perm)
    deallocate(iw)
 
  end subroutine reorder
!----------------------------------------------------------------------


!---------------------------------------------------------------------------
! Multiplies a vector x with the Hamiltonian (M,Mij) and returns the norm of
! the new vetor divided by the norm of the old vector
!---------------------------------------------------------------------------

!!$  REAL ( qp ) FUNCTION test_eigvect_qp(M, Mij, sp_fmt, x)
!!$
!!$    !--------------------------------------------------
!!$    ! IN data
!!$    INTEGER, DIMENSION(:), POINTER :: Mij
!!$    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
!!$    COMPLEX ( qp ), DIMENSION(:), POINTER :: x
!!$    CHARACTER(1) :: sp_fmt
!!$
!!$    ! LOCAL
!!$    COMPLEX ( qp ), DIMENSION(:), ALLOCATABLE, TARGET :: y
!!$    COMPLEX ( qp ), DIMENSION(:), POINTER :: yp
!!$    INTEGER :: err
!!$
!!$    ! avoids using the stack
!!$    ALLOCATE( y(size(x)), STAT=err )
!!$    IF ( err .NE. 0 ) CALL alloc_error( 'sparce_numrec', 'test_eigenvect', 'y' )
!!$
!!$    yp => y
!!$
!!$    call sprs_ax_qp(M, Mij, sp_fmt, x, yp)
!!$    test_eigvect_qp = DOT_PRODUCT(x,yp)/DOT_PRODUCT(x,x)
!!$
!!$    DEALLOCATE(y)
!!$
!!$  END FUNCTION test_eigvect_qp
!**************************************************************************
!*************************************************************************



!**************************************************************************
!**************************************************************************

END MODULE sparse_numrec

  
