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
!---------------------------------------
!Common sorting algorithms
!---------------------------------------
!qsort_index
! swap_i
!
MODULE sort

use precision
use exceptions
IMPLICIT NONE
PRIVATE

PUBLIC :: qsort_index, osort_index
PUBLIC :: qsort
PUBLIC :: sort_csr 


  INTERFACE partition
     module procedure partitionA     
     module procedure partitionAB     
  END INTERFACE

  INTERFACE swap 
     module procedure swap_i     
     module procedure swap_r     
     module procedure swap_z     
  END INTERFACE

CONTAINS
  SUBROUTINE qsort_index(n, array, indx)
    !IN/OUT
    INTEGER n
    REAL ( dp ), DIMENSION( : ) :: array
    INTEGER , DIMENSION( : ) :: indx

    !LOCAL
    INTEGER i, indxt, ir, itemp , j, k, l, jstack, M
    INTEGER, DIMENSION(50) :: istack
    REAL ( dp ) :: a
    M = 7
    ir = size(array)

    !Initiate index array
    DO j = 1, n
       indx(j) = j
    END DO

    DO 
       IF (ir-l .GT. M) THEN
          DO j = l+1, ir
             indxt = indx(j)
             a = array(indxt)
             DO i = j-1, l, -1
                IF (array(indx(i)) .LE. a) EXIT
                indx(i+1) = indx(i)
             END DO
             indx(i+1) = indxt
          END DO
          if (jstack .EQ. 0) EXIT
          jstack = jstack -1
          ir = istack(jstack)
          jstack = jstack -1
          l = istack(jstack)
       ELSE
          k = (l+ir )/2
          CALL swap_i(k,l+1,indx)
          IF (array(indx(l)) .GT. array(indx(ir))) THEN
             CALL swap_i(l,ir,indx)
          END IF
          IF (array(indx(l+1)) .GT. array(indx(ir))) THEN
             CALL swap_i(l+1,ir,indx)
          END IF
          IF (array(indx(l)) .GT. array(indx(l+1))) THEN
             CALL swap_i(l,l+1,indx)
          END IF
          i = l+1
          j = ir
          indxt = indx(l+1)
          a = array(indxt)
          DO
             DO 
                i = i + 1
                IF (array(indx(i)) .LT. a) EXIT
             END DO
             DO 
                j = j - 1
                IF (array(indx(j)) .GT. a) EXIT
             END DO
             IF (j .LT. i) EXIT
             CALL swap_i(i,j,indx)
          END DO
          indx(l+1) = indx(j)
          indx(j) = indxt
          jstack = jstack + 2
          IF (jstack .GT. 50) THEN
             WRITE (*,*) 'SIZE of istack is too small in qsort_index'
             call throw_init_exception(ERR_GENERAL)
          END IF
          IF (ir-i+1 .GE. j-1) THEN
             istack(jstack) = ir
             istack(jstack-1) = i
             ir = j - 1
          ELSE
             istack(jstack) = j-1
             istack(jstack-1) = l
             l = i
          END IF
       END IF
    END DO
  END SUBROUTINE qsort_index

! ------------------------------------------------------
! Sorting routine in ascending order based 
! S. Lacey and R. Box, "A fast easy sort", Byte, 16, 315 (1991)
! As implemented by Erik Oosterwal in O-sort version 3
! ------------------------------------------------------
SUBROUTINE osort_index(length, array, ind)
  
  real(dp), DIMENSION(:), POINTER :: array
  integer,  DIMENSION(:), POINTER :: ind  
  integer :: length

  real(dp), DIMENSION(length) :: larray
  integer :: i,step, tmp
  real(dp) :: sngPhi, sngFib, tmp_d
  
  do i=1,length
     ind(i) = i
  enddo

  ! local copy of dist vector
  larray(1:length) = array(1:length)

  sngPhi = 0.78
  sngFib = length * sngPhi
  step = nint(sngFib)

  do while (step .gt. 0) 

     do i = 1, length - step

        if ( larray(i).gt.larray(i+step) ) then 

           call swap_r(i,i+step,larray)
           call swap_i(i,i+step,ind)
   
        end if

     enddo

     sngFib = sngFib * sngPhi
     step = nint(sngFib)

  end do

END SUBROUTINE osort_index
!*******************************************************
!*******************************************************


  SUBROUTINE swap_i(i,j,indx)
    INTEGER :: i, j, temp
    INTEGER, DIMENSION( : ) :: indx
    temp = indx(i)
    indx(i) = indx(j)
    indx(j) = temp
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(i,j,array)
    INTEGER :: i, j
    REAL (dp), DIMENSION( : ) :: array
    REAL (dp) temp
    temp = array(i)
    array(i) = array(j)
    array(j) = temp
  END SUBROUTINE swap_r

  SUBROUTINE swap_z(i,j,array)
    INTEGER :: i, j
    COMPLEX(dp), DIMENSION( : ) :: array
    COMPLEX(dp) temp
    temp = array(i)
    array(i) = array(j)
    array(j) = temp
  END SUBROUTINE swap_z

  !*********************************************************************************
  ! n     = the row dimension of the matrix
  ! a     = the matrix A in compressed sparse row format.
  ! ja    = the array of column indices of the elements in array a.
  ! ia    = the array of pointers to the rows. 
  ! iwork = integer work array of length max ( n+1, 2*nnz ) 
  !         where nnz = (ia(n+1)-ia(1))  ) .
  ! values= logical indicating whether or not the real values a(*) must 
  !         also be permuted. if (.not. values) then the array a is not
  !         touched by csort and can be a dummy array. 
  SUBROUTINE sort_csr(n,a,ja,ia,values) 
    logical values
    integer n, ja(*), ia(n+1)
    complex(dp) a(*) 
    
    integer, dimension(:), allocatable ::  iwork 
    complex(dp), dimension(:), allocatable ::  zwork 
    integer :: i, k, ln, m, err
  


    do i=1, n
       ln = ia(i+1) - ia(i)
       allocate(iwork(ln), stat=err) 
       allocate(zwork(ln), stat=err) 
       if (err.ne.0) then 
          call throw_init_exception(ERR_GENERAL)
       endif


       ! copy colind on iwork
       do k = 1, ln
          m = ia(i) + k -1
          if (m.lt.1) cycle
          iwork(k) = ja(m)
          zwork(k) = a(m)
       enddo

       ! sort iwork in place and values 
       call qsort(iwork,zwork)

       ! copy back sorted iwork on colind 
       do k = 1, ln
          m = ia(i) + k -1 
          if (m.lt.1) cycle
          ja(m) = iwork(k)
          a(m) = zwork(k)
       enddo

       deallocate(iwork,zwork)
    enddo
     
  END SUBROUTINE sort_csr
  !*********************************************************************************

  !*********************************************************************************
  ! n     = the row dimension of the matrix
  ! a     = the matrix A in compressed sparse row format.
  ! ja    = the array of column indices of the elements in array a.
  ! ia    = the array of pointers to the rows. 
  ! iwork = integer work array of length max ( n+1, 2*nnz ) 
  !         where nnz = (ia(n+1)-ia(1))  ) .
  ! values= logical indicating whether or not the real values a(*) must 
  !         also be permuted. if (.not. values) then the array a is not
  !         touched by csort and can be a dummy array. 
  SUBROUTINE sort_csr_skit(n,a,ja,ia,values) 
    logical values
    integer n, ja(*), ia(n+1)
    complex(dp) a(*) 
    
    integer, dimension(:), allocatable ::  iwork 
    integer i, k, ko, j, ifirst, nnz, next, irow

    nnz = ia(n+1)-1
    allocate(iwork(2*nnz))
    iwork = 0

    do i=1, n
       do k=ia(i), ia(i+1)-1 
          j = ja(k)+1
          iwork(j) = iwork(j)+1
       enddo
    enddo
    !
    ! compute pointers from lengths. 
    !
    iwork(1) = 1
    do i=1,n
       iwork(i+1) = iwork(i) + iwork(i+1)
    enddo
    ! 
    ! get the positions of the nonzero elements in order of columns.
    !
    ifirst = ia(1) 
    nnz = ia(n+1)-ifirst
    do i=1,n
       do k=ia(i),ia(i+1)-1 
          j = ja(k) 
          next = iwork(j) 
          iwork(nnz+next) = k
          iwork(j) = next+1
       enddo
    enddo
    !
    ! convert to coordinate format
    ! 
    do i=1, n
       do k=ia(i), ia(i+1)-1 
          iwork(k) = i
       enddo
    enddo
    !
    ! loop to find permutation: for each element find the correct 
    ! position in (sorted) arrays a, ja. Record this in iwork. 
    ! 
    do k=1, nnz
       ko = iwork(nnz+k) 
       irow = iwork(ko)
       next = ia(irow)
       !
       ! the current element should go in next position in row. iwork
       ! records this position. 
       ! 
       iwork(ko) = next
       ia(irow)  = next+1
    enddo
    !
    ! perform an in-place permutation of the  arrays.
    ! 
    call iperm(nnz, ja(ifirst), iwork) 
    if (values) call zperm(nnz, a(ifirst), iwork) 
    !
    ! reshift the pointers of the original matrix back.
    ! 
    do i=n,1,-1
       ia(i+1) = ia(i)
    enddo
    ia(1) = ifirst 
    deallocate(iwork)

  END SUBROUTINE sort_csr_skit
  !*********************************************************************************

  subroutine iperm(n, x, perm) 
    integer n, perm(n) 
    integer x(n)
    
    integer tmp, tmp1
    integer j,k,ii,next, init
    
    init      = 1
    tmp= x(init)
    ii        = perm(init)
    perm(init)= -perm(init)
    k         = 0

    do while (k .le. n)

       k = k+1
       
       tmp1   = x(ii) 
       x(ii)  = tmp
       next   = perm(ii) 

       if (next .lt. 0 ) then

          do while (perm(init) .lt. 0)    
             init = init+1
             if (init .gt. n) goto 101 
          end do
          tmp= x(init)
          ii= perm(init)
          perm(init)=-perm(init)
          cycle

       endif

       if (k.gt.n) goto 101
       
       tmp       = tmp1
       perm(ii)  = - perm(ii)
       ii        = next 
         
    end do

101 continue    
    do j=1, n
       perm(j) = -perm(j)
    enddo
    
  end subroutine iperm
  !*********************************************************************************

  subroutine zperm(n, x, perm) 
    integer n, perm(n) 
    complex(dp) x(n)
    
    complex(dp) tmp, tmp1
    integer j,k,ii,next, init
    
    init = 1
    tmp = x(init)
    ii = perm(init)
    perm(init)= -perm(init)
    k = 0

    do while (k .le. n)

       k = k+1
       
       tmp1   = x(ii) 
       x(ii)  = tmp
       next   = perm(ii) 

       if (next .lt. 0 ) then

          do while (perm(init) .lt. 0)    
             init = init+1
             if (init .gt. n) goto 101 
          end do
          tmp= x(init)
          ii= perm(init)
          perm(init)=-perm(init)
          cycle

       endif

       if (k.gt.n) goto 101
       
       tmp       = tmp1
       perm(ii)  = - perm(ii)
       ii        = next 
         
    end do

101 continue    
    do j=1, n
       perm(j) = -perm(j)
    enddo
    
  end subroutine zperm
 
 RECURSIVE SUBROUTINE Qsort(a,b)
    INTEGER, INTENT(IN OUT) :: a(:)
    COMPLEX(dp), INTENT(IN OUT), OPTIONAL :: b(:)
    INTEGER :: split
  
    IF (present(b)) THEN
      IF(size(a) > 1) THEN
         CALL PartitionAB(a,b, split)
         CALL Qsort(a(:split-1),b(:split-1))
         CALL Qsort(a(split:),b(split:))
      END IF
    ELSE
      IF(size(a) > 1) THEN
         CALL PartitionA(a, split)
         CALL Qsort(a(:split-1))
         CALL Qsort(a(split:))
      END IF
    END IF

 END SUBROUTINE Qsort
   
 SUBROUTINE PartitionA(a, marker)
    INTEGER, INTENT(IN OUT) :: a(:)
    INTEGER, INTENT(OUT) :: marker
    INTEGER :: left, right, pivot, temp
   
    pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic 
    left = 0                         ! behavior with sorted or reverse sorted data
    right = size(a) + 1
   
    DO WHILE (left < right)
       right = right - 1
       DO WHILE (a(right) > pivot)
          right = right-1
       END DO
       left = left + 1
       DO WHILE (a(left) < pivot)
          left = left + 1
       END DO
       IF (left < right) THEN 
          temp = a(left)
          a(left) = a(right)
          a(right) = temp
       END IF
    END DO
   
    IF (left == right) THEN
       marker = left + 1
    ELSE
       marker = left
    END IF
   
 END SUBROUTINE PartitionA
 
 SUBROUTINE PartitionAB(a, b, marker)
    INTEGER, INTENT(IN OUT) :: a(:)
    COMPLEX(dp), INTENT(IN OUT) :: b(:)
    INTEGER, INTENT(OUT) :: marker
    INTEGER :: left, right, pivot, temp
    COMPLEX(dp) :: tmpz

    pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic 
    left = 0                         ! behavior with sorted or reverse sorted data
    right = size(a) + 1
   
    DO WHILE (left < right)
       right = right - 1
       DO WHILE (a(right) > pivot)
          right = right-1
       END DO
       left = left + 1
       DO WHILE (a(left) < pivot)
          left = left + 1
       END DO
       IF (left < right) THEN 
          temp = a(left)
          a(left) = a(right)
          a(right) = temp
          tmpz = b(left)
          b(left) = b(right)
          b(right) = tmpz
       END IF
    END DO
   
    IF (left == right) THEN
       marker = left + 1
    ELSE
       marker = left
    END IF
   
 END SUBROUTINE PartitionAB

END MODULE sort


