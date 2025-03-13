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
MODULE list_types

  use precision
  implicit none
 
  private

  public :: cplxs, ints, dprec, names
  public :: matrix, v_pointer

  ! 
  ! definition of derived types to be used in conjunction with linked lists
  !

  TYPE cplxs
     
     COMPLEX ( dp )             :: nbr
     TYPE (cplxs),    POINTER   :: next

  END TYPE cplxs

  !___________________________________________________________________________

  TYPE ints

     INTEGER, DIMENSION( : ), POINTER :: nbr
     TYPE (ints),             POINTER :: next
     
  END TYPE ints
  
  !___________________________________________________________________________

  TYPE dprec

     REAL ( dp ), DIMENSION( : ), POINTER :: nbr
     TYPE (dprec),                POINTER :: next
     
  END TYPE dprec
  
  !___________________________________________________________________________


  TYPE names
     
     CHARACTER ( LEN = 2 ), DIMENSION( : ), POINTER :: name
     TYPE (names),                          POINTER :: next
     
  END TYPE names

  !___________________________________________________________________________

  TYPE matrix

     COMPLEX ( dp ), DIMENSION( :, : ), POINTER :: mat
     TYPE(matrix) ,                     POINTER :: next
     
  END TYPE matrix

  !___________________________________________________________________________


  TYPE v_pointer

     COMPLEX ( dp ), DIMENSION( : ), POINTER :: p
     TYPE(v_pointer),                POINTER :: next

  END TYPE v_pointer

  !___________________________________________________________________________





end MODULE list_types
