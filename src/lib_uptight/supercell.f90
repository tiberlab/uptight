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
MODULE supercell

  use precision
  use vectors

  implicit none
  private

  public :: lattice_vectors
  

  INTEGER, PUBLIC, PARAMETER :: latt_max = 1
  INTEGER, PUBLIC, PARAMETER :: latt_dim = ( 2 * latt_max + 1 )**3
  !===========================================================================
  ! 
  ! Subroutine "lattice_vectors" 
  ! calculate all trial lattice vector combinations in the range of one cell
  ! around the central cell. The vectors are sorted by norm.
  ! 
  !=========================================================================
  
CONTAINS

  SUBROUTINE lattice_vectors(latt_vector)
    
    INTEGER, DIMENSION( 3, latt_dim ) :: latt_vector 

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------
    
    INTEGER :: i_latt, l, m, n
    
    !=========================================================================
    
    i_latt = 0
    
    DO l = -latt_max, latt_max
       
       DO m = -latt_max, latt_max
        
          DO n = -latt_max, latt_max
             
             !----------------------------------------------------------------
             
             i_latt = i_latt +1
             latt_vector( :, i_latt ) = (/ l, m, n /)
             
             !----------------------------------------------------------------
             
          END DO
          
       END DO
     
    END DO
    
    latt_vector = sort_by_norm( latt_vector )
    
    !=========================================================================
    
  END SUBROUTINE lattice_vectors




END MODULE supercell
