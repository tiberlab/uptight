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
module constants
  
  use precision

  IMPLICIT NONE
  PRIVATE

  real(kind=dp), parameter,public :: eovh = (1.05420882d-3)   ! A/H
  real(kind=dp), parameter,public :: pi =  3.14159265358979323844_dp ! Greek p real
  real(kind=dp), parameter,public :: HAR = 27.2113845_dp         ! H/eV
  real(kind=dp), parameter,public :: ATU = 0.529177249_dp        ! a.u./Ang
  real(kind=dp), PARAMETER,public :: Kb = (3.166830814d-6)    ! H/K
  
  COMPLEX(kind=dp), PARAMETER,public :: j = (0.d0,1.d0)  ! CMPX unity
  INTEGER, PUBLIC, PARAMETER :: max_ind_D = 100000
  
end module constants

