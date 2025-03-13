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
!!* Contains the type definitions and constants needed by the API routines.
module uptightAPICommon

  use upt_param
  implicit none
  private

  public :: DAC_handlerSize, UPTPointers

  !!* Contains a pointer to a TUPTIn and an OUPT instance
  type UPTPointers
     type(OUPT), pointer :: pUPT
  end type UPTPointers
  
  integer, parameter :: DAC_handlerSize = 4

  
end module uptightAPICommon



