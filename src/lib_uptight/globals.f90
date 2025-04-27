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
MODULE globals

  USE precision

  IMPLICIT NONE
  PRIVATE
  
  INTEGER, PUBLIC, PARAMETER :: LST = 500
  INTEGER, PUBLIC, PARAMETER :: MST = 200
  INTEGER, PUBLIC, PARAMETER :: SST = 100

  
  INTEGER, PUBLIC, PARAMETER :: s   = 1
  INTEGER, PUBLIC, PARAMETER :: px  = 2
  INTEGER, PUBLIC, PARAMETER :: py  = 3
  INTEGER, PUBLIC, PARAMETER :: pz  = 4 
  INTEGER, PUBLIC, PARAMETER :: se  = 5
  INTEGER, PUBLIC, PARAMETER :: dxy = 6
  INTEGER, PUBLIC, PARAMETER :: dyz = 7
  INTEGER, PUBLIC, PARAMETER :: dzx = 8
  INTEGER, PUBLIC, PARAMETER :: dx2y2 = 9
  INTEGER, PUBLIC, PARAMETER :: dz2r2 = 10
  
  INTEGER, PUBLIC, PARAMETER :: n_ref_states = 10


  INTEGER, PUBLIC, PARAMETER :: sss  = 1
  INTEGER, PUBLIC, PARAMETER :: sps  = 2
  INTEGER, PUBLIC, PARAMETER :: pss  = 3
  INTEGER, PUBLIC, PARAMETER :: pps  = 4
  INTEGER, PUBLIC, PARAMETER :: ppp  = 5
  INTEGER, PUBLIC, PARAMETER :: seses= 6  !10
  INTEGER, PUBLIC, PARAMETER :: sess = 7  !7
  INTEGER, PUBLIC, PARAMETER :: sses = 8  !6 
  INTEGER, PUBLIC, PARAMETER :: seps = 9  !9
  INTEGER, PUBLIC, PARAMETER :: pses = 10 !8 
  INTEGER, PUBLIC, PARAMETER :: sds  = 11
  INTEGER, PUBLIC, PARAMETER :: dss  = 12
  INTEGER, PUBLIC, PARAMETER :: pds  = 13
  INTEGER, PUBLIC, PARAMETER :: dps  = 14
  INTEGER, PUBLIC, PARAMETER :: pdp  = 15
  INTEGER, PUBLIC, PARAMETER :: dpp  = 16
  INTEGER, PUBLIC, PARAMETER :: seds = 17
  INTEGER, PUBLIC, PARAMETER :: dses = 18
  INTEGER, PUBLIC, PARAMETER :: dds  = 19
  INTEGER, PUBLIC, PARAMETER :: ddp  = 20
  INTEGER, PUBLIC, PARAMETER :: ddd  = 21

  INTEGER, PUBLIC, PARAMETER :: n_ref_couplings = 21

  INTEGER, PUBLIC, PARAMETER :: n_intra = 7 !Tan: intra couplings
  
! --- Magnetic Field Parameters ---
  INTEGER, PUBLIC, PARAMETER :: GAUGE_NONE = 0
  INTEGER, PUBLIC, PARAMETER :: GAUGE_LANDAU_Z = 1
  INTEGER, PUBLIC, PARAMETER :: GAUGE_SYMMETRIC_Z = 2

  LOGICAL, PUBLIC :: use_magnetic_field = .FALSE.
  REAL(dp), PUBLIC, DIMENSION(3) :: magnetic_field_vector = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  INTEGER, PUBLIC :: gauge_choice = GAUGE_NONE
  ! --- End Magnetic Field Parameters ---
END MODULE globals
