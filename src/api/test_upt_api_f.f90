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
program test_upt_api

integer :: handler(4), ready, original_dim, reduced_dim, n_blocks
real(8) :: cut_fraction

write(*,*) 'initializing UPT ...'
call upt_initsession(handler)

write(*,*) 'handler recieved',handler

! Public API smoke test.  This remains disabled so it is independent of an
! optional METIS installation.
call upt_set_coarse_graining(handler, 0, 1, -1.0d0, 1.0d0, 0.03d0)
call upt_get_coarse_graining_info(handler, ready, original_dim, reduced_dim, n_blocks, cut_fraction)
if (ready /= 0 .or. original_dim /= 0 .or. reduced_dim /= 0 .or. n_blocks /= 1) stop 1

write(*,*) 'kill UPT ...'
call upt_destructsession(handler)

write(*,*) '...done'

end program test_upt_api
