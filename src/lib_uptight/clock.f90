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
module clock

  implicit none
  private
  
 integer, PARAMETER :: NMAXCLKS=5
 integer(4) :: t1(NMAXCLKS),t2(NMAXCLKS), cr, cm
 integer :: nclks=0
 real(4) :: rt1(NMAXCLKS), rt2(NMAXCLKS), ntims=0
 real(4), PARAMETER :: maxtime=3.402823E38

 public :: message_clock, set_clock, write_clock
 public :: set_time, write_time, get_clock, get_sclock
  

contains

 subroutine message_clock(message)

   character(*) :: message
   integer :: l_mess
   character(2) :: str_mess, str_dots

   l_mess=len(message)
   write(str_mess,'(I2)') l_mess
   write(str_dots,'(I2)') 54-l_mess

   write(6,FMT='(A'//str_mess//','//str_dots//'("."))',ADVANCE='NO') message 
   
   call flush(6)

   call set_clock

 end subroutine message_clock

 subroutine set_clock

   if (nclks.lt.NMAXCLKS) then
      nclks = nclks + 1
      call SYSTEM_CLOCK(t1(nclks),cr,cm)
   endif

 end subroutine set_clock

 function get_sclock() result(tm)
   real(4) :: tm

   tm = 0.0
   if (nclks.gt.0) then
      call SYSTEM_CLOCK(t2(nclks),cr,cm)
      
      if (t2(nclks)-t1(nclks).gt.0) then
         tm = (t2(nclks)-t1(nclks))*1.0/cr
      else
         tm = (t2(nclks)-t1(nclks)+cm)*1.0/cr
      endif
      nclks=nclks-1
   else
      write(*,*) 'no clocks'
   end if
 
 end function get_sclock

 subroutine get_clock(tm)
   real(4) :: tm

   tm = 0.0

   if (nclks.gt.0) then
      call SYSTEM_CLOCK(t2(nclks),cr,cm)      
      if (t2(nclks)-t1(nclks).gt.0) then
         tm = (t2(nclks)-t1(nclks))*1.0/cr
      else
         tm = (t2(nclks)-t1(nclks)+cm)*1.0/cr
      endif
      nclks=nclks-1
   else
      write(*,*) 'no clocks'
   endif
   
 end subroutine get_clock


 subroutine write_clock(message)
   character(*), optional :: message
  
   if (present(message)) then
      write(6,FMT='(A)', ADVANCE='NO') message 
   endif   

   if (nclks.gt.0) then
      call SYSTEM_CLOCK(t2(nclks),cr,cm)
      if (t2(nclks)-t1(nclks).gt.0) then
         write(6,*) (t2(nclks)-t1(nclks))*1.d0/cr,"sec"
      else
         write(6,*) (t2(nclks)-t1(nclks)+cm)*1.d0/cr,"sec"
      endif
      nclks=nclks-1
   else 
      write(*,*) 'no clocks'
   endif

 end subroutine write_clock

 subroutine set_time

   if (nclks.lt.NMAXCLKS) then
      ntims = ntims + 1
      call cpu_time(rt1(nclks))
   endif

 end subroutine set_time


 subroutine write_time

   if (nclks.gt.0) then
      call cpu_time(rt2(nclks))
      if (rt2(nclks)-rt1(nclks).gt.0) then
         write(6,*) rt2(nclks)-rt1(nclks),"sec"
      else
         write(6,*) rt2(nclks)-rt1(nclks)+maxtime,"sec"
      endif
      ntims=ntims-1
   endif

 end subroutine write_time


end module clock
