
module allocation

  use precision

  implicit none
  private
  
  type mem_log
     integer :: iolog
     integer(8) :: alloc_mem, peak_mem
  end type mem_log

  public log_allocate, log_deallocate, log_allocatep, log_deallocatep
  public writeMemInfo, writePeakInfo
  public resetMemLog, openMemLog, writeMemLog, closeMemLog, mem_log

  interface log_allocatep
     module procedure allocate_pc
     module procedure allocate_pd, allocate_pi, allocate_pz     
     module procedure allocate_pd2, allocate_pi2, allocate_pz2
  end interface

  interface log_deallocatep
     module procedure deallocate_pc
     module procedure deallocate_pd, deallocate_pi, deallocate_pz
     module procedure deallocate_pd2, deallocate_pi2, deallocate_pz2
  end interface


  interface log_allocate
     module procedure allocate_l, allocate_c
     module procedure allocate_d, allocate_i, allocate_z
     module procedure allocate_d2, allocate_i2, allocate_z2
     module procedure allocate_d3, allocate_i3, allocate_z3
     module procedure allocate_d4
  end interface

  interface log_deallocate
     module procedure deallocate_l, deallocate_c
     module procedure deallocate_d, deallocate_i, deallocate_z
     module procedure deallocate_d2, deallocate_i2, deallocate_z2
     module procedure deallocate_d3, deallocate_i3, deallocate_z3
     module procedure deallocate_d4
  end interface

  !---------------------------------------------------------------
  !---------------------------------------------------------------
contains
  !---------------------------------------------------------------
  subroutine allocate_pc(array,length,mlog)
    character(*), DIMENSION(:), POINTER :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog

    integer :: ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*len(array(1))    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pc


  !---------------------------------------------------------------
  subroutine allocate_pi(array,length,mlog)
    integer(4), DIMENSION(:), POINTER :: array
    integer(4), intent(IN)  :: length
    type(mem_log), intent(INOUT)  :: mlog

    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*4    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pi
  !---------------------------------------------------------------

  subroutine allocate_pd(array,length,mlog)
    real(dp), DIMENSION(:), POINTER :: array
    integer(4), intent(IN)  :: length
    type(mem_log), intent(INOUT)  :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pd
  !---------------------------------------------------------------
  subroutine allocate_pz(array,length,mlog)
    complex(dp), DIMENSION(:), POINTER :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog    
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*2*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pz
  !--------------------------------------------------------------- 
  subroutine allocate_pi2(array,row,col,mlog)
    integer(4), DIMENSION(:,:), POINTER :: array
    integer(4), intent(IN) :: row,col
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*4    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pi2
  !--------------------------------------------------------------- 
  subroutine allocate_pd2(array,row,col,mlog)
    real(dp), DIMENSION(:,:), POINTER :: array
    integer(4), intent(IN) :: row,col
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pd2
  !--------------------------------------------------------------- 
  subroutine allocate_pz2(array,row,col,mlog)
    complex(dp), DIMENSION(:,:), POINTER :: array
    integer(4), intent(IN) :: row,col
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*2*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_pz2
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine allocate_c(array,length,mlog)
    character(*), DIMENSION(:), ALLOCATABLE :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog

    integer :: ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: pointer is already associated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*len(array(1))    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_c
  !---------------------------------------------------------------
  subroutine allocate_i(array,length,mlog)
    integer(4), DIMENSION(:), ALLOCATABLE :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*4    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_i
  !---------------------------------------------------------------
  subroutine allocate_d(array,length,mlog)
    real(kind=dp), DIMENSION(:), ALLOCATABLE :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (ALLOCATED(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_d
  !---------------------------------------------------------------
  subroutine allocate_z(array,length,mlog)
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (ALLOCATED(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*2*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_z
  !---------------------------------------------------------------   
 
  subroutine allocate_i2(array,row,col,mlog)
    integer(4), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*4    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_i2
  !---------------------------------------------------------------
  subroutine allocate_d2(array,row,col,mlog)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_d2
  !---------------------------------------------------------------
  subroutine allocate_z2(array,row,col,mlog)
    complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*2*dp    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_z2
  !---------------------------------------------------------------

  subroutine allocate_i3(array,row,col,dep,mlog)
    integer(4), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col,dep
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else 
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*4    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_i3
  !---------------------------------------------------------------
  subroutine allocate_d3(array,row,col,dep,mlog)
    real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col,dep
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*dp   
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_d3
  !---------------------------------------------------------------
  subroutine allocate_z3(array,row,col,dep,mlog)
    complex(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col,dep
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*2*dp   
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_z3
  !---------------------------------------------------------------
  subroutine allocate_d4(array,row,col,dep,qep,mlog)
    real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: array
    integer(4), intent(IN) :: row,col,dep,qep
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(row,col,dep,qep),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*dp   
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_d4

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine allocate_l(array,length,mlog)
    logical, DIMENSION(:), ALLOCATABLE :: array
    integer(4), intent(IN) :: length
    type(mem_log), intent(INOUT) :: mlog
    integer :: ierr
    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    else
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(*,*) "ALLOCATION ERROR"; STOP
       else
          mlog%alloc_mem= mlog%alloc_mem + size(array)*4    
          if (mlog%alloc_mem.gt.mlog%peak_mem) then
             mlog%peak_mem = mlog%alloc_mem 
          endif
       endif
    endif
  end subroutine allocate_l
  
  !---------------------------------------------------------------
  !--------------------------------------------------------------- 
  !---------------------------------------------------------------
  subroutine deallocate_pc(array,mlog)
    character(*), DIMENSION(:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*len(array(1))
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pc
  !---------------------------------------------------------------
  subroutine deallocate_pi(array,mlog)
    integer(4), DIMENSION(:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*4
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pi
  !---------------------------------------------------------------
  subroutine deallocate_pd(array,mlog)
    real(kind=dp), DIMENSION(:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*dp    
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pd
  !---------------------------------------------------------------
  subroutine deallocate_pz(array,mlog)
    complex(kind=dp), DIMENSION(:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*2*dp
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pz
   !---------------------------------------------------------------
  subroutine deallocate_pi2(array,mlog)
    integer(4), DIMENSION(:,:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*4
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pi2
  !---------------------------------------------------------------
  subroutine deallocate_pd2(array,mlog)
    real(kind=dp), DIMENSION(:,:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*dp
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pd2
  !---------------------------------------------------------------
  subroutine deallocate_pz2(array,mlog)
    complex(kind=dp), DIMENSION(:,:), POINTER :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*2*dp
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_pz2
  !---------------------------------------------------------------
  subroutine deallocate_l(array,mlog)
    logical, DIMENSION(:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*4
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_l
  !---------------------------------------------------------------
  subroutine deallocate_c(array,mlog)
    character(*), DIMENSION(:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (associated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*len(array(1))
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_c
  !---------------------------------------------------------------
  subroutine deallocate_i(array,mlog)
    integer(4), DIMENSION(:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*4
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_i
  !---------------------------------------------------------------
  subroutine deallocate_d(array,mlog)
    real(kind=dp), DIMENSION(:), allocatable :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*dp   
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d
  !---------------------------------------------------------------
  subroutine deallocate_z(array,mlog)
    complex(kind=dp), DIMENSION(:), allocatable :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*2*dp
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_z

  ! ------------------------------------------------------------
  subroutine deallocate_i2(array,mlog)
    integer(4), DIMENSION(:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*4
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_i2
  !---------------------------------------------------------------
  subroutine deallocate_d2(array,mlog)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*dp
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d2
  !---------------------------------------------------------------
  subroutine deallocate_z2(array,mlog)
    complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*2*dp
       deallocate(array)

    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_z2
  ! ------------------------------------------------------------

  subroutine deallocate_i3(array,mlog)
    integer(4), DIMENSION(:,:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*4
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_i3
  !---------------------------------------------------------------
  subroutine deallocate_d3(array,mlog)
    real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*dp
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d3
  !---------------------------------------------------------------
  subroutine deallocate_z3(array,mlog)
    complex(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*2*dp
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_z3
  !---------------------------------------------------------------
  subroutine deallocate_d4(array,mlog)
    real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: array
    type(mem_log), intent(INOUT) :: mlog  

    if (allocated(array)) then
       mlog%alloc_mem= mlog%alloc_mem - size(array)*dp
       deallocate(array)
    else 
       write(*,*) 'Warning in deallocation: array is not allocated' 
    endif
  end subroutine deallocate_d4

  ! ------------------------------------------------------------
  subroutine resetMemLog(mlog)

    type(mem_log) :: mlog       
    mlog%alloc_mem=0
    mlog%peak_mem=0
  end subroutine resetMemLog
  ! ------------------------------------------------------------
  subroutine writeMemInfo(mlog)

    type(mem_log) :: mlog
    character(3) :: str
    integer :: dec

    call memstr(mlog%alloc_mem,dec,str)
    write(mlog%iolog,'(A26,F8.2,A3)') 'current memory allocated: ',mlog%alloc_mem*1.0/dec,str

  end subroutine writeMemInfo
  ! ------------------------------------------------------------
  subroutine writePeakInfo(mlog)

    type(mem_log) :: mlog   
    character(3) :: str
    integer :: dec

    call memstr(mlog%peak_mem,dec,str)
    write(mlog%iolog,'(A26,F8.2,A3)') 'peak memory allocated: ',mlog%peak_mem*1.0/dec,str

  end subroutine writePeakInfo

  ! ------------------------------------------------------------
  subroutine openMemLog(iofile,mlog)
    integer iofile,err
    type(mem_log) :: mlog

    if(iofile.ne.6) then
       open(iofile,file='memory.log',iostat=err)
       if (err.ne.0) then
          write(*,*) 'Cannot open memory log-file'
          stop
       endif
    endif
    mlog%iolog=iofile	

  end subroutine openMemLog

  ! ------------------------------------------------------------
  subroutine writeMemLog(mlog)
    type(mem_log) :: mlog    
    
    call writeMemInfo(mlog)
  end subroutine writeMemLog

  ! ------------------------------------------------------------
  subroutine closeMemLog(mlog)
    type(mem_log) :: mlog    
    
    call writePeakInfo(mlog)
    if(mlog%iolog.ne.6) close(mlog%iolog)
  end subroutine closeMemLog

  ! ------------------------------------------------------------
  subroutine memstr(mem,dec,str)
 
    character(3) :: str
    integer(8) :: mem
    integer :: dec
    
    if(mem.lt.1000) then
       str=' bt'; dec=1
       return      
    endif
    
    if(mem.lt.1e7) then
       str=' kb'; dec=1000
       return      
    endif
    
    if(mem.ge.1e7) then
       str=' Mb'; dec=1000000
       return      
    endif
    
  end subroutine memstr

end module allocation
