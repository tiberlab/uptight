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
module mpi_globals

  implicit none

#ifdef UPT_MPI
  include 'mpif.h'
#endif

  ! MPI COMMON VARIABLES
  INTEGER :: upt_comm 
  INTEGER :: num_procs = 1
  INTEGER :: id = 0
  INTEGER :: upt_rank = 0
  INTEGER :: ierr

  LOGICAL :: id0 = .true.
  CHARACTER(2) :: id_s = " 0"

  INTEGER :: shift_init, shift_end

  INTEGER, DIMENSION(:), ALLOCATABLE :: shift_init_Mi, shift_end_Mi, shift_init_M, shift_end_M

contains

#ifdef UPT_MPI

  subroutine upt_mpi_init(comm)
    integer :: comm

    integer :: upt_group, numchar, is_init=0
    CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname


    call MPI_Initialized(is_init, ierr)
    if (is_init .eq. 0) then
       call MPI_INIT(ierr)
    end if

    if (comm .eq. 0) then
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)

      call MPI_COMM_GROUP(MPI_COMM_WORLD, upt_group, ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD, upt_group, upt_comm, ierr)
    else
      call MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN, ierr)

      call MPI_COMM_GROUP( comm, upt_group, ierr)
      call MPI_COMM_CREATE( comm, upt_group, upt_comm, ierr)
    end if

    call MPI_COMM_RANK( upt_comm, upt_rank, ierr)
    call MPI_COMM_SIZE( upt_comm, num_procs, ierr)
    call MPI_GET_PROCESSOR_NAME(hostname, numchar, ierr)
   
    id = upt_rank
 
    if (upt_rank .eq. 0) then
         id0 = .true.
    else
         id0 = .false.
    endif
    write(id_s,'(i2.2)') id 
 
    if (id0) write(*,*) "Starting UPT PARALLEL Version on #ranks", num_procs
    write(*,*) "@rank "//id_s//" alive "//trim(hostname)

  end subroutine upt_mpi_init


  subroutine upt_mpi_end()
    integer :: ierr
    call MPI_FINALIZE(ierr)
    write(*,*) '(mpi_globals) Finalizing rank= '//id_s
  end subroutine upt_mpi_end


  subroutine init_shifts(nrows, offset)
    INTEGER, intent(in) :: nrows, offset
    
    INTEGER :: ierr
    INTEGER :: row_offset(num_procs)

    if (allocated(shift_init_Mi) .or. allocated(shift_end_Mi) ) then
       write(*,*) "ERROR: shifts already allocated!"
       stop
    end if

    allocate(shift_init_Mi(0:num_procs-1), STAT=ierr)
    allocate(shift_end_Mi(0:num_procs-1), STAT=ierr)
    if (ierr.ne.0) STOP 'init_shifts allocation error'

    if (allocated(shift_init_M) .or. allocated(shift_end_M) ) then
       write(*,*) "ERROR: shifts already allocated!"
       stop
    end if

    allocate(shift_init_M(0:num_procs-1), STAT=ierr)
    allocate(shift_end_M(0:num_procs-1), STAT=ierr)
    if (ierr.ne.0) STOP 'init_shifts allocation error'

    !gather and bcast offset
    call mpi_barrier(upt_comm,ierr)

    call mpi_gather(offset, 1, MPI_INTEGER, row_offset, 1, MPI_INTEGER, 0, upt_comm, ierr) 
    call mpi_barrier(upt_comm,ierr)

    call MPI_Bcast(row_offset,num_procs, MPI_INTEGER, 0, upt_comm,ierr);
    call mpi_barrier(upt_comm,ierr)


    if ( id.eq. 0 .and. id .eq. num_procs-1 ) then 
       shift_init_M(id) = 1 
       shift_end_M(id)  = nrows
    elseif ( id .eq. 0) then
       shift_init_M(id) = 1
       shift_end_M(id) =  row_offset(id+2)
    elseif (  id .ne. num_procs-1 ) then
       shift_init_M(id) = row_offset(id+1)+1
       shift_end_M(id) =  row_offset(id+2)
    elseif (  id .eq. num_procs-1 ) then
       shift_init_M(id) = row_offset(id+1)+1        !for remain data
       shift_end_M(id) = nrows
    end if
    
    !write(*,*) "shift_init_Mi["//id_s//"]",shift_init_M(id)
    !write(*,*) "shift_end_Mi["//id_s//"]",shift_end_M(id)
        
    call mpi_gather(shift_init_M(id), 1, MPI_INTEGER, shift_init_Mi, 1, MPI_INTEGER, 0, upt_comm, ierr) 
    call mpi_barrier(upt_comm,ierr)
    call mpi_gather(shift_end_M(id), 1, MPI_INTEGER, shift_end_Mi, 1, MPI_INTEGER, 0, upt_comm, ierr) 
    call mpi_barrier(upt_comm,ierr)
    
    
    call MPI_Bcast(shift_init_Mi,num_procs, MPI_INTEGER, 0, upt_comm,ierr);
    call MPI_Bcast(shift_end_Mi,num_procs, MPI_INTEGER, 0, upt_comm,ierr);
    call mpi_barrier(upt_comm,ierr)

    shift_init = shift_init_Mi(id)
    shift_end = shift_end_Mi(id)

  end subroutine init_shifts

#else

  subroutine upt_mpi_init(comm)
    integer :: comm

    write(*,*) "Starting UPT Serial Version"

    num_procs = 1
    id = 0
    upt_rank = 0

    id0 = .true.
    id_s = " 0"

  end subroutine upt_mpi_init

  subroutine upt_mpi_end()
  end subroutine upt_mpi_end

  subroutine init_shifts(nrows, offset)
    INTEGER, intent(in) :: nrows, offset
    
    !allocate(shift_init_Mi(0:num_procs-1))
    !allocate(shift_end_Mi(0:num_procs-1))

    !shift_init_Mi(id) = 1 
    !shift_end_Mi(id)  = nrows
    
    shift_init = 1
    shift_end = nrows

  end subroutine init_shifts
    

#endif

  subroutine kill_shifts

    if (allocated(shift_init_Mi)) deallocate(shift_init_Mi, STAT=ierr)
    if (allocated(shift_end_Mi)) deallocate(shift_end_Mi, STAT=ierr)
    if (allocated(shift_init_M)) deallocate(shift_init_M, STAT=ierr)
    if (allocated(shift_end_M)) deallocate(shift_end_M, STAT=ierr)

    if (ierr.ne.0) STOP 'kill_shifts deallocation error'

  end subroutine kill_shifts

end module mpi_globals
