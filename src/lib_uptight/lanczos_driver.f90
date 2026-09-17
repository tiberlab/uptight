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
MODULE lanczos_driver
  USE mpi_globals
  USE precision
  USE upt_param
   USE sparse_matrix, only : CSR
  USE input_output
  USE errors
  USE lanczos_diag
  USE savemofile, only : append_eigenstate
   USE coarse_grain, only : cg_active, icg_active, icgn_active, cg_get_active, cg_lift_active

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lanczos

  contains

    subroutine lanczos(upt)
      
      type(oupt), target :: upt
   
      integer :: num_ev, n_ham, err, file_num, i, k, nv, nc, end_cb, end_vb
      integer :: num_cb, num_vb
      REAL ( dp ),   DIMENSION( : ),  POINTER     :: p_eigen_values
      COMPLEX ( dp ),  DIMENSION( :, : ), POINTER :: p_eigen_vectors

      REAL ( dp ),   DIMENSION( : ), ALLOCATABLE, TARGET   :: eigen_values2
      COMPLEX ( dp ),  DIMENSION( :, : ), ALLOCATABLE, TARGET :: eigen_vectors2
      INTEGER,  DIMENSION( : ), ALLOCATABLE, TARGET :: particles2      

      REAL ( dp ) :: shift
      LOGICAL :: spin_deg
      INTEGER :: len, verbose
 
      CHARACTER(LEN=:), ALLOCATABLE :: states_file

      if (cg_active(upt) .or. icg_active(upt) .or. icgn_active(upt)) then
         call lanczos_coarse(upt)
         return
      end if
      ! clear the states files, but only after all processes have read the states
      ! NOTE: maybe it would be more elegant to broadcast the vectors?
      if (id0) then
        states_file = trim(upt%state_file)
        CALL open_file( states_file, file_num, operation = "write", &
         replace_flag = .TRUE., output_flag = .FALSE. )
        close(file_num)
      end if


      num_ev = upt%num_vb + upt%num_cb
      num_cb = upt%num_cb - upt%start_cb + 1
      num_vb = upt%num_vb - upt%start_vb + 1
      n_ham = upt%ham%nrow
      verbose = upt%verbose

      if (id0 .and. associated(upt%eigen_vectors)) then
        ! write out already present states, but only on rank 0
        do i=1,num_ev

          if (upt%particles(i) .ne. 0) then
            CALL append_eigenstate(states_file, upt%eigen_vectors(:,i), &
                                   upt%eigen_values(i), upt%particles(i))

          end if
        end do
      endif

      err = 0
      nv = 0
      nc = 0

      spin_deg = .true.
      if ( .not.all(equiv(upt%k_point,0.d0,1.0d-13,.false.)) ) then
          spin_deg = .false.
      endif
        if (upt%n_spin == 1) spin_deg = .false.

      ! -------------------------------------------------------------------
      !  ALLOCATIONS
      ! -------------------------------------------------------------------
      if (associated(upt%eigen_values)) then
         if( size(upt%eigen_values).lt.num_ev ) then
            len = size(upt%eigen_values)
            allocate(eigen_values2(len), STAT = err)
            eigen_values2 = upt%eigen_values
            deallocate(upt%eigen_values)
            allocate(upt%eigen_values(num_ev), STAT = err)
            upt%eigen_values = 0.D0
            upt%eigen_values(1:len)=eigen_values2
            deallocate(eigen_values2)
         end if
      else
         allocate(upt%eigen_values(num_ev), STAT = err)
      end if

      IF (err.NE.0) CALL alloc_error('lanczos driver','main','eigen_values')

      if (associated(upt%particles)) then
         if( size(upt%particles).lt.num_ev ) then
            len = size(upt%particles)
            allocate(particles2(len), STAT = err)
            particles2 = upt%particles
            deallocate(upt%particles)
            allocate(upt%particles(num_ev), STAT = err)
            upt%particles = 0
            upt%particles(1:len)=particles2
            deallocate(particles2)
         end if
      else
         allocate(upt%particles(num_ev), STAT = err)
      end if

      IF (err.NE.0) CALL alloc_error('lanczos driver','main','particles')


      if (associated(upt%eigen_vectors)) then
         if( size(upt%eigen_vectors,2).lt.num_ev ) then
            len = size(upt%eigen_vectors,2)
            allocate(eigen_vectors2(n_ham,len), STAT = err)
            eigen_vectors2 = upt%eigen_vectors

            deallocate(upt%eigen_vectors)
            allocate(upt%eigen_vectors(n_ham,num_ev), STAT = err)
            upt%eigen_vectors = ( 0.0D0, 0.0D0 )
            upt%eigen_vectors(:,1:len)=eigen_vectors2
            deallocate(eigen_vectors2)
         end if
      else
         allocate(upt%eigen_vectors(n_ham,num_ev), STAT = err)
      end if
 
      IF (err.NE.0) CALL alloc_error('planczos driver','main','eigen_vectors')
      
      ! ---------------------------------------------------------------------

      if (num_cb .gt. 0) then

         if (verbose.gt.0) write(*,*) '(lanczos) number of conductions ',num_cb     
         if (verbose.gt.0) write(*,*) '(lanczos) start conduction ',upt%start_cb
         if (verbose.gt.0) write(*,*) '(lanczos) tolerance ',upt%long_tol
         if (verbose.gt.0) write(*,*) '(lanczos) dynamic search ',upt%dynamic

         p_eigen_values => upt%eigen_values(upt%num_vb+1:num_ev)
         p_eigen_vectors => upt%eigen_vectors(:,upt%num_vb+1:num_ev)
         upt%particles(upt%num_vb+1:num_ev) = 1
         
         end_cb = upt%start_cb + num_cb - 1
         
         CALL LANCZOS_EV(upt%ham, upt%U, upt%n_spin, upt%min_iter, upt%long_iter, &
                         upt%max_iter, p_eigen_values, p_eigen_vectors, &
                         upt%start_cb, end_cb, n_ham, upt%lambda_cb, &
                         upt%solver_flag, upt%fast_tol, upt%long_tol, upt%ort_tol, &
                         1, upt%dynamic, upt%bitoff, spin_deg, verbose, states_file )
         
      end if !end conductions
      
      ! -------------------------------------------------------------------------
      if (num_vb .gt. 0) then

         if (verbose.gt.0) write(*,*) '(lanczos) number of valence ',num_vb     
         if (verbose.gt.0) write(*,*) '(lanczos) start valence ',upt%start_vb
         if (verbose.gt.0) write(*,*) '(lanczos) tolerance ',upt%long_tol
         if (verbose.gt.0) write(*,*) '(lanczos) dynamic search ',upt%dynamic
         
  
         p_eigen_values => upt%eigen_values(1:upt%num_vb)
         p_eigen_vectors => upt%eigen_vectors(:,1:upt%num_vb)
         upt%particles(1:upt%num_vb) = -1

         end_vb = upt%start_vb + num_vb - 1

         CALL LANCZOS_EV(upt%ham,upt%U, upt%n_spin, upt%min_iter, upt%long_iter, &
              upt%max_iter, p_eigen_values, p_eigen_vectors, &
              upt%start_vb, end_vb, n_ham, upt%lambda_vb, & 
              upt%solver_flag, upt%fast_tol, upt%long_tol, &
              upt%ort_tol, -1, upt%dynamic, upt%bitoff, spin_deg, verbose, states_file )

      end if !end valence
      ! -------------------------------------------------------------------------

    end subroutine lanczos

    subroutine lanczos_coarse(upt)
      type(OUPT), target :: upt
      type(CSR), pointer :: active_ham, active_u
      logical :: active
      integer :: nfull, nred, num_ev, num_cb, num_vb, end_cb, end_vb, err
      integer :: old_shift_init, old_shift_end
      integer :: old_shift_init_mi, old_shift_end_mi
      real(dp), allocatable, target :: raw_values(:)
      complex(dp), allocatable, target :: raw_vectors(:,:)
      real(dp), pointer :: raw_values_slice(:)
      complex(dp), pointer :: raw_vectors_slice(:,:)
      complex(dp), allocatable :: reduced_vectors(:,:), lifted(:,:)

      call cg_get_active(upt, active_ham, active_u, active)
      if (.not.active) return
      nred = active_ham%nrow
      nfull = upt%ham%nrow
       if (upt%verbose > 0 .and. id0) write(*,'(a,i0,a,i0,a,a1,a,i0,a,i0)') &
          '(cg lanczos) active rows ', nred, ', nnz ', active_ham%nnz, &
          ', format ', active_ham%sparse_fmt, ', Mi(1) ', active_ham%Mi(1), &
          ', Mi(end) ', active_ham%Mi(nred+1)
      num_ev = upt%num_vb + upt%num_cb
      num_cb = upt%num_cb - upt%start_cb + 1
      num_vb = upt%num_vb - upt%start_vb + 1
      allocate(raw_values(num_ev), raw_vectors(nred,num_ev), stat=err)
      if (err /= 0) call alloc_error('Lanczos coarse grain','allocate','vectors')
      raw_values = 0.0_dp
      raw_vectors = (0.0_dp, 0.0_dp)
      old_shift_init = shift_init
      old_shift_end = shift_end
      old_shift_init_mi = shift_init_Mi(id)
      old_shift_end_mi = shift_end_Mi(id)
      shift_init = 1
      shift_end = nred
      shift_init_Mi(id) = 1
      shift_end_Mi(id) = nred

      if (num_cb > 0) then
         raw_values_slice => raw_values(upt%num_vb+1:num_ev)
         raw_vectors_slice => raw_vectors(:,upt%num_vb+1:num_ev)
         end_cb = upt%start_cb + num_cb - 1
         call LANCZOS_EV(active_ham, active_u, 1, upt%min_iter, upt%long_iter, &
              upt%max_iter, raw_values_slice, raw_vectors_slice, upt%start_cb, &
              end_cb, nred, upt%lambda_cb, upt%solver_flag, upt%fast_tol, &
              upt%long_tol, upt%ort_tol, 1, upt%dynamic, upt%bitoff, .false., &
              upt%verbose)
      end if
      if (num_vb > 0) then
         raw_values_slice => raw_values(1:upt%num_vb)
         raw_vectors_slice => raw_vectors(:,1:upt%num_vb)
         end_vb = upt%start_vb + num_vb - 1
         call LANCZOS_EV(active_ham, active_u, 1, upt%min_iter, upt%long_iter, &
              upt%max_iter, raw_values_slice, raw_vectors_slice, upt%start_vb, &
              end_vb, nred, upt%lambda_vb, upt%solver_flag, upt%fast_tol, &
              upt%long_tol, upt%ort_tol, -1, upt%dynamic, upt%bitoff, .false., &
              upt%verbose)
      end if

      if (associated(upt%eigen_values)) deallocate(upt%eigen_values)
      if (associated(upt%eigen_vectors)) deallocate(upt%eigen_vectors)
      if (associated(upt%particles)) deallocate(upt%particles)
      allocate(upt%eigen_values(num_ev), upt%eigen_vectors(nfull,num_ev), &
           upt%particles(num_ev), stat=err)
      if (err /= 0) call alloc_error('Lanczos coarse grain','allocate','states')
      upt%eigen_values = raw_values
      upt%particles = 0
      allocate(lifted(nfull,num_ev), stat=err)
      if (err /= 0) call alloc_error('Lanczos coarse grain','allocate','lifted')
      call cg_lift_active(upt, raw_vectors, lifted)
      shift_init = old_shift_init
      shift_end = old_shift_end
      shift_init_Mi(id) = old_shift_init_mi
      shift_end_Mi(id) = old_shift_end_mi
      upt%eigen_vectors = lifted
      deallocate(raw_values, raw_vectors, lifted)
    end subroutine lanczos_coarse


END MODULE lanczos_driver



! This was an old strategy and was removed
! DOUBLE FOLDING TRICK:
!
!       |      |
!       |      |
!  v1 ++++++++++++
!  C2 ------------
!  v1'++++++++++++
!  C1 ------------
!       |      |
!       ________
!
!  s  . . . . . . .
!  s' . . . . . . .
!
!
!       ________
!       |      |
!  V1 ------------
!       |      |
!  V2 ------------
!       |      |
!
!
!   C1 = (C1 - s) + s   = C1
!   C1'= (C1'- s + ds) + s - ds = C1
!
!   v1 = (V1 - s) + s  = v1
!   v1'= (V1 - s - ds) + s - ds = v1 - 2 ds
!
! => Folded eigenvalues differ by 2*ds
!
! This is true only if guess is below conduction and above valence
