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
MODULE jd_driver
  USE mpi_globals
  USE precision
  USE upt_param
  USE input_output
  USE errors
  USE jd_diag
  USE savemofile, only : append_eigenstate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jd

  contains

    subroutine jd(upt)
      
      type(oupt) :: upt
   
      integer :: num_ev, n_ham, err, file_num, i, k, nv, nc, end_cb, end_vb
      integer :: num_cb, num_vb, band_type
      REAL ( dp ),   DIMENSION( : ),  POINTER     :: p_eigen_values
      COMPLEX ( dp ),  DIMENSION( :, : ), POINTER :: p_eigen_vectors

      REAL ( dp ),   DIMENSION( : ), ALLOCATABLE, TARGET   :: eigen_values2
      COMPLEX ( dp ),  DIMENSION( :, : ), ALLOCATABLE, TARGET :: eigen_vectors2
      INTEGER,  DIMENSION( : ), ALLOCATABLE, TARGET :: particles2      

      CHARACTER(LEN=:), ALLOCATABLE :: statesfile

      REAL ( dp ) :: shift
      LOGICAL :: spin_deg
      INTEGER :: len, verbose


      num_ev = upt%num_vb + upt%num_cb
      num_cb = upt%num_cb - upt%start_cb + 1
      num_vb = upt%num_vb - upt%start_vb + 1
      n_ham = upt%ham%nrow
      verbose = upt%verbose
      err = 0
      nv = 0
      nc = 0

      spin_deg = .true.
      if ( .not.all(equiv(upt%k_point,0.d0,1.0d-13,.false.)) ) then
          spin_deg = .false.
      endif

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

      IF (err.NE.0) CALL alloc_error('JD driver','main','eigen_values')

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

      IF (err.NE.0) CALL alloc_error('JD driver','main','particles')


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

         band_type = 1

         if (verbose.gt.0) write(*,*) '(JD) number of conductions ',num_cb     
         if (verbose.gt.0) write(*,*) '(JD) start conduction ',upt%start_cb
         if (verbose.gt.0) write(*,*) '(JD) tolerance ',upt%long_tol
         !if (verbose.gt.0) write(*,*) '(lanczos) dynamic search ',upt%dynamic

         p_eigen_values => upt%eigen_values(upt%num_vb+1:num_ev)
         p_eigen_vectors => upt%eigen_vectors(:,upt%num_vb+1:num_ev)
         upt%particles(upt%num_vb+1:num_ev) = 1
         
         end_cb = upt%start_cb + num_cb - 1
         
         CALL JD_EV(upt%ham, upt%U, upt%n_spin, upt%min_iter, upt%long_iter, &
                         upt%max_iter, p_eigen_values, p_eigen_vectors, &
                         upt%start_cb, end_cb, n_ham, upt%lambda_cb, &
                         upt%solver_flag, upt%fast_tol, upt%long_tol, upt%ort_tol, &
                         1, upt%dynamic, spin_deg, verbose, band_type )
         
      end if !end conductions
      
      ! -------------------------------------------------------------------------
    
      if (num_vb .gt. 0) then

         band_type = 2

         if (verbose.gt.0) write(*,*) '(JD) number of valence ',num_vb     
         if (verbose.gt.0) write(*,*) '(JD) start conduction ',upt%start_vb
         if (verbose.gt.0) write(*,*) '(JD) tolerance ',upt%long_tol
         !if (verbose.gt.0) write(*,*) '(lanczos) dynamic search ',upt%dynamic
         
  
         p_eigen_values => upt%eigen_values(1:upt%num_vb)
         p_eigen_vectors => upt%eigen_vectors(:,1:upt%num_vb)
         upt%particles(1:upt%num_vb) = -1

         end_vb = upt%start_vb + num_vb - 1

         CALL JD_EV(upt%ham, upt%U, upt%n_spin, upt%min_iter, upt%long_iter, &
                        upt%max_iter, p_eigen_values, p_eigen_vectors, &
                        upt%start_vb, end_vb, n_ham, upt%lambda_vb, &
                        upt%solver_flag, upt%fast_tol, upt%long_tol, upt%ort_tol, &
                        -1, upt%dynamic, spin_deg, verbose, band_type )

      end if !end valence
! -------------------------------------------------------------------------

      ! write out states to file
      if (id0) then
        ! delete the file
        statesfile = trim(upt%state_file)
        CALL open_file( statesfile, file_num, operation = "write", &
          replace_flag = .TRUE., output_flag = .FALSE. )
        close(file_num)

        do i = 1,num_ev
         CALL append_eigenstate(statesfile, upt%eigen_vectors(:,i), &
           upt%eigen_values(i), upt%particles(i))
        enddo
      endif


    end subroutine jd


END MODULE jd_driver
