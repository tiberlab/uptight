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
!=============================================================================
!
!                              LANCZOS SOLVER
!
!=============================================================================
!
! Alessandro Pecchia
! Walter Rodrigues
! 
! Dipartimento di Ingegneria Elettronica
! Universita` di Roma "Tor Vergata"
! Tel. +39-6-72597781
! E-mail: dicarlo@ing.uniroma2.it
!         pecchia@ing.uniroma2.it
!     
! 28-05-2013
!
!=============================================================================
!
! This subroutine finds the lowest eigenvalue of the operator:
!
!                       ( H - lambda * I )^2 ,
!
! where lambda is a given real number and H is an Hermitian matrix.
!
! There are three main steps in the iterative solution:
!
!   1) a fast Lanczos loop ( generally a 2x2 T matrix )
!   2) a slow Lanczos loop ( with bigger T matrix )
!   3) a final standard Lanczos iteration
!
! The first two steps are repeated several times. The last step (3)
! should never be used and it is included here just for completeness.
!
!_____________________________________________________________________________
!
! INPUT :
!
! => The H matrix is stored in a sparse format as follows:
!
!    ->  an integer vector "a" contains the row and column indeces for
!        the non-zero elements of H :
!
!          -row (1)
!           col (1)   [ element { row(1), col(1) } ]
!           col (2)   [ element { row(1), col(2) } ]
!             .
!             .
!          -row (2)
!           col (1)   [ element { row(2), col(1) } ]
!             .
!             .
!
!        the rows are marked by a minus sign, thus -r is the r-th row.
!
!     -> the matrix elements are given in a double complex vector "ab".
!        Each "ab" entry is the non-zero element of the H matrix in a row-col
!        order accordig the "a" vector.
!
!_____________________________________________________________________________
!
! OUTPUT :
!
! => lower eigenvalues are stored at each step in the file 'result.dat',
!      which has the following format :
!
!           iteration    first eig.    second eig. etc.
!
!_____________________________________________________________________________
!
! TODO :
!
! -> for memory matrix, calculate ( H - lambda * I )^2 at the beginning.
!
! -> some optimizations (t = b -> loop with temp = b(i), etc.)
!
!=============================================================================
!
! The sparse matrix is defined by two vectors :
!
!  -> "sparse_ind" which is the indeces vector
!  -> "sparse_mat" which contains lower-diagonal values
!      of the matrix corresponding to the indeces found in "sparse_ind".
!
! There are also other parameters :
!
!  -> lambda   : the off-set parameter
!  -> n_ham    : dimension of actual hamiltonian matrix
!
! Since thay are all defined elsewhere, we use them as global input arguments,
! so their type does not appear in this module.
!
!=============================================================================

MODULE LANCZOS_DIAG

!===========================================================================
  USE mpi_globals
  USE precision
  USE errors
  USE exceptions
  USE input_output
  USE savemofile, only : append_eigenstate
  USE sparse_matrix
  USE sparse_numrec
  USE sort               !Collection of sorting routines
  USE clock
  USE omp_lib 
!===========================================================================

  IMPLICIT NONE
  PRIVATE

!===========================================================================
  
  PUBLIC lanczos_ev 

!===========================================================================

  INTERFACE project_out
     module procedure project_out_r
     module procedure project_out_z
     module procedure project_out_sp
     module procedure project_out_z_p
  END INTERFACE
    
CONTAINS
  
!===========================================================================
!
! Subroutines - modified by Walter in 2013 for MPI/OpenMP/CUDA
!
! Driver for Lanczos algorithm
!
!===========================================================================
SUBROUTINE LANCZOS_EV(H, U, n_spin, min_step, long_step, max_step, &
                        energies, eigen_vectors, start_ev, num_ev, n_ham, &
                        lambda, solver_flag, fast_tol, long_tol, ort_tol, &
                        sign, dynamic, bitoff, spin_deg, verbose, statesfile )

    !=========================================================================
    ! Input arguments
    !=========================================================================
    TYPE(CSR)                     :: H, U
    INTEGER,         INTENT( IN ) :: n_spin, min_step, long_step, max_step
    INTEGER,         INTENT( IN ) :: start_ev, num_ev, n_ham
    REAL ( dp ),     INTENT( IN ) :: lambda
    INTEGER,         INTENT( IN ) :: sign
    INTEGER,         INTENT( IN ) :: solver_flag
    REAL ( dp ),     INTENT( IN ) :: fast_tol, long_tol, ort_tol
    LOGICAL,         INTENT( IN ) :: dynamic
    REAL ( dp ),     INTENT( IN ) :: bitoff
    LOGICAL,         INTENT( IN ) :: spin_deg
    INTEGER,         INTENT( IN ) :: verbose
    CHARACTER(*), INTENT(IN), optional :: statesfile 
    !=========================================================================
    ! Output arguments
    !=========================================================================
    REAL ( dp ),   DIMENSION( : ),    POINTER     :: energies
    COMPLEX ( dp ), DIMENSION( :,: ),    POINTER  :: eigen_vectors
    !=========================================================================
    ! Local variables
    !=========================================================================
    INTEGER :: file_num_res
    LOGICAL :: res_flag
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mi !ROWPNT
    INTEGER,        DIMENSION(:), POINTER :: Mj  !COLIND

    COMPLEX ( dp ), DIMENSION( : ), ALLOCATABLE  :: eigen_seed
    COMPLEX ( dp ), DIMENSION( :,: ), ALLOCATABLE :: eigvec
    COMPLEX ( dp ), DIMENSION( : ), ALLOCATABLE  :: p_vec, p_temp    

    REAL ( dp )    :: energy, deltaE, shift, shift_tot, test_E
    REAL ( dp )    :: local_temp_1, local_temp_2, temp_dot_1, temp_dot_2 
    INTEGER :: nr_eigv, counter, err, row, col, tid
    INTEGER :: i,j,k, rank, block_size, colind_l, colind_h

    INTEGER :: colind_low(0:num_procs-1), colind_high(0:num_procs-1)

    TYPE(CSR_real)   :: H_real
    TYPE(CSR_real)   :: H_imag


    ALLOCATE( eigen_seed( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'eigen_seed' )
    ALLOCATE( p_vec( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'eigen_seed' )  
    ALLOCATE( p_temp( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'p_temp' )

    res_flag = .false.
    shift = lambda
    shift_tot = 0.d0 

    !=========================================================================
    ! Initialize counters
    !=========================================================================
    !=========================================================================
    ! LOOP OVER NUMBER OF EIGENVECTORS TO FIND
    !=========================================================================
    
    SELECT CASE (solver_flag)
    CASE(-1)
      if (id0 .and. verbose.gt.0) THEN
            write(*,*) ' '
            write(*,*) 'Serial version of Lanczos solver'
      end if       
    CASE(0)
      if(H%sparse_fmt.eq.'F') then
         if (id0 .and. verbose.gt.0) THEN
            write(*,*) ' '
            write(*,*) 'MPI parallelized on ',num_procs,'processes'
#ifdef __OMP
            write(*,*) 'OpenMP parallelized on',omp_get_max_threads(),'threads'
#endif
         end if
      endif
#ifdef __CUDA
    CASE(1)
       write(*,*) 'GPU accelerated'
    CASE(2)  
       write(*,*) 'GPU accelerated. Splitting matrix'
       call set_clock()
       call split_matrix(H, H_real, H_imag)
       call write_clock
       !print*,'real matrix structure'
       !call analyze(H_real%Mi)
       !print*,'imag matrix structure'
       !call analyze(H_img%Mi)
#endif
     CASE DEFAULT
       WRITE(*,*) 'ERROR Undefined solver flag',solver_flag
       call throw_solve_exception(ERR_LANCZ_DIAG)
    END SELECT 
    
    
    !====================================================================
    ! PARRAVICINI-GROSSO RESTART LANCZOS
    !====================================================================
    if (.not.dynamic) THEN
       call sprs_shift( H%M, H%Mj, H%Mi, -shift )
       if (solver_flag .eq. 2) CALL sprs_shift( H_real%M, H%Mj, H%Mi, -shift)
       shift_tot = shift_tot + shift 
    endif


#ifdef UPT_MPI   
    !////////////////////////////////////////////////////////////////////// 
    ! COMPUTE OVERLAP between nodes
    !//////////////////////////////////////////////////////////////////////

    call compute_ovr(H%M, H%Mj, H%Mi, colind_l, colind_h)
    !write(*,*) '# Rank ', id, colind_l, colind_h

    call mpi_gather(colind_l, 1, MPI_INTEGER, colind_low, 1, &
         MPI_INTEGER, 0, upt_comm, ierr) 
    
    call mpi_gather(colind_h, 1, MPI_INTEGER, colind_high, 1, &
         MPI_INTEGER, 0, upt_comm, ierr) 
    
    call MPI_Bcast(colind_low, num_procs, MPI_INTEGER, 0, upt_comm,ierr);
    
    call MPI_Bcast(colind_high, num_procs, MPI_INTEGER, 0, upt_comm,ierr);
    !=========================================================================
    ! Check that only neighbour nodes need communications
    ! Master does check.
    !=========================================================================
    do  j = 0,num_procs-2
       if(colind_low(j+1) .LT. shift_init_Mi(j) ) then
          if (id0) write(*,*) 
          if (id0) write(*,*) 'The matrix is distributed on too many MPI nodes'
          call throw_solve_exception(ERR_LANCZ_DIAG) 
       endif
    enddo

#endif

    !=========================================================================
    ! LOOP OVER NUMBER OF EIGENVECTORS TO FIND 
    !=========================================================================  
    nr_eigv = start_ev

    DO WHILE (nr_eigv.le.num_ev)  
       
       if (verbose.gt.0) WRITE ( *, * ) 
       if (verbose.gt.0) WRITE ( *, * ) 'shift [eV]  = ', shift
       
       if (dynamic) then
          CALL sprs_shift(  H%M, H%Mj, H%Mi, -shift )
          if (solver_flag .eq. 2) CALL sprs_shift( H_real%M, H%Mj, H%Mi, -shift)
          shift_tot = shift_tot + shift 
       endif
       if (verbose.gt.0) WRITE ( *, * ) 'total shift [eV]  = ', shift_tot

       counter = 0
       
       if (verbose.gt.0) call set_clock()

       SELECT CASE (solver_flag)
         CASE(-1)

            CALL fast_lanczos_ev_old( H%M, H%Mi, H%Mj,  H%sparse_fmt, &
                 min_step, long_step, max_step, eigen_seed, &
                 n_ham, fast_tol, long_tol, ort_tol, counter, & 
                 res_flag, eigen_vectors, nr_eigv, energy, deltaE, verbose)

         CASE(0)

#ifdef UPT_MPI
            CALL fast_lanczos_ev_pold_opt( H%M, H%Mi, H%Mj,  H%sparse_fmt, &
               min_step, long_step, max_step, eigen_seed, &
               n_ham, fast_tol, long_tol, ort_tol, counter, & 
               res_flag, eigen_vectors, num_ev, nr_eigv, energy, deltaE, &
               colind_low, colind_high, verbose)

            call mpi_barrier(upt_comm,ierr)
#else
            CALL fast_lanczos_ev_old( H%M, H%Mi, H%Mj,  H%sparse_fmt, &
                 min_step, long_step, max_step, eigen_seed, &
                 n_ham, fast_tol, long_tol, ort_tol, counter, & 
                 res_flag, eigen_vectors, nr_eigv, energy, deltaE, verbose)
#endif

#ifdef __CUDA
         CASE(1)
            CALL fast_lanczos_ev_cusparse(size(H%M), H%M, H%Mi, H%Mj, &
                 H%sparse_fmt, &
                 min_step, long_step, max_step, eigen_seed, &
                 n_ham, fast_tol, long_tol, ort_tol, counter, & 
                 res_flag, eigen_vectors, nr_eigv, energy, deltaE)
         CASE(2)
            CALL fast_lanczos_ev_cusparse_split(size(H_real%M), size(H_imag%M), & 
                 H_real%M, H_real%Mi, H_real%Mj, H_imag%M, H_imag%Mi, H_imag%Mj, &
                 H%sparse_fmt, &
                 min_step, long_step, max_step, eigen_seed, &
                 n_ham, fast_tol, long_tol, ort_tol, counter, & 
                 res_flag, eigen_vectors, nr_eigv, energy, deltaE)
#endif
         CASE DEFAULT
            WRITE(*,*) 'ERROR Undefined solver flag',solver_flag
            call throw_solve_exception(ERR_LANCZ_DIAG)
       END SELECT

       if (verbose.gt.0) call write_clock()
       
       ! CHECK IF FOLDED STATE.
       IF (energy*sign .lt. 0) THEN
          
          write(*,*) 'WARNING: FOLDED state at ',energy+shift_tot
          IF (dynamic) THEN
             write(*,*) '         Auto change shift and restart'
             ! Signs:   H - s = H + (-s) 
             ! shift opposite to the folded energy is safe by construction 
             shift = -energy
          ELSE
             write(*,*) 'WARNING: CANNOT shift'     
             write(*,*) '         Please use option: dynamic=true'     
             energies(nr_eigv) = energy + shift_tot
             eigen_vectors(:,nr_eigv) = eigen_seed
             nr_eigv = nr_eigv + 1
          ENDIF
          
       ELSE     
          
          if (verbose.gt.0) then
             WRITE(*,*) 
             WRITE(*,*) 'Eigenvalue number',nr_eigv    
             WRITE(*,*) 'Eigenvalue = ', energy+shift_tot, 'Error = ', deltaE
          end if
          
          energies(nr_eigv) = energy + shift_tot
          eigen_vectors(:,nr_eigv) = eigen_seed 
          
          if (id0 .and. present(statesfile)) &
            CALL append_eigenstate(statesfile, eigen_vectors(:,nr_eigv), energies(nr_eigv), sign)
            
          nr_eigv = nr_eigv + 1
          
          if (spin_deg .and. n_spin.eq.2 .and. nr_eigv.le.num_ev) then
             
             if (verbose.gt.0) WRITE(*,*) 
             if (verbose.gt.0) WRITE(*,*) 'Build spin degenerate state', nr_eigv

             call sprs_ax(U%M, U%Mj, U%Mi, U%sparse_fmt, eigen_seed, p_vec)

            
             if (num_procs.gt.1) then
#ifdef UPT_MPI 

                call gather_on_master(p_vec)

                call MPI_Bcast(p_vec, n_ham, MPI_DOUBLE_COMPLEX, 0, upt_comm,ierr)
 
                call mpi_barrier(upt_comm,ierr)

                p_vec = conjg(p_vec)
              

                !test_E = test_eigvect(H%M, H%Mj, H%Mi, H%sparse_fmt, p_vec)

                call sprs_ax( H%M, H%Mj, H%Mi,  H%sparse_fmt, p_vec, p_temp)
                call mpi_barrier(upt_comm,ierr)

                local_temp_1 = DOT_PRODUCT(p_vec(shift_init:shift_end), &
                           p_temp(shift_init:shift_end) )

                call mpi_barrier(upt_comm,ierr)
                call mpi_allreduce(local_temp_1, temp_dot_1, 1 , MPI_double_precision, &
                     MPI_SUM, upt_comm,ierr)
 
                local_temp_2 = DOT_PRODUCT(p_vec(shift_init:shift_end), &
                           p_vec(shift_init:shift_end) )

                call mpi_barrier(upt_comm,ierr)
                call mpi_allreduce(local_temp_2, temp_dot_2, 1 , MPI_double_precision, &
                     MPI_SUM, upt_comm,ierr)

                test_E = temp_dot_1/temp_dot_2


#endif
             else
                p_vec = conjg(p_vec)
                test_E = test_eigvect(H%M, H%Mj, H%Mi, H%sparse_fmt, p_vec)
             endif!end numprocs gt 1


             if (verbose.gt.0) WRITE(*,*) 'Eigenvalue = ', test_E + shift_tot
             
             if (abs(dot_product(eigen_seed,p_vec)).gt.1.0d-13) then
                write(*,*) 'ERROR: up and down are not orthogonal'
                call throw_solve_exception(ERR_LANCZ_DIAG)
             endif
             
             energies(nr_eigv) = test_E + shift_tot
             eigen_vectors(:,nr_eigv) = p_vec(:)  
          
             if (present(statesfile)) &
               CALL append_eigenstate(statesfile, eigen_vectors(:,nr_eigv), energies(nr_eigv), sign)
               
             nr_eigv = nr_eigv + 1            
              
          end if
          
          if (dynamic) then
             shift = energy - sign*bitoff
          end if
          
       END IF

    END DO 


    if (verbose.gt.0) WRITE(*,*) '========================================='
    if (verbose.gt.0) WRITE(*,*)
    ! Restore original H 
    CALL sprs_shift( H%M, H%Mj, H%Mi, shift_tot )
    if (solver_flag.eq.2) then
       call destroy_matrix(H_real)
       call destroy_matrix(H_imag)
    endif


    DEALLOCATE( eigen_seed, STAT = err )
    IF ( err .NE. 0 ) CALL dealloc_error( 'lanczos', 'lanczos', 'eigen_seed' )

    DEALLOCATE( p_vec, STAT = err )
    IF ( err .NE. 0 ) CALL dealloc_error( 'lanczos', 'lanczos', 'p_vec' )

  END SUBROUTINE LANCZOS_EV

  !=========================================================================
  ! END LANCZOS
  !=========================================================================
#ifdef UPT_MPI

SUBROUTINE fast_lanczos_ev_pold_opt( M, Mi, Mj,  sparse_fmt, &
        min_step, long_step, max_step, eigen_seed, n_ham, &
        fast_tol, long_tol, ort_tol_in, counter, res_flag, eigen_v, num_ev, nr_eigv, &
        energy, deltaE, colind_low, colind_high, verbose)


    !=========================================================================
    ! Input arguments
    !=========================================================================
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mi  !ROWPNT
    INTEGER,        DIMENSION(:), POINTER :: Mj  !COLIND
    CHARACTER(1)                  :: sparse_fmt

    INTEGER, INTENT( IN ) :: min_step   ! # of iterations of fast lanczos
    INTEGER, INTENT( IN ) :: long_step  ! # of iterations of long lanczos
    INTEGER, INTENT( IN ) :: max_step   ! # max number of iterations
    LOGICAL, INTENT( IN ) :: res_flag
    INTEGER, INTENT( IN ) :: nr_eigv    ! eigenvalue number
    INTEGER, INTENT( IN ) :: n_ham, num_ev
    REAL ( dp ), INTENT( IN ) :: fast_tol, long_tol, ort_tol_in
    INTEGER, DIMENSION(:) :: colind_low, colind_high 
    INTEGER, INTENT( IN ) :: verbose
    !=========================================================================
    ! Input - output arguments
    !=========================================================================

    COMPLEX( dp ), DIMENSION( : ), INTENT( INOUT ) :: eigen_seed
    COMPLEX( dp ), DIMENSION( : ), ALLOCATABLE, TARGET :: local_eigen_seed

    INTEGER                      , INTENT( INOUT ) :: counter
    COMPLEX( dp ), DIMENSION( :,: ), INTENT( IN )  :: eigen_v
    !=========================================================================
    ! Output arguments
    !=========================================================================

    REAL ( dp ), INTENT( OUT ) :: energy
    REAL ( dp ), INTENT( OUT ) :: deltaE

    !=========================================================================
    ! Local variables
    !=========================================================================
    REAL ( dp ), ALLOCATABLE, DIMENSION( : )    :: temp_v
    REAL ( dp ), POINTER    , DIMENSION( :, : ) :: eigen_vec
    REAL ( dp ), POINTER    , DIMENSION( : )    :: eigen_val

    REAL ( dp ), ALLOCATABLE, DIMENSION( : )    :: ad1,ad2
    REAL ( dp ), ALLOCATABLE, DIMENSION( : )    :: work
    INTEGER, ALLOCATABLE, DIMENSION( : )    :: iwork, ifail

    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE, TARGET    :: j_current
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE,TARGET    :: j_previous
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE, TARGET    :: j_next

    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE, TARGET    :: j_aux
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE,TARGET    :: j_aux2

    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_previous
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_next

    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_current
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_aux
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_aux2

    real( dp ), DIMENSION(10) :: energy_test

    REAL ( dp ) :: alpha, beta, test_E, ort_chk, ort_tol
    REAL ( dp ) :: local_alpha, local_beta, local_ort_chk

    COMPLEX ( dp ) :: eigen_guess

    REAL ( dp )    :: tol, errors, ABSTOL
    REAL ( dp )    :: dev_en, mean_en, tempval

    INTEGER :: j, nstep, info, ldz, VL, VU, IL, IU, M_out, err,N, i, loop, p
    CHARACTER :: JOBZ, RANGE
    CHARACTER(3) :: str_nr
    LOGICAL :: test_fast, converged

    INTEGER :: rank, overlapR, overlapL


    !======================================================================
    ! ALLOCATE LOCAL VECTORS
    !======================================================================
   
    ALLOCATE( j_aux( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_aux' )
    ALLOCATE( j_aux2( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_aux2' )
    ALLOCATE( j_current( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_current' )
    ALLOCATE( j_next( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_next' )
    ALLOCATE( j_previous( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_previous' )
    ALLOCATE( local_eigen_seed( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'local_eigen_seed' )
        
    p_current => j_current
    p_aux => j_aux
    p_aux2 => j_aux2
    p_next => j_next
    p_previous => j_previous

    call mpi_barrier(upt_comm,ierr)

    !======================================================================
    ! Allocate work arrays for eigen values and vectors
    !====================================================================== 
    !----------------------------------------------------------------------
    ! Only node 0 form and diagonalize tri-diagonal T matrix
    !----------------------------------------------------------------------

    if (id0) then

       JOBZ = 'V'
       RANGE = 'I'
       IL = 1; VL = 1
       IU = 1; VU = 1
       ABSTOL = 1e-12
       ldz = max_step
       ort_tol = ort_tol_in
       ort_chk = 0.0D0

       ALLOCATE( eigen_val( 2 ), STAT = err )
       IF ( err .NE. 0 ) &
            CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'eigen_val' )

       ALLOCATE( work( 5*ldz ), STAT = err )
       IF ( err .NE. 0 ) &
            CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'work' )
 
       ALLOCATE( iwork( 5 * ldz ), STAT = err )
       IF ( err .NE. 0 ) &
            CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'iwork' )
       
       ALLOCATE( ifail( ldz ), STAT = err )
       IF ( err .NE. 0 ) &
            CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'ifail' )
    endif

    call MPI_Bcast(ldz,1, MPI_INTEGER, 0, upt_comm,ierr)
    call MPI_Bcast(ort_tol,1, MPI_double_precision, 0, upt_comm,ierr)
    call MPI_Bcast(ort_chk,1, MPI_double_precision, 0, upt_comm,ierr)
    call mpi_barrier(upt_comm,ierr)

    ALLOCATE( eigen_vec( ldz, 2 ), STAT = err )
    IF ( err .NE. 0 ) &
         CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'eigen_vec' )

    ALLOCATE( ad1( ldz ), STAT = err )
    ALLOCATE( ad2( ldz ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos', 'ad' )
    
    ad1 =  0.0D0; ad2= 0.0D0

    call mpi_barrier(upt_comm,ierr)

    !=========================================================================
    if (verbose.gt.0) then
       write(*,*) 
       write(*,*) 'Memory required for vectors:', n_ham*5*16/1e6, 'MBytes'
       write(*,*) 'Memory required for matrix :', ldz*13*8/1e6, 'MBytes'
    endif   
    !=========================================================================

    !=========================================================================
    ! START WITH A RANDOM VECTOR
    !=========================================================================
    converged = .FALSE.
    nstep = min_step
    ALLOCATE( temp_v(SIZE(eigen_seed)), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos', 'lanczos', 'temp_v' )
    
    CALL RANDOM_NUMBER(temp_v)                          
    eigen_seed = (1.0D0,0.0D0)*temp_v                   
    CALL RANDOM_NUMBER(temp_v)                          
    eigen_seed = eigen_seed + (0.0D0,1.0D0)*temp_v      
    eigen_seed = eigen_seed/SQRT(DOT_PRODUCT(eigen_seed,eigen_seed))
    DEALLOCATE(temp_v)
    !=========================================================================

    call mpi_barrier(upt_comm,ierr)
    
    IF (nr_eigv .NE. 1) CALL project_out(eigen_seed, eigen_v, nr_eigv-1)
    
    call mpi_barrier(upt_comm,ierr)
    
    local_eigen_seed = eigen_seed
    call mpi_barrier(upt_comm,ierr)   

    !//////////////////////////////////////////////////////////////////////
    ! DEFINE MAXIMUM OVERLAP FOR ALL NODES
    !//////////////////////////////////////////////////////////////////////
    overlapR = 0; overlapL = 0
    do  j = 0,num_procs-1
       if (overlapR .lt. colind_high(j+1) - shift_end_Mi(j)) then
          overlapR = colind_high(j+1) - shift_end_Mi(j)
       endif

       !ovelap_low is j+1 bcoz its define with dynamic size here and it starts from 1 
       !but actually it starts from zero in the parent function
       if (overlapL .lt. shift_init_Mi(j) - colind_low(j+1)) then
          overlapL = shift_init_Mi(j) - colind_low(j+1) 
       end if 
    enddo
    !/////////////////////////////////////////////////////////////////////////


    call mpi_barrier(upt_comm,ierr)
    CALL MPI_Bcast(nstep, 1, MPI_INTEGER, 0, upt_comm,ierr)
    call mpi_barrier(upt_comm,ierr)

    if (verbose.gt.0) then
       WRITE(str_nr,'(I3.3)') nr_eigv     
       WRITE(*,*) 
       WRITE(*, '(a,2x,a,3x,a,5x,a,8x,a,11x,a)' ) &
            'type','total','#it','eigenvalue ('//str_nr//')','error','ort. error'
       WRITE(*,'(73("-"))')
    end if

    DO WHILE ( .NOT. converged )

       !======================================================================
       !
       ! Allocate work arrays :
       !
       ! The tridiagonal matrix T built on the vectors | j >
       ! has a nstep dimension, where nstep is the number of Lanczos loops.
       !
       !======================================================================
       !======================================================================
       ! Reset counters and arrays
       !======================================================================

       p_previous = ( 0.0D0, 0.0D0 )
       alpha = 0.0D0
       local_alpha = 0.0D0
       beta = 0.0D0
       local_beta = 0.0D0
       p_next = ( 0.0D0, 0.0D0 )
       p_aux = ( 0.0D0, 0.0D0 )
       p_aux2 = ( 0.0D0, 0.0D0 )
       p_current = local_eigen_seed
       test_fast = .FALSE.
       call mpi_barrier(upt_comm,ierr)
       
       !======================================================================
       ! Lanczos loop
       !======================================================================
       
       j = 0
       DO WHILE (j .LE. nstep) !nstep1  !small devices
          j = j + 1

          call sprs_ax( M, Mj, Mi, sparse_fmt, p_current, p_aux)
          
          local_alpha = DOT_PRODUCT( p_aux(shift_init:shift_end), &
               p_aux(shift_init:shift_end) ) 
          
          call mpi_allreduce(local_alpha, alpha, 1 , MPI_double_precision, &
               MPI_SUM, upt_comm,ierr)

          call exchange_vector(p_aux, overlapR, overlapL)

          call sprs_ax(M, Mj, Mi, sparse_fmt, p_aux, p_aux2)

          p_next(shift_init:shift_end) = &
               p_next(shift_init:shift_end) + &
               p_aux2(shift_init:shift_end)

          p_next(shift_init:shift_end) = &
               p_next(shift_init:shift_end) - &
               alpha * p_current(shift_init:shift_end)

          IF (nr_eigv .NE. 1) THEN
            CALL project_out(p_next, eigen_v, nr_eigv-1, 1)
          END IF

          local_beta = DOT_PRODUCT( p_next(shift_init:shift_end), &
               p_next(shift_init:shift_end) )

          call mpi_allreduce(local_beta, beta, 1 , MPI_double_precision, &
               MPI_SUM, upt_comm,ierr)

          beta = SQRT(beta)

          ! CHECK ORTHOGONALITY CONDITION
          IF ( j.EQ.nstep .AND. nstep.GT.2 .AND. nstep.LT.ldz - 10) THEN
             local_ort_chk = DOT_PRODUCT(p_next(shift_init:shift_end), &
                  local_eigen_seed(shift_init:shift_end))
             call mpi_allreduce(local_ort_chk, ort_chk, 1 , &
                  MPI_double_precision, MPI_SUM, upt_comm, ierr)
             ort_chk=abs(ort_chk/beta)
             IF (ort_chk.LT.ort_tol)   nstep = nstep + 10
          END IF


          if (id0) then
             ad2(j) = alpha
             ad1(j) = beta
          endif

          p_previous(shift_init:shift_end) = p_current(shift_init:shift_end)

          call exchange_vector(p_next, overlapR, overlapL)

          p_current = p_next/beta

          p_next(shift_init:shift_end) = -beta*(p_previous(shift_init:shift_end))

       END DO !END LOOP 1
          
       counter = counter + nstep

       !======================================================================
       ! Diagonalize the T matrix (ONLY NODE 0)
       !======================================================================

       if (id0) then
          
          CALL DSTEVX( JOBZ, RANGE, nstep, ad2, ad1, &
               VL, VU, IL, IL, ABSTOL, M_out, eigen_val, eigen_vec, ldz, &
               work, iwork, ifail, info )
          
          IF ( info .NE. 0 ) call throw_solve_exception(ERR_LANCZ_DIAG)
          
          !======================================================================
          ! End of loop test
          !======================================================================
          energy =  SQRT( ABS(eigen_val( 1 )) )
          energy_test( 1 : 9 ) = energy_test( 2 : 10 )
          energy_test( 10 ) = energy
          mean_en = SUM( energy_test ) / 10.0D0
          dev_en = SQRT( SUM( ( energy_test - mean_en )**2 ) / 10.0D0 )
          
          IF ( nstep .EQ. min_step ) THEN             
             IF ( dev_en / mean_en .LT. fast_tol ) test_fast = .TRUE.             
          END IF
       endif

       call mpi_barrier(upt_comm,ierr)
       !======================================================================
       ! Convert ground state vector to the full basis
       !======================================================================
       
       call mpi_bcast(test_fast,1,MPI_LOGICAL, 0, upt_comm,ierr)
       call mpi_bcast(eigen_vec(1,1),ldz,MPI_DOUBLE_COMPLEX, 0, upt_comm,ierr)
       
       alpha = 0.0D0
       beta = 0.0D0
       local_alpha =  0.0D0
       local_beta = 0.0D0
       p_previous = ( 0.0D0, 0.0D0 )
       p_next = ( 0.0D0, 0.0D0 )
       p_aux = ( 0.0D0, 0.0D0 )
       p_aux2 = ( 0.0D0, 0.0D0 )
       p_current = local_eigen_seed
       local_eigen_seed = ( 0.0D0, 0.0D0 )
       
       DO j = 1, nstep !nstep2
          
          call sprs_ax( M, Mj, Mi, sparse_fmt, p_current, p_aux)
          
          local_alpha = DOT_PRODUCT( p_aux(shift_init:shift_end), &
               p_aux(shift_init:shift_end))
          
          call mpi_allreduce(local_alpha, alpha, 1 , MPI_double_precision, &
               MPI_SUM,upt_comm,ierr)

          call exchange_vector(p_aux, overlapR, overlapL)
        
          call sprs_ax(M, Mj, Mi, sparse_fmt, p_aux, p_aux2)
        
          p_next(shift_init:shift_end) = p_next(shift_init:shift_end) + &
               p_aux2(shift_init:shift_end)
          
          p_next(shift_init:shift_end) = &
               p_next(shift_init:shift_end) - &
               alpha * p_current(shift_init:shift_end)          

          IF (nr_eigv .NE. 1) THEN
            CALL project_out(p_next, eigen_v, nr_eigv-1, 1)
          END IF

          local_beta = DOT_PRODUCT( p_next(shift_init:shift_end), &
               p_next(shift_init:shift_end) )

          call mpi_allreduce(local_beta, beta, 1 , MPI_double_precision, MPI_SUM, upt_comm,ierr)

          beta = SQRT(beta)

          local_eigen_seed(shift_init:shift_end) = &
               local_eigen_seed(shift_init:shift_end) &
               + (eigen_vec( j, 1 ) * p_current(shift_init:shift_end))
        
          p_previous(shift_init:shift_end) = p_current(shift_init:shift_end)

          call exchange_vector(p_next, overlapR, overlapL)

          p_current = p_next/beta
        
          p_next(shift_init:shift_end) = -beta* (p_previous(shift_init:shift_end))

       END DO !END LOOP 2

       ! Project out previous eigen vectors
       IF (nr_eigv .NE. 1) THEN
         CALL project_out(local_eigen_seed, eigen_v, nr_eigv-1, 1)
       END IF
     
       !======================================================================
       ! Normalize eigenvector
       !======================================================================
       local_beta = 0.0D0
       beta = 0.0D0
       test_E = 0.0D0
       
       local_beta = DOT_PRODUCT( local_eigen_seed(shift_init:shift_end), &
            local_eigen_seed(shift_init:shift_end) )
       
       call mpi_allreduce(local_beta, beta, 1 , MPI_double_precision, &
            MPI_SUM,upt_comm,ierr)
 
       local_eigen_seed(shift_init:shift_end) = &
            local_eigen_seed(shift_init:shift_end) / SQRT(beta)

       call exchange_vector(local_eigen_seed, overlapR, overlapL)

       !======================================================================
       ! Test vector < psi | H | psi > / <psi | psi>    <u|H|u> = E
       !======================================================================

       call sprs_ax( M, Mj, Mi, sparse_fmt, local_eigen_seed, p_aux)

       local_beta = DOT_PRODUCT(local_eigen_seed(shift_init:shift_end), &
            p_aux(shift_init:shift_end) )

       call mpi_allreduce(local_beta, test_E, 1 , MPI_double_precision, &
            MPI_SUM,upt_comm,ierr)

       call mpi_bcast(energy,1,MPI_double_precision,0,upt_comm,ierr)

       deltaE =  ABS( energy - ABS(test_E))  
    
       !======================================================================
       ! perform some output
       !======================================================================       
       IF (verbose.gt.0) THEN
          
          IF ( nstep .EQ. min_step ) THEN
             WRITE ( *, '( "fast", I7, I6, 2g20.10, 1g14.5 )' ) &
                  counter, nstep, test_E, deltaE, ort_chk
          ELSE
             WRITE ( *, '( "stan", I7, I6, 2g20.10, 1g14.5 )' ) &
                  counter, nstep, test_E, deltaE, ort_chk
          END IF
          
       END IF
   
       !======================================================================
       ! CONVERGENCE TEST 
       !======================================================================
       IF ( deltaE .LT. long_tol) THEN
          converged = .TRUE. ! Eigenvalue has converged to tolerance
          energy = test_E
       ELSE
          converged = .FALSE.
       END IF
      
       !======================================================================
       ! BREAK CONDITION ON STUCKED LOOPS 
       !======================================================================
       IF (counter > max_step) THEN 
         energy = test_E
         converged = .TRUE.
       END IF

       !======================================================================
       ! FAST/SLOW ITERATION SWITCH TEST 
       !======================================================================
       IF ( nstep .EQ. min_step ) THEN
          IF ( test_fast ) THEN
             test_fast = .FALSE.
             nstep = long_step
          END IF
       ELSE
          nstep = min_step
       END IF

    END DO ! LOOP TILL CONVERGENCE

    !======================================================================
    ! Deallocate work arrays
    !======================================================================
    if (id0) then
       
       DEALLOCATE( iwork, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'lanczos_diag', 'lanczos', 'iwork' )
       
       DEALLOCATE( ifail, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'lanczos_diag', 'lanczos', 'ifail' )
       
       DEALLOCATE( work, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'lanczos_diag', 'lanczos', 'work' )
       
       DEALLOCATE( eigen_val, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'lanczos_diag', 'lanczos', 'eigen_val' )
       
       DEALLOCATE( ad1, STAT = err )
       DEALLOCATE( ad2, STAT = err )
       IF ( err .NE. 0 ) &
            CALL dealloc_error( 'lanczos_diag', 'lanczos', 'ad' )
       
    endif
         
    !======================================================================
    ! Finally gather eigenvector on master and bcast to all
    !======================================================================  
    call gather_on_master(local_eigen_seed)
  
    if (id0) eigen_seed = local_eigen_seed

    call mpi_barrier(upt_comm,ierr)

    call MPI_Bcast(eigen_seed, n_ham, MPI_DOUBLE_COMPLEX, 0, upt_comm, ierr)
    !====================================================================== 
    
    DEALLOCATE( eigen_vec, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'eigen_vec' )
    
    DEALLOCATE( local_eigen_seed, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'local_eigen_seed' )
    
    DEALLOCATE( j_aux, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'j_aux' )

    DEALLOCATE( j_aux2, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'j_aux2' )

    DEALLOCATE( j_current, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'j_current' )
    
    DEALLOCATE( j_next, STAT = err )
    IF ( err .NE. 0 ) &
        CALL dealloc_error( 'lanczos_diag', 'lanczos', 'j_next' )

    DEALLOCATE( j_previous, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'j_previous' )
    
    call mpi_barrier(upt_comm,ierr)
    
  END SUBROUTINE fast_lanczos_ev_pold_opt

#endif  

! ==========================================================================

#ifdef UPT_MPI

  SUBROUTINE project_out_z_p(eigen_seed, eigen_vectors, nr_ev, parallel)
    
    !=========================================================================
    ! Input arguments
    !=========================================================================
    COMPLEX ( dp ), DIMENSION(:) :: eigen_seed
    COMPLEX ( dp ), DIMENSION(:,:) :: eigen_vectors
    INTEGER :: nr_ev, parallel
    
    !=========================================================================
    ! Local variables
    !=========================================================================
    INTEGER i
    COMPLEX (dp) local_dot_p, dot_p
    
    DO i = 1, nr_ev

       call mpi_barrier(upt_comm,ierr)
       local_dot_p = DOT_PRODUCT( eigen_vectors(shift_init:shift_end,i), &
       eigen_seed(shift_init:shift_end) )
       call mpi_barrier(upt_comm,ierr)
       call mpi_allreduce(local_dot_p, dot_p, 1 , MPI_double_complex,MPI_SUM,upt_comm,ierr)
       call mpi_barrier(upt_comm,ierr)
       
       eigen_seed = eigen_seed - dot_p * eigen_vectors(:,i)
       call mpi_barrier(upt_comm,ierr)

    END DO
    
  END SUBROUTINE project_out_z_p

#else

  SUBROUTINE project_out_z_p(eigen_seed, eigen_vectors, nr_ev, parallel)
    COMPLEX ( dp ), DIMENSION(:) :: eigen_seed
    COMPLEX ( dp ), DIMENSION(:,:) :: eigen_vectors
    INTEGER :: nr_ev, parallel
  END SUBROUTINE project_out_z_p

#endif
! ============================================================================
! OLD (ORIGINAL) LANCZOS ROUTINE BY JEROME/MARTIN  (readapted)
!
!============================================================================= 
  SUBROUTINE fast_lanczos_ev_old( M, Mi, Mj,  sparse_fmt, &
        min_step, long_step, max_step, eigen_seed, n_ham, &
        fast_tol, long_tol, ort_tol_in, counter, res_flag, eigen_v, nr_eigv, &
        energy, deltaE, verbose)


    !=========================================================================
    ! Input arguments
    !=========================================================================
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mi  !ROWPNT
    INTEGER,        DIMENSION(:), POINTER :: Mj  !COLIND
    CHARACTER(1)                  :: sparse_fmt

    INTEGER, INTENT( IN ) :: min_step   ! # of iterations of fast lanczos
    INTEGER, INTENT( IN ) :: long_step  ! # of iterations of long lanczos  
    INTEGER, INTENT( IN ) :: max_step   ! # max number of iterations 
    LOGICAL, INTENT( IN ) :: res_flag
    INTEGER, INTENT( IN ) :: nr_eigv    ! eigenvalue number

    INTEGER, INTENT( IN ) :: n_ham
    REAL ( dp ), INTENT( IN ) :: fast_tol, long_tol, ort_tol_in
    INTEGER, INTENT( IN ) :: verbose
    !=========================================================================
    ! Input - output arguments
    !=========================================================================

    COMPLEX( dp ), DIMENSION( : ), INTENT( INOUT ) :: eigen_seed
    INTEGER                      , INTENT( INOUT ) :: counter
    COMPLEX( dp ), DIMENSION( :,: ), INTENT( IN )  :: eigen_v    
    !=========================================================================
    ! Output arguments
    !=========================================================================

    REAL ( dp ), INTENT( OUT ) :: energy
    REAL ( dp ), INTENT( OUT ) :: deltaE
    
    !=========================================================================
    ! Local variables
    !=========================================================================
    REAL ( dp ), ALLOCATABLE, DIMENSION( : )    :: temp_v
    REAL ( dp ), POINTER    , DIMENSION( :, : ) :: eigen_vec
    REAL ( dp ), POINTER    , DIMENSION( : )    :: eigen_val

    REAL ( dp ), ALLOCATABLE, DIMENSION( : )    :: ad1,ad2
    REAL ( dp ), ALLOCATABLE, DIMENSION( : )    :: work
    INTEGER, ALLOCATABLE, DIMENSION( : )    :: iwork, ifail

    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE, TARGET    :: j_current
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE,TARGET    :: j_previous
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE, TARGET    :: j_next
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE, TARGET    :: j_aux
    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE,TARGET    :: j_aux2

    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_current
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_previous
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_next
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_aux
    COMPLEX ( dp ), DIMENSION( : ), POINTER    :: p_aux2
    real( dp ), DIMENSION(10) :: energy_test

    REAL ( dp ) :: alpha, beta, test_E, ort_chk, ort_tol
    COMPLEX ( dp ) :: eigen_guess

    REAL ( dp )    :: tol, errors, ABSTOL
    REAL ( dp )    :: dev_en, mean_en

    INTEGER :: j, nstep, info, ldz, VL, VU, IL, IU, M_out, err
    CHARACTER :: JOBZ, RANGE
    LOGICAL :: test_fast, converged 

    !===========================================================================

    ALLOCATE( j_current( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_current' )
    ALLOCATE( j_previous( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_previous' )
    ALLOCATE( j_next( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_next' )
    ALLOCATE( j_aux( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_aux' )
    ALLOCATE( j_aux2( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'j_aux2' )

    p_current => j_current
    p_previous => j_previous
    p_next => j_next
    p_aux => j_aux
    p_aux2 => j_aux2

    JOBZ = 'V'
    RANGE = 'I'
    IL = 1
    IU = 1
    ABSTOL = 1e-12
    ldz = max_step
    ort_tol = ort_tol_in
    ort_chk = 0.d0
    energy_test = 0.d0

    !======================================================================
    ! Allocate work arrays for eigen values and vectors
    !======================================================================
    ALLOCATE( eigen_val( 2 ), STAT = err )
    IF ( err .NE. 0 ) &
         CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'eigen_val' )

    ! eigen_val = 0.0D0

    ALLOCATE( eigen_vec( ldz, 2 ), STAT = err )
    IF ( err .NE. 0 ) &
         CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'eigen_vec' )

    ! eigen_vec = ( 0.0D0, 0.0D0 )

    ALLOCATE( work( 5*ldz ), STAT = err )
    IF ( err .NE. 0 ) &
         CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'work' )

    ! work = ( 0.0D0, 0.0D0 )

    ALLOCATE( iwork( 5 * ldz ), STAT = err )
    IF ( err .NE. 0 ) &
         CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'iwork' )

    ! rwork = 0.0D0

    ALLOCATE( ifail( ldz ), STAT = err )
    IF ( err .NE. 0 ) &
         CALL alloc_error( 'lanczos_diag', 'lanczos_ev', 'ifail' )

    !----------------------------------------------------------------------
    ! Band form of T matrix 
    !----------------------------------------------------------------------
    ALLOCATE( ad1( ldz ), STAT = err )
    ALLOCATE( ad2( ldz ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'lanczos_diag', 'lanczos', 'ad' )

    ad1 =  0.0D0; ad2= 0.0D0 

    !=========================================================================
    write(*,*) 
    write(*,*) 'Memory required for lanczos vectors:', n_ham*5*16, 'bytes' 
    write(*,*) 'Memory required for lanczos matrix :',   ldz*13*8, 'bytes' 
    !=========================================================================
    converged = .FALSE.

    nstep = min_step

    ! CREATE A STARTING RANDOM SEED
    ALLOCATE( temp_v(SIZE(eigen_seed)), STAT = err )
       IF ( err .NE. 0 ) CALL alloc_error( 'lanczos', 'lanczos', 'temp_v' )

    CALL RANDOM_NUMBER(temp_v)
    eigen_seed = (1.0D0,0.0D0)*temp_v
    CALL RANDOM_NUMBER(temp_v)
    eigen_seed = eigen_seed + (0.0D0,1.0D0)*temp_v
    eigen_seed = eigen_seed/SQRT(DOT_PRODUCT(eigen_seed,eigen_seed))

    DEALLOCATE(temp_v)

    IF (nr_eigv .NE. 1) CALL project_out(eigen_seed, eigen_v, nr_eigv-1)

    WRITE(*,*) 
    WRITE(*, '(a,2x,a,3x,a,5x,a,13x,a,11x,a)' ) &
          'type','total','#it','eigenvalue','error','ort. error'
    WRITE(*,'(73("-"))')      
        
    DO WHILE ( .NOT. converged )

       !======================================================================
       !
       ! Allocate work arrays :
       !
       ! The tridiagonal matrix T built on the vectors | j >
       ! has a nstep dimension, where nstep is the number of Lanczos loops.
       !
       !======================================================================

       !======================================================================
       ! Reset counters and arrays
       !======================================================================

       alpha = 0.0D0
       beta = 0.0D0

       p_previous = ( 0.0D0, 0.0D0 )
       p_next = ( 0.0D0, 0.0D0 )
       p_aux = ( 0.0D0, 0.0D0 )

       p_current = eigen_seed

       test_fast = .FALSE.

       !======================================================================
       ! Lanczos loop
       !======================================================================
       j = 0
       DO WHILE (j .LE. nstep )
          j = j + 1

          !-------------------------------------------------------------------
          ! update | J + 1 > to  | J + 1 > + A * A | j >
          !-------------------------------------------------------------------
          !call message_clock("M*v")
          call sprs_ax( M, Mj, Mi, sparse_fmt, p_current, p_aux )
          call sprs_ax( M, Mj, Mi, sparse_fmt, p_aux, p_aux2 )
          !call write_clock()

          p_next = p_next + p_aux2

          !-------------------------------------------------------------------
          ! calculate  < j | A * A | j >
          !-------------------------------------------------------------------
          alpha = ( DOT_PRODUCT( p_aux, p_aux ) )

          !-------------------------------------------------------------------
          ! update | J + 1 > to | J + 1 > - < j | A * A | j > | j >
          !-------------------------------------------------------------------

          p_next = p_next - alpha * p_current

          ! Project out previous eigen vectors 
          IF (nr_eigv .NE. 1)  CALL project_out(p_next, eigen_v, nr_eigv-1)
          !-------------------------------------------------------------------
          ! calculate  SQRT( < J + 1 | J + 1 > )
          !-------------------------------------------------------------------

          beta = SQRT( DOT_PRODUCT( p_next, p_next ) )

          !-------------------------------------------------------------------
          ! Check orthogonality of p_next against eigen_seed and increase 
          !    n_step still below treach hold
          !-------------------------------------------------------------------
          IF ( MOD(j,min_step).eq.0 ) THEN 
             
             ort_chk = ABS(DOT_PRODUCT(p_next,eigen_seed)/beta)
             IF (ort_chk.LT.ort_tol) nstep = nstep + min_step
             IF (nstep.GE.max_step) nstep = max_step
                
          END IF
          !-------------------------------------------------------------------
          ! Store T matrix
          !-------------------------------------------------------------------

          ad2( j ) = alpha
          ad1( j ) = beta

          !----------------------------------------------------------------
          ! update | j - 1 > to | j >
          !----------------------------------------------------------------
          p_previous = p_current
          
          !----------------------------------------------------------------
          ! update | j > to | J + 1 > / SQRT( < J + 1 | J + 1 > )
          !----------------------------------------------------------------
          p_current = p_next / beta
          
          !----------------------------------------------------------------
          ! initialize | J + 1 >   to   - < J | J > | j - 1 > 
          !----------------------------------------------------------------
          p_next = - p_previous * beta

       END DO

       counter = counter + nstep
       !======================================================================
       ! Diagonalize the T matrix
       !======================================================================
       CALL DSTEVX( JOBZ, RANGE, nstep, ad2, ad1, &
              VL, VU, IL, IL, ABSTOL, M_out, eigen_val, eigen_vec, ldz, &
              work, iwork, ifail, info )

       IF ( info .NE. 0 ) call throw_solve_exception(ERR_LANCZ_DIAG)

       !======================================================================
       ! Compute rmsq deviations of the last 10 fast iterations 
       !======================================================================
       energy =  SQRT( ABS(eigen_val( 1 )) )
       energy_test( 1 : 9 ) = energy_test( 2 : 10 )
       energy_test( 10 ) = energy
       mean_en = SUM( energy_test ) / 10.0D0
       dev_en = SQRT(SUM(( energy_test - mean_en )**2 )/ 10.0D0) 
        
       IF ( nstep .EQ. min_step ) THEN
          IF ( dev_en/mean_en .LT. fast_tol ) test_fast = .TRUE.
       END IF
       
       !dev_en =  ABS(energy_test(10)-energy_test(9)) 
       !======================================================================
       ! Convert ground state vector to the full basis
       !======================================================================

       p_current = eigen_seed

       p_previous = ( 0.0D0, 0.0D0 )
       p_next = ( 0.0D0, 0.0D0 )
       p_aux = ( 0.0D0, 0.0D0 )

       eigen_seed = ( 0.0D0, 0.0D0 )

       alpha = ( 0.0D0, 0.0D0 )
       beta = ( 0.0D0, 0.0D0 )

       DO j = 1, nstep

          !-------------------------------------------------------------------
          ! update | J + 1 > to  | J + 1 > + A | j >
          !-------------------------------------------------------------------
          call sprs_ax( M, Mj, Mi, sparse_fmt, p_current, p_aux ) 
          call sprs_ax( M, Mj, Mi, sparse_fmt, p_aux, p_aux2 )  

          p_next = p_next + p_aux2

          !-------------------------------------------------------------------
          ! calculate  < j | A | j >
          !-------------------------------------------------------------------
          alpha = DOT_PRODUCT( p_aux, p_aux )

          !-------------------------------------------------------------------
          ! update | J + 1 > to | J + 1 > - < j | A | j > | j >
          !-------------------------------------------------------------------

          p_next = p_next - alpha * j_current
          ! Project out previous eigen vectors 
          IF (nr_eigv .NE. 1)  CALL project_out(p_next, eigen_v, nr_eigv-1)
          !-------------------------------------------------------------------
          ! calculate  SQRT( < J + 1 | J + 1 > )
          !-------------------------------------------------------------------

          beta = SQRT( DOT_PRODUCT( p_next, p_next ) )

          eigen_seed = eigen_seed + eigen_vec( j, 1 ) * p_current

          !----------------------------------------------------------------
          ! update | j - 1 > to | j >
          !----------------------------------------------------------------
          p_previous = p_current
          
          !----------------------------------------------------------------
          ! update | j > to | J + 1 > / SQRT( < J + 1 | J + 1 > )
          !----------------------------------------------------------------
          p_current = p_next / beta
          
          !----------------------------------------------------------------
          ! initialize | J + 1 >   to   - < J | J > | j - 1 > 
          !----------------------------------------------------------------
          p_next = - p_previous * beta


       END DO
       ! Project out previous eigen vectors 
       IF (nr_eigv .NE. 1)  CALL project_out(eigen_seed, eigen_v, nr_eigv-1) 

       !======================================================================
       ! Normalize eigenvector
       !======================================================================
       eigen_seed = eigen_seed / SQRT( DOT_PRODUCT( eigen_seed, eigen_seed ) )
       !======================================================================
       ! Computes the expectation value of the eigenvector: < x | A | x >  
       !======================================================================
       test_E = test_eigvect(M, Mj, Mi, sparse_fmt, eigen_seed)
       !======================================================================
       deltaE =  ABS(energy - ABS(test_E))
       !======================================================================
       ! perform some output
       !======================================================================

       IF ( nstep .EQ. min_step ) THEN

          WRITE ( *, '( "fast", I7, I6, 2g20.10, 1g14.5 )' ) &
               counter, nstep, test_E, deltaE, ort_chk

       ELSE

          WRITE ( *, '( "stan", I7, I6, 2g20.10, 1g14.5 )' ) &
               counter, nstep, test_E, deltaE, ort_chk

       END IF

       IF ( deltaE .LT. long_tol .OR. dev_en .LT. long_tol) THEN
 
          converged = .TRUE. ! Eigenvalue has converged to tolerance
          energy = test_E
          !WRITE (*,*) 'Eig value = ', energy, 'Error = ', deltaE
          
       ELSE
          
          converged = .FALSE.
          
       END IF
       !======================================================================
       ! End of loop condition
       !======================================================================
       
       IF ( nstep .EQ. min_step ) THEN

          IF ( test_fast ) THEN

             test_fast = .FALSE.
             nstep = long_step

          END IF

       ELSE

          nstep = min_step

       END IF

    END DO
    !======================================================================
    ! Deallocate work arrays
    !======================================================================

    DEALLOCATE( iwork, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'iwork' )

    DEALLOCATE( ifail, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'ifail' )

    DEALLOCATE( work, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'work' )

    DEALLOCATE( eigen_vec, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'eigen_vec' )

    DEALLOCATE( eigen_val, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'eigen_val' )
       
    DEALLOCATE( ad1, STAT = err )
    DEALLOCATE( ad2, STAT = err )
    IF ( err .NE. 0 ) &
         CALL dealloc_error( 'lanczos_diag', 'lanczos', 'ad' )  
    
    !=========================================================================

  END SUBROUTINE fast_lanczos_ev_old

! ==========================================================================
  
  SUBROUTINE project_out_z(eigen_seed, eigen_vectors, nr_ev)
    
!=========================================================================
! Input arguments
!=========================================================================
    COMPLEX ( dp ), DIMENSION(:) :: eigen_seed
    COMPLEX ( dp ), DIMENSION(:,:) :: eigen_vectors
    INTEGER :: nr_ev
    
!=========================================================================
! Local variables
!=========================================================================
    INTEGER i
    
    DO i = 1, nr_ev
       eigen_seed = eigen_seed - DOT_PRODUCT(eigen_vectors(:,i),eigen_seed)* &
                                                      eigen_vectors(:,i)
    END DO
    
  END SUBROUTINE project_out_z










  SUBROUTINE project_out_r(eigen_seed, eigen_vectors, nr_ev)
    
!=========================================================================
! Input arguments
!=========================================================================
    REAL ( dp ), DIMENSION(:) :: eigen_seed
    REAL ( dp ), DIMENSION(:,:) :: eigen_vectors
    INTEGER :: nr_ev
    
!=========================================================================
! Local variables
!=========================================================================
    INTEGER i
    
    DO i = 1, nr_ev
       eigen_seed = eigen_seed - DOT_PRODUCT(eigen_vectors(:,i),eigen_seed)* &
                                                      eigen_vectors(:,i)
    END DO
    
  END SUBROUTINE project_out_r


  SUBROUTINE project_out_sp(eigen_seed, eigen_vectors, nr_ev)
    
!=========================================================================
! Input arguments
!=========================================================================
    COMPLEX ( sp ), DIMENSION(:) :: eigen_seed
    COMPLEX ( sp ), DIMENSION(:,:) :: eigen_vectors
    INTEGER :: nr_ev
    
!=========================================================================
! Local variables
!=========================================================================
    INTEGER i
    
    DO i = 1, nr_ev
       eigen_seed = eigen_seed - DOT_PRODUCT(eigen_vectors(:,i),eigen_seed)* &
                                                      eigen_vectors(:,i)
    END DO
    
  END SUBROUTINE project_out_sp



#ifdef UPT_MPI

  SUBROUTINE exchange_vector(vect, overlapR, overlapL)
    
    complex(dp), DIMENSION(:) :: vect
    integer ::  overlapR, overlapL
    
    integer :: reqs(num_procs+1)
    integer :: rank
    integer :: status(MPI_STATUS_SIZE)

    if (num_procs.gt.1) then
       do rank =0, num_procs-2
          if(id .eq. rank+1) then
             call MPI_irecv(vect(shift_end_Mi(rank)-overlapL+1), &
                  overlapL,MPI_DOUBLE_COMPLEX, &
                  rank,877+rank, upt_comm, reqs(rank+1), ierr)
             call MPI_wait(reqs(rank+1), status, ierr)
          else if(id .eq. rank) then
             call MPI_isend(vect(shift_end_Mi(rank)-overlapL+1), &
                  overlapL,MPI_DOUBLE_COMPLEX,&
                  rank+1,877+rank, upt_comm, reqs(rank+1), ierr)
          end if
       enddo
       
       call mpi_barrier(upt_comm,ierr)
       
       do rank =0, num_procs-2                
          if(id .eq. rank) then
             call MPI_irecv(vect(shift_init_Mi(rank+1)), overlapR, &
                  MPI_DOUBLE_COMPLEX,&
                  rank+1,978+rank, upt_comm, reqs(rank+1), ierr)
             call MPI_wait(reqs(rank+1), status, ierr)
          else if(id .eq. rank+1) then
             call MPI_isend(vect(shift_init_Mi(rank+1)), overlapR, &
                  MPI_DOUBLE_COMPLEX,&
                  rank,978+rank, upt_comm, reqs(rank+1), ierr)
          end if
       enddo
       
       call mpi_barrier(upt_comm,ierr)
    endif
    
  END SUBROUTINE exchange_vector

  
  SUBROUTINE gather_on_master(vect)

    complex(dp), DIMENSION(:) :: vect

    integer :: reqs(num_procs+1)
    integer :: rank
    integer :: status(MPI_STATUS_SIZE)

    if (num_procs.gt.1) then
       do rank =1, num_procs-1
          !if ( id.lt.num_procs-1) then
          if(id .eq. 0) then
             call MPI_irecv(vect(shift_init_Mi(rank)),&
                  shift_end_Mi(rank)-shift_init_Mi(rank)+1, & 
                  MPI_DOUBLE_COMPLEX, rank, 111+rank, upt_comm, reqs(rank+1), ierr)
          end if
          if(id .ne. 0 .and. id .eq. rank) then
             call MPI_isend(vect(shift_init_Mi(rank)), &
                  shift_end_Mi(rank)-shift_init_Mi(rank)+1, &
                  MPI_DOUBLE_COMPLEX, 0, 111+rank, upt_comm, reqs(rank+1), ierr)
             call MPI_wait(reqs(rank+1), status, ierr)
          end if
       enddo

       call mpi_barrier(upt_comm,ierr)

    endif

  END SUBROUTINE gather_on_master

#endif

END MODULE LANCZOS_DIAG
