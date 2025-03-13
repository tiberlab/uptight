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
!                              JD SOLVER
!
!=============================================================================
!
! Walter Rodrigues
! 
! Dipartimento di Ingegneria Elettronica
! Universita` di Roma "Tor Vergata"
! Tel. +39-6-72597781
! E-mail: dicarlo@ing.uniroma2.it
!         pecchia@ing.uniroma2.it
!     
! 13-04-2013
!
!=============================================================================

!=============================================================================

MODULE JD_DIAG

!===========================================================================
  USE mpi_globals
  USE precision
  USE errors
  USE exceptions
  USE input_output
  USE sparse_matrix
  USE sparse_numrec
  USE sort               !Collection of sorting routines
  USE clock
  USE omp_lib 
!===========================================================================

  IMPLICIT NONE
  PRIVATE

!===========================================================================
  
  PUBLIC JD_EV

CONTAINS

!===========================================================================
!
! Subroutines - by Walter in 2015 for MPI/OpenMP/CUDA
!
! Driver for JD algorithm
!
!===========================================================================
SUBROUTINE JD_EV(H, U, n_spin, min_step, long_step, max_step, &
                 energies, eigen_vectors, start_ev, num_ev, n_ham, &
                 lambda, solver_flag_jd, fast_tol, long_tol, ort_tol, &
                 sign, dynamic, spin_deg, verbose,  band_type )

    !=========================================================================
    ! Input arguments
    !=========================================================================
    TYPE(CSR)                     :: H, U
    INTEGER                       :: n_spin
    INTEGER,         INTENT( IN ) :: min_step, long_step, max_step
    INTEGER,         INTENT( IN ) :: start_ev, num_ev, n_ham, band_type
    REAL ( dp ),     INTENT( IN ) :: lambda
    INTEGER,         INTENT( IN ) :: sign
    INTEGER,         INTENT( IN ) :: solver_flag_jd
    REAL ( dp ),     INTENT( IN ) :: fast_tol, long_tol, ort_tol
    LOGICAL,         INTENT( IN ) :: dynamic
    LOGICAL,         INTENT( IN ) :: spin_deg
    INTEGER,         INTENT( IN ) :: verbose
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

    TYPE(CSR_real_sp_mx_prec)   :: H_real
    TYPE(CSR_real_sp_mx_prec)   :: H_imag

    REAL ( dp )    ::  jd_tol
    INTEGER :: jd_min, jd_max, ls_restart, ls_max
    REAL (dp) :: ls_tol
    !Complex (dp), DIMENSION(:), POINTER :: eigen_vec_out
    !Real (dp), DIMENSION(:), POINTER :: eigen_val_out
    !INTEGER, DIMENSION( : ), POINTER :: row_offset

    ALLOCATE( eigen_seed( n_ham ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'jd_diag', 'jd_ev', 'eigen_seed' )


    jd_min = num_ev+4
    jd_max = num_ev+10
    ls_tol = 0.1
    ls_max = 4
    ls_restart = 10
    jd_tol = long_tol

    res_flag = .false.
    shift = lambda
    !fast_tol_sp = fast_tol
    !long_tol_sp = long_tol
    !ort_tol_sp = ort_tol

#ifdef UPT_MPI 
    if (id .eq. 0) write(*,*) " "
    !if (id .eq. 0) write(*,*) '(pLanczos) Tolerance: ',long_tol_sp
    if (id .eq. 0) write(*,*) " JD tol =", jd_tol
#else
    write(*,*) " "
    !write(*,*) '(pLanczos) Tolerance: ',long_tol_sp
    write(*,*) " JD tol =", jd_tol
#endif

    SELECT CASE (solver_flag_jd)
    CASE(0)
      if(H%sparse_fmt.eq.'F') then
         if (verbose.gt.0) THEN
            write(*,*) ' '
            write(*,*) 'MPI parallelized on ',num_procs,'processes'
            write(*,*) 'OpenMP parallelized on',omp_get_max_threads(),'threads'
         end if
       !call set_clock()
       call split_matrix_sp_mx_prec(H, H_real, H_imag)
       !call write_clock
      endif
#ifdef __CUDA
    CASE(2)
       write(*,*) 'GPU accelerated'
       call split_matrix_sp_mx_prec(H, H_real, H_imag)
    CASE(1)
       if(id == 0) write(*,*) 'GPU accelerated. Splitting matrix parallel'
       !call set_clock()
       call split_matrix_sp_mx_prec(H, H_real, H_imag)
       !call write_clock

#endif
     CASE DEFAULT
       WRITE(*,*) 'ERROR Undefined solver flag',solver_flag_jd
       call throw_solve_exception(ERR_JD_DIAG)
    END SELECT 

#ifdef UPT_MPI   
    !////////////////////////////////////////////////////////////////////// 
    ! COMPUTE OVERLAP between nodes
    !//////////////////////////////////////////////////////////////////////

    call compute_ovr(H%M, H%Mj, H%Mi, colind_l, colind_h)

    call mpi_gather(colind_l, 1, MPI_INTEGER, colind_low, 1, &
         MPI_INTEGER, 0, upt_comm, ierr) 
    
    call mpi_gather(colind_h, 1, MPI_INTEGER, colind_high, 1, &
         MPI_INTEGER, 0, upt_comm, ierr) 
    
    call MPI_Bcast(colind_low, num_procs, MPI_INTEGER, 0, upt_comm,ierr);
    
    call MPI_Bcast(colind_high, num_procs, MPI_INTEGER, 0, upt_comm,ierr);


#endif

    if (verbose.gt.0) call set_clock()

    SELECT CASE (solver_flag_jd)
    CASE(0)
    
#ifdef UPT_MPI

    CALL jd_cpu_no_pc_split_mxprec_pal(n_ham, size(H%M), H_imag%nnz, H_real%M, H_real%Mi,&
                                       H_real%Mj, H_imag%M, H_imag%Mi, H_imag%Mj, H%sparse_fmt, &
                                       band_type, jd_tol, shift, jd_min, jd_max, num_ev, &
                                       energies, eigen_vectors, ls_tol, ls_restart, ls_max, &
                                       colind_low, colind_high, shift_init, shift_end, num_procs, id,  upt_comm)

    !call mpi_barrier(upt_comm,ierr)

#endif


       
      CASE(1)

#if ( defined __CUDA && defined UPT_MPI )
       
      CALL jd_single_gpu_no_pc_split_mxprec_pal(n_ham, size(H%M), H_imag%nnz, H_real%M, H_real%Mi,&
                                            H_real%Mj, H_imag%M, H_imag%Mi, H_imag%Mj, H%sparse_fmt, &
                                            band_type, jd_tol, shift, jd_min, jd_max, num_ev, &
                                            energies, eigen_vectors, ls_tol, ls_restart, ls_max, &
                                            colind_low, colind_high, shift_init, shift_end, num_procs, id,  upt_comm)
#endif

#ifdef __CUDA

      CASE(2)

      CALL jd_single_gpu_no_pc_split_mxprec(n_ham, size(H%M), size(H_imag%M), H_real%M, H_real%Mi,&
                                            H_real%Mj, H_imag%M, H_imag%Mi, H_imag%Mj, H%sparse_fmt, &
                                            band_type, jd_tol, shift, jd_min, jd_max, num_ev, &
                                            energies, eigen_vectors, ls_tol, ls_restart, ls_max)

#endif
      CASE DEFAULT
            WRITE(*,*) 'ERROR Undefined solver flag',solver_flag_jd
            call throw_solve_exception(ERR_JD_DIAG)
      END SELECT

      if (verbose.gt.0) call write_clock()

      call destroy_matrix(H_real)
      call destroy_matrix(H_imag)
      !DEALLOCATE(eigen_vec_out, STAT = err)
      DEALLOCATE(eigen_seed, STAT = err)
      !DEALLOCATE(eigen_val_out, STAT = err)

END SUBROUTINE JD_EV

END MODULE JD_DIAG
