! cpa_multi_solver.f90
!
! N-component Soven self-consistency and spectral function for multi-component alloy.
!
! Generalizes cpa_solver.f90 to:
!   - TWO self-energies: sigma1(n_ref_states) for sl1, sigma2 for sl2
!   - n_pseudo1 = M*N pseudo-elements per sublattice (general N-component CPA)
!   - The BZ-sum F_alpha is computed separately for each sublattice
!
! Soven equation per orbital alpha, sublattice s (=1 or 2):
!
!   F_alpha^s(z) = (1/N_k) SUM_k [ z*I - HK(k) ]^{-1}_{sl_s,alpha,alpha}
!                  (spin-averaged: mean of up and down diagonal elements)
!
!   DO i = 1, n_pseudo_s
!      v_i = pseudo_s(i)%onsite(alpha) - sigma_s_old(alpha)
!      t_i = v_i / (1 - F_alpha^s * v_i)
!   END DO
!
!   t_avg^s = SUM_i conc_i * t_i
!
!   sigma_s_new(alpha) = sigma_s_old(alpha) + t_avg^s / (1 + F_alpha^s * t_avg^s)
!
! Both sublattices are updated simultaneously within each Soven iteration.
!
! Spectral function:
!   A(k,E) = -(1/pi) * Im Tr[ G(k,E+i*eta) ]   (trace over all 4*n_ref_states)
!
MODULE cpa_multi_solver

  USE precision,       ONLY : dp
  USE globals,         ONLY : n_ref_states
  USE constants,       ONLY : pi
  USE cpa_multi_types, ONLY : cpa_multi_params
  USE cpa_multi_hk,    ONLY : build_multi_HK
  USE cpa_linalg,      ONLY : invert_complex_matrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cpa_multi_solve_sigma_at_z
  PUBLIC :: cpa_multi_compute_spectral_function

  ! Pre-computed constants
  INTEGER, PARAMETER :: NDIM = 4*n_ref_states   ! = 40

CONTAINS

  ! Inline inversion using pre-allocated work — avoids repeated alloc in hot loop
  SUBROUTINE invert_40x40( A, ipiv, work, lwork, ierr )
    INTEGER,                     INTENT(IN)    :: lwork
    COMPLEX(dp), DIMENSION(NDIM,NDIM), INTENT(INOUT) :: A
    INTEGER,     DIMENSION(NDIM),      INTENT(INOUT) :: ipiv
    COMPLEX(dp), DIMENSION(lwork),     INTENT(INOUT) :: work
    INTEGER,                           INTENT(OUT)   :: ierr
    CALL ZGETRF( NDIM, NDIM, A, NDIM, ipiv, ierr )
    IF ( ierr /= 0 ) RETURN
    CALL ZGETRI( NDIM, A, NDIM, ipiv, work, lwork, ierr )
  END SUBROUTINE invert_40x40

  !===========================================================================
  ! Subroutine cpa_multi_solve_sigma_at_z
  !
  ! Solve sigma1(z) and sigma2(z) via Soven self-consistency for z=E+i*eta.
  !
  ! INPUT:
  !   E, eta        : energy and broadening (eta > 0)
  !   params        : cpa_multi_params
  !   kpts(3,n_k)   : BZ mesh (fractional coordinates)
  !   n_k           : number of BZ k-points
  !   tol           : convergence threshold |sigma_new - sigma_old|_max
  !   max_iter      : maximum number of Soven iterations
  !   mixing_alpha  : linear mixing factor (1 = no mixing)
  !
  ! OUTPUT:
  !   sigma1(n_ref_states) : converged CPA self-energy for sublattice 1
  !   sigma2(n_ref_states) : converged CPA self-energy for sublattice 2
  !   converged            : .TRUE. if converged before max_iter
  !   n_iter_used          : actual number of iterations
  !===========================================================================

  SUBROUTINE cpa_multi_solve_sigma_at_z( E, eta, params, kpts, n_k, &
                                          tol, max_iter, mixing_alpha, &
                                          sigma1, sigma2, converged, n_iter_used )

    REAL(dp),               INTENT(IN)  :: E, eta
    TYPE(cpa_multi_params), INTENT(IN)  :: params
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: kpts
    INTEGER,                INTENT(IN)  :: n_k
    REAL(dp),               INTENT(IN)  :: tol
    INTEGER,                INTENT(IN)  :: max_iter
    REAL(dp),               INTENT(IN)  :: mixing_alpha

    COMPLEX(dp), DIMENSION(n_ref_states), INTENT(OUT) :: sigma1, sigma2
    LOGICAL,                              INTENT(OUT) :: converged
    INTEGER,                              INTENT(OUT) :: n_iter_used

    !=========================================================================
    ! Local variables
    !=========================================================================

    COMPLEX(dp), DIMENSION(n_ref_states) :: sigma1_old, sigma1_new
    COMPLEX(dp), DIMENSION(n_ref_states) :: sigma2_old, sigma2_new
    COMPLEX(dp), DIMENSION(n_ref_states) :: F1_alpha, F2_alpha

    COMPLEX(dp) :: HK(NDIM,NDIM), Gk(NDIM,NDIM)
    INTEGER      :: ipiv(NDIM)
    COMPLEX(dp), ALLOCATABLE :: work_arr(:)
    INTEGER :: lwork

    COMPLEX(dp) :: z, v_i, t_i, t_avg
    REAL(dp) :: max_diff

    INTEGER :: i_iter, i_k, i_orb, ierr, idx
    INTEGER :: off1_up, off1_dn, off2_up, off2_dn

    lwork   = NDIM * NDIM
    ALLOCATE(work_arr(lwork))
    off1_up = 0
    off1_dn = n_ref_states
    off2_up = 2*n_ref_states
    off2_dn = 3*n_ref_states

    z = CMPLX(E, eta, dp)

    !=========================================================================
    ! Initialize sigma from VCA (real, imaginary part = 0)
    !   sigma_s(alpha) = SUM_i conc_i * pseudo_s(i)%onsite(alpha)
    !=========================================================================

    sigma1_old = CMPLX(0.0_dp, 0.0_dp, dp)
    sigma2_old = CMPLX(0.0_dp, 0.0_dp, dp)

    DO idx = 1, params%n_pseudo1
      DO i_orb = 1, n_ref_states
        sigma1_old(i_orb) = sigma1_old(i_orb) + &
          params%pseudo1(idx)%conc * params%pseudo1(idx)%onsite(i_orb)
        sigma2_old(i_orb) = sigma2_old(i_orb) + &
          params%pseudo2(idx)%conc * params%pseudo2(idx)%onsite(i_orb)
      END DO
    END DO

    converged = .FALSE.

    !=========================================================================
    ! Soven loop
    !=========================================================================

    soven_loop: DO i_iter = 1, max_iter

      !------------------------------------------------------------------------
      ! Step 1: BZ sum for F1_alpha and F2_alpha
      !------------------------------------------------------------------------

      F1_alpha = CMPLX(0.0_dp, 0.0_dp, dp)
      F2_alpha = CMPLX(0.0_dp, 0.0_dp, dp)

      DO i_k = 1, n_k

        CALL build_multi_HK( kpts(:,i_k), params, sigma1_old, sigma2_old, HK )

        ! G(k,z) = (z*I - HK)^{-1}
        Gk = -HK
        DO i_orb = 1, ndim
          Gk(i_orb, i_orb) = Gk(i_orb, i_orb) + z
        END DO

        CALL invert_40x40( Gk, ipiv, work_arr, lwork, ierr )

        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): matrix inversion failed at i_k=', &
                     i_k, '  ierr=', ierr
          STOP 1
        END IF

        ! Accumulate spin-averaged diagonal for each sublattice
        DO i_orb = 1, n_ref_states
          F1_alpha(i_orb) = F1_alpha(i_orb) + &
            0.5_dp * ( Gk(off1_up + i_orb, off1_up + i_orb) + &
                       Gk(off1_dn + i_orb, off1_dn + i_orb) )
          F2_alpha(i_orb) = F2_alpha(i_orb) + &
            0.5_dp * ( Gk(off2_up + i_orb, off2_up + i_orb) + &
                       Gk(off2_dn + i_orb, off2_dn + i_orb) )
        END DO

      END DO

      F1_alpha = F1_alpha / REAL(n_k, dp)
      F2_alpha = F2_alpha / REAL(n_k, dp)

      !------------------------------------------------------------------------
      ! Step 2: Soven update for sigma1 and sigma2, per orbital
      !------------------------------------------------------------------------

      sigma1_new = CMPLX(0.0_dp, 0.0_dp, dp)
      sigma2_new = CMPLX(0.0_dp, 0.0_dp, dp)

      DO i_orb = 1, n_ref_states

        !-- sublattice 1 --
        t_avg = CMPLX(0.0_dp, 0.0_dp, dp)
        DO idx = 1, params%n_pseudo1
          v_i = CMPLX(params%pseudo1(idx)%onsite(i_orb), 0.0_dp, dp) - sigma1_old(i_orb)
          t_i = v_i / ( CMPLX(1.0_dp,0.0_dp,dp) - F1_alpha(i_orb)*v_i )
          t_avg = t_avg + params%pseudo1(idx)%conc * t_i
        END DO
        sigma1_new(i_orb) = sigma1_old(i_orb) + &
          t_avg / ( CMPLX(1.0_dp,0.0_dp,dp) + F1_alpha(i_orb)*t_avg )

        !-- sublattice 2 --
        t_avg = CMPLX(0.0_dp, 0.0_dp, dp)
        DO idx = 1, params%n_pseudo2
          v_i = CMPLX(params%pseudo2(idx)%onsite(i_orb), 0.0_dp, dp) - sigma2_old(i_orb)
          t_i = v_i / ( CMPLX(1.0_dp,0.0_dp,dp) - F2_alpha(i_orb)*v_i )
          t_avg = t_avg + params%pseudo2(idx)%conc * t_i
        END DO
        sigma2_new(i_orb) = sigma2_old(i_orb) + &
          t_avg / ( CMPLX(1.0_dp,0.0_dp,dp) + F2_alpha(i_orb)*t_avg )

      END DO

      !------------------------------------------------------------------------
      ! Step 3: Convergence check (max over both sublattices)
      !------------------------------------------------------------------------

      max_diff = MAX( &
        MAXVAL( ABS(sigma1_new - sigma1_old) ), &
        MAXVAL( ABS(sigma2_new - sigma2_old) ) )

      IF ( max_diff < tol ) THEN
        sigma1_old   = sigma1_new
        sigma2_old   = sigma2_new
        converged    = .TRUE.
        n_iter_used  = i_iter
        EXIT soven_loop
      END IF

      !------------------------------------------------------------------------
      ! Step 4: Linear mixing
      !------------------------------------------------------------------------

      sigma1_old = mixing_alpha * sigma1_new + (1.0_dp - mixing_alpha) * sigma1_old
      sigma2_old = mixing_alpha * sigma2_new + (1.0_dp - mixing_alpha) * sigma2_old

      n_iter_used = i_iter

    END DO soven_loop

    sigma1 = sigma1_old
    sigma2 = sigma2_old

    IF ( .NOT. converged ) THEN
      WRITE(*,'(a,f10.4,a,i6,a,es12.4)') &
        ' WARNING (cpa_multi_solver): not converged at E=', E, &
        '  n_iter=', n_iter_used, '  max_diff=', max_diff
    END IF

    DEALLOCATE(work_arr)

  END SUBROUTINE cpa_multi_solve_sigma_at_z


  !===========================================================================
  ! Subroutine cpa_multi_compute_spectral_function
  !
  ! Compute A(k,E) along L-Gamma-X k-path over an energy range.
  !
  ! For each energy E:
  !   1. Solve sigma1(E+i*eta), sigma2(E+i*eta) via BZ-mesh Soven loop
  !   2. For each k on the path: build HK(k), invert, compute trace
  !      A(k,E) = -(1/pi) * Im Tr[ G(k,E+i*eta) ]
  !
  ! OUTPUT file format: gnuplot-compatible block data
  !   k_dist   E(eV)   A(k,E)
  !   (blank line after each E-block)
  !===========================================================================

  SUBROUTINE cpa_multi_compute_spectral_function( E_min, E_max, n_E, eta, &
                                                   params, &
                                                   kpts_bz, n_k_bz, &
                                                   kpts_path, k_dist, n_k_path, &
                                                   tol, max_iter, mixing_alpha, &
                                                   output_filename )

    REAL(dp),               INTENT(IN) :: E_min, E_max
    INTEGER,                INTENT(IN) :: n_E
    REAL(dp),               INTENT(IN) :: eta
    TYPE(cpa_multi_params), INTENT(IN) :: params
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: kpts_bz
    INTEGER,                INTENT(IN) :: n_k_bz
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: kpts_path
    REAL(dp), DIMENSION(:),   INTENT(IN) :: k_dist
    INTEGER,                INTENT(IN) :: n_k_path
    REAL(dp),               INTENT(IN) :: tol
    INTEGER,                INTENT(IN) :: max_iter
    REAL(dp),               INTENT(IN) :: mixing_alpha
    CHARACTER(LEN=*),       INTENT(IN) :: output_filename

    !=========================================================================
    ! Local variables
    !=========================================================================

    COMPLEX(dp), DIMENSION(n_ref_states) :: sigma1, sigma2
    COMPLEX(dp) :: HK(NDIM,NDIM), Gk(NDIM,NDIM)
    INTEGER      :: ipiv(NDIM)
    COMPLEX(dp), ALLOCATABLE :: work_arr(:)
    INTEGER :: lwork

    COMPLEX(dp) :: z
    REAL(dp)    :: E, A_kE
    LOGICAL     :: converged
    INTEGER     :: n_iter_used

    INTEGER :: i_E, i_kp, i_orb, ierr, file_num
    
    ! Peak detection at Gamma
    INTEGER :: i_gamma
    REAL(dp), ALLOCATABLE :: E_array(:), A_gamma(:)
    REAL(dp), ALLOCATABLE :: peak_energies(:)
    INTEGER :: n_peaks

    lwork = NDIM * NDIM
    ALLOCATE(work_arr(lwork))
    ALLOCATE(E_array(n_E), A_gamma(n_E))

    !=========================================================================

    ! Find Gamma point index (k = [0,0,0])
    i_gamma = -1
    DO i_kp = 1, n_k_path
      IF ( ABS(kpts_path(1,i_kp)) < 1.0d-6 .AND. &
           ABS(kpts_path(2,i_kp)) < 1.0d-6 .AND. &
           ABS(kpts_path(3,i_kp)) < 1.0d-6 ) THEN
        i_gamma = i_kp
        EXIT
      END IF
    END DO

    OPEN( NEWUNIT=file_num, FILE=TRIM(output_filename), STATUS='REPLACE', ACTION='WRITE' )

    WRITE(file_num, '(a)') '# k_dist   E(eV)   A(k,E)'

    WRITE(*,'(a)') ' ============================================'
    WRITE(*,'(a)') ' cpa_multi_compute_spectral_function: starting'
    WRITE(*,'(a,f10.4,a,f10.4,a,i6)') '   E_min=', E_min, '  E_max=', E_max, &
                                        '  n_E=', n_E
    WRITE(*,'(a,i6,a,i6)') '   n_k_bz=', n_k_bz, '  n_k_path=', n_k_path
    WRITE(*,'(a)') ' ============================================'

    energy_loop: DO i_E = 1, n_E

      IF ( n_E > 1 ) THEN
        E = E_min + (E_max - E_min) * REAL(i_E - 1, dp) / REAL(n_E - 1, dp)
      ELSE
        E = E_min
      END IF

      !------------------------------------------------------------------------
      ! Step 1: Soven self-consistency for sigma1, sigma2 at z = E + i*eta
      !------------------------------------------------------------------------

      CALL cpa_multi_solve_sigma_at_z( E, eta, params, kpts_bz, n_k_bz, &
                                        tol, max_iter, mixing_alpha, &
                                        sigma1, sigma2, converged, n_iter_used )

      WRITE(*,'(a,f10.4,a,i5,a,l2)') &
        '   E=', E, '   n_iter=', n_iter_used, '   converged=', converged

      z = CMPLX(E, eta, dp)

      !------------------------------------------------------------------------
      ! Step 2: Compute A(k,E) along k-path using converged sigma
      !------------------------------------------------------------------------

      ! Store energy value
      E_array(i_E) = E

      DO i_kp = 1, n_k_path

        CALL build_multi_HK( kpts_path(:,i_kp), params, sigma1, sigma2, HK )

        Gk = -HK
        DO i_orb = 1, ndim
          Gk(i_orb,i_orb) = Gk(i_orb,i_orb) + z
        END DO

        CALL invert_40x40( Gk, ipiv, work_arr, lwork, ierr )

        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): inversion failed at E=', E, &
                     '  k_path idx=', i_kp
          STOP 1
        END IF

        ! A(k,E) = -(1/pi) Im Tr[Gk]  (full trace over all 4*n_ref_states)
        A_kE = 0.0_dp
        DO i_orb = 1, ndim
          A_kE = A_kE - AIMAG(Gk(i_orb,i_orb)) / pi
        END DO

        WRITE(file_num, '(3(es16.8))') k_dist(i_kp), E, A_kE
        
        ! Store A at Gamma for peak detection
        IF ( i_kp == i_gamma .AND. i_gamma > 0 ) THEN
          A_gamma(i_E) = A_kE
        END IF

      END DO

      ! Blank line separates E-blocks (gnuplot block-data format)
      WRITE(file_num, '(a)') ''

    END DO energy_loop

    CLOSE(file_num)

    !=========================================================================
    ! Peak detection and band gap analysis at Gamma
    !=========================================================================

    IF ( i_gamma > 0 ) THEN
      CALL detect_peaks_and_report( E_array, A_gamma, n_E )
    ELSE
      WRITE(*,'(a)') ' WARNING: Gamma point not found in k-path, skipping peak analysis'
    END IF

    WRITE(*,'(a)') ' ============================================'
    WRITE(*,'(a,a)') ' Output written to: ', TRIM(output_filename)
    WRITE(*,'(a)') ' Plot with gnuplot:'
    WRITE(*,'(a,a,a)') '   set pm3d map; splot "', TRIM(output_filename), '" u 1:2:3'
    WRITE(*,'(a)') ' ============================================'

    DEALLOCATE(work_arr, E_array, A_gamma)

  END SUBROUTINE cpa_multi_compute_spectral_function


  !===========================================================================
  ! Subroutine detect_peaks_and_report
  !
  ! Find local maxima (peaks) in A(Gamma, E) and report VBM, CBM, band gap.
  !
  ! Peak criterion: A(i) > A(i-1) AND A(i) > A(i+1) AND A(i) > threshold
  ! where threshold = 0.05 * max(A_gamma) to filter noise
  !===========================================================================

  SUBROUTINE detect_peaks_and_report( E_array, A_gamma, n_E )
    REAL(dp), DIMENSION(:), INTENT(IN) :: E_array, A_gamma
    INTEGER,                INTENT(IN) :: n_E

    REAL(dp), ALLOCATABLE :: peak_E(:)
    INTEGER :: i, n_peaks, i_peak
    REAL(dp) :: A_max, threshold, VBM, CBM, E_gap
    LOGICAL :: found_VBM, found_CBM

    ! Threshold for peak detection (5% of max spectral function)
    A_max = MAXVAL(A_gamma)
    threshold = 0.05_dp * A_max

    ! Count peaks
    n_peaks = 0
    DO i = 2, n_E - 1
      IF ( A_gamma(i) > A_gamma(i-1) .AND. &
           A_gamma(i) > A_gamma(i+1) .AND. &
           A_gamma(i) > threshold ) THEN
        n_peaks = n_peaks + 1
      END IF
    END DO

    IF ( n_peaks == 0 ) THEN
      WRITE(*,'(a)') ' ============================================'
      WRITE(*,'(a)') ' Peak Analysis at Gamma: No peaks detected'
      WRITE(*,'(a,es12.4)') '   (threshold = 5% of A_max = ', threshold, ')'
      WRITE(*,'(a)') ' ============================================'
      RETURN
    END IF

    ! Extract peak energies
    ALLOCATE( peak_E(n_peaks) )
    i_peak = 0
    DO i = 2, n_E - 1
      IF ( A_gamma(i) > A_gamma(i-1) .AND. &
           A_gamma(i) > A_gamma(i+1) .AND. &
           A_gamma(i) > threshold ) THEN
        i_peak = i_peak + 1
        peak_E(i_peak) = E_array(i)
      END IF
    END DO

    ! Sort peak energies (bubble sort - fine for small n_peaks)
    CALL sort_real_array( peak_E, n_peaks )

    ! Identify VBM and CBM
    ! VBM = highest peak with E < 0
    ! CBM = lowest peak with E > 0
    found_VBM = .FALSE.
    found_CBM = .FALSE.
    VBM = -999.0_dp
    CBM =  999.0_dp

    DO i = 1, n_peaks
      IF ( peak_E(i) < 0.0_dp ) THEN
        VBM = peak_E(i)  ! Keep updating to get the highest (last one < 0)
        found_VBM = .TRUE.
      ELSE IF ( peak_E(i) > 0.0_dp .AND. .NOT. found_CBM ) THEN
        CBM = peak_E(i)  ! First one > 0
        found_CBM = .TRUE.
      END IF
    END DO

    ! Report
    WRITE(*,'(a)') ' ============================================'
    WRITE(*,'(a)') ' Peak Analysis at Gamma Point'
    WRITE(*,'(a)') ' ============================================'
    WRITE(*,'(a,i4)') ' Number of peaks detected: ', n_peaks
    WRITE(*,'(a)') ' Peak energies (eV), sorted:'
    DO i = 1, n_peaks
      WRITE(*,'(a,i4,a,f10.5)') '   Peak ', i, ':  E = ', peak_E(i)
    END DO
    WRITE(*,'(a)') ' --------------------------------------------'

    IF ( found_VBM ) THEN
      WRITE(*,'(a,f10.5,a)') ' VBM (Valence Band Maximum):  ', VBM, ' eV'
    ELSE
      WRITE(*,'(a)') ' VBM: Not found (no peak with E < 0)'
    END IF

    IF ( found_CBM ) THEN
      WRITE(*,'(a,f10.5,a)') ' CBM (Conduction Band Minimum):', CBM, ' eV'
    ELSE
      WRITE(*,'(a)') ' CBM: Not found (no peak with E > 0)'
    END IF

    IF ( found_VBM .AND. found_CBM ) THEN
      E_gap = CBM - VBM
      WRITE(*,'(a,f10.5,a)') ' Band Gap (CBM - VBM):         ', E_gap, ' eV'
    ELSE
      WRITE(*,'(a)') ' Band Gap: Cannot determine'
    END IF

    WRITE(*,'(a)') ' ============================================'

    DEALLOCATE(peak_E)

  END SUBROUTINE detect_peaks_and_report


  !===========================================================================
  ! Simple bubble sort for real array (ascending order)
  !===========================================================================

  SUBROUTINE sort_real_array( arr, n )
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: arr
    INTEGER,                INTENT(IN)    :: n
    INTEGER :: i, j
    REAL(dp) :: tmp

    DO i = 1, n - 1
      DO j = i + 1, n
        IF ( arr(j) < arr(i) ) THEN
          tmp    = arr(i)
          arr(i) = arr(j)
          arr(j) = tmp
        END IF
      END DO
    END DO

  END SUBROUTINE sort_real_array

END MODULE cpa_multi_solver
