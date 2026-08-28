! cpa_multi_solver.f90
!
! N-component Soven self-consistency and spectral function for the
! multi-component alloy, under the VCA-before-CPA scheme (see
! implementation_plan.md, "Critical Revision").
!
! REVISED: the Soven loop now runs over the REAL ELEMENTS of each
! sublattice (elem_sl1(:), elem_sl2(:), sizes n_real_elem_sl1 /
! n_real_elem_sl2), not over the M*N binaries. If a sublattice has only
! ONE real element (no disorder there), its self-energy is fixed at the
! (real, since imaginary part = 0) VCA onsite value of that element and
! is NEVER updated by the Soven iteration -- there is no disorder to
! self-consistently resum on that sublattice. If BOTH sublattices have
! only one real element (pure binary case), no Soven iteration is
! performed at all.
!
! SOC is never part of the Soven loop (see cpa_multi_types.f90 /
! cpa_multi_build.f90): it always enters build_multi_HK via the
! precomputed so_p_sl1_vca / so_p_sl2_vca (2nd VCA pass).
!
! Soven equation per orbital alpha, sublattice s (=1 or 2), only for
! sublattices with disorder (n_real_elem_s > 1):
!
!   F_alpha^s(z) = (1/N_k) SUM_k [ z*I - HK(k) ]^{-1}_{sl_s,alpha,alpha}
!                  (spin-averaged: mean of up and down diagonal elements)
!
!   DO i = 1, n_real_elem_s
!      v_i = elem_s(i)%onsite(alpha) - sigma_s_old(alpha)
!      t_i = v_i / (1 - F_alpha^s * v_i)
!   END DO
!
!   t_avg^s = SUM_i elem_s(i)%conc * t_i
!
!   sigma_s_new(alpha) = sigma_s_old(alpha) + t_avg^s / (1 + F_alpha^s * t_avg^s)
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
  ! If a sublattice has only one real element (n_real_elem_s == 1), its
  ! self-energy is fixed to that element's VCA onsite (no disorder, no
  ! Soven update). If BOTH sublattices have n_real_elem == 1, this routine
  ! returns immediately with converged=.TRUE., n_iter_used=0, and no BZ
  ! sum is ever performed.
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
  !   sigma1(n_ref_states) : converged (or fixed) self-energy for sublattice 1
  !   sigma2(n_ref_states) : converged (or fixed) self-energy for sublattice 2
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

    LOGICAL :: disorder_sl1, disorder_sl2

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

    !=========================================================================
    ! Determine which sublattices carry disorder (n_real_elem > 1).
    ! A sublattice with a single real element has NO CPA disorder: its
    ! self-energy is fixed at that element's VCA onsite value (real,
    ! imaginary part exactly zero) for the whole run.
    !=========================================================================

    disorder_sl1 = ( params%n_real_elem_sl1 > 1 )
    disorder_sl2 = ( params%n_real_elem_sl2 > 1 )

    sigma1_old = CMPLX(0.0_dp, 0.0_dp, dp)
    sigma2_old = CMPLX(0.0_dp, 0.0_dp, dp)

    DO idx = 1, params%n_real_elem_sl1
      DO i_orb = 1, n_ref_states
        sigma1_old(i_orb) = sigma1_old(i_orb) + &
          params%elem_sl1(idx)%conc * params%elem_sl1(idx)%onsite(i_orb)
      END DO
    END DO
    DO idx = 1, params%n_real_elem_sl2
      DO i_orb = 1, n_ref_states
        sigma2_old(i_orb) = sigma2_old(i_orb) + &
          params%elem_sl2(idx)%conc * params%elem_sl2(idx)%onsite(i_orb)
      END DO
    END DO

    !=========================================================================
    ! Pure-binary shortcut: no disorder anywhere -> sigma = VCA onsite,
    ! no Soven iteration, no BZ sum ever performed.
    !=========================================================================

    IF ( .NOT. disorder_sl1 .AND. .NOT. disorder_sl2 ) THEN
      sigma1      = sigma1_old
      sigma2      = sigma2_old
      converged   = .TRUE.
      n_iter_used = 0
      RETURN
    END IF

    lwork   = NDIM * NDIM
    ALLOCATE(work_arr(lwork))
    off1_up = 0
    off1_dn = n_ref_states
    off2_up = 2*n_ref_states
    off2_dn = 3*n_ref_states

    z = CMPLX(E, eta, dp)

    converged = .FALSE.

    !=========================================================================
    ! Soven loop
    !=========================================================================

    soven_loop: DO i_iter = 1, max_iter

      !------------------------------------------------------------------------
      ! Step 1: BZ sum for F1_alpha and F2_alpha
      ! (still needed for both sublattices even if only one has disorder,
      !  since HK couples sl1 and sl2 via hopping)
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
      ! Step 2: Soven update for sigma1 and sigma2, per orbital.
      ! Only sublattices with disorder are updated; the other keeps its
      ! fixed VCA onsite value (sigma_old carried through unchanged).
      !------------------------------------------------------------------------

      sigma1_new = sigma1_old
      sigma2_new = sigma2_old

      DO i_orb = 1, n_ref_states

        IF ( disorder_sl1 ) THEN
          t_avg = CMPLX(0.0_dp, 0.0_dp, dp)
          DO idx = 1, params%n_real_elem_sl1
            v_i = CMPLX(params%elem_sl1(idx)%onsite(i_orb), 0.0_dp, dp) - sigma1_old(i_orb)
            t_i = v_i / ( CMPLX(1.0_dp,0.0_dp,dp) - F1_alpha(i_orb)*v_i )
            t_avg = t_avg + params%elem_sl1(idx)%conc * t_i
          END DO
          sigma1_new(i_orb) = sigma1_old(i_orb) + &
            t_avg / ( CMPLX(1.0_dp,0.0_dp,dp) + F1_alpha(i_orb)*t_avg )
        END IF

        IF ( disorder_sl2 ) THEN
          t_avg = CMPLX(0.0_dp, 0.0_dp, dp)
          DO idx = 1, params%n_real_elem_sl2
            v_i = CMPLX(params%elem_sl2(idx)%onsite(i_orb), 0.0_dp, dp) - sigma2_old(i_orb)
            t_i = v_i / ( CMPLX(1.0_dp,0.0_dp,dp) - F2_alpha(i_orb)*v_i )
            t_avg = t_avg + params%elem_sl2(idx)%conc * t_i
          END DO
          sigma2_new(i_orb) = sigma2_old(i_orb) + &
            t_avg / ( CMPLX(1.0_dp,0.0_dp,dp) + F2_alpha(i_orb)*t_avg )
        END IF

      END DO

      !------------------------------------------------------------------------
      ! Step 3: Convergence check (max over the sublattice(s) with disorder)
      !------------------------------------------------------------------------

      max_diff = 0.0_dp
      IF ( disorder_sl1 ) max_diff = MAX( max_diff, MAXVAL( ABS(sigma1_new - sigma1_old) ) )
      IF ( disorder_sl2 ) max_diff = MAX( max_diff, MAXVAL( ABS(sigma2_new - sigma2_old) ) )

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

      IF ( disorder_sl1 ) &
        sigma1_old = mixing_alpha * sigma1_new + (1.0_dp - mixing_alpha) * sigma1_old
      IF ( disorder_sl2 ) &
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
  !      (skipped internally if neither sublattice has disorder)
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
    
    ! VBM/CBM detection across all (k,E)
    REAL(dp), ALLOCATABLE :: A_matrix(:,:)
    REAL(dp) :: VBM_E, VBM_k, CBM_E, CBM_k, band_gap
    LOGICAL :: found_VBM, found_CBM

    lwork = NDIM * NDIM
    ALLOCATE(work_arr(lwork))
    ALLOCATE(A_matrix(n_k_path, n_E))

    !=========================================================================

    OPEN( NEWUNIT=file_num, FILE=TRIM(output_filename), STATUS='REPLACE', ACTION='WRITE' )

    WRITE(file_num, '(a)') '# k_dist   E(eV)   A(k,E)'

    WRITE(*,'(a)') ' ============================================'
    WRITE(*,'(a)') ' cpa_multi_compute_spectral_function: starting'
    WRITE(*,'(a,f10.4,a,f10.4,a,i6)') '   E_min=', E_min, '  E_max=', E_max, &
                                        '  n_E=', n_E
    WRITE(*,'(a,i6,a,i6)') '   n_k_bz=', n_k_bz, '  n_k_path=', n_k_path
    WRITE(*,'(a,i3,a,i3)') '   n_real_elem_sl1=', params%n_real_elem_sl1, &
                            '  n_real_elem_sl2=', params%n_real_elem_sl2
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

        A_matrix(i_kp, i_E) = A_kE
        WRITE(file_num, '(3(es16.8))') k_dist(i_kp), E, A_kE

      END DO

      ! Blank line separates E-blocks (gnuplot block-data format)
      WRITE(file_num, '(a)') ''

    END DO energy_loop

    CLOSE(file_num)

    !=========================================================================
    ! VBM/CBM detection and band gap analysis
    !=========================================================================

    CALL detect_VBM_CBM_bandgap( A_matrix, k_dist, E_min, E_max, n_k_path, n_E, &
                                  VBM_E, VBM_k, CBM_E, CBM_k, &
                                  found_VBM, found_CBM )

    IF ( found_VBM .AND. found_CBM ) THEN
      band_gap = CBM_E - VBM_E
      
      WRITE(*,'(a)') ' ============================================'
      WRITE(*,'(a)') ' Band Structure Analysis'
      WRITE(*,'(a)') ' ============================================'
      WRITE(*,'(a,f10.5,a,f10.5)') ' VBM:  E = ', VBM_E, ' eV,  k-dist = ', VBM_k
      WRITE(*,'(a,f10.5,a,f10.5)') ' CBM:  E = ', CBM_E, ' eV,  k-dist = ', CBM_k
      WRITE(*,'(a,f10.5,a)') ' Band Gap: ', band_gap, ' eV'
      
      ! Check if direct or indirect
      IF ( ABS(VBM_k - CBM_k) < 1.0d-6 ) THEN
        WRITE(*,'(a)') ' Band Gap Type: DIRECT'
      ELSE
        WRITE(*,'(a)') ' Band Gap Type: INDIRECT'
      END IF
      WRITE(*,'(a)') ' ============================================'
    ELSE
      WRITE(*,'(a)') ' ============================================'
      WRITE(*,'(a)') ' Band Structure Analysis'
      WRITE(*,'(a)') ' ============================================'
      IF ( .NOT. found_VBM ) WRITE(*,'(a)') ' VBM: Not found (no peak with E < 0)'
      IF ( .NOT. found_CBM ) WRITE(*,'(a)') ' CBM: Not found (no peak with E > 0)'
      WRITE(*,'(a)') ' Band Gap: Cannot determine'
      WRITE(*,'(a)') ' ============================================'
    END IF

    WRITE(*,'(a)') ' ============================================'
    WRITE(*,'(a,a)') ' Output written to: ', TRIM(output_filename)
    WRITE(*,'(a)') ' Plot with gnuplot:'
    WRITE(*,'(a,a,a)') '   set pm3d map; splot "', TRIM(output_filename), '" u 1:2:3'
    WRITE(*,'(a)') ' ============================================'

    DEALLOCATE(work_arr, A_matrix)

  END SUBROUTINE cpa_multi_compute_spectral_function


  !===========================================================================
  ! Subroutine detect_VBM_CBM_bandgap
  !
  ! Find VBM and CBM across all (k,E) points in the spectral function.
  !
  ! Strategy:
  !   Step 1: For EACH k-point, identify two key peaks:
  !           - VB_peak(k): highest-energy peak with E < 0
  !           - CB_peak(k): lowest-energy peak with E > 0
  !   Step 2: Among all k-points:
  !           - VBM = highest VB_peak(k) across all k
  !           - CBM = lowest CB_peak(k) across all k
  !
  ! Peak criterion: A(k,i_E) > A(k,i_E-1) AND A(k,i_E) > A(k,i_E+1) 
  !                 AND A(k,i_E) > threshold
  ! where threshold = 0.05 * max(A_matrix) to filter noise
  !===========================================================================

  SUBROUTINE detect_VBM_CBM_bandgap( A_matrix, k_dist, E_min, E_max, &
                                      n_k_path, n_E, &
                                      VBM_E, VBM_k, CBM_E, CBM_k, &
                                      found_VBM, found_CBM )

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: A_matrix  ! (n_k_path, n_E)
    REAL(dp), DIMENSION(:),   INTENT(IN) :: k_dist    ! (n_k_path)
    REAL(dp),                 INTENT(IN) :: E_min, E_max
    INTEGER,                  INTENT(IN) :: n_k_path, n_E

    REAL(dp), INTENT(OUT) :: VBM_E, VBM_k, CBM_E, CBM_k
    LOGICAL,  INTENT(OUT) :: found_VBM, found_CBM

    REAL(dp) :: E, A_max, threshold
    REAL(dp) :: VB_peak_E, CB_peak_E  ! Peaks at current k-point
    LOGICAL  :: found_VB_at_k, found_CB_at_k
    INTEGER  :: i_k, i_E

    ! Initialize global VBM/CBM
    VBM_E = -999.0_dp
    VBM_k = 0.0_dp
    CBM_E = 999.0_dp
    CBM_k = 0.0_dp
    found_VBM = .FALSE.
    found_CBM = .FALSE.

    ! Threshold for peak detection (5% of max spectral function)
    A_max = MAXVAL(A_matrix)
    threshold = 0.05_dp * A_max

    ! Loop over all k-points
    DO i_k = 1, n_k_path

      ! Initialize peaks for this k-point
      VB_peak_E = -999.0_dp
      CB_peak_E = 999.0_dp
      found_VB_at_k = .FALSE.
      found_CB_at_k = .FALSE.

      ! Step 1: Find VB and CB peaks at this k-point
      DO i_E = 2, n_E - 1

        ! Check if this is a local maximum
        IF ( A_matrix(i_k, i_E) > A_matrix(i_k, i_E-1) .AND. &
             A_matrix(i_k, i_E) > A_matrix(i_k, i_E+1) .AND. &
             A_matrix(i_k, i_E) > threshold ) THEN

          ! Calculate energy
          IF ( n_E > 1 ) THEN
            E = E_min + (E_max - E_min) * REAL(i_E - 1, dp) / REAL(n_E - 1, dp)
          ELSE
            E = E_min
          END IF

          ! Update VB peak: highest energy peak with E < 0
          IF ( E < 0.0_dp .AND. E > VB_peak_E ) THEN
            VB_peak_E = E
            found_VB_at_k = .TRUE.
          END IF

          ! Update CB peak: lowest energy peak with E > 0
          IF ( E > 0.0_dp .AND. E < CB_peak_E ) THEN
            CB_peak_E = E
            found_CB_at_k = .TRUE.
          END IF

        END IF

      END DO

      ! Step 2: Update global VBM (highest among all VB_peak_E)
      IF ( found_VB_at_k .AND. VB_peak_E > VBM_E ) THEN
        VBM_E = VB_peak_E
        VBM_k = k_dist(i_k)
        found_VBM = .TRUE.
      END IF

      ! Step 2: Update global CBM (lowest among all CB_peak_E)
      IF ( found_CB_at_k .AND. CB_peak_E < CBM_E ) THEN
        CBM_E = CB_peak_E
        CBM_k = k_dist(i_k)
        found_CBM = .TRUE.
      END IF

    END DO

  END SUBROUTINE detect_VBM_CBM_bandgap

END MODULE cpa_multi_solver