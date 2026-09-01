! cpa_multi_solver.f90
!
! N-component Soven self-consistency and spectral function for the
! multi-component alloy, under the VCA-before-CPA scheme (see
! implementation_plan.md, "Critical Revision").
!
! Two self-energy structures are used, depending on params%odd_hopping:
!
! (A) odd_hopping = .FALSE. (VCA hopping, diagonal disorder only):
!     ONE self-energy SHARED by every species of a sublattice (unchanged
!     from the original implementation). This is exactly the "standard
!     CPA limit" of Blackman, Esterling, Berk, PRB 4, 2412 (1971), Sec.
!     IV: shared self-energy is correct ONLY when the hopping integral is
!     non-random (same value regardless of which species sit at the bond
!     ends). Soven equation per orbital alpha, sublattice s:
!
!       F_alpha^s(z) = (1/N_k) SUM_k [ z*I - HK(k) ]^{-1}_{sl_s,alpha,alpha}
!       t_i = v_i / (1 - F_alpha^s v_i),  v_i = elem_s(i)%onsite(alpha) - sigma_s(alpha)
!       t_avg = SUM_i c_i t_i
!       sigma_s_new(alpha) = sigma_s(alpha) + t_avg / (1 + F_alpha^s t_avg)
!
! (B) odd_hopping = .TRUE. (BEB, random/off-diagonal hopping):
!     a SEPARATE self-energy PER ACTIVE SPECIES of a sublattice. A shared
!     self-energy is explicitly NOT valid here -- Blackman Sec. III, Eq.
!     (3.14)-(3.22) (verified directly from the original scanned paper,
!     not from OCR/paraphrase) shows the renormalized locators gamma^A,
!     gamma^B require DIFFERENT self-energies U_1 != U_2 whenever hopping
!     is random; only in the nonrandom-hopping limit does U_1=U_2 (path A
!     above). Quote: "the inherent nonequivalence of atomic type, even on
!     averaging, is reflected in the self-energy corrections of the
!     renormalized locators... U1 and U2 are generally unequal."
!     Papaconstantopoulos, Gonis, Laufer, PRB 40, 12196 (1989), Eq. (2.8)
!     shows the identical structure (self-energy matrix diag(sigma^A,
!     sigma^B), generally sigma^A != sigma^B) for the multiband
!     tight-binding case. Koepernik, Velicky, Hayn, Eschrig, PRB 55, 5717
!     (1997), Eq. (65), gives the general multi-sublattice single-site
!     CPA condition as a concentration-weighted SUM over the possible
!     occupying species of a site.
!
!     Our lattice is bipartite (zinc-blende: sublattice 1 and sublattice
!     2, bonds exist ONLY between the two sublattices, never within one),
!     unlike Blackman's original single-sublattice A/B alloy (where bonds
!     can be AA, BB, or AB, forcing an extra "interactor" U_3 between A
!     and B because A and B directly bond to each other). Because no two
!     species of the SAME sublattice are ever directly bonded here, the
!     self-energy can be taken diagonal in species (no same-sublattice
!     interactor term is generated at the single-site level used by all
!     three references above) while still keeping a DISTINCT diagonal
!     value per species -- this is the direct generalization of
!     Blackman's Eq. (3.14)-(3.15) to M/N components per sublattice.
!
!     Each active species Q of sublattice s occupies its OWN block in the
!     augmented (BEB) Hamiltonian (see cpa_multi_hk.f90::build_multi_HK_odd).
!     The bare hopping placed in that block's off-diagonal (cross-sublattice)
!     coupling is TRANSLATIONALLY INVARIANT and NEVER scaled by
!     concentration -- this is the defining trick of the augmented-space
!     (BEB) method: species Q's block hops with the SAME bare t_{QQ'} to
!     every neighbor regardless of how likely species Q actually is.
!     Concentration enters ENTIRELY through the single-site scattering
!     condition below, not through the Hamiltonian matrix elements.
!
!     For a given realization where the physical atom actually occupying
!     this site is species q, Koepernik Eq. (63)-(64) fixes the single-
!     site scattering potential b^{(q)QQ'} = delta_{QQ',Qq} x (bare atomic
!     inverse propagator of species q): i.e. b^{(q)} is EXACTLY ZERO for
!     every block Q != q. This is NOT "species Q has bare onsite energy
!     zero" (which is what an earlier, INCORRECT version of this code
!     used) -- it is the v -> infinity limit of the ordinary single-site
!     T-matrix v/(1-Fv): the block simply has no atom in it for this
!     realization, i.e. is decoupled with infinite scattering strength,
!     giving t_vacancy = lim_{v->inf} v/(1-Fv) = -1/F_Q, a UNIVERSAL value
!     that does not depend on sigma_Q at all (only on the current F_Q).
!     Averaging Koepernik Eq. (65) over which species q actually sits at
!     the site (species Q itself, weight c_Q; anything else, weight
!     1-c_Q) gives the single-site CPA condition PER SPECIES Q:
!
!       F_Q,alpha(z) = (1/N_k) SUM_k [ z*I - HK(k) ]^{-1}_{Q,alpha,alpha}
!                      (spin-averaged; HK built with the CURRENT guess of
!                      EVERY active species' own self-energy, both
!                      sublattices, since all species are coupled through
!                      the bare cross-sublattice hopping)
!
!       v_own = elem%onsite(alpha) - sigma_Q(alpha)
!       t_own = v_own / (1 - F_Q,alpha v_own)
!       t_vac = -1 / F_Q,alpha
!       0 = c_Q t_own + (1-c_Q) t_vac        <- single-site CPA condition
!
!     which is solved iteratively via the same robust Soven-type update
!     step used by the standard (shared-sigma) path:
!
!       t_avg = c_Q t_own + (1-c_Q) t_vac
!       sigma_Q_new(alpha) = sigma_Q(alpha) + t_avg / (1 + F_Q,alpha t_avg)
!
!     Reduction check: with a single active species (c_Q=1), the vacancy
!     term drops out (weight 1-c_Q=0) and sigma_Q=onsite(alpha) is the
!     trivial fixed point (v_own=0), matching the existing "no disorder"
!     shortcut. A species with c_Q=0 must still be excluded from the
!     augmented basis entirely (see cpa_multi_types.f90) rather than
!     relying on this formula alone, because even a correctly-vanishing
!     weight on the "own" term would leave its BARE hopping channel to
!     the other sublattice fully active in HK(k) for every OTHER
!     species' realization q -- i.e. c_Q=0 is a singular limit of the
!     equation above (it forces t_vac=0 i.e. F_Q -> infinity, a pole
!     that is not approached smoothly), so it is handled by removing the
!     block outright rather than by this self-consistency condition.
!
! SOC is never part of the Soven loop (see cpa_multi_types.f90 /
! cpa_multi_build.f90): it always enters build_multi_HK/build_multi_HK_odd
! via the precomputed so_p_sl1_vca / so_p_sl2_vca (2nd VCA pass).
!
! Spectral function:
!   A(k,E) = -(1/pi) * Im Tr[ G(k,E+i*eta) ]   (trace over the full HK dimension)
!
MODULE cpa_multi_solver

  USE precision,       ONLY : dp
  USE globals,         ONLY : n_ref_states
  USE constants,       ONLY : pi
  USE cpa_multi_types, ONLY : cpa_multi_params
  USE cpa_multi_hk,    ONLY : build_multi_HK, build_multi_HK_odd, get_multi_HK_ndim
  USE cpa_linalg,      ONLY : invert_complex_matrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cpa_multi_solve_sigma_at_z
  PUBLIC :: cpa_multi_compute_spectral_function

  ! Fixed dimension used when params%odd_hopping = .FALSE. (VCA hopping,
  ! unchanged from before). When params%odd_hopping = .TRUE., the actual
  ! Bloch Hamiltonian dimension is obtained at runtime via
  ! get_multi_HK_ndim(params) and all work arrays are allocated accordingly.
  INTEGER, PARAMETER :: NDIM = 4*n_ref_states   ! = 40

CONTAINS

  ! Inline inversion using pre-allocated work — avoids repeated alloc in hot loop.
  ! ndim_loc is the actual matrix dimension (== NDIM for VCA hopping, or the
  ! enlarged BEB dimension from get_multi_HK_ndim(params) for ODD hopping).
  SUBROUTINE invert_nxn( ndim_loc, A, ipiv, work, lwork, ierr )
    INTEGER,                       INTENT(IN)    :: ndim_loc
    INTEGER,                       INTENT(IN)    :: lwork
    COMPLEX(dp), DIMENSION(:,:),   INTENT(INOUT) :: A
    INTEGER,     DIMENSION(:),     INTENT(INOUT) :: ipiv
    COMPLEX(dp), DIMENSION(:),     INTENT(INOUT) :: work
    INTEGER,                       INTENT(OUT)   :: ierr
    CALL ZGETRF( ndim_loc, ndim_loc, A, ndim_loc, ipiv, ierr )
    IF ( ierr /= 0 ) RETURN
    CALL ZGETRI( ndim_loc, A, ndim_loc, ipiv, work, lwork, ierr )
  END SUBROUTINE invert_nxn

  !===========================================================================
  ! Subroutine cpa_multi_solve_sigma_at_z
  !
  ! Solve sigma1(z) and sigma2(z) via Soven self-consistency for z=E+i*eta.
  !
  ! sigma1/sigma2 are returned as ALLOCATABLE(:,:), shape (n_species,n_ref_states):
  !   - odd_hopping = .FALSE. : n_species = 1 (single value shared by the
  !     whole sublattice; use sigma1(1,:) when building the Hamiltonian).
  !   - odd_hopping = .TRUE.  : n_species = n_active_sl1**2
  !     (resp. n_active_sl2**2), with row (Q-1)*n_active+Qp holding
  !     Sigma^{Q,Qp} for the BEB chemical-space coherent potential.
  !
  ! If a sublattice has only one ACTIVE species, its self-energy is fixed
  ! at that element's VCA onsite (no disorder, no Soven update). If BOTH
  ! sublattices have only one active species, this routine returns
  ! immediately with converged=.TRUE., n_iter_used=0, and no BZ sum is
  ! ever performed.
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
  !   sigma1, sigma2 : converged (or fixed) self-energy, see shape above
  !   converged      : .TRUE. if converged before max_iter
  !   n_iter_used    : actual number of iterations
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

    COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: sigma1, sigma2
    LOGICAL,                              INTENT(OUT) :: converged
    INTEGER,                              INTENT(OUT) :: n_iter_used

    IF ( params%odd_hopping ) THEN
      CALL solve_sigma_odd( E, eta, params, kpts, n_k, tol, max_iter, &
                             mixing_alpha, sigma1, sigma2, converged, n_iter_used )
    ELSE
      CALL solve_sigma_shared( E, eta, params, kpts, n_k, tol, max_iter, &
                                mixing_alpha, sigma1, sigma2, converged, n_iter_used )
    END IF

  END SUBROUTINE cpa_multi_solve_sigma_at_z


  !===========================================================================
  ! Subroutine solve_sigma_shared
  !
  ! odd_hopping = .FALSE. path: ONE self-energy shared by the whole
  ! sublattice (unchanged logic from the original implementation).
  ! Returns sigma1(1,:), sigma2(1,:).
  !===========================================================================

  SUBROUTINE solve_sigma_shared( E, eta, params, kpts, n_k, &
                                  tol, max_iter, mixing_alpha, &
                                  sigma1, sigma2, converged, n_iter_used )

    REAL(dp),               INTENT(IN)  :: E, eta
    TYPE(cpa_multi_params), INTENT(IN)  :: params
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: kpts
    INTEGER,                INTENT(IN)  :: n_k
    REAL(dp),               INTENT(IN)  :: tol
    INTEGER,                INTENT(IN)  :: max_iter
    REAL(dp),               INTENT(IN)  :: mixing_alpha

    COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: sigma1, sigma2
    LOGICAL,                              INTENT(OUT) :: converged
    INTEGER,                              INTENT(OUT) :: n_iter_used

    LOGICAL :: disorder_sl1, disorder_sl2

    COMPLEX(dp), DIMENSION(n_ref_states) :: sigma1_old, sigma1_new
    COMPLEX(dp), DIMENSION(n_ref_states) :: sigma2_old, sigma2_new
    COMPLEX(dp), DIMENSION(n_ref_states) :: F1_alpha, F2_alpha

    INTEGER :: ndim_loc
    COMPLEX(dp), ALLOCATABLE :: HK(:,:), Gk(:,:)
    INTEGER,     ALLOCATABLE :: ipiv(:)
    COMPLEX(dp), ALLOCATABLE :: work_arr(:)
    INTEGER :: lwork

    COMPLEX(dp) :: z, v_i, t_i, t_avg
    REAL(dp) :: max_diff

    INTEGER :: i_iter, i_k, i_orb, ierr, idx
    INTEGER :: off1_up, off1_dn, off2_up, off2_dn

    ALLOCATE( sigma1(1,n_ref_states), sigma2(1,n_ref_states) )

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

    IF ( .NOT. disorder_sl1 .AND. .NOT. disorder_sl2 ) THEN
      sigma1(1,:) = sigma1_old
      sigma2(1,:) = sigma2_old
      converged   = .TRUE.
      n_iter_used = 0
      RETURN
    END IF

    ndim_loc = get_multi_HK_ndim( params )
    lwork    = ndim_loc * ndim_loc
    ALLOCATE( HK(ndim_loc,ndim_loc), Gk(ndim_loc,ndim_loc), ipiv(ndim_loc) )
    ALLOCATE(work_arr(lwork))

    off1_up = 0
    off1_dn = n_ref_states
    off2_up = 2*n_ref_states
    off2_dn = 3*n_ref_states

    z = CMPLX(E, eta, dp)

    converged = .FALSE.

    soven_loop: DO i_iter = 1, max_iter

      F1_alpha = CMPLX(0.0_dp, 0.0_dp, dp)
      F2_alpha = CMPLX(0.0_dp, 0.0_dp, dp)

      DO i_k = 1, n_k

        CALL build_multi_HK( kpts(:,i_k), params, sigma1_old, sigma2_old, HK )

        Gk = -HK
        DO i_orb = 1, ndim_loc
          Gk(i_orb, i_orb) = Gk(i_orb, i_orb) + z
        END DO

        CALL invert_nxn( ndim_loc, Gk, ipiv, work_arr, lwork, ierr )

        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): matrix inversion failed at i_k=', &
                     i_k, '  ierr=', ierr
          STOP 1
        END IF

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

      IF ( disorder_sl1 ) &
        sigma1_old = mixing_alpha * sigma1_new + (1.0_dp - mixing_alpha) * sigma1_old
      IF ( disorder_sl2 ) &
        sigma2_old = mixing_alpha * sigma2_new + (1.0_dp - mixing_alpha) * sigma2_old

      n_iter_used = i_iter

    END DO soven_loop

    sigma1(1,:) = sigma1_old
    sigma2(1,:) = sigma2_old

    IF ( .NOT. converged ) THEN
      WRITE(*,'(a,f10.4,a,i6,a,es12.4)') &
        ' WARNING (cpa_multi_solver): not converged at E=', E, &
        '  n_iter=', n_iter_used, '  max_diff=', max_diff
    END IF

    DEALLOCATE(work_arr)
    DEALLOCATE(HK, Gk, ipiv)

  END SUBROUTINE solve_sigma_shared


  !=========================================================================
  ! Matrix BEB update for one sublattice.  The chemical-space convention is
  ! row (Q-1)*n_species+Qp.  This implements Koepernik et al. Eq. (65):
  ! <T> = sum_q c_q [b^(q)-Gamma]^-1 = 0, with
  ! delta Sigma = [I + <T> Gamma]^-1 <T>.
  !===========================================================================

  SUBROUTINE update_beb_sigma_block( n_species, onsite, concentration, gamma, &
                                     sigma_old, sigma_new, ierr )

    INTEGER, INTENT(IN) :: n_species
    REAL(dp), DIMENSION(n_species,n_ref_states), INTENT(IN) :: onsite
    REAL(dp), DIMENSION(n_species), INTENT(IN) :: concentration
    COMPLEX(dp), DIMENSION(n_species,n_species,n_ref_states), INTENT(IN) :: gamma
    COMPLEX(dp), DIMENSION(n_species*n_species,n_ref_states), INTENT(IN) :: sigma_old
    COMPLEX(dp), DIMENSION(n_species*n_species,n_ref_states), INTENT(OUT) :: sigma_new
    INTEGER, INTENT(OUT) :: ierr

    COMPLEX(dp), ALLOCATABLE :: b_q(:,:), t_q(:,:), t_avg(:,:), lhs(:,:), delta(:,:)
    COMPLEX(dp) :: v_q
    INTEGER :: i_orb, q, qp, info
    REAL(dp), PARAMETER :: INV_POTENTIAL_CUTOFF = 1.0d12
    REAL(dp), PARAMETER :: V_TOL = 1.0d-12

    ALLOCATE( b_q(n_species,n_species), t_q(n_species,n_species), &
              t_avg(n_species,n_species), lhs(n_species,n_species), &
              delta(n_species,n_species) )

    ierr = 0
    sigma_new = sigma_old

    DO i_orb = 1, n_ref_states
      t_avg = CMPLX(0.0_dp, 0.0_dp, dp)

      DO q = 1, n_species
        ! Only the actually occupied q channel is finite in b^(q).  The
        ! zero entries represent absent chemical states, not zero onsite.
        b_q = CMPLX(0.0_dp, 0.0_dp, dp)
        v_q = CMPLX(onsite(q,i_orb), 0.0_dp, dp) - &
              sigma_old((q-1)*n_species+q,i_orb)
        IF ( ABS(v_q) > V_TOL ) THEN
          b_q(q,q) = CMPLX(1.0_dp, 0.0_dp, dp) / v_q
        ELSE
          ! v_q=0 is the zero-scattering limit; represent its b->infinity
          ! limit without risking a floating-point division by zero.
          b_q(q,q) = CMPLX(INV_POTENTIAL_CUTOFF, 0.0_dp, dp)
        END IF

        t_q = b_q - gamma(:,:,i_orb)
        CALL invert_complex_matrix( t_q, n_species, info )
        IF ( info /= 0 ) THEN
          ierr = info
          EXIT
        END IF
        t_avg = t_avg + concentration(q) * t_q
      END DO
      IF ( ierr /= 0 ) EXIT

      lhs = MATMUL(t_avg, gamma(:,:,i_orb))
      DO q = 1, n_species
        lhs(q,q) = lhs(q,q) + CMPLX(1.0_dp, 0.0_dp, dp)
      END DO
      CALL invert_complex_matrix( lhs, n_species, info )
      IF ( info /= 0 ) THEN
        ierr = info
        EXIT
      END IF

      delta = MATMUL(lhs, t_avg)
      DO q = 1, n_species
        DO qp = 1, n_species
          sigma_new((q-1)*n_species+qp,i_orb) = &
            sigma_old((q-1)*n_species+qp,i_orb) + delta(q,qp)
        END DO
      END DO
    END DO

    DEALLOCATE( b_q, t_q, t_avg, lhs, delta )

  END SUBROUTINE update_beb_sigma_block


  !=========================================================================
  ! odd_hopping = .TRUE. path: matrix BEB CPA. sigma1 and sigma2 are
  ! flattened chemical-space matrices, shape (n_active_sl* ** 2,n_ref_states).
  !=========================================================================

  SUBROUTINE solve_sigma_odd( E, eta, params, kpts, n_k, &
                               tol, max_iter, mixing_alpha, &
                               sigma1, sigma2, converged, n_iter_used )

    REAL(dp),               INTENT(IN)  :: E, eta
    TYPE(cpa_multi_params), INTENT(IN)  :: params
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: kpts
    INTEGER,                INTENT(IN)  :: n_k
    REAL(dp),               INTENT(IN)  :: tol
    INTEGER,                INTENT(IN)  :: max_iter
    REAL(dp),               INTENT(IN)  :: mixing_alpha

    COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: sigma1, sigma2
    LOGICAL,                              INTENT(OUT) :: converged
    INTEGER,                              INTENT(OUT) :: n_iter_used

    LOGICAL :: disorder_sl1, disorder_sl2
    INTEGER :: n1, n2

    COMPLEX(dp), ALLOCATABLE :: sigma1_old(:,:), sigma1_new(:,:)
    COMPLEX(dp), ALLOCATABLE :: sigma2_old(:,:), sigma2_new(:,:)
    COMPLEX(dp), ALLOCATABLE :: gamma1(:,:,:), gamma2(:,:,:)
    REAL(dp), ALLOCATABLE :: onsite1(:,:), onsite2(:,:), concentration1(:), concentration2(:)

    INTEGER :: ndim_loc
    COMPLEX(dp), ALLOCATABLE :: HK(:,:), Gk(:,:)
    INTEGER,     ALLOCATABLE :: ipiv(:)
    COMPLEX(dp), ALLOCATABLE :: work_arr(:)
    INTEGER :: lwork

    COMPLEX(dp) :: z
    REAL(dp) :: max_diff
    REAL(dp), DIMENSION(n_ref_states) :: vca_onsite1, vca_onsite2

    INTEGER :: i_iter, i_k, i_orb, ierr, q1, q2, qp1, qp2
    INTEGER, ALLOCATABLE :: off_up_sl1(:), off_dn_sl1(:)
    INTEGER, ALLOCATABLE :: off_up_sl2(:), off_dn_sl2(:)

    n1 = params%n_active_sl1
    n2 = params%n_active_sl2

    ALLOCATE( sigma1(n1*n1,n_ref_states), sigma2(n2*n2,n_ref_states) )
    ALLOCATE( sigma1_old(n1*n1,n_ref_states), sigma1_new(n1*n1,n_ref_states) )
    ALLOCATE( sigma2_old(n2*n2,n_ref_states), sigma2_new(n2*n2,n_ref_states) )
    ALLOCATE( gamma1(n1,n1,n_ref_states), gamma2(n2,n2,n_ref_states) )
    ALLOCATE( onsite1(n1,n_ref_states), onsite2(n2,n_ref_states) )
    ALLOCATE( concentration1(n1), concentration2(n2) )

    DO q1 = 1, n1
      onsite1(q1,:) = params%elem_sl1(params%active_sl1(q1))%onsite(:)
      concentration1(q1) = params%elem_sl1(params%active_sl1(q1))%conc
    END DO
    DO q2 = 1, n2
      onsite2(q2,:) = params%elem_sl2(params%active_sl2(q2))%onsite(:)
      concentration2(q2) = params%elem_sl2(params%active_sl2(q2))%conc
    END DO

    disorder_sl1 = ( n1 > 1 )
    disorder_sl2 = ( n2 > 1 )

    sigma1_old = CMPLX(0.0_dp, 0.0_dp, dp)
    sigma2_old = CMPLX(0.0_dp, 0.0_dp, dp)
    vca_onsite1 = MATMUL(concentration1, onsite1)
    vca_onsite2 = MATMUL(concentration2, onsite2)
    DO q1 = 1, n1
      IF ( disorder_sl1 ) THEN
        sigma1_old((q1-1)*n1+q1,:) = CMPLX(vca_onsite1, -eta, dp)
      ELSE
        sigma1_old((q1-1)*n1+q1,:) = CMPLX(onsite1(q1,:), 0.0_dp, dp)
      END IF
    END DO
    DO q2 = 1, n2
      IF ( disorder_sl2 ) THEN
        sigma2_old((q2-1)*n2+q2,:) = CMPLX(vca_onsite2, -eta, dp)
      ELSE
        sigma2_old((q2-1)*n2+q2,:) = CMPLX(onsite2(q2,:), 0.0_dp, dp)
      END IF
    END DO

    IF ( .NOT. disorder_sl1 .AND. .NOT. disorder_sl2 ) THEN
      sigma1      = sigma1_old
      sigma2      = sigma2_old
      converged   = .TRUE.
      n_iter_used = 0
      DEALLOCATE( sigma1_old, sigma1_new, sigma2_old, sigma2_new, gamma1, gamma2 )
      DEALLOCATE( onsite1, onsite2, concentration1, concentration2 )
      RETURN
    END IF

    ndim_loc = get_multi_HK_ndim( params )
    lwork    = ndim_loc * ndim_loc
    ALLOCATE( HK(ndim_loc,ndim_loc), Gk(ndim_loc,ndim_loc), ipiv(ndim_loc) )
    ALLOCATE(work_arr(lwork))

    ! Species-block offsets, matching build_multi_HK_odd's layout exactly:
    ! ACTIVE sl1 species blocks first, then ACTIVE sl2 species blocks.
    ALLOCATE( off_up_sl1(n1), off_dn_sl1(n1) )
    ALLOCATE( off_up_sl2(n2), off_dn_sl2(n2) )
    DO q1 = 1, n1
      off_up_sl1(q1) = 2*n_ref_states*(q1-1)
      off_dn_sl1(q1) = off_up_sl1(q1) + n_ref_states
    END DO
    DO q2 = 1, n2
      off_up_sl2(q2) = 2*n_ref_states*n1 + 2*n_ref_states*(q2-1)
      off_dn_sl2(q2) = off_up_sl2(q2) + n_ref_states
    END DO

    z = CMPLX(E, eta, dp)

    converged = .FALSE.

    soven_loop_odd: DO i_iter = 1, max_iter

      ! Step 1: on-site coherent Green matrix Gamma_s^{Q,Qp}. The Q /= Qp
      ! entries are required by the matrix BEB condition.

      gamma1 = CMPLX(0.0_dp, 0.0_dp, dp)
      gamma2 = CMPLX(0.0_dp, 0.0_dp, dp)

      DO i_k = 1, n_k

        CALL build_multi_HK_odd( kpts(:,i_k), params, sigma1_old, sigma2_old, HK )

        Gk = -HK
        DO i_orb = 1, ndim_loc
          Gk(i_orb, i_orb) = Gk(i_orb, i_orb) + z
        END DO

        CALL invert_nxn( ndim_loc, Gk, ipiv, work_arr, lwork, ierr )

        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): matrix inversion failed at i_k=', &
                     i_k, '  ierr=', ierr
          STOP 1
        END IF

        DO i_orb = 1, n_ref_states
          DO q1 = 1, n1
            DO qp1 = 1, n1
              gamma1(q1,qp1,i_orb) = gamma1(q1,qp1,i_orb) + &
                0.5_dp * ( Gk(off_up_sl1(q1)+i_orb, off_up_sl1(qp1)+i_orb) + &
                           Gk(off_dn_sl1(q1)+i_orb, off_dn_sl1(qp1)+i_orb) )
            END DO
          END DO
          DO q2 = 1, n2
            DO qp2 = 1, n2
              gamma2(q2,qp2,i_orb) = gamma2(q2,qp2,i_orb) + &
                0.5_dp * ( Gk(off_up_sl2(q2)+i_orb, off_up_sl2(qp2)+i_orb) + &
                           Gk(off_dn_sl2(q2)+i_orb, off_dn_sl2(qp2)+i_orb) )
            END DO
          END DO
        END DO

      END DO

      gamma1 = gamma1 / REAL(n_k, dp)
      gamma2 = gamma2 / REAL(n_k, dp)

      ! Step 2: concentration-weighted matrix BEB update.

      sigma1_new = sigma1_old
      sigma2_new = sigma2_old

      IF ( disorder_sl1 ) THEN
        CALL update_beb_sigma_block(n1, onsite1, concentration1, gamma1, &
                                    sigma1_old, sigma1_new, ierr)
        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): BEB update failed on sublattice 1, ierr=', ierr
          STOP 1
        END IF
      END IF

      IF ( disorder_sl2 ) THEN
        CALL update_beb_sigma_block(n2, onsite2, concentration2, gamma2, &
                                    sigma2_old, sigma2_new, ierr)
        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): BEB update failed on sublattice 2, ierr=', ierr
          STOP 1
        END IF
      END IF

      !------------------------------------------------------------------------
      ! Step 3: Convergence check (max over ALL active species, both
      ! sublattices, that carry disorder)
      !------------------------------------------------------------------------

      max_diff = 0.0_dp
      IF ( disorder_sl1 ) max_diff = MAX( max_diff, MAXVAL( ABS(sigma1_new - sigma1_old) ) )
      IF ( disorder_sl2 ) max_diff = MAX( max_diff, MAXVAL( ABS(sigma2_new - sigma2_old) ) )

      IF ( max_diff < tol ) THEN
        sigma1_old   = sigma1_new
        sigma2_old   = sigma2_new
        converged    = .TRUE.
        n_iter_used  = i_iter
        EXIT soven_loop_odd
      END IF

      !------------------------------------------------------------------------
      ! Step 4: Linear mixing
      !------------------------------------------------------------------------

      IF ( disorder_sl1 ) &
        sigma1_old = mixing_alpha * sigma1_new + (1.0_dp - mixing_alpha) * sigma1_old
      IF ( disorder_sl2 ) &
        sigma2_old = mixing_alpha * sigma2_new + (1.0_dp - mixing_alpha) * sigma2_old

      n_iter_used = i_iter

    END DO soven_loop_odd

    sigma1 = sigma1_old
    sigma2 = sigma2_old

    IF ( .NOT. converged ) THEN
      WRITE(*,'(a,f10.4,a,i6,a,es12.4)') &
        ' WARNING (cpa_multi_solver): not converged at E=', E, &
        '  n_iter=', n_iter_used, '  max_diff=', max_diff
    END IF

    DEALLOCATE(work_arr)
    DEALLOCATE(HK, Gk, ipiv)
    DEALLOCATE( off_up_sl1, off_dn_sl1, off_up_sl2, off_dn_sl2 )
    DEALLOCATE( sigma1_old, sigma1_new, sigma2_old, sigma2_new, gamma1, gamma2 )
    DEALLOCATE( onsite1, onsite2, concentration1, concentration2 )

  END SUBROUTINE solve_sigma_odd


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

    ! Coherent potential: shape (1,n_ref_states) when odd_hopping = .FALSE.,
    ! and (n_active_sl* ** 2,n_ref_states) in the matrix BEB path. Row
    ! (Q-1)*n_active+Qp is Sigma^{Q,Qp}.
    COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: sigma1, sigma2
    INTEGER :: ndim_loc
    COMPLEX(dp), ALLOCATABLE :: HK(:,:), Gk(:,:)
    INTEGER,     ALLOCATABLE :: ipiv(:)
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

    ndim_loc = get_multi_HK_ndim( params )
    lwork = ndim_loc * ndim_loc
    ALLOCATE( HK(ndim_loc,ndim_loc), Gk(ndim_loc,ndim_loc), ipiv(ndim_loc) )
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
    WRITE(*,'(a,i3,a,i3)') '   n_active_sl1=', params%n_active_sl1, &
                            '  n_active_sl2=', params%n_active_sl2
    WRITE(*,'(a,l2,a,i5)') '   odd_hopping=', params%odd_hopping, &
                            '  HK dimension=', ndim_loc
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

      IF ( ALLOCATED(sigma1) ) DEALLOCATE(sigma1)
      IF ( ALLOCATED(sigma2) ) DEALLOCATE(sigma2)

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

        IF ( params%odd_hopping ) THEN
          CALL build_multi_HK_odd( kpts_path(:,i_kp), params, sigma1, sigma2, HK )
        ELSE
          CALL build_multi_HK( kpts_path(:,i_kp), params, sigma1(1,:), sigma2(1,:), HK )
        END IF

        Gk = -HK
        DO i_orb = 1, ndim_loc
          Gk(i_orb,i_orb) = Gk(i_orb,i_orb) + z
        END DO

        CALL invert_nxn( ndim_loc, Gk, ipiv, work_arr, lwork, ierr )

        IF ( ierr /= 0 ) THEN
          WRITE(*,*) 'ERROR (cpa_multi_solver): inversion failed at E=', E, &
                     '  k_path idx=', i_kp
          STOP 1
        END IF

        ! A(k,E) = -(1/pi) Im Tr[Gk]  (full trace over the Bloch Hamiltonian
        ! dimension, which is 4*n_ref_states for VCA hopping or the
        ! enlarged BEB dimension for ODD hopping)
        A_kE = 0.0_dp
        DO i_orb = 1, ndim_loc
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
      IF ( ABS(VBM_k - CBM_k) < 1.0d-2 ) THEN
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
    DEALLOCATE(HK, Gk, ipiv)
    IF ( ALLOCATED(sigma1) ) DEALLOCATE(sigma1)
    IF ( ALLOCATED(sigma2) ) DEALLOCATE(sigma2)

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