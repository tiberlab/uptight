! cpa_multi_hk.f90
!
! Bloch Hamiltonian H_K(k) builder for the multi-component CPA alloy.
!
! Generalization of build_hk.f90::build_HK to use:
!   - sigma1(alpha) : CPA self-energy for sublattice 1 (complex, n_ref_states)
!   - sigma2(alpha) : CPA self-energy for sublattice 2
!   - hopping_vca   : VCA-averaged hopping (real, n_ref_couplings)
!   - so_p_sl1_vca, so_p_sl2_vca : VCA-averaged SOC per sublattice
!
! Basis layout (identical to old build_hk.f90):
!   [sl1_up (n_ref_states) | sl1_dn (n_ref_states) |
!    sl2_up (n_ref_states) | sl2_dn (n_ref_states)]
!   Total: 4*n_ref_states = 40
!
! K-vector convention: fractional (same as in build_hk.f90 and sparse_ham).
! Phase = 2*pi * DOT(k_frac, bond_frac)
!
! The subroutines get_1nn_bond_vectors and add_soc_p_block are ported
! directly from the prototype build_hk.f90 (logic verified correct there).
!
MODULE cpa_multi_hk

  USE precision,      ONLY : dp
  USE globals,        ONLY : n_ref_states, px, py, pz
  USE constants,      ONLY : pi
  USE cpa_multi_types, ONLY : cpa_multi_params
  USE TB_ham_min,     ONLY : koster_slater_min

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: build_multi_HK
  PUBLIC :: build_multi_HK_odd
  PUBLIC :: get_multi_HK_ndim
  PUBLIC :: get_1nn_bond_vectors
  PUBLIC :: add_soc_p_block

CONTAINS

  !===========================================================================
  ! Subroutine get_1nn_bond_vectors
  !
  ! 4 zinc-blende 1NN bond vectors in fractional units of the cubic
  ! lattice constant a (from sl1 atom at origin toward 4 sl2 1NN atoms).
  !
  ! OUTPUT: bonds(3,4)  units of a (lattice constant a = 1)
  !===========================================================================

  SUBROUTINE get_1nn_bond_vectors( bonds )

    REAL(dp), DIMENSION(3,4), INTENT(OUT) :: bonds

    bonds(:,1) = (/  1.0_dp,  1.0_dp,  1.0_dp /) * 0.25_dp
    bonds(:,2) = (/  1.0_dp, -1.0_dp, -1.0_dp /) * 0.25_dp
    bonds(:,3) = (/ -1.0_dp,  1.0_dp, -1.0_dp /) * 0.25_dp
    bonds(:,4) = (/ -1.0_dp, -1.0_dp,  1.0_dp /) * 0.25_dp

  END SUBROUTINE get_1nn_bond_vectors


  !===========================================================================
  ! Subroutine build_multi_HK
  !
  ! Build the 40×40 Bloch Hamiltonian for the multi-component CPA alloy.
  ! NOTE: this routine assumes VCA hopping (params%odd_hopping = .FALSE.);
  ! it must not be called when params%odd_hopping = .TRUE. -- callers must
  ! dispatch to build_multi_HK_odd instead in that case (see
  ! get_multi_HK_ndim to size the enlarged matrix). This routine's fixed
  ! 4*n_ref_states signature is left completely untouched so existing
  ! VCA-hopping behavior and performance are unaffected.
  !
  ! Build the 40×40 Bloch Hamiltonian for the multi-component CPA alloy.
  !
  ! INPUT:
  !   kvec(3)    : k-vector in fractional coordinates (units of 2*pi/a)
  !   params     : cpa_multi_params (hopping_vca, so_p VCA values)
  !   sigma1(n_ref_states) : current CPA self-energy for sublattice 1 (COMPLEX)
  !   sigma2(n_ref_states) : current CPA self-energy for sublattice 2 (COMPLEX)
  !
  ! OUTPUT:
  !   HK(4*n_ref_states, 4*n_ref_states) : Bloch Hamiltonian, COMPLEX
  !===========================================================================

  SUBROUTINE build_multi_HK( kvec, params, sigma1, sigma2, HK )

    REAL(dp), DIMENSION(3),     INTENT(IN) :: kvec
    TYPE(cpa_multi_params),     INTENT(IN) :: params
    COMPLEX(dp), DIMENSION(n_ref_states), INTENT(IN) :: sigma1, sigma2

    COMPLEX(dp), DIMENSION(4*n_ref_states, 4*n_ref_states), INTENT(OUT) :: HK

    !=========================================================================
    ! Local variables
    !=========================================================================

    REAL(dp), DIMENSION(3,4) :: bonds_frac
    REAL(dp), DIMENSION(3)   :: cos_latt
    REAL(dp)                 :: phase

    REAL(dp), DIMENSION(n_ref_states, n_ref_states) :: tb_bloc

    COMPLEX(dp) :: exp_fac

    INTEGER :: i_bond, i_orb, j_orb
    INTEGER :: off1_up, off1_dn, off2_up, off2_dn

    !=========================================================================
    ! Block offsets (row/col base indices, 0-based):
    !   sl1_up:  rows 1..n_ref,          off1_up = 0
    !   sl1_dn:  rows n_ref+1..2*n_ref,  off1_dn = n_ref_states
    !   sl2_up:  rows 2*n_ref+1..3*n_ref, off2_up = 2*n_ref_states
    !   sl2_dn:  rows 3*n_ref+1..4*n_ref, off2_dn = 3*n_ref_states
    !=========================================================================

    HK = (0.0_dp, 0.0_dp)

    off1_up = 0
    off1_dn = n_ref_states
    off2_up = 2*n_ref_states
    off2_dn = 3*n_ref_states

    !=========================================================================
    ! Onsite block sublattice 1: sigma1(alpha) (same for both spins)
    !=========================================================================

    DO i_orb = 1, n_ref_states
      HK( off1_up + i_orb, off1_up + i_orb ) = sigma1(i_orb)
      HK( off1_dn + i_orb, off1_dn + i_orb ) = sigma1(i_orb)
    END DO

    !=========================================================================
    ! Onsite block sublattice 2: sigma2(alpha)
    !=========================================================================

    DO i_orb = 1, n_ref_states
      HK( off2_up + i_orb, off2_up + i_orb ) = sigma2(i_orb)
      HK( off2_dn + i_orb, off2_dn + i_orb ) = sigma2(i_orb)
    END DO

    !=========================================================================
    ! Off-diagonal block: hopping sl1->sl2, sum over 4 1NN bonds with Bloch phase
    !=========================================================================

    CALL get_1nn_bond_vectors( bonds_frac )

    DO i_bond = 1, 4

      cos_latt = bonds_frac(:,i_bond) / SQRT( SUM(bonds_frac(:,i_bond)**2) )

      ! Phase: 2*pi * k_frac . bond_frac (same as sparse_ham convention)
      phase = 2.0_dp * pi * DOT_PRODUCT( kvec, bonds_frac(:,i_bond) )
      exp_fac = EXP( CMPLX(0.0_dp, phase, dp) )

      CALL koster_slater_min( cos_latt, params%hopping_vca, tb_bloc )

      DO i_orb = 1, n_ref_states
        DO j_orb = 1, n_ref_states
          ! sl1_up -> sl2_up
          HK( off1_up + i_orb, off2_up + j_orb ) = &
            HK( off1_up + i_orb, off2_up + j_orb ) + &
            exp_fac * CMPLX(tb_bloc(i_orb, j_orb), 0.0_dp, dp)
          ! sl1_dn -> sl2_dn
          HK( off1_dn + i_orb, off2_dn + j_orb ) = &
            HK( off1_dn + i_orb, off2_dn + j_orb ) + &
            exp_fac * CMPLX(tb_bloc(i_orb, j_orb), 0.0_dp, dp)
        END DO
      END DO

    END DO

    !=========================================================================
    ! Hermitian conjugate blocks (sl2 -> sl1)
    !=========================================================================

    DO i_orb = 1, n_ref_states
      DO j_orb = 1, n_ref_states
        HK( off2_up + j_orb, off1_up + i_orb ) = &
          CONJG( HK( off1_up + i_orb, off2_up + j_orb ) )
        HK( off2_dn + j_orb, off1_dn + i_orb ) = &
          CONJG( HK( off1_dn + i_orb, off2_dn + j_orb ) )
      END DO
    END DO

    !=========================================================================
    ! p-orbital SOC blocks (VCA): added to both sublattice 1 and 2
    ! SOC is intra-atomic (no coupling between sl1 and sl2 in SOC term)
    !=========================================================================

    CALL add_soc_p_block( params%so_p_sl1_vca, off1_up, off1_dn, HK )
    CALL add_soc_p_block( params%so_p_sl2_vca, off2_up, off2_dn, HK )

  END SUBROUTINE build_multi_HK


  !===========================================================================
  ! Subroutine add_soc_p_block
  !
  ! Add p-orbital spin-orbit coupling to the Hamiltonian block of one atom.
  ! Ported from build_hk.f90::add_soc_p_block (verified correct).
  !
  ! INPUT:
  !   so_p           : effective p-SOC parameter (Delta_p)
  !   off_up, off_dn : row/col offset for spin-up / spin-down blocks
  !
  ! INPUT/OUTPUT:
  !   HK : Hamiltonian matrix (modified in-place)
  !===========================================================================

  SUBROUTINE add_soc_p_block( so_p, off_up, off_dn, HK )

    REAL(dp),    INTENT(IN)    :: so_p
    INTEGER,     INTENT(IN)    :: off_up, off_dn
    COMPLEX(dp), DIMENSION(:,:), INTENT(INOUT) :: HK

    COMPLEX(dp) :: ci
    ci = CMPLX(0.0_dp, 1.0_dp, dp)

    !------- up-up -------
    HK(off_up+px, off_up+py) = HK(off_up+px, off_up+py) - ci*so_p
    HK(off_up+py, off_up+px) = HK(off_up+py, off_up+px) + ci*so_p

    !------- down-down -------
    HK(off_dn+px, off_dn+py) = HK(off_dn+px, off_dn+py) + ci*so_p
    HK(off_dn+py, off_dn+px) = HK(off_dn+py, off_dn+px) - ci*so_p

    !------- up-down -------
    HK(off_up+px, off_dn+pz) = HK(off_up+px, off_dn+pz) - so_p
    HK(off_up+py, off_dn+pz) = HK(off_up+py, off_dn+pz) + ci*so_p
    HK(off_up+pz, off_dn+px) = HK(off_up+pz, off_dn+px) + so_p
    HK(off_up+pz, off_dn+py) = HK(off_up+pz, off_dn+py) - ci*so_p

    !------- down-up -------
    HK(off_dn+px, off_up+pz) = HK(off_dn+px, off_up+pz) + so_p
    HK(off_dn+py, off_up+pz) = HK(off_dn+py, off_up+pz) + ci*so_p
    HK(off_dn+pz, off_up+px) = HK(off_dn+pz, off_up+px) - so_p
    HK(off_dn+pz, off_up+py) = HK(off_dn+pz, off_up+py) - ci*so_p

  END SUBROUTINE add_soc_p_block


  !===========================================================================
  ! Function get_multi_HK_ndim
  !
  ! Return the Bloch Hamiltonian dimension for the current params. When
  ! odd_hopping = .FALSE., this is the fixed 4*n_ref_states (2 sublattices
  ! x spin), matching build_multi_HK. When odd_hopping = .TRUE., the basis
  ! is expanded to one block per ACTIVE species per sublattice (BEB scheme):
  !   ndim = n_spin * n_ref_states * (n_active_sl1 + n_active_sl2)
  ! Species with zero concentration are excluded from the augmented basis
  ! entirely (see cpa_multi_types.f90), so n_active_sl*, not
  ! n_real_elem_sl*, sets the enlarged dimension.
  !===========================================================================

  FUNCTION get_multi_HK_ndim( params ) RESULT( ndim )

    TYPE(cpa_multi_params), INTENT(IN) :: params
    INTEGER :: ndim

    IF ( params%odd_hopping ) THEN
      ndim = 2 * n_ref_states * &
             ( params%n_active_sl1 + params%n_active_sl2 )
    ELSE
      ndim = 4 * n_ref_states
    END IF

  END FUNCTION get_multi_HK_ndim


  !===========================================================================
  ! Subroutine build_multi_HK_odd
  !
  ! Build the enlarged Bloch Hamiltonian under off-diagonal disorder (ODD)
  ! for hopping, following the BEB extended-Hilbert-space scheme (Blackman
  ! 1971; Papaconstantopoulos 1989 multiband generalization; Koepernik
  ! 1997/1998 multi-sublattice generalization).
  !
  ! Basis layout: one spin-up/spin-down block of size n_ref_states per
  ! ACTIVE (nonzero-concentration) REAL ELEMENT per sublattice, ordered as
  !   [ sl1_elem_1_up, sl1_elem_1_dn, sl1_elem_2_up, sl1_elem_2_dn, ...,
  !     sl2_elem_1_up, sl2_elem_1_dn, sl2_elem_2_up, sl2_elem_2_dn, ... ]
  ! Total dimension = get_multi_HK_ndim(params).
  !
  ! Onsite blocks: the BEB coherent potential is a FULL local matrix in
  ! chemical-species space. sigma1((Q-1)*n1+Qp,alpha) is Sigma_1^{Q,Qp}
  ! (and analogously for sublattice 2), so Q /= Qp couples chemical blocks
  ! on the SAME physical site. These entries are the BEB interactor, not
  ! physical same-sublattice hopping; see Blackman Eq. (3.14)-(3.22),
  ! Papaconstantopoulos Eq. (2.8), and Koepernik et al. Eq. (62)-(65).
  !
  ! Off-diagonal (hopping) blocks: for each ACTIVE species pair (Q on sl1,
  ! Q' on sl2), the block h^{Q,Q'}(k) is built from THAT PAIR's OWN
  ! hopping parameters (params%pairs(:)%t), summed over the 4 zinc-blende
  ! 1NN bonds with Bloch phase -- no averaging across species pairs.
  !
  ! SOC: unchanged, VCA (so_p_sl1_vca / so_p_sl2_vca), applied identically
  ! to every species block of the corresponding sublattice (SOC never
  ! enters the Soven loop, see cpa_multi_types.f90).
  !
  ! INPUT:
  !   kvec(3) : k-vector in fractional coordinates
  !   params  : cpa_multi_params with odd_hopping = .TRUE. and pairs(:) set
  !   sigma1(n_active_sl1**2,n_ref_states) : BEB chemical-space potential,
  !                                           sublattice 1
  !   sigma2(n_active_sl2**2,n_ref_states) : BEB chemical-space potential,
  !                                           sublattice 2
  !
  ! OUTPUT:
  !   HK(ndim,ndim) : Bloch Hamiltonian, COMPLEX, ndim = get_multi_HK_ndim(params)
  !===========================================================================

  SUBROUTINE build_multi_HK_odd( kvec, params, sigma1, sigma2, HK )

    REAL(dp), DIMENSION(3),     INTENT(IN) :: kvec
    TYPE(cpa_multi_params),     INTENT(IN) :: params
    COMPLEX(dp), DIMENSION(:,:), INTENT(IN) :: sigma1, sigma2

    COMPLEX(dp), DIMENSION(:,:), INTENT(OUT) :: HK

    !=========================================================================
    ! Local variables
    !=========================================================================

    REAL(dp), DIMENSION(3,4) :: bonds_frac
    REAL(dp), DIMENSION(3)   :: cos_latt
    REAL(dp)                 :: phase

    REAL(dp), DIMENSION(n_ref_states, n_ref_states) :: tb_bloc

    COMPLEX(dp) :: exp_fac

    INTEGER :: i_bond, i_orb, j_orb
    INTEGER :: n1, n2, i_pair, q1, q2, qp1, qp2
    INTEGER :: off1_up, off1_dn, off2_up, off2_dn
    INTEGER, ALLOCATABLE :: off_up_sl1(:), off_dn_sl1(:)
    INTEGER, ALLOCATABLE :: off_up_sl2(:), off_dn_sl2(:)

    !=========================================================================
    ! Block offsets: 2 blocks (up,dn) x n_ref_states per ACTIVE species,
    ! sl1 species first, then sl2 species. Offsets are 0-based row/col
    ! base indices. Zero-concentration species have no block at all.
    !=========================================================================

    n1 = params%n_active_sl1
    n2 = params%n_active_sl2

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

    HK = (0.0_dp, 0.0_dp)

    !=========================================================================
    ! Local BEB coherent potential. Preserve Q /= Qp elements: the
    ! species-space interactor is local, although zinc-blende has no
    ! physical same-sublattice hopping.
    !=========================================================================

    DO q1 = 1, n1
      DO qp1 = 1, n1
        DO i_orb = 1, n_ref_states
          HK( off_up_sl1(q1)+i_orb, off_up_sl1(qp1)+i_orb ) = &
            sigma1((q1-1)*n1+qp1,i_orb)
          HK( off_dn_sl1(q1)+i_orb, off_dn_sl1(qp1)+i_orb ) = &
            sigma1((q1-1)*n1+qp1,i_orb)
        END DO
      END DO
    END DO

    DO q2 = 1, n2
      DO qp2 = 1, n2
        DO i_orb = 1, n_ref_states
          HK( off_up_sl2(q2)+i_orb, off_up_sl2(qp2)+i_orb ) = &
            sigma2((q2-1)*n2+qp2,i_orb)
          HK( off_dn_sl2(q2)+i_orb, off_dn_sl2(qp2)+i_orb ) = &
            sigma2((q2-1)*n2+qp2,i_orb)
        END DO
      END DO
    END DO

    !=========================================================================
    ! Off-diagonal (hopping) blocks: one per species pair, own parameters
    !=========================================================================

    CALL get_1nn_bond_vectors( bonds_frac )

    DO i_pair = 1, SIZE(params%pairs)

      q1 = params%pairs(i_pair)%i_elem1
      q2 = params%pairs(i_pair)%i_elem2

      off1_up = off_up_sl1(q1);  off1_dn = off_dn_sl1(q1)
      off2_up = off_up_sl2(q2);  off2_dn = off_dn_sl2(q2)

      DO i_bond = 1, 4

        cos_latt = bonds_frac(:,i_bond) / SQRT( SUM(bonds_frac(:,i_bond)**2) )

        phase = 2.0_dp * pi * DOT_PRODUCT( kvec, bonds_frac(:,i_bond) )
        exp_fac = EXP( CMPLX(0.0_dp, phase, dp) )

        CALL koster_slater_min( cos_latt, params%pairs(i_pair)%t, tb_bloc )

        DO i_orb = 1, n_ref_states
          DO j_orb = 1, n_ref_states
            HK( off1_up + i_orb, off2_up + j_orb ) = &
              HK( off1_up + i_orb, off2_up + j_orb ) + &
              exp_fac * CMPLX(tb_bloc(i_orb, j_orb), 0.0_dp, dp)
            HK( off1_dn + i_orb, off2_dn + j_orb ) = &
              HK( off1_dn + i_orb, off2_dn + j_orb ) + &
              exp_fac * CMPLX(tb_bloc(i_orb, j_orb), 0.0_dp, dp)
          END DO
        END DO

      END DO

      DO i_orb = 1, n_ref_states
        DO j_orb = 1, n_ref_states
          HK( off2_up + j_orb, off1_up + i_orb ) = &
            CONJG( HK( off1_up + i_orb, off2_up + j_orb ) )
          HK( off2_dn + j_orb, off1_dn + i_orb ) = &
            CONJG( HK( off1_dn + i_orb, off2_dn + j_orb ) )
        END DO
      END DO

    END DO

    !=========================================================================
    ! SOC (VCA, unchanged): applied identically to every species block
    !=========================================================================

    DO q1 = 1, n1
      CALL add_soc_p_block( params%so_p_sl1_vca, off_up_sl1(q1), off_dn_sl1(q1), HK )
    END DO
    DO q2 = 1, n2
      CALL add_soc_p_block( params%so_p_sl2_vca, off_up_sl2(q2), off_dn_sl2(q2), HK )
    END DO

    DEALLOCATE( off_up_sl1, off_dn_sl1, off_up_sl2, off_dn_sl2 )

  END SUBROUTINE build_multi_HK_odd

END MODULE cpa_multi_hk
