! cpa_multi_build.f90
!
! Build fully-populated cpa_multi_params from CPA config.
!
!   Onsite at a_vegard: ion%offset + ion%energy(alpha)  [bare]
!     + Tan onsite correction (fac_I, fac_O) over 4 1NN bonds [Tan only]
!   ion%offset = etb_file_offset + E_v  (E_v read from INIshell .dat)
!
!   SOC: ion%so_energy(1)  [bare]
!     + 4 * so_corr_per_bond (fac_so) [Tan only]
!
!   Hopping at dist_vegard:
!     Jancu: energy * (dist_ref/dist_vegard)^scaling
!     Tan:   energy * EXP(-scaling * (dist_vegard - dist_ref + delta_d))
!
MODULE cpa_multi_build

  USE precision,   ONLY : dp
  USE globals,     ONLY : LST, n_ref_states, n_ref_couplings
  USE type_defs,   ONLY : material_data, parent_data, ion_orbit, pair_coupling, &
                          destroy_material
  USE input_data,  ONLY : read_data
  USE states_and_couplings, ONLY : ref_states_and_couplings, sort_states
  USE cpa_config,  ONLY : cpa_config_t
  USE cpa_multi_types, ONLY : cpa_multi_params, pseudo_element
  USE globals,     ONLY : sps, pss, sses, sess, seps, pses, sds, dss, &
                          seds, dses, pds, dps, pdp, dpp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: build_multi_cpa_params

CONTAINS

  !===========================================================================
  SUBROUTINE build_multi_cpa_params( cfg, params )

    TYPE(cpa_config_t),     INTENT(IN)  :: cfg
    TYPE(cpa_multi_params), INTENT(OUT) :: params

    INTEGER :: i, j, idx, n_pseudo
    REAL(dp) :: conc_ij
    CHARACTER(LEN=LST) :: work_path
    TYPE(material_data) :: mat_ij

    CHARACTER(LEN=5), DIMENSION(:), POINTER :: ref_states    => NULL()
    CHARACTER(LEN=5), DIMENSION(:), POINTER :: ref_couplings => NULL()
    INTEGER :: n_ref_st, n_ref_cpl

    REAL(dp) :: a_ij
    REAL(dp) :: onsite1(n_ref_states), onsite2(n_ref_states)
    REAL(dp) :: so_p1, so_p2
    REAL(dp) :: t_ij(n_ref_couplings)
    INTEGER  :: i_ion_sl1, i_ion_sl2, n_st_ion1, n_st_ion2, n_so, alpha

    CALL ref_states_and_couplings( ref_states, n_ref_st, ref_couplings, n_ref_cpl )

    params%n_sl1     = cfg%n_sl1
    params%n_sl2     = cfg%n_sl2
    params%n_pseudo1 = cfg%n_sl1 * cfg%n_sl2
    params%n_pseudo2 = cfg%n_sl1 * cfg%n_sl2
    params%scheme    = TRIM(cfg%scheme)

    n_pseudo = params%n_pseudo1
    ALLOCATE( params%pseudo1(n_pseudo), params%pseudo2(n_pseudo) )
    ALLOCATE( params%hopping_vca(n_ref_couplings) )

    DO idx = 1, n_pseudo
      ALLOCATE( params%pseudo1(idx)%onsite(n_ref_states) )
      ALLOCATE( params%pseudo2(idx)%onsite(n_ref_states) )
      params%pseudo1(idx)%onsite = 0.0_dp
      params%pseudo2(idx)%onsite = 0.0_dp
      params%pseudo1(idx)%so_p   = 0.0_dp
      params%pseudo2(idx)%so_p   = 0.0_dp
      params%pseudo1(idx)%conc   = 0.0_dp
      params%pseudo2(idx)%conc   = 0.0_dp
    END DO
    params%hopping_vca  = 0.0_dp
    params%so_p_sl1_vca = 0.0_dp
    params%so_p_sl2_vca = 0.0_dp

    !--- Vegard's law lattice constant ---
    params%a_vegard = 0.0_dp
    work_path = './'
    DO i = 1, cfg%n_sl1
      DO j = 1, cfg%n_sl2
        conc_ij = cfg%x1(i) * cfg%x2(j)
        a_ij    = read_lattice_param( cfg, i, j )
        params%a_vegard = params%a_vegard + conc_ij * a_ij
      END DO
    END DO
    params%dist_vegard = params%a_vegard * SQRT(3.0_dp) / 4.0_dp
    WRITE(*,'(a,f12.6)') ' (cpa_multi_build) a_vegard  (Ang) = ', params%a_vegard
    WRITE(*,'(a,f12.6)') ' (cpa_multi_build) d_vegard  (Ang) = ', params%dist_vegard

    !--- Per-binary loop ---
    DO i = 1, cfg%n_sl1
      DO j = 1, cfg%n_sl2

        idx     = (i-1)*cfg%n_sl2 + j
        conc_ij = cfg%x1(i) * cfg%x2(j)

        params%pseudo1(idx)%idx_self  = i
        params%pseudo1(idx)%idx_other = j
        params%pseudo1(idx)%conc      = conc_ij
        params%pseudo2(idx)%idx_self  = j
        params%pseudo2(idx)%idx_other = i
        params%pseudo2(idx)%conc      = conc_ij

        ! Init mat_ij: .etb filename + E_v from .dat
        CALL init_mat_data_for_binary( cfg, i, j, mat_ij )

        ! read_data searches etb_path for the .etb file
        CALL read_data( mat_ij, work_path, cfg%etb_path )
        CALL sort_states( mat_ij, ref_states, ref_couplings )

        CALL find_ion_indices( mat_ij, cfg%sl1(i), cfg%sl2(j), i_ion_sl1, i_ion_sl2 )

        n_st_ion1 = SIZE(mat_ij%parent(1)%ion(i_ion_sl1)%energy)
        n_st_ion2 = SIZE(mat_ij%parent(1)%ion(i_ion_sl2)%energy)

        ! Bare onsite = offset + energy
        ! ion%offset already includes E_v (added by read_data: offset += e_v)
        onsite1 = 0.0_dp
        onsite2 = 0.0_dp
        DO alpha = 1, n_st_ion1
          onsite1(alpha) = mat_ij%parent(1)%ion(i_ion_sl1)%offset + &
                           mat_ij%parent(1)%ion(i_ion_sl1)%energy(alpha)
        END DO
        DO alpha = 1, n_st_ion2
          onsite2(alpha) = mat_ij%parent(1)%ion(i_ion_sl2)%offset + &
                           mat_ij%parent(1)%ion(i_ion_sl2)%energy(alpha)
        END DO

        ! Tan: add 4-bond onsite correction
        IF ( TRIM(cfg%scheme) == 'tan' ) THEN
          CALL add_tan_onsite_correction( mat_ij, i_ion_sl1, &
               params%dist_vegard, n_st_ion1, onsite1 )
          CALL add_tan_onsite_correction( mat_ij, i_ion_sl2, &
               params%dist_vegard, n_st_ion2, onsite2 )
        END IF

        params%pseudo1(idx)%onsite = onsite1
        params%pseudo2(idx)%onsite = onsite2

        ! SOC: bare + Tan correction
        n_so = SIZE(mat_ij%parent(1)%ion(i_ion_sl1)%so_energy)
        so_p1 = MERGE(mat_ij%parent(1)%ion(i_ion_sl1)%so_energy(1), 0.0_dp, n_so >= 1)
        n_so = SIZE(mat_ij%parent(1)%ion(i_ion_sl2)%so_energy)
        so_p2 = MERGE(mat_ij%parent(1)%ion(i_ion_sl2)%so_energy(1), 0.0_dp, n_so >= 1)

        IF ( TRIM(cfg%scheme) == 'tan' ) THEN
          CALL add_tan_so_correction( mat_ij, cfg%sl1(i), so_p1 )
          CALL add_tan_so_correction( mat_ij, cfg%sl2(j), so_p2 )
        END IF

        params%pseudo1(idx)%so_p = so_p1
        params%pseudo2(idx)%so_p = so_p2

        ! Hopping VCA accumulation
        CALL compute_hopping_at_vegard( mat_ij, cfg%scheme, i_ion_sl1, i_ion_sl2, &
                                        params%dist_vegard, t_ij )
        params%hopping_vca  = params%hopping_vca  + conc_ij * t_ij
        params%so_p_sl1_vca = params%so_p_sl1_vca + conc_ij * so_p1
        params%so_p_sl2_vca = params%so_p_sl2_vca + conc_ij * so_p2

        WRITE(*,'(a,a,a,f8.4,a,f8.5,a,f8.5)') &
          ' (cpa_multi_build) ', TRIM(cfg%sl1(i))//TRIM(cfg%sl2(j)), &
          '  conc=', conc_ij, &
          '  so_p(sl1)=', so_p1, '  so_p(sl2)=', so_p2

        CALL destroy_material( mat_ij )

      END DO
    END DO

    WRITE(*,'(a,f10.6)') ' (cpa_multi_build) so_p_sl1_vca = ', params%so_p_sl1_vca
    WRITE(*,'(a,f10.6)') ' (cpa_multi_build) so_p_sl2_vca = ', params%so_p_sl2_vca
    WRITE(*,'(a,5f10.4)') ' (cpa_multi_build) hopping(1..5) = ', params%hopping_vca(1:5)

    DEALLOCATE( ref_states, ref_couplings )

  END SUBROUTINE build_multi_cpa_params


  !===========================================================================
  ! Read lattice constant (nm -> Angstrom) from INIshell .dat file
  !===========================================================================

  FUNCTION read_lattice_param( cfg, i, j ) RESULT(a_lat)
    TYPE(cpa_config_t), INTENT(IN) :: cfg
    INTEGER,            INTENT(IN) :: i, j
    REAL(dp) :: a_lat

    CHARACTER(LEN=LST) :: fname, line, lhs, rhs
    INTEGER :: fnum, ios, eq_pos
    LOGICAL :: in_lattice_section
    REAL(dp) :: a_nm

    fname = TRIM(cfg%dat_path)//TRIM(cfg%sl1(i))//TRIM(cfg%sl2(j))//'.dat'
    OPEN(NEWUNIT=fnum, FILE=TRIM(fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF ( ios /= 0 ) THEN
      WRITE(*,*) 'ERROR (cpa_multi_build): cannot open ', TRIM(fname); STOP 1
    END IF

    a_lat = 0.0_dp; in_lattice_section = .FALSE.
    DO
      READ(fnum, '(A)', IOSTAT=ios) line;  IF ( ios /= 0 ) EXIT
      eq_pos = INDEX(line,'#');  IF ( eq_pos > 0 ) line = line(1:eq_pos-1)
      line = ADJUSTL(line);  IF ( LEN_TRIM(line) == 0 ) CYCLE
      IF ( line(1:1) == '[' ) THEN
        in_lattice_section = ( INDEX(line,'[lattice]') > 0 );  CYCLE
      END IF
      IF ( in_lattice_section ) THEN
        eq_pos = INDEX(line,'=');  IF ( eq_pos < 2 ) CYCLE
        lhs = ADJUSTL(line(1:eq_pos-1));  rhs = ADJUSTL(line(eq_pos+1:))
        IF ( TRIM(lhs) == 'a' ) THEN
          READ(rhs,*,IOSTAT=ios) a_nm
          IF ( ios /= 0 ) THEN
            WRITE(*,*) 'ERROR: cannot read lattice constant from ',TRIM(fname); STOP 1
          END IF
          a_lat = a_nm * 10.0_dp   ! nm -> Angstrom
          CLOSE(fnum);  RETURN
        END IF
      END IF
    END DO
    CLOSE(fnum)
    IF ( ABS(a_lat) < 1.0d-10 ) THEN
      WRITE(*,*) 'ERROR: lattice constant not found in ',TRIM(fname);  STOP 1
    END IF
  END FUNCTION read_lattice_param


  !===========================================================================
  ! Read E_v (eV) from INIshell .dat file ('E_v = ...' anywhere in file)
  !===========================================================================

  FUNCTION read_ev_from_dat( cfg, i, j ) RESULT(ev)
    TYPE(cpa_config_t), INTENT(IN) :: cfg
    INTEGER,            INTENT(IN) :: i, j
    REAL(dp) :: ev

    CHARACTER(LEN=LST) :: fname, line, lhs, rhs
    INTEGER :: fnum, ios, eq_pos

    fname = TRIM(cfg%dat_path)//TRIM(cfg%sl1(i))//TRIM(cfg%sl2(j))//'.dat'
    OPEN(NEWUNIT=fnum, FILE=TRIM(fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF ( ios /= 0 ) THEN
      WRITE(*,*) 'ERROR (cpa_multi_build): cannot open ', TRIM(fname); STOP 1
    END IF

    ev = 0.0_dp
    DO
      READ(fnum,'(A)',IOSTAT=ios) line;  IF ( ios /= 0 ) EXIT
      eq_pos = INDEX(line,'#');  IF ( eq_pos > 0 ) line = line(1:eq_pos-1)
      line = ADJUSTL(line);  IF ( LEN_TRIM(line) == 0 ) CYCLE
      eq_pos = INDEX(line,'=');  IF ( eq_pos < 2 ) CYCLE
      lhs = ADJUSTL(line(1:eq_pos-1));  rhs = ADJUSTL(line(eq_pos+1:))
      IF ( TRIM(lhs) == 'E_v' ) THEN
        READ(rhs,*,IOSTAT=ios) ev
        IF ( ios /= 0 ) ev = 0.0_dp
        CLOSE(fnum);  RETURN
      END IF
    END DO
    CLOSE(fnum)
  END FUNCTION read_ev_from_dat


  !===========================================================================
  ! Initialize material_data for one binary:
  !   - .etb filename set so read_data finds it in etb_path
  !   - e_v read from .dat and set so read_data applies: offset += e_v
  !===========================================================================

  SUBROUTINE init_mat_data_for_binary( cfg, i, j, mat )
    TYPE(cpa_config_t),  INTENT(IN)  :: cfg
    INTEGER,             INTENT(IN)  :: i, j
    TYPE(material_data), INTENT(OUT) :: mat

    NULLIFY(mat%nr_parents, mat%parent, mat%ion, mat%bowing)
    ALLOCATE( mat%nr_parents(1) )
    mat%nr_parents(1) = 1
    ALLOCATE( mat%parent(1) )
    CALL nullify_ion_array_dummy( mat%parent(1) )

    mat%parent(1)%data_file = TRIM(cfg%sl1(i))//TRIM(cfg%sl2(j))//'.etb'
    mat%parent(1)%content   = 1.0_dp
    ! E_v from .dat file: read_data adds e_v to ion%offset (offset += e_v)
    mat%parent(1)%e_v       = read_ev_from_dat( cfg, i, j )
    WRITE(*,'(a,a,a,f8.4)') ' (cpa_multi_build) E_v for ', &
      TRIM(cfg%sl1(i))//TRIM(cfg%sl2(j)), ' = ', mat%parent(1)%e_v

  END SUBROUTINE init_mat_data_for_binary


  SUBROUTINE nullify_ion_array_dummy( parent )
    USE type_defs, ONLY : parent_data
    TYPE(parent_data), INTENT(INOUT) :: parent
    NULLIFY(parent%ion, parent%pair)
    parent%n_ion   = 0
    parent%n_pair  = 0
    parent%e_v     = 0.0_dp
    parent%content = 1.0_dp
  END SUBROUTINE nullify_ion_array_dummy


  !===========================================================================
  ! Find ion indices for sl1_name and sl2_name within the parent
  !===========================================================================

  SUBROUTINE find_ion_indices( mat, name1, name2, i_ion1, i_ion2 )
    TYPE(material_data), INTENT(IN)  :: mat
    CHARACTER(LEN=*),    INTENT(IN)  :: name1, name2
    INTEGER,             INTENT(OUT) :: i_ion1, i_ion2
    INTEGER :: k
    i_ion1 = 0;  i_ion2 = 0
    DO k = 1, mat%parent(1)%n_ion
      IF ( TRIM(mat%parent(1)%ion(k)%name) == TRIM(name1) ) i_ion1 = k
      IF ( TRIM(mat%parent(1)%ion(k)%name) == TRIM(name2) ) i_ion2 = k
    END DO
    IF ( i_ion1 == 0 ) THEN
      WRITE(*,*) 'ERROR: ion "'//TRIM(name1)//'" not found in binary ETB'; STOP 1
    END IF
    IF ( i_ion2 == 0 ) THEN
      WRITE(*,*) 'ERROR: ion "'//TRIM(name2)//'" not found in binary ETB'; STOP 1
    END IF
  END SUBROUTINE find_ion_indices


  !===========================================================================
  ! Tan onsite correction — exact replication of TB_ham::get_ons_corr
  ! accumulated over 4 equivalent 1NN bonds.
  !===========================================================================

  SUBROUTINE add_tan_onsite_correction( mat, i_ion, dist_vegard, n_st, onsite )
    TYPE(material_data), INTENT(IN)    :: mat
    INTEGER,             INTENT(IN)    :: i_ion
    REAL(dp),            INTENT(IN)    :: dist_vegard
    INTEGER,             INTENT(IN)    :: n_st
    REAL(dp),            INTENT(INOUT) :: onsite(n_ref_states)

    TYPE(pair_coupling), POINTER :: pair
    TYPE(ion_orbit),     POINTER :: ion
    INTEGER  :: ip, n_orb_b, i_state, ii, iii, i_offset
    REAL(dp) :: dist_phase, ons_corr(n_ref_states)

    DO ip = 1, mat%parent(1)%n_pair
      pair => mat%parent(1)%pair(ip)
      IF ( .NOT. ASSOCIATED(pair%fac_I) ) CYCLE
      IF ( .NOT. ASSOCIATED(pair%l_I)  ) CYCLE

      ion     => mat%parent(1)%ion(i_ion)
      n_orb_b  = SIZE(ion%deg)

      ! i_offset: 0 if atom is pair%name(2); n_orb_b if atom is pair%name(1)
      IF ( TRIM(ion%name) /= TRIM(pair%name(1)) ) THEN
        i_offset = 0
      ELSE
        i_offset = n_orb_b
      END IF

      dist_phase = dist_vegard - pair%dist_ref + pair%delta_d

      ons_corr = 0.0_dp
      i_state  = 0
      DO ii = 1, n_orb_b
        DO iii = 1, ion%deg(ii)
          i_state = i_state + 1
          IF ( i_state > n_st ) EXIT
          ons_corr(i_state) = pair%fac_I(i_offset+ii) * &
                              EXP( -pair%l_I(i_offset+ii) * dist_phase )
        END DO
        IF ( i_state > n_st ) EXIT
      END DO

      DO i_state = 1, n_st
        ons_corr(i_state) = ons_corr(i_state) + pair%fac_O * EXP( -pair%l_O * dist_phase )
      END DO

      ! 4 equivalent bonds
      onsite(1:n_st) = onsite(1:n_st) + 4.0_dp * ons_corr(1:n_st)

    END DO
  END SUBROUTINE add_tan_onsite_correction


  !===========================================================================
  ! Tan SOC correction — exact replication of TB_ham::get_so_corr
  ! accumulated over 4 equivalent 1NN bonds.
  !===========================================================================

  SUBROUTINE add_tan_so_correction( mat, atom_name, so_p )
    TYPE(material_data), INTENT(IN)    :: mat
    CHARACTER(LEN=*),    INTENT(IN)    :: atom_name
    REAL(dp),            INTENT(INOUT) :: so_p

    TYPE(pair_coupling), POINTER :: pair
    INTEGER :: ip, n_so

    DO ip = 1, mat%parent(1)%n_pair
      pair => mat%parent(1)%pair(ip)
      IF ( .NOT. ASSOCIATED(pair%fac_so) ) CYCLE
      n_so = SIZE(pair%fac_so) / 2
      IF ( n_so < 1 ) CYCLE
      ! name(1) -> fac_so(n_so+1); name(2) -> fac_so(1)
      IF ( TRIM(atom_name) == TRIM(pair%name(1)) ) THEN
        so_p = so_p + 4.0_dp * pair%fac_so(n_so + 1)
      ELSE
        so_p = so_p + 4.0_dp * pair%fac_so(1)
      END IF
    END DO
  END SUBROUTINE add_tan_so_correction


  !===========================================================================
  ! Compute hopping at dist_vegard for one binary
  !===========================================================================

  SUBROUTINE compute_hopping_at_vegard( mat, scheme, i_ion_sl1, i_ion_sl2, &
                                        dist_vegard, t_out )
    TYPE(material_data), INTENT(IN)  :: mat
    CHARACTER(LEN=*),    INTENT(IN)  :: scheme
    INTEGER,             INTENT(IN)  :: i_ion_sl1, i_ion_sl2
    REAL(dp),            INTENT(IN)  :: dist_vegard
    REAL(dp),            INTENT(OUT) :: t_out(n_ref_couplings)

    TYPE(pair_coupling), POINTER :: pair
    INTEGER  :: n_cpl, k, ref_pair(n_ref_couplings)
    REAL(dp) :: dist_ref, dist_phase

    t_out = 0.0_dp
    IF ( mat%parent(1)%n_pair < 1 ) RETURN

    pair     => mat%parent(1)%pair(1)
    n_cpl    =  SIZE(pair%energy)
    dist_ref =  pair%dist_ref

    DO k = 1, n_ref_couplings;  ref_pair(k) = k;  END DO

    ! Swap if pair%name(1) is sl2 atom (not sl1)
    IF ( TRIM(pair%name(1)) /= TRIM(mat%parent(1)%ion(i_ion_sl1)%name) ) THEN
      CALL apply_swap_permutation( ref_pair )
    END IF

    IF ( TRIM(scheme) == 'jancu' ) THEN
      DO k = 1, n_cpl
        t_out(ref_pair(k)) = pair%energy(k) * (dist_ref/dist_vegard)**pair%scaling(k)
      END DO
    ELSE IF ( TRIM(scheme) == 'tan' ) THEN
      dist_phase = dist_vegard - dist_ref + pair%delta_d
      DO k = 1, n_cpl
        t_out(ref_pair(k)) = pair%energy(k) * EXP( -pair%scaling(k) * dist_phase )
      END DO
    END IF

  END SUBROUTINE compute_hopping_at_vegard


  !===========================================================================
  ! Antisymmetric-pair swap permutation (mirrors check_ref_pair in checks.f90)
  !===========================================================================

  SUBROUTINE apply_swap_permutation( ref_pair )
    INTEGER, INTENT(INOUT) :: ref_pair(n_ref_couplings)
    INTEGER :: tmp
    tmp=ref_pair(sps) ; ref_pair(sps) =ref_pair(pss) ; ref_pair(pss) =tmp
    tmp=ref_pair(sses); ref_pair(sses)=ref_pair(sess); ref_pair(sess)=tmp
    tmp=ref_pair(seps); ref_pair(seps)=ref_pair(pses); ref_pair(pses)=tmp
    tmp=ref_pair(sds) ; ref_pair(sds) =ref_pair(dss) ; ref_pair(dss) =tmp
    tmp=ref_pair(seds); ref_pair(seds)=ref_pair(dses); ref_pair(dses)=tmp
    tmp=ref_pair(pds) ; ref_pair(pds) =ref_pair(dps) ; ref_pair(dps) =tmp
    tmp=ref_pair(pdp) ; ref_pair(pdp) =ref_pair(dpp) ; ref_pair(dpp) =tmp
  END SUBROUTINE apply_swap_permutation

END MODULE cpa_multi_build
