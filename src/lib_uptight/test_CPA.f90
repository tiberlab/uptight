! test_CPA.f90
!
! Driver program for multi-component diagonal CPA on zinc-blende alloys,
! under the VCA-before-CPA scheme (see implementation_plan.md,
! "Critical Revision").
!
! Usage:
!   ./test_CPA config.txt
!
! Config file format (see cpa_config.f90 for full documentation):
!   etb_path  = './Tan/'
!   dat_path  = './datafiles/'
!   sl1 = Al Ga In
!   x1  = 0.3 0.8 0.1
!   sl2 = N P As Sb
!   x2  = 0.1 0.2 0.3 0.4
!   E_min = -5.0
!   E_max =  5.0
!   n_E = 201
!   eta = 0.02
!   nk_bz_x = 8
!   nk_bz_y = 8
!   nk_bz_z = 8
!   n_per_segment = 100
!   tol = 1.0e-6
!   max_iter = 200
!   mixing_alpha = 0.5
!   output_filename = 'spectral_function.dat'   ! optional; auto-named if absent
!
PROGRAM test_CPA

  USE precision,            ONLY : dp, set_machine_acc
  USE globals,              ONLY : n_ref_states
  USE cpa_config,           ONLY : cpa_config_t, read_and_validate_config, &
                                    build_output_filename
  USE cpa_multi_types,      ONLY : cpa_multi_params, destroy_cpa_multi_params
  USE cpa_multi_build,      ONLY : build_multi_cpa_params
  USE kmesh_gen,            ONLY : generate_mp_mesh, generate_kpath_LGX
  USE cpa_multi_solver,     ONLY : cpa_multi_compute_spectral_function

  IMPLICIT NONE

  !=============================================================================
  ! Variables
  !=============================================================================

  CHARACTER(LEN=500) :: config_file

  TYPE(cpa_config_t)    :: cfg
  TYPE(cpa_multi_params) :: params

  ! Reciprocal lattice: use identity (k in fractional, consistent with
  ! build_multi_HK / build_hk.f90 convention)
  REAL(dp), DIMENSION(3,3) :: rec_latt

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: kpts_bz
  INTEGER :: n_k_bz
  INTEGER, DIMENSION(3) :: nk_bz

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: kpts_path
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: k_dist
  INTEGER, DIMENSION(3) :: label_pos
  INTEGER :: n_k_path

  CHARACTER(LEN=500) :: output_filename
  INTEGER :: n_args, i, j

  !=============================================================================
  ! Read config filename from command line
  !=============================================================================

  CALL set_machine_acc()

  n_args = COMMAND_ARGUMENT_COUNT()
  IF ( n_args < 1 ) THEN
    WRITE(*,*) 'Usage: ./test_CPA <config_file>'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, config_file)

  WRITE(*,'(a)') ' ============================================='
  WRITE(*,'(a)') ' test_CPA: General N-component diagonal CPA'
  WRITE(*,'(a)') '           for zinc-blende multi-component alloy'
  WRITE(*,'(a)') '           (VCA-before-CPA scheme)'
  WRITE(*,'(a)') ' ============================================='
  WRITE(*,'(a,a)') ' Config file: ', TRIM(config_file)

  !=============================================================================
  ! Step 1: Read and validate configuration
  !=============================================================================

  CALL read_and_validate_config( TRIM(config_file), cfg )

  !=============================================================================
  ! Step 2: Build CPA parameters (Vegard lattice, real elements via VCA,
  !         hopping VCA), see cpa_multi_build.f90
  !=============================================================================

  CALL build_multi_cpa_params( cfg, params )

  !=============================================================================
  ! Print summary of the real-element structure (VCA-before-CPA)
  !=============================================================================

  WRITE(*,'(a)') ''
  WRITE(*,'(a,i4,a,i4)') ' Real elements: n_real_elem_sl1=', params%n_real_elem_sl1, &
    '  n_real_elem_sl2=', params%n_real_elem_sl2
  WRITE(*,'(a,l2)') ' odd_hopping (BEB off-diagonal disorder for hopping) = ', &
    params%odd_hopping
  WRITE(*,'(a)') ' Sublattice 1 real elements (VCA-averaged onsite/SOC):'
  DO i = 1, params%n_real_elem_sl1
    WRITE(*,'(a,a2,a,f8.4)') &
      '   ', params%elem_sl1(i)%name, '  conc=', params%elem_sl1(i)%conc
  END DO
  WRITE(*,'(a)') ' Sublattice 2 real elements (VCA-averaged onsite/SOC):'
  DO j = 1, params%n_real_elem_sl2
    WRITE(*,'(a,a2,a,f8.4)') &
      '   ', params%elem_sl2(j)%name, '  conc=', params%elem_sl2(j)%conc
  END DO

  IF ( params%n_real_elem_sl1 == 1 ) THEN
    WRITE(*,'(a)') ' NOTE: sublattice 1 has no disorder (single real element) -- &
                    &sigma1 is fixed at VCA onsite, no Soven update on sl1.'
  END IF
  IF ( params%n_real_elem_sl2 == 1 ) THEN
    WRITE(*,'(a)') ' NOTE: sublattice 2 has no disorder (single real element) -- &
                    &sigma2 is fixed at VCA onsite, no Soven update on sl2.'
  END IF
  IF ( params%n_real_elem_sl1 == 1 .AND. params%n_real_elem_sl2 == 1 ) THEN
    WRITE(*,'(a)') ' NOTE: pure binary case -- no CPA disorder at all, Soven &
                    &loop is skipped entirely.'
  END IF

  !=============================================================================
  ! Step 3: Generate BZ mesh and k-path
  !
  ! NOTE: rec_latt = Identity  (k-vectors in fractional units, consistent
  ! with build_multi_HK / sparse_ham convention; see build_hk.f90 comments)
  !=============================================================================

  rec_latt = 0.0_dp
  rec_latt(1,1) = 1.0_dp
  rec_latt(2,2) = 1.0_dp
  rec_latt(3,3) = 1.0_dp

  nk_bz = (/ cfg%nk_bz_x, cfg%nk_bz_y, cfg%nk_bz_z /)
  CALL generate_mp_mesh( rec_latt, nk_bz, kpts_bz, n_k_bz )

  CALL generate_kpath_LGX( rec_latt, cfg%n_per_segment, &
                            kpts_path, k_dist, n_k_path, label_pos )

  WRITE(*,'(a,3(i4,1x))') ' BZ mesh points (Nx,Ny,Nz): ', nk_bz
  WRITE(*,'(a,i6)')  ' Total BZ mesh k-points:   ', n_k_bz
  WRITE(*,'(a,i6)')  ' k-path total points:       ', n_k_path
  WRITE(*,'(a,3(i6,1x))') ' L, Gamma, X positions:     ', label_pos
  WRITE(*,'(a,3(f8.5,1x))') ' k_dist at L, G, X:        ', &
    k_dist(label_pos(1)), k_dist(label_pos(2)), k_dist(label_pos(3))

  !=============================================================================
  ! Step 4: Determine output filename
  !=============================================================================

  IF ( LEN_TRIM(cfg%output_filename) == 0 ) THEN
    output_filename = build_output_filename( cfg )
    WRITE(*,'(a,a)') ' Auto-generated output filename: ', TRIM(output_filename)
  ELSE
    output_filename = TRIM(cfg%output_filename)
    WRITE(*,'(a,a)') ' Output filename (from config): ', TRIM(output_filename)
  END IF

  !=============================================================================
  ! Step 5: Run CPA Soven self-consistency + spectral function
  !=============================================================================

  CALL cpa_multi_compute_spectral_function( &
    cfg%E_min, cfg%E_max, cfg%n_E, cfg%eta, &
    params, &
    kpts_bz, n_k_bz, &
    kpts_path, k_dist, n_k_path, &
    cfg%tol, cfg%max_iter, cfg%mixing_alpha, &
    TRIM(output_filename) )

  !=============================================================================
  ! Cleanup
  !=============================================================================

  CALL destroy_cpa_multi_params( params )
  IF (ALLOCATED(kpts_bz))   DEALLOCATE(kpts_bz)
  IF (ALLOCATED(kpts_path)) DEALLOCATE(kpts_path)
  IF (ALLOCATED(k_dist))    DEALLOCATE(k_dist)
  IF (ALLOCATED(cfg%sl1))   DEALLOCATE(cfg%sl1)
  IF (ALLOCATED(cfg%sl2))   DEALLOCATE(cfg%sl2)
  IF (ALLOCATED(cfg%x1))    DEALLOCATE(cfg%x1)
  IF (ALLOCATED(cfg%x2))    DEALLOCATE(cfg%x2)

  WRITE(*,'(a)') ''
  WRITE(*,'(a)') ' test_CPA: DONE.'
  WRITE(*,'(a,a)') ' Results: ', TRIM(output_filename)

END PROGRAM test_CPA