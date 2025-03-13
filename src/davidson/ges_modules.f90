!=============================================================================
!
! Contains any modules that are commonly used by 1D/2D/3D version.
!
!_____________________________________________________________________________
!
!  1. MODULE dimension
!
!  2. MODULE limits
!
!  3. MODULE poisson_control
!
!  4. MODULE math_constants
!
!  5. MODULE sparse_matrix
!
!  6. MODULE sparse_matrix_real
!
!  7. MODULE sparse_matrix_real_nonsym
!
!  8. MODULE derived_constants
!
!  9. MODULE lapack_mat
!
!=============================================================================

MODULE dimension

  !===========================================================================
  !
  ! dimensionC = '1D' , '2D' , '3D'
  !
  !===========================================================================
  
  CHARACTER( LEN = 2 ) :: dimensionC
  
  !===========================================================================

END MODULE dimension

!=============================================================================

MODULE limits

  !===========================================================================
  !
  ! temp_limit --> Temperature cannot get below that value.
  !   
  ! lastbandL = .TRUE. --> Values in E_c1DM for bands which are not defined
  !                        are set to value of the latest valid band edge.
  !             .FALSE.--> are set to E_cb / E_vb       (in e0 !)
  !
  ! cond_band_max   Maximum value for conduction bands (in e0 !)
  !  val_band_max   Minimum value for valence    bands (in e0 !)
  !
  !___________________________________________________________________________
  !
  ! min_dens ...  Minimum density that can appear in current integration.
  !               Valid for electrons and holes.
  !               Units  [particles/cm^3]
  !
  !___________________________________________________________________________
  !
  ! min_precision_sg              ...  Minimum precision for effective-mass sg
  ! min_precision_kp                     ...  Minimum precision for kp
  ! min_iter_itere_kp,min_iter_itsub_kp  ...  for kp and Jacek Majewski
  ! min_iter_itere_sg,min_iter_itsub_sg  ...  for sg and Jacek Majewski
  !
  !===========================================================================

  REAL ( KIND(0.0D0) ) :: temp_limit = 1.0D0
  REAL ( KIND(0.0D0) ) :: E_cb = 100.0D0, E_vb = -100.0D0
  REAL ( KIND(0.0D0) ) :: cond_band_max = 100.0D0, val_band_max = -100.0D0
  REAL ( KIND(0.0D0) ) :: min_dens = 1.0D-10
  
  ! set to lower value to allow adjustment of schroedinger-kp-residual
  ! via input file ($numeric-control)

  ! REAL ( KIND(0.0D0) ) :: min_precision_sg = 1.0D-10, &
  ! min_precision_kp = 1.0D-10

  REAL ( KIND(0.0D0) ) :: min_precision_sg = 1.0D-3 , min_precision_kp = 1.0D-3

  INTEGER              :: min_iter_itere_kp = 10, min_iter_itsub_kp = 20
  INTEGER              :: min_iter_itere_sg = 10, min_iter_itsub_sg = 20

  LOGICAL              :: lastbandL  = .TRUE.

  !===========================================================================
  !
  ! Comment G. Zandler:
  !
  ! Make sure that in routine ges_modules.f90 the limits
  ! for E_vb, E_cb, E...max are set to +-100 (instead of +-10).
  !
  ! Otherwise you may get problems for wide-gap materials.
  ! 
  ! ALSO: set min_dens to 1d-100
  ! 
  ! Comment S. Hackenbuchner:
  !         
  ! set min_dens to 1d-10 is fine
  !
  !===========================================================================

END MODULE limits

!=============================================================================

MODULE poisson_control

  !===========================================================================
  !
  ! built_inL = .TRUE.  --> The routine poisson_problem calculates the built-
  !                         in potential.
  !             .FALSE. --> The routine poisson_problem calculates
  !                         the potential;
  !                         the built_in potential is assumed to exist in 
  !                         built_in_potXDV in module potentialsXD.
  !
  ! classL    = .TRUE.  --> The routines poisson_problem and
  !                         solve_current_problem1D
  !                         calculate the potentials assuming
  !                         a pure classical treatment.
  !             .FALSE. --> The routines poisson_problem and
  !                         solve_current_problem1D
  !                         calculate the potentials as specified
  !                         in input parser.
  !
  !===========================================================================
  
  LOGICAL :: built_inL = .TRUE. , classL = .TRUE. 
  
  !===========================================================================

END MODULE poisson_control

!=============================================================================

MODULE math_constants

  !===========================================================================
  !
  ! Contains mathematical constants : pi
  !
  !===========================================================================

  !  REAL(KIND(0.0D0))           :: pi = 3.1415926535897932385d0 
  REAL ( KIND(0.0D0) ), PARAMETER :: pi = 3.1415926535897932385d0 
  
  !===========================================================================

END MODULE math_constants

!=============================================================================

MODULE sparse_matrix

  !===========================================================================
  !
  ! Matrix stored in row-indexed sparse storage mode.
  ! Numerical recipes p.71
  ! Note: complex version
  !
  !___________________________________________________________________________
  !
  ! max_num_off_diag  ...  Maximum number offdiagonal elements
  !
  !___________________________________________________________________________
  !
  ! sa, ija   (1..N+1+number_of_off_diagonal_elements)
  !
  !===========================================================================
  ! Additional infos 
  INTEGER :: nrow, ncol
  CHARACTER(1) :: sparse_fmt     ! L(ower),U(pper),F(ull)       
  !  
  COMPLEX ( KIND(0.0D0) ), DIMENSION( : ), POINTER :: sa
  INTEGER,                 DIMENSION( : ), POINTER :: ija

  !===========================================================================

END MODULE sparse_matrix

!=============================================================================

MODULE sparse_matrix_real

  !===========================================================================
  !
  ! Matrix stored in row-indexed sparse storage mode.
  ! Numerical recipes p.71
  ! Note: real version
  !
  !___________________________________________________________________________
  !
  ! max_num_off_diag  ...  Maximum number offdiagonal elements
  !
  !___________________________________________________________________________
  !
  ! sa, ija   (1..N+1+number_of_off_diagonal_elements)
  !
  !===========================================================================
 
  REAL     ( KIND(0.0D0) ), DIMENSION( : ), POINTER :: sa
  INTEGER,                  DIMENSION( : ), POINTER :: ija

  !===========================================================================

END MODULE sparse_matrix_real

!=============================================================================

MODULE sparse_matrix_real_nonsym

  !===========================================================================
  !
  ! Matrix stored in row-indexed sparse storage mode.
  ! Numerical recipes p.71
  ! Note: real version
  !
  !___________________________________________________________________________
  !
  ! max_num_off_diag  ...  Maximum number offdiagonal elements
  !
  !___________________________________________________________________________
  !
  ! sa, ija   (1..N+1+number_of_off_diagonal_elements)
  !
  !===========================================================================
 
  REAL     ( KIND(0.0D0) ), DIMENSION( : ), POINTER :: sa
  INTEGER,                  DIMENSION( : ), POINTER :: ija
  
  !===========================================================================

END MODULE sparse_matrix_real_nonsym

!=============================================================================

MODULE derived_constants

  !===========================================================================
  !
  ! Contains physical constants derived from values in database.
  !
  ! h2b2m_evAA2  = reduce_planck**2/(2*ELECTRON_MASS)
  !                /ABS(electron_charge)  * (1d10)**2  = 3.81 eV AA**2
  !
  ! h2b2m_Jm2 = 3.81*ABS(electron_charge) * (1d-10)**2 = 
  !                         reduce_planck**2/(2*ELECTRON_MASS)
  !
  !===========================================================================

  REAL ( KIND(0.0D0) ) :: h2b2m_evAA2, h2b2m_Jm2

  !===========================================================================

END MODULE derived_constants

!=============================================================================

MODULE lapack_mat

  !===========================================================================
  !
  ! kp_lapackM  -->   k.p matrix for LAPACK (dim_kp,dim_kp)
  !
  !===========================================================================

  COMPLEX ( KIND(0.0D0) ), DIMENSION( :, : ), POINTER :: kp_lapackM
  INTEGER                                             :: dim_kp
 
  !===========================================================================

END MODULE lapack_mat

!=============================================================================
