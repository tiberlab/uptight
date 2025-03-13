!=============================================================================
!
!  1. Module normalization_factors_sg1D
!
!=============================================================================

MODULE normalization_factors_sg1D

  !===========================================================================
  !
  ! Contains vectors used to calculate the normalization
  ! and backtransformation of k.p / Schroedinger states.
  !
  !===========================================================================

  REAL( KIND(0.0D0) ), DIMENSION( : ), POINTER :: diag_kpV, norm_kpV

  !===========================================================================

END MODULE normalization_factors_sg1D
