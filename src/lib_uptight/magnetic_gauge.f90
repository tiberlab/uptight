MODULE magnetic_gauge
  USE precision, only: dp
  USE globals, only: GAUGE_NONE, GAUGE_LANDAU_Z, GAUGE_SYMMETRIC_Z, use_magnetic_field

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calculate_peierls_phase

  ! --- Physical Constants ---
  ! Using approximate values for now. Need to verify consistency with uptight's units.
  ! Elementary charge in Coulombs
  REAL(dp), PARAMETER :: E_CHARGE = 1.602176634E-19_dp
  ! Reduced Planck constant in J*s
  REAL(dp), PARAMETER :: HBAR = 1.054571817E-34_dp
  ! Factor e/hbar in SI units (C / (J*s) = 1 / V) -> Units of 1 / (T * m^2) = 1 / Weber
  REAL(dp), PARAMETER :: E_OVER_HBAR = E_CHARGE / HBAR

CONTAINS

  !---------------------------------------------------------------------
  FUNCTION get_vector_potential(r, B_vector, gauge_choice) RESULT(A_vector)
    ! Calculates the vector potential A at position r for a given B_vector and gauge choice.
    ! Currently only supports uniform B along z-axis for Landau and Symmetric gauges.
    ! Assumes r is in Angstroms and B_vector is in Tesla. Result A_vector is in T*Angstrom.
    REAL(dp), INTENT(IN) :: r(3)              ! Position vector (Angstroms)
    REAL(dp), INTENT(IN) :: B_vector(3)       ! Magnetic field vector (Tesla)
    INTEGER, INTENT(IN)  :: gauge_choice      ! Selected gauge (GAUGE_NONE, GAUGE_LANDAU_Z, etc.)
    REAL(dp)             :: A_vector(3)       ! Resulting vector potential (T*Angstrom)

    REAL(dp), PARAMETER :: TOL = 1.0E-9_dp ! Tolerance for checking if B is along z

    A_vector = 0.0_dp ! Default to zero

    SELECT CASE (gauge_choice)
      CASE (GAUGE_LANDAU_Z)
        ! Check if B is along z within tolerance
        IF (ABS(B_vector(1)) < TOL .AND. ABS(B_vector(2)) < TOL) THEN
          ! A = (-B_z * y, 0, 0)
          A_vector(1) = -B_vector(3) * r(2)
          A_vector(2) = 0.0_dp
          A_vector(3) = 0.0_dp
        ELSE
          ! Placeholder: Add warning or error for non-z field in Landau Z gauge
          WRITE(*,*) 'WARNING: Landau Z gauge currently only supports B along z-axis.'
        END IF

      CASE (GAUGE_SYMMETRIC_Z)
        ! Check if B is along z within tolerance
        IF (ABS(B_vector(1)) < TOL .AND. ABS(B_vector(2)) < TOL) THEN
          ! A = 0.5 * (-B_z * y, B_z * x, 0)
          A_vector(1) = -0.5_dp * B_vector(3) * r(2)
          A_vector(2) =  0.5_dp * B_vector(3) * r(1)
          A_vector(3) = 0.0_dp
        ELSE
          ! Placeholder: Add warning or error for non-z field in Symmetric Z gauge
          WRITE(*,*) 'WARNING: Symmetric Z gauge currently only supports B along z-axis.'
        END IF

      CASE (GAUGE_NONE)
        ! A_vector remains zero

      CASE DEFAULT
        WRITE(*,*) 'WARNING: Unknown gauge choice in get_vector_potential. Setting A=0.'
        ! A_vector remains zero

    END SELECT

  END FUNCTION get_vector_potential
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  FUNCTION calculate_peierls_phase(R_i, R_j, B_vector, gauge_choice) RESULT(phase_factor)
    ! Calculates the Peierls phase factor exp(i * e/hbar * integral(A.dl))
    ! using the midpoint approximation for the integral.
    ! Assumes R_i, R_j are in Angstroms, B_vector in Tesla.
    REAL(dp), INTENT(IN) :: R_i(3)            ! Position of site i (Angstroms)
    REAL(dp), INTENT(IN) :: R_j(3)            ! Position of site j (Angstroms)
    REAL(dp), INTENT(IN) :: B_vector(3)       ! Magnetic field vector (Tesla)
    INTEGER, INTENT(IN)  :: gauge_choice      ! Selected gauge
    COMPLEX(dp)          :: phase_factor      ! Resulting complex phase factor (dimensionless)

    REAL(dp) :: R_mid(3)
    REAL(dp) :: dR(3)
    REAL(dp) :: A_mid(3)        ! Units: T*Angstrom
    REAL(dp) :: integral_approx ! Units: T*Angstrom^2
    REAL(dp) :: phase_arg       ! Dimensionless
    REAL(dp), PARAMETER :: ANGSTROM_TO_METER_SQUARED = 1.0E-20_dp

    ! If magnetic field is not used globally, return phase factor of 1.
    IF (.not. use_magnetic_field) THEN
       phase_factor = CMPLX(1.0_dp, 0.0_dp, kind=dp)
       RETURN
    END IF

    ! Midpoint approximation for the line integral path
    R_mid = 0.5_dp * (R_i + R_j)
    dR = R_j - R_i

    ! Get vector potential at the midpoint (Units: T*Angstrom)
    A_mid = get_vector_potential(R_mid, B_vector, gauge_choice)

    ! Calculate the line integral approximation: A(midpoint) . (Rj - Ri)
    ! Units: (T*Angstrom) . (Angstrom) = T*Angstrom^2
    integral_approx = dot_product(A_mid, dR)

    ! Convert integral to SI units (T*m^2 = Weber) for use with SI constants e/hbar
    integral_approx = integral_approx * ANGSTROM_TO_METER_SQUARED

    ! Calculate the phase argument (dimensionless)
    ! (e/hbar [1/Weber]) * (integral [Weber])
    phase_arg = E_OVER_HBAR * integral_approx

    ! Calculate the complex phase factor exp(i * phase_arg)
    phase_factor = CMPLX(COS(phase_arg), SIN(phase_arg), kind=dp)

  END FUNCTION calculate_peierls_phase
  !---------------------------------------------------------------------

END MODULE magnetic_gauge