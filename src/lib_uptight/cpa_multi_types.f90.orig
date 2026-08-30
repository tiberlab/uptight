! cpa_multi_types.f90
!
! Data structures for VCA-before-CPA general N-component diagonal CPA.
!
! REVISED (see implementation_plan.md, "Critical Revision"):
!
!   The previous "pseudo-element" scheme treated e.g. Ga in GaAs and Ga in
!   GaSb as two independent disorder species and applied CPA directly on
!   the M*N binaries. This is WRONG: Ga-As and Ga-Sb bonds are perfectly
!   correlated (Ga always sits with either As or Sb, never randomly), so
!   the pseudo-elements are not an uncorrelated-disorder ensemble as CPA
!   requires.
!
!   VCA-before-CPA fixes this:
!     1) Onsite and SOC are first VCA-averaged over binaries BY REAL
!        ELEMENT (e.g. all Ga-containing binaries are merged into one
!        "Ga" real_element).
!     2) CPA (Soven) is then only run over the REAL disorder that remains
!        in each sublattice, i.e. over the distinct real_element(:) list,
!        not over the M*N binaries.
!     3) Hopping is VCA-averaged directly over all M*N binaries (never
!        grouped by element) because a pair (bond) has a unique identity
!        and cannot be conflated across binaries.
!     4) SOC undergoes a SECOND VCA pass: after real_element%so_p is
!        obtained (VCA by element), it is VCA-averaged AGAIN across all
!        real elements in the sublattice to give a single effective
!        so_p_sl*_vca used in the off-diagonal SOC block. SOC never enters
!        the Soven loop.
!
MODULE cpa_multi_types

  USE precision, ONLY : dp
  USE globals,   ONLY : n_ref_states, n_ref_couplings

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: real_element, cpa_multi_params
  PUBLIC :: destroy_cpa_multi_params

  !===========================================================================
  ! TYPE real_element
  !
  ! One physical chemical element on a sublattice (e.g. "Ga", "As"), after
  ! VCA-averaging its onsite/SOC parameters over every binary that contains
  ! it. This is the object that participates in the Soven CPA loop (for
  ! onsite only).
  !===========================================================================

  TYPE real_element

    CHARACTER(LEN=2) :: name    ! physical element name, e.g. 'Ga', 'As'
    REAL(dp) :: conc            ! total concentration of this element in its sublattice

    ! Onsite energies at Gamma, VCA-averaged (by element) over all binaries
    ! containing this element, in reference order (size n_ref_states).
    REAL(dp), ALLOCATABLE :: onsite(:)

    ! SOC p-channel parameter, VCA-averaged (by element, "1st VCA pass")
    ! over all binaries containing this element.
    REAL(dp) :: so_p

  END TYPE real_element


  !===========================================================================
  ! TYPE cpa_multi_params
  !
  ! All parameters needed for the multi-component CPA Hamiltonian build
  ! and Soven self-consistency, under the VCA-before-CPA scheme.
  !===========================================================================

  TYPE cpa_multi_params

    INTEGER :: n_sl1         ! M: number of elements defining sl1 binaries
    INTEGER :: n_sl2         ! N: number of elements defining sl2 binaries

    INTEGER :: n_real_elem_sl1   ! number of DISTINCT real elements on sl1 (<= M)
    INTEGER :: n_real_elem_sl2   ! number of DISTINCT real elements on sl2 (<= N)

    TYPE(real_element), ALLOCATABLE :: elem_sl1(:)   ! size n_real_elem_sl1
    TYPE(real_element), ALLOCATABLE :: elem_sl2(:)   ! size n_real_elem_sl2

    ! Hopping: single VCA average over all M*N binaries (size n_ref_couplings).
    ! NEVER grouped by real element -- a bond/pair has a unique identity.
    REAL(dp), ALLOCATABLE :: hopping_vca(:)

    ! SOC: 2nd VCA pass -- average of elem_sl*(:)%so_p weighted by
    ! elem_sl*(:)%conc, giving one effective SOC parameter per sublattice.
    ! Used only in the off-diagonal SOC block; never enters the Soven loop.
    REAL(dp) :: so_p_sl1_vca
    REAL(dp) :: so_p_sl2_vca

    ! Vegard's law lattice constant and nearest-neighbour distance
    REAL(dp) :: a_vegard
    REAL(dp) :: dist_vegard     ! = a_vegard * sqrt(3) / 4

    ! ETB scheme, validated identical across all binaries: 'jancu' or 'tan'
    CHARACTER(LEN=10) :: scheme

  END TYPE cpa_multi_params


CONTAINS

  !===========================================================================
  ! Subroutine destroy_cpa_multi_params
  !
  ! Deallocate all allocatable components.
  !===========================================================================

  SUBROUTINE destroy_cpa_multi_params( params )

    TYPE(cpa_multi_params), INTENT(INOUT) :: params

    INTEGER :: idx

    IF ( ALLOCATED(params%elem_sl1) ) THEN
      DO idx = 1, SIZE(params%elem_sl1)
        IF ( ALLOCATED(params%elem_sl1(idx)%onsite) ) &
          DEALLOCATE(params%elem_sl1(idx)%onsite)
      END DO
      DEALLOCATE(params%elem_sl1)
    END IF

    IF ( ALLOCATED(params%elem_sl2) ) THEN
      DO idx = 1, SIZE(params%elem_sl2)
        IF ( ALLOCATED(params%elem_sl2(idx)%onsite) ) &
          DEALLOCATE(params%elem_sl2(idx)%onsite)
      END DO
      DEALLOCATE(params%elem_sl2)
    END IF

    IF ( ALLOCATED(params%hopping_vca) ) DEALLOCATE(params%hopping_vca)

  END SUBROUTINE destroy_cpa_multi_params

END MODULE cpa_multi_types