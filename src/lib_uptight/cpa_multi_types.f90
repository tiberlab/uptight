! cpa_multi_types.f90
!
! Data structures for general N-component diagonal CPA.
! Supports arbitrary M elements on sublattice 1, N elements on sublattice 2.
! Total pseudo-elements per sublattice = M*N (general-CPA treatment).
!
MODULE cpa_multi_types

  USE precision, ONLY : dp
  USE globals,   ONLY : n_ref_states, n_ref_couplings

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pseudo_element, cpa_multi_params
  PUBLIC :: destroy_cpa_multi_params

  !===========================================================================
  ! TYPE pseudo_element
  !
  ! Represents one "pseudo-element" on a sublattice in the general-CPA sense:
  ! the atom sl1(idx_self) as it appears in binary sl1(idx_self)+sl2(idx_other)
  ! evaluated at a_vegard.
  !
  ! Index mapping: for M elements on sl1 and N on sl2,
  !   flat index idx = (i-1)*N + j   where i in 1..M, j in 1..N
  !===========================================================================

  TYPE pseudo_element

    INTEGER  :: idx_self   ! index into sl1(:) or sl2(:) (own sublattice)
    INTEGER  :: idx_other  ! index into the OTHER sublattice (defines binary)
    REAL(dp) :: conc       ! = x1(i) * x2(j)  for pseudo1; same for pseudo2

    ! Onsite energies at Gamma, in reference order (size n_ref_states).
    ! For Jancu: = raw offset + energy  (summed over 4 bonds divided by 4)
    ! For Tan:   = same but with Tan onsite correction also accumulated
    REAL(dp), ALLOCATABLE :: onsite(:)   ! size n_ref_states

    ! Effective p-channel SOC parameter (final, after Tan so_corr if applicable)
    REAL(dp) :: so_p

  END TYPE pseudo_element


  !===========================================================================
  ! TYPE cpa_multi_params
  !
  ! All parameters needed for the multi-component CPA Hamiltonian build
  ! and Soven self-consistency.
  !===========================================================================

  TYPE cpa_multi_params

    INTEGER :: n_sl1         ! M: number of elements on sublattice 1
    INTEGER :: n_sl2         ! N: number of elements on sublattice 2
    INTEGER :: n_pseudo1     ! = n_sl1 * n_sl2  pseudo-elements on sl1
    INTEGER :: n_pseudo2     ! = n_sl1 * n_sl2  pseudo-elements on sl2

    ! Flat arrays of pseudo-elements (size n_pseudo1 and n_pseudo2)
    ! Index idx = (i-1)*n_sl2 + j  maps to binary sl1(i)+sl2(j)
    TYPE(pseudo_element), ALLOCATABLE :: pseudo1(:)
    TYPE(pseudo_element), ALLOCATABLE :: pseudo2(:)

    ! VCA-averaged hopping (size n_ref_couplings), Harrison/Tan-scaled to a_vegard.
    ! hopping_vca(k) = SUM_{i,j} conc(i,j) * t_ij(k)
    REAL(dp), ALLOCATABLE :: hopping_vca(:)

    ! VCA-averaged SOC per sublattice (for off-diagonal SOC block)
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

    IF ( ALLOCATED(params%pseudo1) ) THEN
      DO idx = 1, SIZE(params%pseudo1)
        IF ( ALLOCATED(params%pseudo1(idx)%onsite) ) &
          DEALLOCATE(params%pseudo1(idx)%onsite)
      END DO
      DEALLOCATE(params%pseudo1)
    END IF

    IF ( ALLOCATED(params%pseudo2) ) THEN
      DO idx = 1, SIZE(params%pseudo2)
        IF ( ALLOCATED(params%pseudo2(idx)%onsite) ) &
          DEALLOCATE(params%pseudo2(idx)%onsite)
      END DO
      DEALLOCATE(params%pseudo2)
    END IF

    IF ( ALLOCATED(params%hopping_vca) ) DEALLOCATE(params%hopping_vca)

  END SUBROUTINE destroy_cpa_multi_params

END MODULE cpa_multi_types
