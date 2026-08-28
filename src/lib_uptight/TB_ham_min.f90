! TB_ham_min.f90
!
! Ban rut gon, doc lap cua koster_slater tu TB_ham.f90 (uptight goc).
!
! NOI DUNG DUOC COPY Y NGUYEN cong thuc Slater-Koster tu TB_ham.f90::koster_slater
! (khong sua doi bat ky he so nao) - chi doi ten module va bo cac dependency
! khong can thiet (mpi_globals, sparse_matrix, neighbours, v.v.) vi CPA
! khong dung cau truc sparse/supercell cua uptight.
!
! Dung cho build_hk.f90 de tinh hopping block giua 2 atom (cation-anion)
! trong primitive cell zinc-blende.
!
MODULE TB_ham_min

  USE precision, ONLY : dp
  USE globals,   ONLY : n_ref_states, n_ref_couplings, &
                         s, px, py, pz, se, dxy, dyz, dzx, dx2y2, dz2r2, &
                         sss, sps, pss, pps, ppp, seses, sess, sses, seps, &
                         pses, sds, dss, pds, dps, pdp, dpp, seds, dses, &
                         dds, ddp, ddd

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: koster_slater_min

CONTAINS

  !===========================================================================
  !
  ! Subroutine "koster_slater_min" :
  !
  ! COPY Y NGUYEN cong thuc tu TB_ham.f90::koster_slater (khong sua doi).
  !
  ! Tinh cac Koster-Slater two-center integrals cho atomic orbitals.
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => cos_latt( 3 )  - direction cosines cua bond vector.
  ! => t( n_ref_couplings ) - coupling matrix elements doc theo bond axis.
  !
  ! OUTPUT :
  !
  ! => v( n_ref_states, n_ref_states ) - matrix elements theo cartesian axes.
  !
  ! NOTE : thu tu orbital giong state.dat cua uptight :
  !
  !    1  2   3   4   5    6    7    8     9     10
  !    s  px  py  pz  s*  dxy  dyz  dzx  dx2y2  dz2r2
  !
  !===========================================================================

  SUBROUTINE koster_slater_min( cos_latt, t, v )

    !=========================================================================

    REAL ( dp ), DIMENSION( : ), INTENT( IN ) :: t
    REAL ( dp ), DIMENSION( 3 ), INTENT( IN ) :: cos_latt

    REAL ( dp ), DIMENSION(:,:) :: v

    !=========================================================================
    ! Local variables
    !=========================================================================

    REAL ( dp ) :: al, am, an, al2, am2, an2, dsq3

    !=========================================================================

    al = cos_latt( 1 )
    am = cos_latt( 2 )
    an = cos_latt( 3 )

    al2 = al**2
    am2 = am**2
    an2 = an**2

    dsq3 = SQRT( 3.0d0 )

    v = 0.0_dp

    !=========================================================================
    ! s-s
    v( s, s )   = t( sss )

    ! s-p, p-s
    v( s, px )  =   al * t( sps )
    v( px, s )  = - al * t( pss )

    v( s, py )  =   am * t( sps )
    v( py, s )  = - am * t( pss )

    v( s, pz )  =   an * t( sps )
    v( pz, s )  = - an * t( pss )

    ! p-p
    v( px, px ) = al2 * t( pps ) + ( 1.0d0 - al2 ) * t( ppp )
    v( py, py ) = am2 * t( pps ) + ( 1.0d0 - am2 ) * t( ppp )
    v( pz, pz ) = an2 * t( pps ) + ( 1.0d0 - an2 ) * t( ppp )

    v( px, py ) = al * am * ( t( pps ) - t( ppp ) )
    v( py, px ) = al * am * ( t( pps ) - t( ppp ) )

    v( px, pz ) = al * an * ( t( pps ) - t( ppp ) )
    v( pz, px ) = al * an * ( t( pps ) - t( ppp ) )

    v( py, pz ) = am * an * ( t( pps ) - t( ppp ) )
    v( pz, py ) = am * an * ( t( pps ) - t( ppp ) )

    ! se-s, s-se, se-se
    v( se, se ) = t( seses )
    v( s, se )  = t( sses )
    v( se, s )  = t( sess )

    ! se-p, p-se
    v( se, px ) =   al * t( seps )
    v( px, se ) = - al * t( pses )

    v( se, py ) =   am * t( seps )
    v( py, se ) = - am * t( pses )

    v( se, pz ) =   an * t( seps )
    v( pz, se ) = - an * t( pses )

    !_________________________________________________________________________
    ! s-d, d-s
    v( s, dxy )   = dsq3 * al * am * t( sds )
    v( dxy, s )   = dsq3 * al * am * t( dss )

    v( s, dyz )   = dsq3 * am * an * t( sds )
    v( dyz, s )   = dsq3 * am * an * t( dss )

    v( s, dzx )   = dsq3 * an * al * t( sds )
    v( dzx, s )   = dsq3 * an * al * t( dss )

    v( s, dx2y2 ) = dsq3 * 0.5d0 * ( al2 - am2 ) * t( sds )
    v( dx2y2, s ) = dsq3 * 0.5d0 * ( al2 - am2 ) * t( dss )

    v( s, dz2r2 ) = ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( sds )
    v( dz2r2, s ) = ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dss )

    !-------------------------------------------------------------------------
    ! se-d, d-se
    v( se, dxy ) = dsq3 * al * am * t( seds )
    v( dxy, se ) = dsq3 * al * am * t( dses )

    v( se, dyz ) = dsq3 * am * an * t( seds )
    v( dyz, se ) = dsq3 * am * an * t( dses )

    v( se, dzx ) = dsq3 * an * al * t( seds )
    v( dzx, se ) = dsq3 * an * al * t( dses )

    v( se, dx2y2 ) = dsq3 * 0.5d0 * ( al2 - am2 ) * t( seds )
    v( dx2y2, se ) = dsq3 * 0.5d0 * ( al2 - am2 ) * t( dses )

    v( se, dz2r2 ) = ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( seds )
    v( dz2r2, se ) = ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dses )

    !-------------------------------------------------------------------------
    ! p-d, d-p
    v( px, dxy ) =   am * ( dsq3 * al2 * t( pds ) &
         + ( 1.0d0 - 2.0d0 * al2 ) * t( pdp ) )
    v( dxy, px ) = - am * ( dsq3 * al2 * t( dps ) &
         + ( 1.0d0 - 2.0d0 * al2 ) * t( dpp ) )

    v( px, dyz ) =   al * am * an * ( dsq3 * t( pds ) - 2.0d0 * t( pdp ) )
    v( dyz, px ) = - al * am * an * ( dsq3 * t( dps ) - 2.0d0 * t( dpp ) )

    v( px, dzx ) =   an * ( dsq3 * al2 * t( pds ) &
         + ( 1.0d0 - 2.0d0 * al2 ) * t( pdp ) )
    v( dzx, px ) = - an * ( dsq3 * al2 * t( dps ) &
         + ( 1.0d0 - 2.0d0 * al2 ) * t( dpp ) )

    v( px, dx2y2 ) =   al * ( 0.5d0 * dsq3 * ( al2 - am2 ) * t( pds ) &
         + ( 1.0d0 - al2 + am2 ) * t( pdp ) )
    v( dx2y2, px ) = - al * ( 0.5d0 * dsq3 * ( al2 - am2 ) * t( dps ) &
         + ( 1.0d0 - al2 + am2 ) * t( dpp ) )

    v( px, dz2r2 ) =   al * ( ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( pds ) &
         - dsq3 * an2 * t( pdp ) )
    v( dz2r2, px ) = - al * ( ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dps ) &
         - dsq3 * an2 * t( dpp ) )

    !-------------------------------------------------------------------------

    v( py, dxy ) =   al * ( dsq3 * am2 * t( pds ) &
         + ( 1.0d0 - 2.0d0 * am2 ) * t( pdp ) )
    v( dxy, py ) = - al * ( dsq3 * am2 * t( dps ) &
         + ( 1.0d0 - 2.0d0 * am2 ) * t( dpp ) )

    v( py, dyz ) =   an * ( dsq3 * am2 * t( pds ) &
         + ( 1.0d0 - 2.0d0 * am2 ) * t( pdp ) )
    v( dyz, py ) = - an * ( dsq3 * am2 * t( dps ) &
         + ( 1.0d0 - 2.0d0 * am2 ) * t( dpp ) )

    v( py, dzx ) =   al * am * an * ( dsq3 * t( pds ) - 2.0d0 * t( pdp ) )
    v( dzx, py ) = - al * am * an * ( dsq3 * t( dps ) - 2.0d0 * t( dpp ) )

    v( py, dx2y2 ) =   am * ( 0.5d0 * dsq3 * ( al2 - am2 ) * t( pds ) &
         - ( 1.0d0 - am2 + al2 ) * t( pdp ) )
    v( dx2y2, py ) = - am * ( 0.5d0 * dsq3 * ( al2 - am2 ) * t( dps ) &
         - ( 1.0d0 - am2 + al2 ) * t( dpp ) )

    v( py, dz2r2 ) =   am * ( ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( pds ) &
         - dsq3 * an2 * t( pdp ) )
    v( dz2r2, py ) = - am * ( ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dps ) &
         - dsq3 * an2 * t( dpp ) )

    !-------------------------------------------------------------------------

    v( pz, dxy ) =   al * am * an * ( dsq3 * t( pds ) - 2.0d0 * t( pdp ) )
    v( dxy, pz ) = - al * am * an * ( dsq3 * t( dps ) - 2.0d0 * t( dpp ) )

    v( pz, dyz ) =   am * ( dsq3 * an2 * t( pds ) &
         + ( 1.0d0 - 2.0d0 * an2 ) * t( pdp ) )
    v( dyz, pz ) = - am * ( dsq3 * an2 * t( dps ) &
         + ( 1.0d0 - 2.0d0 * an2 ) * t( dpp ) )

    v( pz, dzx ) =   al * ( dsq3 * an2 * t( pds ) &
         + ( 1.0d0 - 2.0d0 * an2 ) * t( pdp ) )
    v( dzx, pz ) = - al * ( dsq3 * an2 *t( dps ) &
         + ( 1.0d0 - 2.0d0 * an2 ) * t( dpp ) )

    v( pz, dx2y2 ) =   an * ( 0.5d0 * dsq3 * ( al2 - am2 ) * t( pds ) &
         - ( al2 - am2 ) * t( pdp ) )
    v( dx2y2, pz ) = - an * ( 0.5d0 * dsq3 * ( al2 - am2 ) * t( dps ) &
         - ( al2 - am2 ) * t( dpp ) )

    v( pz, dz2r2 ) =   an * ( ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( pds ) &
         + dsq3 * ( al2 + am2 ) * t( pdp ) )
    v( dz2r2, pz ) = - an * ( ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dps ) &
         + dsq3 * ( al2 + am2 ) * t( dpp ) )

    !-------------------------------------------------------------------------
    ! d-d
    v( dxy, dxy ) = 3.0d0 * al2 * am2 * t( dds ) &
         + ( al2 + am2 - 4.d0 * al2 * am2 ) * t( ddp ) &
         + ( an2 + al2 * am2 ) * t( ddd )

    v( dyz, dyz ) = 3.0d0 * am2 * an2 * t( dds ) &
         + ( am2 + an2 - 4.d0 * am2 * an2 ) * t( ddp ) &
         + ( al2 + am2 * an2 ) * t( ddd )

    v( dzx, dzx ) = 3.0d0 * an2 * al2 * t( dds ) &
         + ( an2 + al2 - 4.d0 * an2 * al2 ) * t( ddp ) &
         + ( am2 + an2 * al2 ) * t( ddd )

    v( dx2y2, dx2y2 ) = 0.75d0 * ( al2 - am2 )**2 * t( dds )  &
         + ( al2 + am2 - ( al2 - am2 )**2 ) * t( ddp ) &
         + ( an2 + 0.25d0 * ( al2 - am2 )**2 ) * t( ddd )

    v( dz2r2, dz2r2 ) = ( an2 - 0.5d0 * ( al2 + am2) )**2 * t( dds ) &
         + 3.0d0 * an2 * ( al2 + am2 ) * t( ddp ) &
         + 0.75d0 * ( al2 + am2 )**2 * t( ddd )

    !-------------------------------------------------------------------------

    v( dxy, dyz ) = al * an * ( 3.0d0 * am2 * t( dds ) &
         + ( 1.0d0 - 4.0d0 * am2 ) * t( ddp ) + ( am2 - 1.0d0 ) * t( ddd ) )
    v( dyz, dxy ) = al * an * ( 3.0d0 * am2 * t( dds ) &
         + ( 1.0d0 - 4.0d0 * am2 ) * t( ddp ) + ( am2 - 1.0d0 ) * t( ddd ) )

    v( dxy, dzx ) = am * an * ( 3.0d0 * al2 * t( dds ) &
         + ( 1.0d0 - 4.0d0 * al2 ) * t( ddp ) + ( al2 - 1.0d0 ) * t( ddd ) )
    v( dzx, dxy ) = am * an * ( 3.0d0 * al2 * t( dds ) &
         + ( 1.0d0 - 4.0d0 * al2 ) * t( ddp ) + ( al2 - 1.0d0 ) * t( ddd ) )

    v( dxy, dx2y2 ) = al * am * ( al2 - am2 ) * ( 1.5d0 * t( dds ) &
         - 2.0d0 * t( ddp ) + 0.5d0 * t( ddd ) )
    v( dx2y2, dxy ) = al * am * ( al2 - am2 ) * ( 1.5d0 * t( dds ) &
         - 2.0d0 * t( ddp ) + 0.5d0 * t( ddd ) )

    v( dxy, dz2r2 ) = dsq3 * al * am * ( &
         ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         - 2.0d0 * an2 * t( ddp ) + 0.5d0 * ( 1.0d0 + an2 ) * t( ddd ) )
    v( dz2r2, dxy ) = dsq3 * al * am * ( &
         ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         - 2.0d0 * an2 * t( ddp ) + 0.5d0 * ( 1.0d0 + an2 ) * t( ddd ) )

    v( dyz, dzx ) = al * am * ( 3.0d0 * an2 * t( dds ) &
         + ( 1.0d0 - 4.0d0 * an2 ) * t( ddp ) + ( an2 - 1.0d0 ) * t( ddd ) )
    v( dzx, dyz ) = al * am * ( 3.0d0 * an2 * t( dds ) &
         + ( 1.0d0 - 4.0d0 * an2 ) * t( ddp ) + ( an2 - 1.0d0 ) * t( ddd ) )

    v( dyz, dx2y2 ) = am * an * ( 1.5d0 * ( al2 - am2 ) * t( dds ) &
         - ( 1.0d0 + 2.0d0 * ( al2 - am2 ) ) * t( ddp ) &
         + ( 1.0d0 + 0.5d0 * ( al2 - am2 ) ) * t( ddd ) )
    v( dx2y2, dyz ) = am * an * ( 1.5d0 * ( al2 - am2 ) * t( dds ) &
         - ( 1.0d0 + 2.0d0 * ( al2 - am2 ) ) * t( ddp ) &
         + ( 1.0d0 + 0.5d0 * ( al2 - am2 ) ) * t( ddd ) )

    v( dyz, dz2r2 ) = dsq3 * am * an * ( &
         ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         + ( al2 + am2 - an2 ) * t( ddp ) - 0.5d0 * ( al2 + am2 ) * t( ddd ) )
    v( dz2r2, dyz ) = dsq3 * am * an * ( &
         ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         + ( al2 + am2 - an2 ) * t( ddp ) - 0.5d0 * ( al2 + am2 ) * t( ddd ) )

    v( dzx, dx2y2 ) = al * an * ( 1.5d0 * ( al2 - am2 ) * t( dds )  &
         + ( 1.0d0 - 2.0d0 * ( al2 - am2 ) ) * t( ddp ) &
         - ( 1.0d0 - 0.5d0 * ( al2 - am2 ) ) * t( ddd ) )
    v( dx2y2, dzx ) = al * an * ( 1.5d0 * ( al2 - am2 ) * t( dds ) &
         + ( 1.0d0 - 2.0d0 * ( al2 - am2 ) ) * t( ddp ) &
         - ( 1.0d0 - 0.5d0 * ( al2 - am2 ) ) * t( ddd ) )

    v( dzx, dz2r2 ) = dsq3 * al * an * ( &
         ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         + ( al2 + am2 - an2 ) * t( ddp ) - 0.5d0 * ( al2 + am2 ) * t( ddd ) )
    v( dz2r2, dzx ) = dsq3 * al * an * ( &
         ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         + ( al2 + am2 - an2 ) * t( ddp ) - 0.5d0 * ( al2 + am2 ) * t( ddd ) )

    v( dx2y2, dz2r2 ) = dsq3 * ( al2 - am2 ) * &
         ( 0.5d0 * ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         - an2 * t( ddp ) + 0.25d0 * ( 1.0d0 + an2 ) * t( ddd ) )
    v( dz2r2, dx2y2 ) = dsq3 * ( al2 - am2 ) * &
         ( 0.5d0 * ( an2 - 0.5d0 * ( al2 + am2 ) ) * t( dds ) &
         - an2 * t( ddp ) + 0.25d0 * ( 1.0d0 + an2 ) * t( ddd ) )

    !=========================================================================

  END SUBROUTINE koster_slater_min

END MODULE TB_ham_min
