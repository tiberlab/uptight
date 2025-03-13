! This file is part of uptight.
!
! uptight is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! uptight is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with uptight. If not, see <https://www.gnu.org/licenses/>.
!
!=============================================================================
!
!             Module "latt_fun" - (c) Jerome Gleize - 2002
!
!=============================================================================
!
! => lattice_output
! => D_basis
! => cosine 
!
!=============================================================================

MODULE latt_fun

  !===========================================================================

  USE precision, only : dp, emach
  USE type_defs, only : ion_basis
  USE errors
  USE vectors

  !===========================================================================
  
  IMPLICIT NONE
  PRIVATE
  !===========================================================================

  !PUBLIC :: lattice_output
  PUBLIC :: D_basis, cosine, unit_displacement

CONTAINS


  !===========================================================================
  !
  ! Subroutine "lattice_output" : screen output of lattice data.
  !
  !===========================================================================

!!$  SUBROUTINE lattice_output
!!$
!!$    !=========================================================================
!!$
!!$    ! local variables :
!!$
!!$    TYPE (pair_coupling), POINTER :: new_pair
!!$
!!$    CHARACTER ( LEN = 2 ), DIMENSION( 2 ) :: pair_name, pair_name_swap
!!$    INTEGER,               DIMENSION( 3 ) :: n_latt_vec
!!$
!!$    REAL ( dp ) :: n_dist
!!$
!!$    INTEGER :: i_a, i_b, i_b_test, i_neighbour, n_order
!!$
!!$    LOGICAL :: absent
!!$
!!$    !=========================================================================
!!$
!!$    ! title output :
!!$
!!$    IF ( ioutput_flag .GE. 1 ) THEN
!!$    
!!$       !______________________________________________________________________
!!$
!!$       IF ( ioutput_flag .GE. 2 ) THEN
!!$       
!!$          !-------------------------------------------------------------------
!!$   
!!$          write_format = ' '
!!$          write_format = '( "LATTICE STRUCTURE :", //, 3x, "basis atom pairs&
!!$               & with coordinates in primitive supercell vectors units /",&
!!$               & /, 3x, "equivalent nearest neighbours with pair distance&
!!$               & in Angstrom", / )'
!!$          
!!$          WRITE ( *, write_format )
!!$          
!!$          !-------------------------------------------------------------------
!!$          
!!$       ELSE
!!$
!!$          !-------------------------------------------------------------------
!!$          
!!$          write_format = ' '
!!$          write_format = '( "SHORT LATTICE STRUCTURE :", //, 3x, &
!!$               &"basis atom pairs / nearest neighbours order / &
!!$               &pair distance in Angstrom", / )'
!!$          
!!$          WRITE ( *, write_format )
!!$          
!!$          !-------------------------------------------------------------------
!!$
!!$       END IF
!!$
!!$       !______________________________________________________________________
!!$       
!!$    END IF
!!$    
!!$    !=========================================================================
!!$
!!$    ! basis output :
!!$    
!!$    IF ( ioutput_flag .GE. 1 ) THEN
!!$       
!!$       current_near => start_near
!!$       
!!$       DO i_a = 1, n_basis
!!$
!!$          i_b = 0
!!$          
!!$          DO i_neighbour = 1, SIZE( current_near%ind )
!!$             
!!$             !________________________________________________________________
!!$             
!!$             i_b_test = i_b
!!$             i_b = current_near%ind( i_neighbour )
!!$             
!!$             IF ( i_b .GT. i_a ) CYCLE
!!$             
!!$             n_order = current_near%order( i_neighbour )
!!$             n_dist = current_near%dist( i_neighbour )
!!$             n_latt_vec = current_near%vec( :, i_neighbour )
!!$             
!!$             !________________________________________________________________
!!$                          
!!$             ! atom names and coordinates :
!!$             
!!$             IF ( ioutput_flag .GE. 2 ) THEN
!!$                                
!!$                !-------------------------------------------------------------
!!$
!!$                IF ( i_b .NE. i_b_test ) THEN
!!$                   
!!$                   !..........................................................
!!$                   
!!$                   write_format = ' '
!!$                   write_format = &
!!$                        '( /, x, ' // TRIM( int_format( i_a ) ) // &
!!$                        ', x, a2, " ( ", 3( f8.4, x ), ") <-> ", ' // &
!!$                        TRIM( int_format( i_b ) ) // &
!!$                        ', x, a2, " ( ", 3( f8.4, x ), ")" )'
!!$                   
!!$                   WRITE ( *, write_format ) i_a, basis%name( i_a ), &
!!$                         basis%coord( :, i_a ), i_b, basis%name( i_b ), &
!!$                         basis%coord( :, i_b )
!!$
!!$                   !..........................................................
!!$                   
!!$                END IF
!!$
!!$                !-------------------------------------------------------------
!!$                
!!$                ! equivalent pairs :
!!$                
!!$                write_format = ' ' 
!!$                write_format = '( 8x, "(", 3i3, " )", 1x, f12.4 )'
!!$                
!!$                WRITE ( *, write_format ) n_latt_vec, n_dist
!!$             
!!$                !-------------------------------------------------------------
!!$                
!!$             ELSE
!!$
!!$                !-------------------------------------------------------------
!!$                
!!$                write_format = ' '
!!$                write_format = '( /, 1x, ' // TRIM( int_format( i_a ) ) // &
!!$                     ', 1x, a5, 3x, ' // TRIM( int_format( i_b ) ) // &
!!$                     ', 1x, a5, 3x, ' // TRIM( int_format( n_order ) ) // &
!!$                     ', 1x, f10.5 )'
!!$                
!!$                WRITE ( *, TRIM( write_format ) ) i_a, basis%name( i_a ), &
!!$                     i_b, basis%name( i_b ), n_order, n_dist
!!$                
!!$                !-------------------------------------------------------------
!!$
!!$             END IF
!!$
!!$             !________________________________________________________________
!!$                          
!!$          END DO
!!$          
!!$          current_near => current_near%next
!!$             
!!$       END DO
!!$       
!!$       WRITE ( *, '( /, "***" )' )
!!$       
!!$    END IF
!!$    
!!$    !=========================================================================
!!$    
!!$  END SUBROUTINE lattice_output


  !===========================================================================
  !
  ! Function "D_basis" : calculates the vector :
  !
  !  D_basis = basis( beta ) - basis( alpha )
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => i_alpha, i_beta - integers : indexes of basis atoms.
  !
  ! OUTPUT :
  !
  ! => D_basis( 3 ) - REAL array : vector joining alpha and beta.
  !
  !===========================================================================

  FUNCTION D_basis( basis, i_beta, i_alpha )

    !=========================================================================

    ! input arguments :
    TYPE (ion_basis) :: basis
    INTEGER, INTENT( IN ) :: i_alpha, i_beta

    !_________________________________________________________________________

    ! output result :

    REAL ( dp ), DIMENSION( 3 ) :: D_basis

    !=========================================================================

    D_basis = basis%coord( :,i_beta ) - basis%coord( :, i_alpha )
    D_basis = MATMUL( basis%prim, D_basis )

    !=========================================================================

  END FUNCTION D_basis


  !===========================================================================
  !
  ! Function "cosine" : computes the direction cosine of vector D :
  !
  !  cosine( i ) = cos( D, i-axis )
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => D( 3 ) - REAL array : vector D.
  !
  !===========================================================================

  FUNCTION cosine( D )

    !=========================================================================

    ! input arguments :

    REAL ( dp ), DIMENSION( 3 ), INTENT( IN ) :: D

    !_________________________________________________________________________

    ! output result :

    REAL ( dp ), DIMENSION( 3 ) :: cosine

    !=========================================================================

    cosine = 0.0d0

    IF ( norm( D ) .LE. emach ) RETURN

    cosine = D / norm( D )

    !=========================================================================

  END FUNCTION cosine


  FUNCTION unit_displacement( basis, D_latt, i_origin, i_tip )

    !=========================================================================
    
    ! input
    TYPE (ion_basis) :: basis
    INTEGER, INTENT( IN ) :: i_origin, i_tip
    REAL(dp), DIMENSION(3), INTENT(IN) :: D_latt

    ! output
    REAL(dp), DIMENSION(3) :: unit_displacement

    ! local
    REAL(dp), DIMENSION(3) :: D_tmp

    D_tmp = MATMUL( basis%prim, D_latt) + D_basis( basis, i_tip, i_origin )
    unit_displacement = cosine( D_tmp )

    !=========================================================================

  END FUNCTION unit_displacement

  !===========================================================================

END MODULE latt_fun

!=============================================================================
















