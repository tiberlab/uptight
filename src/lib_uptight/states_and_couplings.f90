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
!            Module "states_and_couplings" - (c) Jerome Gleize - 2002
!
!=============================================================================
!
! contains routines related to the supercell
!
!_____________________________________________________________________________
!
! => ref_states_and_couplings
! => states_sym
! => number_states_atom
!
!=============================================================================


MODULE states_and_couplings

  !===========================================================================

  USE globals
  USE precision
  USE type_defs
  USE checks, only : check_ion, check_ref_order, check_ref
  USE input_output
  USE errors

  !===========================================================================

  IMPLICIT NONE
  PRIVATE
  !===========================================================================

  PUBLIC :: ref_states_and_couplings, sort_couplings, number_states_atom
  PUBLIC :: init_n_st, set_max_order
  PUBLIC :: sort_states 

CONTAINS


  !===========================================================================
  !
  ! Subroutine "ref_states_and_couplings"
  !
  !===========================================================================

  SUBROUTINE ref_states_and_couplings( ref_states, n_ref_st, & 
                                       ref_couplings, n_ref_cpl)

    CHARACTER(CPLLEN),  DIMENSION(:), POINTER :: ref_couplings     
    CHARACTER(STATELEN),  DIMENSION(:), POINTER :: ref_states
    INTEGER, INTENT(OUT) :: n_ref_st, n_ref_cpl
    
    !=========================================================================
    INTEGER :: err, i_ion, i_par

    n_ref_st = 10
    !-------------------------------------------------------------------------
    ! Allocate the reference states array (global array)
    !-------------------------------------------------------------------------

    ALLOCATE( ref_states( n_ref_st ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
         'ref_states_and_couplings', 'ref_states' )

    ref_states = ' '
  
    !-------------------------------------------------------------------------
    ! Loop on reference states
    ! s = 1, px = 2, ...
    !-------------------------------------------------------------------------
    ref_states( s )     = 's'
    ref_states( px )    = 'px'
    ref_states( py )    = 'py'
    ref_states( pz )    = 'pz'
    ref_states( se )    = 's*'
    ref_states( dxy )   = 'dxy'
    ref_states( dyz )   = 'dyz'     
    ref_states( dzx )   = 'dzx'
    ref_states( dx2y2 ) = 'dx2y2'
    ref_states( dz2r2 ) = 'dz2r2'
 
    !-------------------------------------------------------------------------
    ! Get the number of reference couplings
    !-------------------------------------------------------------------------

    n_ref_cpl = 21
    !-------------------------------------------------------------------------
    ! Allocate the reference couplings array (global array)
    !-------------------------------------------------------------------------
    
    ALLOCATE( ref_couplings( n_ref_cpl ), STAT = err )
    IF ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
         'ref_states_and_couplings', 'ref_couplings' )

    ref_couplings = ' '
    
    !-------------------------------------------------------------------------
    ! Loop on the reference couplings
    !-------------------------------------------------------------------------
    ref_couplings( sss ) = 'sss'
    ref_couplings( sps ) = 'sps'
    ref_couplings( pss ) = 'pss'
    ref_couplings( pps ) = 'pps'
    ref_couplings( ppp ) = 'ppp'
    ref_couplings( seses ) = 's*s*s'
    ref_couplings( sess ) = 's*ss'
    ref_couplings( sses ) = 'ss*s'
    ref_couplings( seps ) = 's*ps'
    ref_couplings( pses ) = 'ps*s'
    ref_couplings( sds ) = 'sds'
    ref_couplings( dss ) = 'dss'
    ref_couplings( pds ) = 'pds'
    ref_couplings( dps ) = 'dps'
    ref_couplings( pdp ) = 'pdp'
    ref_couplings( dpp ) = 'dpp'
    ref_couplings( seds ) = 's*ds'
    ref_couplings( dses ) = 'ds*s' 
    ref_couplings( dds ) = 'dds'
    ref_couplings( ddp ) = 'ddp' 
    ref_couplings( ddd ) = 'ddd'
   

  END SUBROUTINE ref_states_and_couplings
  !===========================================================================
  !===========================================================================
  !
  ! Subroutine to initialize basis%n_st (contains the n. states per atom)
  ! NOTE: the field mat_data%ion should be already initialized (init_mat_ion)
  !
  !===========================================================================

  SUBROUTINE init_n_st(basis,mat_data)

    TYPE(ion_basis) :: basis
    TYPE(material_data), DIMENSION(:), POINTER :: mat_data
    
    INTEGER :: i


    DO i=1,basis%n_basis

       basis%n_st(i) = number_states_atom(basis,i,mat_data)
       !!write(*,*) "atom",i,"n. states =",basis%n_st(i)


    END DO


  END SUBROUTINE init_n_st


  SUBROUTINE sort_states(mat_data, ref_states, ref_couplings)
    TYPE(material_data) :: mat_data
    CHARACTER(STATELEN), DIMENSION(:) :: ref_states
    CHARACTER(CPLLEN), DIMENSION(:) :: ref_couplings

    ! LOCALS
    INTEGER :: i_mat, i_par, i_ion, i_pair, i_st, err!, i_mat_row, i_mat_column
    INTEGER :: n_states, n_ref_st, n_cpl, n_ref_cpl
    INTEGER, DIMENSION(:), POINTER :: ind_ref  
    REAL(dp), DIMENSION(:), POINTER :: energy, scaling, P, S, Q
    CHARACTER(CPLLEN), DIMENSION(:), POINTER :: strings
  
    n_ref_cpl = size(ref_couplings)
    n_ref_st = size(ref_states)
    
    !----------------------------------------------------------------------
    ! Check and order the states for all atoms in the material
    !---------------------------------------------------------------------- 

    do i_par = 1, sum(mat_data%nr_parents)  

        do i_ion = 1, mat_data%parent(i_par)%n_ion

            !...............................................................
            ! Place states in reference order
            !
            ! if tab(i) == ref(j) => ind_ref(i) = j 
            !...............................................................
            n_states = size(mat_data%parent(i_par)%ion(i_ion)%state)
            
            allocate(ind_ref(n_states), stat=err)
            allocate(energy(n_states), stat=err)
            allocate(strings(n_states), stat=err)
            if ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
                                 'sort_states', 'energy' )
            
            CALL check_ref_order( mat_data%parent(i_par)%ion(i_ion)%state, &
                ref_states, ind_ref )
           
            energy = mat_data%parent(i_par)%ion(i_ion)%energy
            strings = mat_data%parent(i_par)%ion(i_ion)%state


            do i_st=1, n_states
               mat_data%parent(i_par)%ion(i_ion)%state(ind_ref(i_st)) = &
                                 strings(i_st) !DO NOT UNDERSTAND WHY TO DO THIS
               mat_data%parent(i_par)%ion(i_ion)%energy(ind_ref(i_st)) = &
                                 energy(i_st) 

            end do
            
            CALL check_ref( mat_data%parent(i_par)%ion(i_ion)%state, &
                ref_states, ind_ref )

            mat_data%parent(i_par)%ion(i_ion)%ind_ref = ind_ref
            
            mat_data%parent(i_par)%ion(i_ion)%ind_ref_inv = 0
            
            do i_st=1, n_states
              mat_data%parent(i_par)%ion(i_ion)%ind_ref_inv(ind_ref(i_st)) = i_st
            end do
           
            deallocate(energy, ind_ref)
            deallocate(strings)

        end do 
            
        !...............................................................
        ! Place couplings in reference order
        ! 2019-02-15, M. Auf der Maur : they are now placed at the position of the
        !    reference couplings, because after the code works on the full coupling
        !    matrix.
        !...............................................................
        do i_pair = 1, mat_data%parent(i_par)%n_pair   
        
          n_cpl = size(mat_data%parent(i_par)%pair(i_pair)%cpl)
          allocate(ind_ref(n_cpl), stat=err)
          allocate(energy(n_ref_cpl), stat=err)
          allocate(scaling(n_ref_cpl), stat=err)
          allocate(P(n_ref_cpl), stat=err)
          allocate(S(n_ref_cpl), stat=err)
          allocate(Q(n_ref_cpl), stat=err)
          allocate(strings(n_ref_cpl), stat=err)
          energy = 0.0
          scaling = 0.0
          P = 0.0
          S = 0.0
          Q = 0.0
          strings = ''
          if ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
                                 'sort_states', 'ind_ref' )
          
          CALL check_ref(mat_data%parent(i_par)&
                        %pair(i_pair)%cpl, ref_couplings, ind_ref )
           
          do i_st=1, n_cpl 
            strings(ind_ref(i_st)) = mat_data%parent(i_par)%pair(i_pair)%cpl(i_st)
            energy(ind_ref(i_st)) = mat_data%parent(i_par)%pair(i_pair)%energy(i_st)
            scaling(ind_ref(i_st)) = mat_data%parent(i_par)%pair(i_pair)%scaling(i_st)
            P(ind_ref(i_st)) = mat_data%parent(i_par)%pair(i_pair)%fac_P(i_st)
            S(ind_ref(i_st)) = mat_data%parent(i_par)%pair(i_pair)%fac_S(i_st)
            Q(ind_ref(i_st)) = mat_data%parent(i_par)%pair(i_pair)%fac_Q(i_st)

            !mat_data(i_mat)%parent(i_par)%pair(i_pair)&
            !%cpl(ind_ref(i_st)) = strings(i_st) 
             
            !mat_data(i_mat)%parent(i_par)%pair(i_pair)&
            !%energy(ind_ref(i_st)) = energy(i_st) 
               
            !mat_data(i_mat)%parent(i_par)%pair(i_pair)&
            !%scaling(ind_ref(i_st)) = scaling(i_st) 
          end do
          
          ! the following is a bit quirky
          deallocate(mat_data%parent(i_par)%pair(i_pair)%cpl)
          deallocate(mat_data%parent(i_par)%pair(i_pair)%energy)
          deallocate(mat_data%parent(i_par)%pair(i_pair)%scaling)
          deallocate(mat_data%parent(i_par)%pair(i_pair)%fac_P)
          deallocate(mat_data%parent(i_par)%pair(i_pair)%fac_S)
          deallocate(mat_data%parent(i_par)%pair(i_pair)%fac_Q)
          allocate(mat_data%parent(i_par)%pair(i_pair)%cpl(n_ref_cpl))
          allocate(mat_data%parent(i_par)%pair(i_pair)%energy(n_ref_cpl))
          allocate(mat_data%parent(i_par)%pair(i_pair)%scaling(n_ref_cpl))
          allocate(mat_data%parent(i_par)%pair(i_pair)%fac_P(n_ref_cpl))
          allocate(mat_data%parent(i_par)%pair(i_pair)%fac_S(n_ref_cpl))
          allocate(mat_data%parent(i_par)%pair(i_pair)%fac_Q(n_ref_cpl))

          mat_data%parent(i_par)%pair(i_pair)%cpl = strings 
          mat_data%parent(i_par)%pair(i_pair)%energy = energy 
          mat_data%parent(i_par)%pair(i_pair)%scaling = scaling 
          mat_data%parent(i_par)%pair(i_pair)%fac_P = P
          mat_data%parent(i_par)%pair(i_pair)%fac_S = S
          mat_data%parent(i_par)%pair(i_pair)%fac_Q = Q
              
           
          deallocate(strings)
          deallocate(energy)
          deallocate(scaling)
          deallocate(P)
          deallocate(S)
          deallocate(Q)
          deallocate(ind_ref)
           
        end do 
    end do

  !------------------------------------------
  ! For interfaces between original materials
  !------------------------------------------
  ! do i_mat_row = 1, size(mat_data)
    
  !   do i_mat_column = 1, size(mat_data)
      
  !     do i_par = 1, sum(interface_data(i_mat_row, i_mat_column)%nr_parents)  

  !         do i_ion = 1, interface_data(i_mat_row, i_mat_column)%parent(i_par)%n_ion

  !             !...............................................................
  !             ! Place states in reference order
  !             !
  !             ! if tab(i) == ref(j) => ind_ref(i) = j 
  !             !...............................................................
  !             n_states = size(interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%state)
              
  !             allocate(ind_ref(n_states), stat=err)
  !             allocate(energy(n_states), stat=err)
  !             allocate(strings(n_states), stat=err)
  !             if ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
  !                                  'sort_states', 'energy' )
              
  !             CALL check_ref_order( interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%state, &
  !                 ref_states, ind_ref )
             
  !             energy = interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%energy
  !             strings = interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%state


  !             do i_st=1, n_states
  !                interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%state(ind_ref(i_st)) = &
  !                                  strings(i_st) !DO NOT UNDERSTAND WHY TO DO THIS
  !                interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%energy(ind_ref(i_st)) = &
  !                                  energy(i_st) 

  !             end do
              
  !             CALL check_ref( interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%state, &
  !                 ref_states, ind_ref )

  !             interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%ind_ref = ind_ref
              
  !             interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%ind_ref_inv = 0
              
  !             do i_st=1, n_states
  !               interface_data(i_mat_row, i_mat_column)%parent(i_par)%ion(i_ion)%ind_ref_inv(ind_ref(i_st)) = i_st
  !             end do
             
  !             deallocate(energy, ind_ref)
  !             deallocate(strings)

  !         end do 
              
  !         !...............................................................
  !         ! Place couplings in reference order
  !         ! 2019-02-15, M. Auf der Maur : they are now placed at the position of the
  !         !    reference couplings, because after the code works on the full coupling
  !         !    matrix.
  !         !...............................................................
  !         do i_pair = 1, interface_data(i_mat_row, i_mat_column)%parent(i_par)%n_pair   
          
  !           n_cpl = size(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%cpl)
  !           allocate(ind_ref(n_cpl), stat=err)
  !           allocate(energy(n_ref_cpl), stat=err)
  !           allocate(scaling(n_ref_cpl), stat=err)
  !           allocate(strings(n_ref_cpl), stat=err)
  !           energy = 0.0
  !           scaling = 0.0
  !           strings = ''
  !           if ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
  !                                  'sort_states', 'ind_ref' )
            
  !           CALL check_ref(interface_data(i_mat_row, i_mat_column)%parent(i_par)&
  !                         %pair(i_pair)%cpl, ref_couplings, ind_ref )
             
  !           do i_st=1, n_cpl 
  !             strings(ind_ref(i_st)) = interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%cpl(i_st)
  !             energy(ind_ref(i_st)) = interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%energy(i_st)
  !             scaling(ind_ref(i_st)) = interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%scaling(i_st)

  !             !mat_data(i_mat)%parent(i_par)%pair(i_pair)&
  !             !%cpl(ind_ref(i_st)) = strings(i_st) 
               
  !             !mat_data(i_mat)%parent(i_par)%pair(i_pair)&
  !             !%energy(ind_ref(i_st)) = energy(i_st) 
                 
  !             !mat_data(i_mat)%parent(i_par)%pair(i_pair)&
  !             !%scaling(ind_ref(i_st)) = scaling(i_st) 
  !           end do
            
  !           ! the following is a bit quirky
  !           deallocate(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%cpl)
  !           deallocate(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%energy)
  !           deallocate(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%scaling)
  !           allocate(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%cpl(n_ref_cpl))
  !           allocate(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%energy(n_ref_cpl))
  !           allocate(interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%scaling(n_ref_cpl))

  !           interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%cpl = strings 
                 
  !           interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%energy = energy 
                 
  !           interface_data(i_mat_row, i_mat_column)%parent(i_par)%pair(i_pair)%scaling = scaling 
                
             
  !           deallocate(strings)
  !           deallocate(energy)
  !           deallocate(scaling)
  !           deallocate(ind_ref)
             
  !         end do 
  !     end do

  !   end do
      
  ! end do

  end subroutine sort_states


  !===========================================================================
  !
  ! Function "number_states_atom" : calculates the number of atomic states
  ! included for a given atom in the basis.
  !
  !===========================================================================
  ! 
  ! INPUT : basis, mat_data, i_atom
  !
  ! OUTPUT : number of states
  !
  ! NOTE: Requires the fiels 'mat_data%ion' to be initialized after a call to
  !       "init_mat_ion" in alloys.f90
  !
  !===========================================================================

  FUNCTION number_states_atom( basis, i_atom, mat_data )

    !=========================================================================
    ! input arguments :
    TYPE(ion_basis) :: basis
    TYPE(material_data), DIMENSION(:), POINTER :: mat_data
    INTEGER :: i_atom

    !_________________________________________________________________________

    ! output result :

    INTEGER :: number_states_atom

    !_________________________________________________________________________
    
    ! local variables :

    CHARACTER ( LEN = 3 ) :: name_test

    INTEGER :: i_mat, i_ion

    !=========================================================================
  
    CALL check_ion( basis, i_atom, mat_data, i_mat, i_ion )
    
    number_states_atom = SIZE( mat_data( i_mat )%ion( i_ion )%state)
    
    !=========================================================================
    
  END FUNCTION number_states_atom

  !===========================================================================
  !
  ! Subroutine "sort_couplings" :
  !
  ! sort the TB couplings parameters with respect to the reference list.
  !
  !===========================================================================

  SUBROUTINE sort_couplings( i_parent, i_ref, TB_cpl, TB_pow )
    
    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER,                 INTENT( IN ) :: i_parent
    
    !=========================================================================
    ! input - output arguments
    !=========================================================================

    INTEGER,          DIMENSION( : ),    INTENT( INOUT ) :: i_ref
    DOUBLE PRECISION, DIMENSION( :, : ), INTENT( INOUT ) :: TB_cpl, TB_pow
        
    !=========================================================================
    ! Local variables
    !=========================================================================

    DOUBLE PRECISION, DIMENSION( SIZE( i_ref ) ) :: TB_cpl_aux
    DOUBLE PRECISION, DIMENSION( SIZE( i_ref ) ) :: TB_pow_aux
    DOUBLE PRECISION, DIMENSION( SIZE( i_ref ) ) :: cpl_tmp, pow_tmp

    INTEGER, DIMENSION( SIZE( i_ref ) ) :: ref_aux, ref_tmp

    LOGICAL, DIMENSION( SIZE( i_ref ) ) :: sort_test
    
    INTEGER :: mini, maxi, i_sort, n_test
    
    !========================================================================
    
    ref_aux = i_ref
    TB_cpl_aux = TB_cpl( i_parent, : )
    TB_pow_aux = TB_pow( i_parent, : )

    maxi = MAXVAL( ref_aux )
    
    i_sort = 0
    
    !=========================================================================
    
    DO WHILE ( i_sort .LT. SIZE( i_ref ) )
       
       !______________________________________________________________________
       
       mini = MINVAL( ref_aux )
       
       sort_test = ( ref_aux .EQ. mini )
       
       n_test = COUNT( sort_test )

       !______________________________________________________________________

       ref_tmp( i_sort + 1 : i_sort + n_test ) = PACK( ref_aux, sort_test )
       cpl_tmp( i_sort + 1 : i_sort + n_test ) = PACK( TB_cpl_aux, sort_test )
       pow_tmp( i_sort + 1 : i_sort + n_test ) = PACK( TB_pow_aux, sort_test )
      
       !______________________________________________________________________
       
       i_sort = i_sort + n_test
       
       WHERE( sort_test ) ref_aux = maxi + 1
       
       !______________________________________________________________________
       
    END DO
    
    !=========================================================================

    TB_cpl( i_parent, : ) = cpl_tmp
    TB_pow( i_parent, : ) = pow_tmp

    i_ref = ref_tmp

    !=========================================================================

  END SUBROUTINE  sort_couplings

  !===========================================================================

  SUBROUTINE set_max_order(mat_data)

    TYPE(material_data), DIMENSION(:), POINTER :: mat_data
    !TYPE(material_data), DIMENSION(:,:), POINTER :: interface_data

    INTEGER :: n_mat, i_mat!, i_mat_row, i_mat_column
  

    n_mat = SIZE(mat_data)
    
    DO i_mat = 1, n_mat
    IF (size(mat_data(i_mat)%parent) .gt. 0) THEN 
       mat_data(i_mat)%max_order= &
            MAXVAL( mat_data(i_mat)%parent(1)%pair(:)%order )
    END IF

    END DO

    !-------------------
    ! For the interfaces
    !-------------------
    ! DO i_mat_row = 1, n_mat
    !   DO i_mat_column = 1, n_mat
    !     IF (size(interface_data(i_mat_row, i_mat_column)%parent) .gt. 0) THEN 
    !       interface_data(i_mat_row, i_mat_column)%max_order= &
    !           MAXVAL( interface_data(i_mat_row, i_mat_column)%parent(1)%pair(:)%order )
    !     END IF
    !   END DO
    ! END DO

  END SUBROUTINE set_max_order


END MODULE states_and_couplings

!=============================================================================



