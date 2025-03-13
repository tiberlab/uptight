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
module checks

  USE precision
  USE type_defs
  USE errors
  USE exceptions

  implicit none
  private


  public :: check_ion, check_pair, check_interface_pair, check_uncommon_interface_pair, check_ref_pair
  public :: check_ref, check_ref_order 
  private :: swap_couplings

contains
  !===========================================================================
  !
  ! Subroutine "check_ion" : 
  !
  !===========================================================================
  !
  ! checks if an atom "a" is included in the ion data lists.  
  !
  ! if so, finds which element in the list matches this ion. 
  ! and returns the corresponding material and ion indexes :
  !
  ! atom "a" <--> ion_data( i_mat )%ion( i_ion )
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => i_atom - index of the "a" atom in the super basis.
  !
  !___________________________________________________________________________
  !
  ! OUTPUT :
  !
  ! => i_mat : integer = index of the material containing the atom "a".
  ! => i_ion : integer = index of the ion matching "a" in the data list.
  !
  !===========================================================================
  
  SUBROUTINE check_ion( basis, i_atom, material, i_mat, i_ion )

    !=========================================================================

    ! input arguments :
    TYPE(ion_basis) :: basis
    TYPE(material_data), DIMENSION(:), POINTER :: material
    INTEGER, INTENT( IN ) :: i_atom
 
    !_________________________________________________________________________
    
    ! output arguments :

    INTEGER, INTENT( OUT ) :: i_mat, i_ion
    CHARACTER(TYPELEN) :: str

    !=========================================================================
    i_mat = basis%mat(i_atom)
    str = basis%atomtypes(basis%type(i_atom))

    DO i_ion = 1, SIZE( material(i_mat)%ion )
       IF (material(i_mat)%ion(i_ion)%name .EQ. str ) RETURN
    END DO
    
    WRITE(*,*) 'ERROR: *',TRIM(str),'* not found ', &
                'in material', i_mat
                
    WRITE(*,*) material(i_mat)%ion(:)%name

    !!ERROR - incorrect reference atom list : check program algorithm.'
    call throw_solve_exception(ERR_REF_ATMLIST)    
    !=========================================================================

  END SUBROUTINE check_ion


  !===========================================================================
  !
  ! Subroutine "check_pair" :
  !
  !===========================================================================
  !
  ! checks if an atom pair "a" - "b" belonging to a single material is 
  ! included in the pair data list for the corresponding material index i_mat.
  !
  ! if so, finds which element in the list matches this pair. 
  ! and returns the corresponding parent material and pair index :
  !
  ! pair "a" - "b" <--> pair_data( i_mat )%parent( i_parent )%pair( i_pair )
  !
  ! For ALLOYS this is rather fake since different parents are mixed...
  !===========================================================================
  !
  ! INPUT :
  !
  ! => i_a, i_b : integer - index of the "a" and "b" atoms in the super basis.
  ! => i_mat    : integer = index of the material containing the pair.
  !___________________________________________________________________________
  !
  ! OUTPUT :
  !
  ! => i_parent : integer = index of the parent in the data list.
  ! => i_pair   : integer = index of the pair in the data list.
  !
  !===========================================================================

  SUBROUTINE check_pair(name_a, name_b, parent, i_order, o_parent, o_pair )

    !=========================================================================
    ! input arguments :
    CHARACTER(ATOMLEN) :: name_a, name_b
    TYPE(parent_data), DIMENSION(:) :: parent
    INTEGER, INTENT( IN ) :: i_order

    ! output arguments :
    INTEGER, INTENT( OUT ) :: o_parent, o_pair
    !_________________________________________________________________________

    ! local variables :
    CHARACTER(ATOMLEN), DIMENSION(2) :: name_pair, name_swap, name_test
    TYPE(pair_coupling), DIMENSION(:), POINTER :: pair_p
    INTEGER :: order_test, n_pair, i_parent, i_pair
    !=========================================================================
    ! Get the name of the atoms in the pair to check
    ! Swap them in case the order of the atoms is not that of the reference.
    ! 
    name_pair = (/ name_a, name_b /)
    name_swap = (/ name_b, name_a /)

    DO i_parent = 1, SIZE( parent )
        
       ! this check avoids to find pair with content = 0.0 
       ! LUAN: This IF is a bug when content < EPS but nonzero (e.g. single impurity case)
       !IF (abs(parent(i_parent)%content).lt.EPS) cycle

       ! scan the pair list for the current parent
       n_pair =  parent( i_parent )%n_pair
       pair_p => parent( i_parent )%pair 

       DO i_pair = 1, n_pair

          ! compare the pair name + the pair order to the reference
          name_test = pair_p( i_pair )%name
          order_test = pair_p( i_pair )%order
          
          ! if this matches the current data : return the parent + pair index
          IF ( ( ALL( name_pair .EQ. name_test ) .OR. &
               ALL( name_swap .EQ. name_test ) ) .AND. &
               ( i_order .EQ. order_test )  ) THEN
             o_parent = i_parent
             o_pair = i_pair
             RETURN
          END IF
       END DO
          
    END DO
    
    !=========================================================================
    
    ! if not a proper pair, data are not consistent :
    ! For instance, the program finds that Cation - Cation couplings are
    ! consistent with the order of neares neighbours required, but the
    ! coupling parameters for such pairs are missing in the datafile.
    
    o_pair = 0
    o_parent = 0

    !write_format = ' '
    !write_format = '( "ERROR in routine ''", a, "'' : ", //&
    !     &"no coupling parameters were included in the input&
    !     & files for the pair ''", a, " - ", a, "'', order ", i2,".")'
                                               
    !WRITE ( *, write_format ) &
    !     'check_pair', basis%name( i_a ), basis%name( i_b ), i_order
    
    
    !=========================================================================

  END SUBROUTINE check_pair

  !===========================================================================
  ! FOR AN INTERFACE PAIR OF ATOMS RETURNS THE 
  ! MATERIAL, PARENT AND PAIR 
  !
  SUBROUTINE check_interface_pair( basis, i_a, i_b, material, i_order, &
                                                interface_mat, parent, pair )

    !=========================================================================
    ! input arguments :
    TYPE(ion_basis) :: basis
    TYPE(material_data), DIMENSION(:), POINTER :: material
    INTEGER, INTENT( IN ) :: i_a, i_b, i_order

    !_________________________________________________________________________
    ! output arguments :

    INTEGER, INTENT( OUT ) :: parent, pair!ab_mat, 
    TYPE(material_data), POINTER :: interface_mat

    !_________________________________________________________________________

    ! local variables :
    INTEGER, DIMENSION(2) :: mat
    INTEGER :: order_test, n_pair, mat_a, mat_b
    INTEGER :: i_ion, i_mat, i_parent, i_pair, i_ab_mat
    CHARACTER(ATOMLEN) :: name_a, name_b

    !=========================================================================

    !
    ! For the moment :
    !
    ! the code looks at the pair which causes the problem 
    !
    ! --> atom a in material 1 / atom b in material 2
    !
    ! The code looks for the matching pair in the data basis of each
    ! material :
    !
    ! WE ASSUME HERE THAT THIS IS POSSIBLE --> COMMON ATOM in 1 and 2

    ! (for  Si/SiGe interface it's ok, 
    ! in "SiGe" there is information on Si <-> Ge and Si<-> Ge cpl 
    !
    ! When the code finds a matching pair ( e.g. in 1 ) it "pushes" the
    ! atom b in material 1, temporarily to get the ETB data.
    !

    ! Conventionally checks first in material with lowest content
    ! NOTE:
    ! Old BUG GaAs/In(100)Ga(0)As  and GaAs/InAs treated differently was solved (check_pair)
    !
    ! In general assigning the interface to either material A or B is a problem
    ! The couple As-Ga is taken always from one side but a better average could be made
    ! TODO: work at interface materials 
    ! RND Alloy works because matrix elements are not averged
    ! 
    ! ---------------------------------------------------------------------------------

    mat_a = basis%mat( i_a )
    mat_b = basis%mat( i_b )

    if (material(mat_a)%parent(1)%content .le. material(mat_b)%parent(1)%content) then  
        mat = (/ mat_a, mat_b /)
    else
        mat = (/ mat_b, mat_a /)
    endif

    name_a = basis%atomtypes(basis%type(i_a))
    name_b = basis%atomtypes(basis%type(i_b))
    
    DO i_mat = 1, 2   

       i_ab_mat = mat( i_mat )
       call check_pair( name_a, name_b, material(i_ab_mat)%parent, &
                                                      i_order, parent, pair)
       IF (parent .ne. 0) THEN
          interface_mat => material(i_ab_mat)
          RETURN
       END IF   

    END DO

    !=========================================================================
    
    ! if not a proper pair, data are not consistent :
    ! For instance, the program finds that Cation - Cation couplings are
    ! consistent with the order of neares neighbours required, but the
    ! coupling parameters for such pairs are missing in the datafile.
    pair = 0
    parent = 0
    
    !write_format = ' '
    !write_format = '( "ERROR in routine ''", a, "'' : ", //&
    !     &"no coupling parameters were included in the input&
    !     & files for the pair ''", a, " - ", a, "''." )'
    
    !WRITE ( *, write_format ) &
    !     'check_interface_pair', basis%name( i_a ), basis%name( i_b )
    
    
    !=========================================================================

  END SUBROUTINE check_interface_pair
  
  
  !===========================================================================
  ! FOR AN INTERFACE PAIR OF ATOMS RETURNS THE
  ! MATERIAL, PARENT AND PAIR
  
  SUBROUTINE check_uncommon_interface_pair( basis, i_a, i_b, interfaces, i_order, &
                                                interface_mat, parent, pair )

    !=========================================================================
    ! input arguments :
    TYPE(ion_basis) :: basis
    TYPE(material_data), DIMENSION(:,:), POINTER :: interfaces
    INTEGER, INTENT( IN ) :: i_a, i_b, i_order

    !_________________________________________________________________________
    ! output arguments :

    INTEGER, INTENT( OUT ) :: parent, pair
    TYPE(material_data), POINTER :: interface_mat

    !_________________________________________________________________________

    ! local variables :
    INTEGER, DIMENSION(2) :: mat
    INTEGER :: order_test, n_pair, mat_a, mat_b
    INTEGER :: i_ion, i_pair, i_ab_mat
    CHARACTER(ATOMLEN) :: name_a, name_b



    mat_a = basis%mat( i_a )
    mat_b = basis%mat( i_b )

    name_a = basis%atomtypes(basis%type(i_a))
    name_b = basis%atomtypes(basis%type(i_b))
    !write(*,*) "Check ", name_a, name_b, ' in mat ', mat_a, mat_b

    IF (associated(interfaces(mat_a, mat_b)%nr_parents)) THEN
    !write(*,*) 'n_parents in mat ab: ', size(interfaces(mat_a, mat_b)%nr_parents)

      call check_pair( name_a, name_b, &
                         interfaces(mat_a, mat_b)%parent, &
                         i_order, parent, pair)

      IF (parent .ne. 0) THEN
        interface_mat => interfaces(mat_a, mat_b)
        !write(*,*) "Found parent: ", interfaces(mat_a, mat_b)%parent(parent)%name
        RETURN
      END IF

    END IF

    IF (associated(interfaces(mat_b, mat_a)%nr_parents)) THEN

      !write(*,*) 'n_parents in mat ba: ', size(interfaces(mat_b, mat_a)%nr_parents)
      call check_pair( name_a, name_b, &
                         interfaces(mat_b, mat_a)%parent, &
                         i_order, parent, pair)

      IF (parent .ne. 0) THEN
        interface_mat => interfaces(mat_b, mat_a)
        !write(*,*) "Found parent: ", interfaces(mat_b, mat_a)%parent(parent)%name
        RETURN
      END IF

    END IF
                                                
    !=========================================================================

    ! if not a proper pair, data are not consistent :
    ! For instance, the program finds that Cation - Cation couplings are
    ! consistent with the order of neares neighbours required, but the
    ! coupling parameters for such pairs are missing in the datafile.
    pair = 0
    parent = 0

    !=========================================================================

  END SUBROUTINE check_uncommon_interface_pair

!!$  SUBROUTINE check_interface_pair( basis, i_a, i_b, parent, nr_parents, i_order, &
!!$                                               o_parent, o_pair, o_nr_parents )
!!$
!!$    !=========================================================================
!!$    ! input arguments :
!!$    TYPE(ion_basis) :: basis
!!$    TYPE(parent_data), DIMENSION(:) :: parent
!!$    INTEGER, DIMENSION(:), INTENT(IN) :: nr_parents
!!$    INTEGER, INTENT( IN ) :: i_a, i_b, i_order
!!$
!!$    ! output arguments :
!!$    INTEGER, INTENT( OUT ) :: o_parent, o_pair, o_nr_parents
!!$    !_________________________________________________________________________
!!$
!!$    ! local variables :
!!$    CHARACTER(TYPELEN) :: type_a, type_b
!!$    CHARACTER(ATOMLEN), DIMENSION(3) :: name_a, name_b
!!$    CHARACTER(ATOMLEN), DIMENSION(6) :: name_list
!!$    INTEGER :: nr_a, nr_b, nr, mat_a, mat_b
!!$    TYPE(ion_orbit), DIMENSION(:), POINTER :: p_ion
!!$    !=========================================================================
!!$
!!$    type_a = basis%atomtypes(basis%type(i_a))
!!$    type_b = basis%atomtypes(basis%type(i_b))
!!$  
!!$    ! -------------------------------------------------------------------------
!!$    ! Count the number of ions in type_a and type_b
!!$    ! Assumes: TYPELEN = 6 !!!!!!!!!!!!!!
!!$    IF (trim(name_a(ATOMLEN+1:2*ATOMLEN)).EQ."") THEN 
!!$       name_a(1) = type_a(1:ATOMLEN)
!!$       nr_a = 1
!!$    ELSE IF (trim(name_a(2*ATOMLEN+1:3*ATOMLEN)).EQ."") THEN
!!$       name_a(2) = type_a(ATOMLEN+1:2*ATOMLEN)
!!$       nr_a = 2
!!$    ELSE 
!!$       name_a(3) = type_a(2*ATOMLEN+1:3*ATOMLEN)
!!$       nr_a = 3
!!$    END IF
!!$    
!!$    IF ( trim(name_b(ATOMLEN+1:2*ATOMLEN)).EQ."") THEN 
!!$       name_b(1) = type_b(1:ATOMLEN)
!!$       nr_b = 1
!!$    ELSE IF (trim(name_b(2*ATOMLEN+1:3*ATOMLEN)).EQ."") THEN
!!$       name_b(2) = type_b(ATOMLEN+1:2*ATOMLEN)
!!$       nr_b = 2
!!$    ELSE 
!!$       name_b(3) = type_b(2*ATOMLEN+1:3*ATOMLEN)
!!$       nr_b = 3
!!$    END IF
!!$    ! -------------------------------------------------------------------------
!!$
!!$    nr = nr_a*nr_b
!!$
!!$    IF (nr_parents(1) .ne. nr_parents(2)) THEN
!!$       !
!!$       ! Al(x)GaAs / AlAs(y)P
!!$       !  4 1 AlAs AlP GaAs GaP AlAs xy x(1-y) y(1-x) (1-x)(1-y) 1
!!$       ! (AlGa)-(AsP) => nr_a=2,nr_b=2, o_parent = 1, nr_parents=4
!!$       ! Al-As, nr_a=1, nr_b=1, o_parent=4, nr_parents=1
!!$       
!!$       IF (nr == nr_parents(1)) THEN
!!$
!!$          CALL check_pair(name_a(1), name_b(1), parent(1:nr), &
!!$                                  i_order, i_parent, i_pair)
!!$
!!$          IF (i_parent .ne. 0) THEN
!!$             o_parent = 1
!!$             o_nr_parents = nr
!!$             o_pair = i_pair
!!$             RETURN
!!$          END IF
!!$
!!$       ELSE IF (nr == nr_parents(2)) THEN
!!$
!!$          nr_a = nr_parents(1)
!!$          CALL check_pair(name_a(1), name_b(1), parent(nr_a+1: nr_a+nr), &
!!$                                  i_order, i_parent, i_pair)
!!$
!!$          IF (i_parent .ne. 0) THEN
!!$             o_parent = nr_a + 1
!!$             o_nr_parents = nr
!!$             o_pair = i_pair
!!$             RETURN
!!$          END IF
!!$
!!$       ELSE
!!$          print*,"(check interface) ERROR" 
!!$       END IF
!!$
!!$    ELSE
!!$       !
!!$       ! Al(x)InAs / Al(y)GaP
!!$       !  2 2 AlP InP AlAs GaAs x 1-x y 1-y
!!$       ! AlIn-P => nr_a=2,nr_b=1, o_parent = 1, nr_parents=2
!!$       ! AlGa-As => nr_a=1,nr_b=2, o_parent = 3, nr_parents=2 
!!$       !
!!$       ! Al(x)GaAs(y)P / Al(z)GaAs(w)Sb
!!$       !  4 4 AlAs AlSb GaAs GaSb | AlAs AlP GaAs GaP  
!!$       ! (AlGa)-(AsP) =>   
!!$       ! (AlGa)-(AsSb)
!!$       !
!!$       DO i = 1, nr_a
!!$          name_list( ) = name_a(i)
!!$          DO j = 1, nr_b
!!$             insert(name_b(j))
!!$          ENDDO
!!$       ENDDO
!!$ 
!!$       DO i_parent = 1, nr_parents(1)
!!$          p_ion => parent(i_parent)%ion
!!$          DO i_ion = 1, size(p_ion)
!!$             
!!$
!!$
!!$
!!$
!!$    END IF
!!$
!!$
!!$    CALL check_pair(name_a(i), name_b(j), parent(1:nr), &
!!$         i_order, i_parent, i_pair)
!!$    
!!$    nr_a = nr_parents(1) + 1
!!$
!!$ 
!!$    !=========================================================================
!!$    
!!$    ! if not a proper pair, data are not consistent :
!!$    ! For instance, the program finds that Cation - Cation couplings are
!!$    ! consistent with the order of neares neighbours required, but the
!!$    ! coupling parameters for such pairs are missing in the datafile.
!!$
!!$    o_pair = 0
!!$    o_parent = 0
!!$    
!!$    !write_format = ' '
!!$    !write_format = '( "ERROR in routine ''", a, "'' : ", //&
!!$    !     &"no coupling parameters were included in the input&
!!$    !     & files for the pair ''", a, " - ", a, "''." )'    
!!$    !=========================================================================
!!$
!!$  END SUBROUTINE check_interface_pair
!!$  
  !===========================================================================

  !===========================================================================
  !
  ! Subroutine "check_ref_pair"
  !
  ! Sort the pair couplings according to the reference list.
  ! 
  ! for giver atoms i_a and i_b  and a refernce coupling list  
  ! returns a permutation array such that:
  ! 
  ! mat(i_mat)%parent(i_parent)%pair(i_pair)%cpl(k) = ref_cpl( ref_pair(k) )
  !===========================================================================

  SUBROUTINE check_ref_pair(basis, i_a, i_b, pair, ref_pair )

    !=========================================================================  
    ! input arguments :
    TYPE(ion_basis) :: basis
    TYPE(pair_coupling) :: pair
    INTEGER, INTENT( IN ) :: i_a, i_b

    ! output arguments :

    INTEGER, DIMENSION( : ), POINTER :: ref_pair

    ! local variables :
    
    !CHARACTER(CPLLEN), DIMENSION( : ), ALLOCATABLE :: coupling_swap
    
    CHARACTER(ATOMLEN), DIMENSION( 2 ) :: pair_name, pair_name_test
    CHARACTER(TYPELEN) :: name_a, name_b
    INTEGER :: n_cpl, err
    
    !=========================================================================
    name_a = basis%atomtypes(basis%type(i_a))
    name_b = basis%atomtypes(basis%type(i_b))

    pair_name(1) = name_a(1:ATOMLEN)
    pair_name(2) = name_b(1:ATOMLEN)
    pair_name_test = pair%name
    
    !-------------------------------------------------------------------------
    ! Test on the ion order
    !-------------------------------------------------------------------------

    IF ( ALL( pair_name .EQ. pair_name_test ) ) THEN
   
       !----------------------------------------------------------------------
       ! The ions agree with the data input
       !----------------------------------------------------------------------
       ref_pair = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21] 
       !CALL check_ref( material( i_mat )%parent( i_parent )&
       !     %pair( i_pair )%cpl, ref_couplings, ref_pair )        
      
    ELSE
       !----------------------------------------------------------------------
       ! the ions must be swapped to agree with the data input
       !----------------------------------------------------------------------

       ref_pair = [1,3,2,4,5,6,8,7,10,9,12,11,14,13,16,15,18,17,19,20,21] 

       !ref_pair = [1,3,2,4,5,7,6,9,8,10,12,11,14,13,16,15,18,17,19,20,21]        
       !n_cpl = SIZE( material( i_mat )%parent( i_parent )%pair( i_pair )%cpl )

       !ALLOCATE( coupling_swap( n_cpl ), STAT = err )
       !IF ( err .NE. 0 ) CALL alloc_error( 'states_and_couplings', &
       !     'check_ref_pair', 'coupling_swap' )

       !CALL swap_couplings(material(i_mat)%parent(i_parent)%pair(i_pair)%cpl, & 
       !                       coupling_swap )

       !CALL check_ref( coupling_swap, ref_couplings, ref_pair )
       
       !DEALLOCATE( coupling_swap, STAT = err )
       !IF ( err .NE. 0 ) CALL dealloc_error( 'states_and_couplings', &
       !     'check_ref_pair', 'coupling_swap' )
       
       
    END IF

    !=========================================================================

  END SUBROUTINE check_ref_pair

  !===========================================================================
  !
  ! Subroutine "check_ref" :
  !
  ! check if all elements in the array "tab" are contained in the reference
  ! array "ref_tab". If so, calculate the index array "ind_tab" giving the
  ! mapping of "tab" onto "ref_tab" :
  !
  !                       tab( i ) = ref_tab( ind_tab( i ) )
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => tab - character strings array : contains the names list to chek.
  ! => ref_tab - character strings array : contains a reference name list.
  !
  ! OUTPUT :
  !
  ! => ind_tab - integer array : gives the mapping of "tab" onto "ref_tab".
  !
  !===========================================================================
  
  SUBROUTINE check_ref( tab, ref_tab, ind_tab )

    CHARACTER(CPLLEN), DIMENSION( : ), INTENT( IN ) :: tab, ref_tab

    !_________________________________________________________________________

    INTEGER, DIMENSION( SIZE( tab ) ), INTENT( INOUT ) :: ind_tab
    !_________________________________________________________________________

    ! local variables :
    INTEGER :: i_tab, i_ref

    DO i_tab = 1, SIZE( tab )

       i_ref = 1
       
       DO WHILE ( tab( i_tab ) .NE. ref_tab( i_ref ) )
          
          i_ref = i_ref + 1
          IF ( i_ref .GT. SIZE( ref_tab ) ) then
             CALL throw_init_exception(ERR_REF_NAME) 
          END IF

       END DO

       ind_tab( i_tab ) = i_ref

    END DO

  END SUBROUTINE check_ref

   
  SUBROUTINE check_ref_order( tab, ref_tab, ind_tab )

    CHARACTER(CPLLEN), DIMENSION( : ), INTENT( IN ) :: tab, ref_tab

    !_________________________________________________________________________

    INTEGER, DIMENSION( SIZE( tab ) ), INTENT( INOUT ) :: ind_tab
    !_________________________________________________________________________

    ! local variables :
    INTEGER :: i_tab, i_ref, n_ref

    n_ref = 0
    DO i_ref = 1, SIZE( ref_tab )

      DO i_tab = 1, SIZE( tab )
       
        IF ( tab( i_tab ) .EQ. ref_tab( i_ref ) ) THEN
          
          n_ref = n_ref + 1
          ind_tab( i_tab ) = n_ref

          EXIT

        END IF
      END DO
    END DO

    IF ( n_ref .NE. SIZE( tab ) ) THEN
      CALL throw_init_exception(ERR_REF_NAME) 
    END IF

  END SUBROUTINE check_ref_order


  !===========================================================================
  !
  ! Subroutine "swap_couplings" :
  !
  !==========================================================================

  SUBROUTINE swap_couplings( coupling, coupling_swap )
    
    !=========================================================================
    
    ! input output arguments :
    
    CHARACTER( CPLLEN ),  DIMENSION( : ), POINTER :: coupling
    
    ! local variables :

    CHARACTER( CPLLEN ), DIMENSION( SIZE( coupling ) ), INTENT( OUT ) :: &
         coupling_swap

    !=========================================================================

    coupling_swap = ' '

    coupling_swap = coupling

    WHERE ( coupling .EQ. 'ss*s' ) coupling_swap = 's*ss'
    WHERE ( coupling .EQ. 's*ss' ) coupling_swap = 'ss*s'
    
    WHERE ( coupling .EQ. 'sps' ) coupling_swap = 'pss'
    WHERE ( coupling .EQ. 'pss' ) coupling_swap = 'sps'

    WHERE ( coupling .EQ. 's*ps' ) coupling_swap = 'ps*s'
    WHERE ( coupling .EQ. 'ps*s' ) coupling_swap = 's*ps'

    WHERE ( coupling .EQ. 'sds' ) coupling_swap = 'dss'
    WHERE ( coupling .EQ. 'dss' ) coupling_swap = 'sds'

    WHERE ( coupling .EQ. 's*ds' ) coupling_swap = 'ds*s'
    WHERE ( coupling .EQ. 'ds*s' ) coupling_swap = 's*ds'

    WHERE ( coupling .EQ. 'pds' ) coupling_swap = 'dps'
    WHERE ( coupling .EQ. 'dps' ) coupling_swap = 'pds'

    WHERE ( coupling .EQ. 'pdp' ) coupling_swap = 'dpp'
    WHERE ( coupling .EQ. 'dpp' ) coupling_swap = 'pdp'

    !=========================================================================

  END SUBROUTINE swap_couplings


end module checks
