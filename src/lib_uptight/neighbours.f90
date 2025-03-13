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
!        Module "nearest_neighbours" - Alessandro Pecchia  
!
!=============================================================================
!
! contains routines used to build and refine a neighbour map in a supercell :
!
! Public:
! => make_neighbour_map: build the neighbour list for unstrained supercell
! => refine_neighbour_map: refine a neighbour map of 2nn starting from 1nn
! => is_dg_bond     : function used to test dg_bonds 
!
!
!=============================================================================

MODULE neighbours

  !===========================================================================

  USE precision, only : dp, equiv, accur, fuzzy
  USE type_defs, only : ATOMLEN, ion_basis, TStructure, material_data
  USE checks, only : check_pair, check_interface_pair
  USE supercell, only : latt_max, latt_dim, lattice_vectors
  USE errors, only :  alloc_error, dealloc_error, lattice_errors
  USE vectors
  USE input_output
  USE sort, only : osort_index  
  USE exceptions
  !===================================================================

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: nearest_neighbours
  PUBLIC :: create_neighbours_map
  PUBLIC :: make_neighbours_map
  PUBLIC :: refine_neighbours_map
  PUBLIC :: write_neighbours_map
  PUBLIC :: destroy_neighbours_map

  PUBLIC :: check_input_nn_list
  PUBLIC :: check_nn_map
  PUBLIC :: is_dg_bond
  PUBLIC :: write_am_coo

  PRIVATE :: select_nearest_neighbours, nearest_test 
  PRIVATE :: sort_near_by_ind, sort_near_by_dist
  PRIVATE :: reallocate

  INTEGER, PARAMETER :: TMPDIM = 50

  !________________________________________________________________
  ! NEW WAY USING AN ARRAY
  ! 
  ! ind contains the atom indeces
  ! dist contains the atom distances
  ! vec contains the image cell position of each atom
  ! order contains the atom tb order (1st, 2nd, ...)
  ! n_cpl is the effective n of couplings (not self-counting)
  !       so n_cpl=size(ind)-1
  TYPE nearest_neighbours
     REAL ( dp ), DIMENSION( : ),    POINTER :: dist
     INTEGER,     DIMENSION( :, : ), POINTER :: vec
     INTEGER,     DIMENSION( : ),    POINTER :: ind
     INTEGER,     DIMENSION( : ),    POINTER :: order
     INTEGER                                 :: n_cpl
  END TYPE nearest_neighbours

  !________________________________________________________________
  
  INTERFACE reallocate
     MODULE PROCEDURE reallocate_i, reallocate_d
     MODULE PROCEDURE reallocate_i2
  END INTERFACE
  

  INTEGER, DIMENSION(3, latt_dim) :: latt_vector

CONTAINS

  subroutine create_neighbours_map(nn_map, natoms)
    type(nearest_neighbours), dimension(:), pointer :: nn_map
    integer :: natoms, err

    if (.not.associated(nn_map)) then
      allocate(nn_map(natoms), stat=err)
      if (err .ne. 0 ) then
       write(*,*) "Allocation error"
       call throw_init_exception(ERR_NN_LIST)
      endif
    else
      write(*,*) 'neighbour map already allocated'
      call throw_init_exception(ERR_NN_LIST)
    endif 

  end subroutine create_neighbours_map

  !===========================================================================
  ! make_nearest_neighbours: 
  !
  ! build the neighbour list from a general atomic structure
  ! no 1st n.n. list are already present.
  !===========================================================================
  subroutine make_neighbours_map(str, basis, mat_data, neig_map)

  TYPE( TStructure), INTENT( IN ) :: str 
  TYPE( ion_basis ), INTENT( INOUT ) :: basis  
  TYPE( material_data ), DIMENSION(:), POINTER :: mat_data
  TYPE( nearest_neighbours ), DIMENSION(:), POINTER :: neig_map

  ! ------------------------------------------------------------------------- 
  ! locals: 
  INTEGER :: n_n, H_count, ia, ib, err
  INTEGER :: i,j,k,l, i_mat, i_cell

  INTEGER, DIMENSION(:), ALLOCATABLE :: nn_ind, nn_order
  REAL(dp), DIMENSION(:), ALLOCATABLE :: nn_dist
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nn_vec

  REAL ( dp ), DIMENSION(:), POINTER :: tmp_dist
  REAL ( dp ), DIMENSION(3) :: vec
  INTEGER, DIMENSION(3) :: t_vec
  REAL ( dp ) :: dist, ref_dist

  ! -------------------------------------------------------------------------
  ! Set the lattice vectors using this routine defined in supercell
  ! -------------------------------------------------------------------------

  IF (.not.ASSOCIATED(neig_map)) THEN
    call create_neighbours_map(neig_map, basis%n_basis+basis%n_dg_bond)
  END IF

  DO i = 1, basis%n_basis + basis%n_dg_bond

     ! ------------------------------------------------------------------
     ! Construct neighbour list if not defined 
     ! ------------------------------------------------------------------
     call select_nearest_neighbours(i, basis, mat_data, nn_ind, &
                                                   nn_order, nn_dist, nn_vec )

     n_n = size(nn_ind)

     !-------------------------------------------------------------------------
     ! Allocates arrays
     !-------------------------------------------------------------------------
     ALLOCATE( neig_map(i)%ind( n_n + 1 ), STAT = err )
     ALLOCATE( neig_map(i)%order( n_n + 1 ), STAT = err )
     ALLOCATE( neig_map(i)%dist( n_n + 1 ), STAT = err )
     ALLOCATE( neig_map(i)%vec( 3, n_n + 1 ), STAT = err )
     IF(err.NE.0) CALL alloc_error('nearest_neighbours','make_nn','new_near')

     ! ------------------------------------------------------------------------
     ! Initializations
     ! ------------------------------------------------------------------------
     neig_map(i)%dist = 0.0D0 
     neig_map(i)%vec = 0     
     neig_map(i)%order = 1

     neig_map(i)%order(1) = 0 ! onsite
     neig_map(i)%ind(1) = i   ! onsite

     !-------------------------------------------------------------------------
     ! Set neighbour indeces and count Hydrogens (dg bonds)
     !-------------------------------------------------------------------------
     H_count = 0
     DO j = 1, n_n

       neig_map(i)%ind( j + 1 ) = nn_ind( j )
       neig_map(i)%order( j + 1 ) = nn_order( j )
       neig_map(i)%dist( j + 1 ) = nn_dist( j )
       neig_map(i)%vec(:, j + 1 ) = nn_vec(:, j )

        IF( nn_ind(j) .GT. basis%n_basis ) THEN
           H_count = H_count + 1                         
        END IF 

     ENDDO

     basis%n_dg(i) = H_count    

     deallocate(nn_ind, nn_order, nn_dist, nn_vec)
     !----------------------------------------------------------------------
     ! Sort with respect to atom index
     !----------------------------------------------------------------------     
     CALL sort_near_by_ind( neig_map(i) )
     !----------------------------------------------------------------------
     ! Count number of couplings in folded cell (drops periodic copies)
     !----------------------------------------------------------------------
     neig_map(i)%n_cpl = 0
     DO j = 1, SIZE( neig_map(i)%ind ) - H_count

        k = neig_map(i)%ind(j)
        neig_map(i)%n_cpl = neig_map(i)%n_cpl + 1

        check2: DO l = 1, j-1
           IF (k .EQ. neig_map(i)%ind(l)) THEN 
              neig_map(i)%n_cpl = neig_map(i)%n_cpl - 1              
              EXIT check2
           END IF
        END DO check2

     END DO

  END DO
  
  END subroutine make_neighbours_map

  !===========================================================================
  ! refine_nearest_neighbours: 
  !
  ! Starts from a 1st n.n. list on input and add 2nd or 3rd n.n. if needed. 
  ! 2nd and 3rd nearest_neighbours are searched among the neighbours of neighbours
  !===========================================================================
  SUBROUTINE refine_neighbours_map(str, basis, mat_data, neig_map)

  TYPE( TStructure), INTENT( IN ) :: str 
  TYPE( ion_basis ), INTENT( INOUT ) :: basis  
  TYPE( material_data ), DIMENSION(:), POINTER :: mat_data
  TYPE( nearest_neighbours ), DIMENSION(:), POINTER :: neig_map

  ! ------------------------------------------------------------------------- 
  ! locals: 
  INTEGER :: n_n, H_count, ia, ib, ifnn, n2_n, err
  INTEGER :: i,j,k,l, i_mat, i_cell, i_parent, i_pair, i_ord, n_dg
  INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_ind, tmp_ord
  INTEGER, DIMENSION(:), POINTER :: ind_p
  REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET :: tmp_dist
  REAL(dp), DIMENSION(:), POINTER :: dist_p
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmp_vec
  REAL ( dp ), DIMENSION(3) :: vec
  INTEGER, DIMENSION(3) :: t_vec
  REAL ( dp ) :: dist, dist_ref
  CHARACTER(ATOMLEN) :: name_a, name_b 
  ! -------------------------------------------------------------------------
  ! Set the lattice vectors using this routine defined in supercell
  ! -------------------------------------------------------------------------
  CALL lattice_vectors( latt_vector )

  IF (.not.ASSOCIATED(neig_map)) THEN
    call create_neighbours_map(neig_map, basis%n_basis+basis%n_dg_bond)
  END IF

  DO i = 1, basis%n_basis + basis%n_dg_bond

     ! define original position of atom i using inv_indexa(i) 
     ia = str%inv_indexa(i)

     ! ------------------------------------------------------------------
     ! Construct neighbour list if not defined 
     ! ------------------------------------------------------------------
     if ( ANY( str%nn_list(:,ia) .NE. 0) ) then

        n_n = COUNT( str%nn_list(:, ia) .GT. 0  ) 
   
     else

        call throw_init_exception(ERR_NN_LIST)

     end if

     !-------------------------------------------------------------------------
     ! Allocates arrays
     !-------------------------------------------------------------------------
     ALLOCATE( neig_map(i)%ind(n_n + 1), STAT = err )
     ALLOCATE( neig_map(i)%order( n_n + 1 ), STAT = err )
     ALLOCATE( neig_map(i)%dist( n_n + 1 ), STAT = err )
     ALLOCATE( neig_map(i)%vec( 3, n_n + 1 ), STAT = err )
     IF(err.NE.0) CALL alloc_error('nearest_neighbours','make_nn','new_near')

     !-------------------------------------------------------------------------
     ! Set neighbour indeces and count Hydrogens (dg bonds)
     !-------------------------------------------------------------------------
     neig_map(i)%ind = 0
     neig_map(i)%ind(1) = i
     neig_map(i)%order = 1
     neig_map(i)%order(1) = 0 !onsite
     neig_map(i)%dist = 0.d0
     neig_map(i)%vec = 0
     H_count = 0

     DO j = 1, n_n

        ib = str%indexa( str%nn_list( j, ia) ) 
        neig_map(i)%ind( j + 1 ) = ib

        IF( ib .GT. basis%n_basis ) THEN
           H_count = H_count + 1                         
        END IF 

     ENDDO

     basis%n_dg(i) = H_count    

     ! ------------------------------------------------------------------------
     ! Loop on all nearest_neighbours and find vectors and distances
     ! ------------------------------------------------------------------------
     ! NOTES:
     ! 
     ! near%vec()   :   is in fractional coordinates
     ! near%dist()  :   is in absolute coordinates
     ! near%order() :   follows the distance order
     ! ------------------------------------------------------------------------
     H_count = 0
     l = 2

     IF (basis%periodic) THEN

        ALLOCATE(tmp_ind(latt_dim * n_n + 1), STAT = err )
        ALLOCATE(tmp_dist(latt_dim * n_n + 1), STAT = err )
        ALLOCATE(tmp_vec(3,latt_dim * n_n + 1), STAT = err )
        tmp_ind = 0
        tmp_ind(1) = neig_map(i)%ind(1)
        tmp_dist = 0.d0
        tmp_vec = 0

        ! loop over neighbour cells to find the right translation vector
        DO i_cell = 1, latt_dim

           t_vec(:) = latt_vector(:, i_cell)
                 
           DO j = 2, n_n + 1
              
              ! computes distance vector in absolute coordinates:
              k = neig_map(i)%ind(j)
              
              vec =MATMUL(basis%prim,basis%coord(:,i)-basis%coord(:,k))
              
              dist = norm( vec -  MATMUL(basis%prim, t_vec) )
              
              IF (.not.already_in(l, k, tmp_ind, tmp_vec, t_vec)) THEN

                 !write(*,'(i4,i4,f8.4,i4,i4,i4)') l,k,dist, t_vec
                 tmp_ind( l ) = k
                 tmp_dist( l ) = dist
                 tmp_vec(1:3, l ) = t_vec(1:3)
                 
                 l = l + 1

              END IF

           ENDDO
        ENDDO
          
     ELSE

        ALLOCATE(tmp_ind(n_n+1), STAT = err )
        ALLOCATE(tmp_dist(n_n+1), STAT = err )
        ALLOCATE(tmp_vec(3,n_n+1), STAT = err )
        tmp_ind = 0
        tmp_ind(1) = neig_map(i)%ind(1)
        tmp_dist = 0.d0
        tmp_vec = 0

        DO j = 2, n_n + 1
        
           ! computes distance vector in absolute coordinates:
           k = neig_map(i)%ind(j)
           
           vec =MATMUL(basis%prim,basis%coord(:,i))-&
                MATMUL(basis%prim,basis%coord(:,k))
           
           dist = norm( vec )

           tmp_ind( l ) = k
           tmp_dist( l ) = dist
           tmp_vec(1:3, l ) = 0

           l = l + 1

        ENDDO
        
     END IF

     !-------------------------------------------------------------------
     ! Reorder the nearest_neighbours by distance and find tb order.
     !------------------------------------------------------------------- 
     ALLOCATE( ind_p(l-1) )
     dist_p => tmp_dist
     CALL osort_index(l-1,dist_p,ind_p)

     n_dg = SUM( basis%n_dg( 1:i-1 ) )
     !----------------------------------------------------------------------
     ! Set TB order 1 for neighbour up to the ion valence 
     !----------------------------------------------------------------------
     first:DO j = 2, n_n + 1
        
        neig_map(i)%ind(j) = tmp_ind( ind_p(j) )
        neig_map(i)%dist(j) = tmp_dist( ind_p(j) )
        neig_map(i)%vec(1:3,j) = tmp_vec(1:3, ind_p(j) )
        neig_map(i)%order(j) = 1

        IF (  neig_map(i)%ind( j ) .GT. basis%n_basis ) THEN
              
           H_count = H_count + 1
           basis%dg_coord( :, n_dg + H_count ) = &
                basis%coord(:, neig_map(i)%ind(j)) + neig_map(i)%vec(:,j)
           
           neig_map(i)%vec(:,j) = 0
           
        END IF

     END DO first

     DEALLOCATE(tmp_ind,tmp_dist,tmp_vec,ind_p)

     !----------------------------------------------------------------------
     ! Search for 2nd nearest_neighbours
     !---------------------------------------------------------------------- 
     l = n_n + 1  ! set at the last element of tmp_ind

     i_mat = basis%mat( i )
     i_ord = 2
     
     IF ( mat_data(i_mat)%max_order .GT. 1 ) THEN

        ALLOCATE(tmp_ind(TMPDIM))
        ALLOCATE(tmp_dist(TMPDIM))
        ALLOCATE(tmp_vec(3,TMPDIM))
        ALLOCATE(tmp_ord(TMPDIM))
        tmp_ind(1:n_n+1) = neig_map(i)%ind(1:n_n+1)
        tmp_dist(1:n_n+1) = neig_map(i)%dist(1:n_n+1)
        tmp_vec(:,1:n_n+1) = neig_map(i)%vec(:,1:n_n+1)
        tmp_ord(1:n_n+1) = neig_map(i)%order(1:n_n+1)        

        DO j = 1, n_n

           ifnn = str%inv_indexa( tmp_ind( j + 1 ) )
           n2_n=COUNT( str%nn_list(:, ifnn) .GT. 0  )

           DO k = 1, n2_n  

              ib = str%indexa( str%nn_list(k, ifnn) )

              name_a = basis%atomtypes(basis%type(ia))
              name_b = basis%atomtypes(basis%type(ib))


              CALL check_pair(name_a, name_b, mat_data(i_mat)%parent, i_ord, &
                                                         i_parent, i_pair )

              
              IF (i_parent .NE. 0) THEN
     
                 dist_ref = mat_data(i_mat)%parent(i_parent) &
                                                        %pair(i_pair)%dist_ref

         
                 CALL nearest_test( ia, ib, basis, dist_ref, 2, &
                                      l, tmp_ind, tmp_ord, tmp_dist, tmp_vec) 
                 
                 
              END IF

           END DO
        
        END DO
 
        !-------------------------------------------------------------------
        ! Rellocates arrays
        !-------------------------------------------------------------------
        n_n = l - 1 

        call REALLOCATE( neig_map(i)%ind, n_n + 1)
        call REALLOCATE( neig_map(i)%order, n_n + 1 )
        call REALLOCATE( neig_map(i)%dist, n_n + 1 )
        call REALLOCATE( neig_map(i)%vec, 3, n_n + 1 )

        ! ------------------------------------------------------------------
        ! Copy from temp arrays
        ! ------------------------------------------------------------------
        
        neig_map(i)%order(1:n_n+1) = tmp_ord(1:n_n+1)
        neig_map(i)%ind(1:n_n+1) = tmp_ind(1:n_n+1)
        neig_map(i)%dist(1:n_n+1) =  tmp_dist(1:n_n+1) 
        neig_map(i)%vec(:,1:n_n+1) = tmp_vec(:,1:n_n+1)   

        DEALLOCATE(tmp_ind,tmp_ord,tmp_dist,tmp_vec)
        
     END IF

     !----------------------------------------------------------------------
     ! Sort with respect to atom index
     !----------------------------------------------------------------------     
     CALL sort_near_by_ind( neig_map(i) )
     !----------------------------------------------------------------------
     ! Count number of couplings in folded cell (drops periodic copies)
     !----------------------------------------------------------------------
     neig_map(i)%n_cpl = 0
     DO j = 1, SIZE( neig_map(i)%ind ) - H_count

        k = neig_map(i)%ind(j)
        neig_map(i)%n_cpl = neig_map(i)%n_cpl + 1

        check2: DO l = 1, j-1
           IF (k .EQ. neig_map(i)%ind(l)) THEN 
              neig_map(i)%n_cpl = neig_map(i)%n_cpl - 1              
              EXIT check2
           END IF
        END DO check2

     END DO

  END DO  ! i = 1 , basis%n_basis

  END subroutine refine_neighbours_map

  !===========================================================================
  !
  ! Subroutine select_nearest_neighbours( i_a, basis, mat_data, neighbour_ind)
  !
  !===========================================================================
  !
  ! Select all the basis atoms which lie in the 3D nearest nearest_neighbours' 
  ! range of an atom.
  !
  !===========================================================================
  !
  ! * input arguments :
  !
  !   --> i_a : index of atom "a" in the super basis
  !
  !===========================================================================
  ! Example:
  !
  ! (1) (2)(1) (2)
  !   (3)    (3)
  ! (4) (5)(4) (5)     5:  3 4 2 4
  ! (1) (2)(1) (2)
  !   (3)    (3)
  !
  !---------------------------------------------------------------------------
  SUBROUTINE select_nearest_neighbours( i_a, basis, mat_data, nn_ind, nn_order, &
                                                             nn_dist, nn_vec )

    !=========================================================================
    ! Input arguments
    !=========================================================================

    INTEGER, INTENT( IN ) :: i_a
    TYPE (ion_basis), INTENT(IN) :: basis       
    TYPE (material_data), DIMENSION(:), POINTER :: mat_data
    TYPE(material_data), POINTER :: interface_data
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nn_ind, nn_order
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nn_dist 
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: nn_vec 

    !=========================================================================
    ! Local variables
    !=========================================================================
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_ind, tmp_ord
    REAL(dp), DIMENSION(:), ALLOCATABLE :: tmp_dist
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmp_vec

    INTEGER :: i_n, i_ord, i_parent, i_pair, i_mat, i_a_mat
    INTEGER :: k, n_n, n_tmp, err
    LOGICAL :: test
    REAL(dp) :: dist_ref
    CHARACTER(ATOMLEN) :: name_a, name_b 
    !=========================================================================
    ! Loop on all basis atoms (can be made more efficient by pre-indexing)
    !=========================================================================
    ALLOCATE(tmp_ind(TMPDIM))
    ALLOCATE(tmp_ord(TMPDIM))
    ALLOCATE(tmp_dist(TMPDIM))
    ALLOCATE(tmp_vec(3,TMPDIM))

    n_n = 0
    i_a_mat = basis%mat(i_a)

    DO i_ord = 1, mat_data( i_a_mat )%max_order 

       DO i_n = 1, basis%n_basis + basis%n_dg_bond

          IF (i_n .EQ. i_a) CYCLE

          i_mat = basis%mat(i_n)

          IF (i_mat .EQ. i_a_mat) THEN

             name_a = basis%atomtypes(basis%type(i_a))
             name_b = basis%atomtypes(basis%type(i_n)) 

             CALL check_pair(name_a, name_b, mat_data(i_mat)%parent, i_ord, &
                                                         i_parent, i_pair )

          ELSE
             
             CALL check_interface_pair( basis, i_a, i_n, mat_data, i_ord, &
                                             interface_data, i_parent, i_pair)

          ENDIF

          IF (i_parent .EQ. 0) CYCLE

          dist_ref = interface_data%parent(i_parent)%pair(i_pair)%dist_ref
          
          CALL nearest_test( i_a, i_n, basis, dist_ref, i_ord, &
                                     n_n, tmp_ind, tmp_ord, tmp_dist, tmp_vec) 


       END DO

    END DO

    !=========================================================================
    ! Allocate output array
    !=========================================================================

    ALLOCATE( nn_ind( n_n ), STAT = err )
    ALLOCATE( nn_order( n_n ), STAT = err )
    ALLOCATE( nn_dist( n_n ), STAT = err )
    ALLOCATE( nn_vec( 3, n_n ), STAT = err )

    IF (err.NE.0) CALL alloc_error('nearest_neighbours', 'select_neighbours', 'nn')

    !=========================================================================
    ! Set neighbour indeces
    !=========================================================================
    nn_ind(1:n_n) = tmp_ind(1:n_n)
    nn_order(1:n_n) = tmp_ord(1:n_n)
    nn_dist(1:n_n) = tmp_dist(1:n_n)
    nn_vec(:,1:n_n) = tmp_vec(:,1:n_n)  

    DEALLOCATE(tmp_ind,tmp_ord,tmp_dist,tmp_vec)

  END SUBROUTINE select_nearest_neighbours


  !===========================================================================
  !
  ! "nearest_test"
  !
  !===========================================================================
  !
  ! Checks if a given couple of atoms "a" and "b" belong to a 3D
  ! nearest nearest_neighbours' range.
  !
  ! Atom "a" has a fixed position inside the central supercell.
  !
  !___________________________________________________________________________
  !
  ! INPUT :
  !
  ! i_a, i_b : integer - super - basis index of atoms "a" and "b"
  !
  !___________________________________________________________________________
  !
  ! OUTPUT :
  !
  ! returns the number of nearest_neighbours (repeated in supercells)
  !
  !===========================================================================

  SUBROUTINE nearest_test( i_a, i_b, basis, ref_dist, i_ord, n_n,  &
                                         tmp_ind, tmp_ord, tmp_dist, tmp_vec )

    !-------------------------------------------------------------------------
    ! Input arguments
    !----------------------------------------------------------------------

    INTEGER, INTENT( IN ) :: i_a, i_b
    TYPE ( ion_basis ), INTENT(IN) :: basis    
    REAL ( dp ), INTENT( IN ) :: ref_dist
    INTEGER, INTENT( IN ) :: i_ord

    !-------------------------------------------------------------------------
    ! Output 
    !-------------------------------------------------------------------------
    
    INTEGER, INTENT(INOUT) :: n_n
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_ind, tmp_ord
    REAL(dp), DIMENSION(:), ALLOCATABLE :: tmp_dist
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmp_vec

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------

    INTEGER :: i_cell,err
    REAL ( dp ), DIMENSION(3) :: vec 
    INTEGER, DIMENSION(3) :: t_vec
    REAL ( dp ) :: dist

    !-------------------------------------------------------------------------
    ! Calculate the distance between atoms. 
    ! If distance is less than a cutoff, the two atoms  
    ! are considered interacting nearest_neighbours
    !-------------------------------------------------------------------------
    vec = MATMUL( basis%prim, basis%coord(:,i_a) - basis%coord(:,i_b) )
    
    if (basis%periodic) then
       
       do i_cell = 1, latt_dim
          
          t_vec(:) = latt_vector(:, i_cell )
          dist = norm( vec -  MATMUL(basis%prim, t_vec) )

          IF ( equiv( dist, ref_dist, 0.30d0, .false.) .AND. &
             .NOT. already_in(n_n+1, i_b, tmp_ind, tmp_vec, t_vec) ) THEN

             n_n = n_n + 1
             tmp_ind( n_n ) = i_b
             tmp_ord( n_n ) = i_ord
             tmp_dist( n_n ) = dist
             tmp_vec(:, n_n ) = t_vec(:)
 
          END IF

       end do
       
    else

       dist = norm( vec )

       IF ( equiv( dist, ref_dist, 0.30d0, .false.) ) THEN
          
          n_n = n_n + 1
          tmp_ind( n_n ) = i_b
          tmp_ord( n_n ) = i_ord
          tmp_dist( n_n ) = dist
          tmp_vec(:, n_n ) = t_vec(:)
          
       END IF
      
    end if


  END SUBROUTINE  nearest_test



  !===========================================================================
  !
  ! Subroutine "sort_near_by_ind". 
  ! Resort the nearest_neighbours with increasing atom indeces. 
  ! It must be done for correct computation of the sparse Hamiltonian
  !
  !===========================================================================
  SUBROUTINE sort_near_by_ind( near )

    !-------------------------------------------------------------------------
    ! Input -output argument
    !-------------------------------------------------------------------------

    TYPE (nearest_neighbours), INTENT( INOUT ) :: near

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------

    REAL( dp ), DIMENSION( SIZE( near%ind ) )    :: dist
    INTEGER,    DIMENSION( 3, SIZE( near%ind ) ) :: vec
    INTEGER,    DIMENSION( SIZE( near%ind ) )    :: index
    INTEGER,    DIMENSION( SIZE( near%ind ) )    :: order
    LOGICAL,    DIMENSION( SIZE( near%ind ) )    :: ind_test

    INTEGER :: mini, maxi, i_sort, i_coord, n_test

    !-------------------------------------------------------------------------
    ! Initialization
    !-------------------------------------------------------------------------

    index = near%ind
    order = near%order  
    dist  = near%dist
    vec   = near%vec

    maxi = MAXVAL( index )

    i_sort = 0

    !-------------------------------------------------------------------------
    ! Loop on sorted values
    !-------------------------------------------------------------------------

    DO WHILE ( i_sort .LT. SIZE( near%ind ) )

       !----------------------------------------------------------------------
       ! Get next index
       !----------------------------------------------------------------------

       mini = MINVAL( index )

       ind_test = ( index .EQ. mini )

       n_test = COUNT( ind_test )

       !----------------------------------------------------------------------
       ! Store current index's values
       !----------------------------------------------------------------------

       near%ind( i_sort + 1 : i_sort + n_test ) = PACK( index, ind_test )
       near%order( i_sort + 1 : i_sort + n_test ) = PACK( order, ind_test )
       near%dist( i_sort + 1 : i_sort + n_test ) = PACK( dist, ind_test )

       DO i_coord = 1, 3

          near%vec( i_coord, i_sort + 1 : i_sort + n_test ) = &
               PACK( vec( i_coord, : ), ind_test )

       END DO

       !----------------------------------------------------------------------
       ! Update counters
       !----------------------------------------------------------------------

       i_sort = i_sort + n_test

       WHERE( ind_test ) index = maxi + 1

       !----------------------------------------------------------------------
       ! End loop on sorted values
       !----------------------------------------------------------------------

    END DO

    !=========================================================================

  END SUBROUTINE  sort_near_by_ind


  !===========================================================================
  !
  ! Subroutine "sort_near_by_dist"
  !
  ! PROBABLY NOT WORKING !!! CHECK
  !
  !===========================================================================
  SUBROUTINE sort_near_by_dist(near)
    
    TYPE (nearest_neighbours), INTENT( INOUT ) :: near

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------
    REAL    ( dp ),  DIMENSION( : ),   ALLOCATABLE :: dist_aux
    INTEGER,         DIMENSION( :,:),  ALLOCATABLE :: vec_aux
    INTEGER,         DIMENSION( : ),   ALLOCATABLE :: ind_aux
    LOGICAL,         DIMENSION( : ),   ALLOCATABLE :: dist_test

    REAL ( dp ) :: mini, maxi

    INTEGER :: i_sort, n_test, i_coord, i_test, order, i, n_dist, n_order

    !-------------------------------------------------------------------------
    ! Initialize auxiliary arrays
    !-------------------------------------------------------------------------
    n_dist = size(near%ind)
    allocate( dist_aux(n_dist))
    allocate( vec_aux(3, n_dist ))
    allocate( ind_aux(n_dist))
    allocate( dist_test(n_dist));

    dist_aux = near%dist
    ind_aux = near%ind
    vec_aux = near%vec

    !-------------------------------------------------------------------------
    ! Initialize counters
    !-------------------------------------------------------------------------

    maxi = MAXVAL( dist_aux )

    i_sort = 0
    order = -1
    n_order = -1

    !-------------------------------------------------------------------------
    ! Loop on sorted values
    !-------------------------------------------------------------------------

    DO WHILE ( i_sort .LT. n_dist )

       !----------------------------------------------------------------------
       ! Get next order of nearest distances
       !----------------------------------------------------------------------

       mini = MINVAL( dist_aux )

       DO i_test = 1, SIZE( dist_test )

          dist_test( i_test ) = equiv( dist_aux( i_test ), mini, accur, fuzzy )

       END DO

       n_test = COUNT( dist_test )

       !----------------------------------------------------------------------
       ! Store current order data
       !----------------------------------------------------------------------

       near%dist( i_sort + 1 : i_sort + n_test ) = PACK( dist_aux, dist_test )
       near%ind( i_sort + 1 : i_sort + n_test ) = PACK( ind_aux, dist_test )

       DO i_coord = 1, 3

          near%vec( i_coord, i_sort + 1 : i_sort + n_test ) = &
                                   PACK( vec_aux( i_coord, : ), dist_test )

       END DO

       near%order( i_sort + 1 : i_sort + n_test ) = order + 1

       !----------------------------------------------------------------------
       ! Update counters
       !----------------------------------------------------------------------

       i_sort = i_sort + n_test
       order = order + 1


       !-------------------------------------------------------------------
       !new code 17 Feb 2006 M.Povolotskyi
       DO i =1, n_dist
          IF (dist_test(i)) THEN
             dist_aux(i) = maxi + 1.0D0
          END IF
       END DO
       !-------------------------------------------------------------------

    END DO
    deallocate(dist_aux,vec_aux,ind_aux,dist_test)
    !=========================================================================

  END SUBROUTINE  sort_near_by_dist


  !===========================================================================
  !
  ! Function "is_dangling_bond( i_n )"
  !
  !===========================================================================

  FUNCTION is_dg_bond( i_n, basis, current_near )

    !=========================================================================

    ! input arguments :

    INTEGER, INTENT( IN ) :: i_n
    TYPE ( ion_basis ) :: basis    
    TYPE (nearest_neighbours) :: current_near
    !_________________________________________________________________________

    ! output result :

    LOGICAL :: is_dg_bond
    CHARACTER ( LEN = 500 ) :: write_format
    !=========================================================================

    is_dg_bond = .FALSE.

    IF ( current_near%order( i_n ) .NE. 1 ) THEN

       is_dg_bond = .FALSE.
       RETURN

    END IF

    !IF (struct_building_flag .OR. load_precon_flag) THEN

       IF (current_near%ind(i_n) .GT. basis%n_basis) THEN
          is_dg_bond = .TRUE.
          RETURN
       ELSE
          is_dg_bond = .FALSE.
          RETURN
       END IF

    !END IF
    !=========================================================================

    SELECT CASE ( basis%periodic_BC )

       !______________________________________________________________________

    CASE ( 'yz' )

       IF ( current_near%vec( 1, i_n ) .NE. 0 ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'xz' )

       IF ( current_near%vec( 2, i_n ) .NE. 0 ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'xy' )

       IF ( current_near%vec( 3, i_n  ) .NE. 0 ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'z' )

       IF ( ANY( current_near%vec( 1 : 2, i_n ) .NE. 0 ) ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'x' )

       IF ( ANY( current_near%vec( 2 : 3, i_n ) .NE. 0 ) ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'y' )

       IF ( ( current_near%vec( 1, i_n ) .NE. 0 ) .OR. &
            ( current_near%vec( 3, i_n ) .NE. 0 ) ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'no' )

       IF ( ANY( current_near%vec( :, i_n ) .NE. 0 ) ) is_dg_bond = .TRUE.

       !______________________________________________________________________

    CASE ( 'xyz' ) 

       is_dg_bond = .FALSE.

       !______________________________________________________________________

    CASE DEFAULT

       write_format = ' '
       write_format = '( /, "DATA ERROR in file ''", a, "'' :", /,&
            & "wrong type in list ''", a, "''.", / )'

       WRITE ( *, TRIM( write_format ) ) &
            "data_driver", "periodic boundary conditions"

       call throw_init_exception(ERR_NN_LIST)
       !______________________________________________________________________

    END SELECT

    !=========================================================================

  END FUNCTION is_dg_bond

  !==============================================================================

  FUNCTION already_in( nn, k, ind, vec, vec_t)
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: vec
    INTEGER, DIMENSION(3) :: vec_t
    INTEGER :: k, nn
    
    LOGICAL :: already_in
    ! ---------------------------------------------------------------------
    INTEGER :: i
    
    already_in = .false.
    
    DO i = 1, nn-1
       
       IF( ind(i).EQ.k .AND. ALL(vec(:,i).EQ.vec_t(:),DIM=1) ) THEN
          
          already_in = .true.
          
       END IF
       
    END DO
    
  END function already_in
  !===========================================================================
  
  SUBROUTINE write_neighbours_map(nn_map)

    TYPE(nearest_neighbours), DIMENSION(:), intent(IN) :: nn_map

    INTEGER :: i_n, i_a, i_b


    DO i_a = 1, SIZE(nn_map)
   
       write(*,*) '--------------------------------------------'
       write(*,*) 'ATOM #',i_a
       !write(*,*) 'ncpl=',current_near%n_cpl
 
       DO i_n = 1, SIZE( nn_map(i_a)%ind )
          
          i_b = nn_map(i_a)%ind( i_n )
          
          write(*,'(i7,a,f10.5,a,i2)') i_b,'   d=', &
          nn_map(i_a)%dist( i_n ), '  o=',nn_map(i_a)%order( i_n )
          write(*,'(10x,a,3(i3),a)') 'v=(',nn_map(i_a)%vec( :, i_n ),')'

       END DO

    END DO
    write(*,*) '============================================='

  END SUBROUTINE write_neighbours_map

  !===========================================================================
  SUBROUTINE check_input_nn_list(str)

    TYPE( TStructure), INTENT( IN ) :: str 

    INTEGER, DIMENSION(:,:), POINTER     :: nn_list

    INTEGER :: ia, ib, j, k, n_n_ia, n_n_ib
    LOGICAL :: test
    
    nn_list => str%nn_list


    do ia = 1, str%n_atoms

      n_n_ia = COUNT( nn_list(:, ia) .GT. 0  ) 

      do j = 1, n_n_ia

         test = .false.
         ib = nn_list(j,ia)
         n_n_ib = COUNT( nn_list(:, ib) .GT. 0  )

         do k = 1, n_n_ib 
            
            if(nn_list(k,ib).eq.ia) then
               test = .true.
               exit
            end if

         end do
         
         if (.not.test) then
            write(*,*) 'ERROR: neighbour map not symmetric:',ia,ib
            call throw_init_exception(ERR_NN_LIST)
         end if

      end do

    end do

  END SUBROUTINE check_input_nn_list
  !=========================================================================== 
  
  SUBROUTINE check_nn_map(nn_map)

    TYPE(nearest_neighbours), DIMENSION(:) :: nn_map  
 
    INTEGER :: it, j, k, n_n_ia, n_n_ib, ia, ib    
    LOGICAL :: test
    
    ! Loop on all atoms
    DO ia = 1, SIZE(nn_map)

       ! Loop on nearest_neighbours of each atom
       DO j = 1, SIZE( nn_map(ia)%ind(:) )
          
          ib = nn_map(ia)%ind( j )

          ! check if symmetric neighbour is present
          test = .false.
          DO k = 1, SIZE( nn_map(ib)%ind(:) )  
             
             if( nn_map(ib)%ind(k) .eq. ia ) then
                test = .true.
                exit
             end if
             
          END DO

          if (.not.test) then
            write(*,*) 'ERROR: neighbour map not symmetric:',ia,ib
            call throw_init_exception(ERR_NN_LIST)
         end if

       END DO

    END DO


  END SUBROUTINE check_nn_map

  !=========================================================================== 

  SUBROUTINE destroy_neighbours_map( nn_map )

    TYPE(nearest_neighbours), DIMENSION(:), POINTER :: nn_map  
    INTEGER :: ia

    DO ia = 1, SIZE(nn_map) 
      call delete_near(nn_map(ia))
    END DO
    DEALLOCATE(nn_map)

  END SUBROUTINE destroy_neighbours_map

  !===========================================================================
  SUBROUTINE delete_near( near )
    TYPE (nearest_neighbours) :: near

    DEALLOCATE( near%ind )
    DEALLOCATE( near%order )
    DEALLOCATE( near%dist )
    DEALLOCATE( near%vec )    
    
  END SUBROUTINE delete_near
  !===========================================================================
  

  ! =========================================================================
  SUBROUTINE reallocate_i(array,n)
    
    INTEGER, DIMENSION(:), POINTER :: array
    INTEGER :: n

    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER :: oldsize, err

    oldsize = SIZE(array)
    ALLOCATE( tmp(oldsize), STAT=err )
    tmp(1:oldsize) = array(1:oldsize)
    DEALLOCATE(array)
    ALLOCATE(array(n), STAT=err)
    array(1:oldsize) = tmp(1:oldsize)
    DEALLOCATE(tmp)       

    IF(err.NE.0) CALL alloc_error('nearest_neighbours','make_nn_list','reallocate')

  END SUBROUTINE reallocate_i 

  ! =========================================================================
  SUBROUTINE reallocate_d(array,n)
    
    REAL(dp), DIMENSION(:), POINTER :: array
    INTEGER :: n

    REAL(dp), DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER :: oldsize, err

    oldsize = SIZE(array)
    ALLOCATE( tmp(oldsize), STAT=err )
    tmp = array
    DEALLOCATE(array)
    ALLOCATE(array(n), STAT=err)
    array(1:oldsize) = tmp(1:oldsize)
    DEALLOCATE(tmp)       

    IF(err.NE.0) CALL alloc_error('nearest_neighbours','make_nn_list','reallocate')

  END SUBROUTINE reallocate_d
  
  !=========================================================================

  SUBROUTINE reallocate_i2(array,m,n)

    INTEGER, DIMENSION(:,:), POINTER :: array
    INTEGER :: m,n

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmp
    INTEGER :: oldsize1,oldsize2, err

    oldsize1 = SIZE(array, DIM=1)
    oldsize2 = SIZE(array, DIM=2)
    ALLOCATE( tmp(oldsize1,oldsize2), STAT=err )
    tmp(:, :) = array(:, :)
    DEALLOCATE(array)
    ALLOCATE(array(m,n), STAT=err)
    array(1:oldsize1, 1:oldsize2) = tmp(1:oldsize1, 1:oldsize2)
    DEALLOCATE(tmp)   

    IF(err.NE.0) CALL alloc_error('nearest_neighbours','make_nn_list','reallocate2')

  END SUBROUTINE reallocate_i2


  SUBROUTINE write_am_coo(basis,nn_map)
    TYPE(ion_basis) :: basis
    TYPE(nearest_neighbours), DIMENSION(:), INTENT(IN) :: nn_map

    INTEGER :: i,k
    
    open(1001, file = 'AM_COO.dat')

    DO i = 1, basis%n_basis + basis%n_dg_bond    
       
       do k = 1, SIZE(nn_map(i)%ind)
  
          write(1001,*) i, nn_map(i)%ind(k), 1

       end do

    END DO

    close(1001)

  END SUBROUTINE write_am_coo


END MODULE neighbours
