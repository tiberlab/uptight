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
!              Module "TB_ham" - (c) Jerome Gleize - 2002
!                                (c) Alessandro Pecchia - 2005-2014 
!
!=============================================================================
! 
!
!=============================================================================

MODULE TB_ham

  !===========================================================================
  USE mpi_globals
  USE precision
  USE constants
  USE globals
  USE allocation, only : memstr
  USE type_defs
  USE list_types, only : matrix
  USE checks
  USE latt_fun
  USE neighbours, only : nearest_neighbours, is_dg_bond
  USE upt_param
  USE sparse_matrix
  USE sparse_numrec
  USE errors
  USE exceptions
  USE input_output
  USE crystal_field_d
  USE sort
  USE magnetic_gauge, only : calculate_peierls_phase
  USE globals, only : magnetic_field_vector, gauge_choice ! Added for Peierls
  !===========================================================================

  IMPLICIT NONE
  PRIVATE

  !===========================================================================

  PUBLIC sparse_ham
  PUBLIC check_if_hermitian
  PUBLIC check_if_antihermitian
  PUBLIC check_if_real
  PUBLIC hermitianize
  PUBLIC get_ion_block_size
  PUBLIC get_ion_orbitals

  PRIVATE koster_slater
  PRIVATE spin_sparse
  PRIVATE rotate_d, rotate_p, crystal_field_zb
  PRIVATE hydrogen_coupling
  PRIVATE store_ind
  PRIVATE store_val

  INTERFACE store_ind
     MODULE PROCEDURE store_ind_file
     MODULE PROCEDURE store_ind_mem     
  END INTERFACE

  INTERFACE store_val
     MODULE PROCEDURE store_val_file
     MODULE PROCEDURE store_val_mem     
  END INTERFACE


  INTERFACE hermitianize
     MODULE PROCEDURE hermitianize_csr
     MODULE PROCEDURE hermitianize_ex  
  END INTERFACE

  INTERFACE check_if_hermitian
     MODULE PROCEDURE check_if_hermitian_csr
     MODULE PROCEDURE check_if_hermitian_ex  
  END INTERFACE

  !===========================================================================
CONTAINS

  !===========================================================================
  !
  ! Subroutine sparse_ham
  !
  !===========================================================================
  !
  ! Compute the ETB Hamiltonian matrix in a sparse format
  !
  ! The format is the following :
  !
  ! --> sparse index vector :
  !
  !     - row_index / column index ...
  !
  !     i.e. the index of a non zero element are stored row by row,
  !          the row index being a negative value.
  !
  ! --> sparse values vector :
  !
  !     they are stored according to the column index sequence
  !
  ! This format can be converted to Compressed Sparse Row format.
  !
  !===========================================================================
  
  SUBROUTINE sparse_ham(upt)

    TYPE(OUPT), TARGET :: upt

    !=========================================================================
    !
    ! Local variables
    !    
    !=========================================================================

    !=========================================================================
    !
    ! Hamiltonian blocs :
    !
    ! --> cpl : matrix structure containing all blocs describing the coupling
    !           of an atom in the super basis to its nearest neighbours.
    !           ( including itself = onsite bloc )
    !
    ! --> cpl_ind : array containing the index of the nearest neighbours
    !               of the current atom in the super basis.
    ! 
    ! NOTE : These two structures are allocated for each atom in the basis
    !        according to the number of its nearest neighbours to consider.
    !
    !.........................................................................
    !
    ! EXAMPLE : 4 couplings + 1 onsite / 4 atoms
    !
    !   atom 1 is coupled to atoms 1 ( onsite ), 2, and 4 ( three times )
    ! 
    !   => cpl_ind =     1       2         4           4           4          
    !  
    !      [cpl]   =  [ 1/1 ] [ 1/2 ] [ 1/4 (1)]  [ 1/4 (2) ]  [ 1/4 (3) ]  
    !
    ! NOTE : couplings of an atom to equivalent atoms in the super basis 
    !        (-> they form a "star" ) are treated separately for the 
    !        calcultion. In the final Hamiltonian, however, they must
    !        be added to a single block :
    !
    !      column :      1       2      3    4
    !
    ! row 1 :        [ 1/1 ]  [ 1/2 ]   0   [ 1/4 ]
    !
    ! with : [ 1/4 ] =  [ 1/4 (1)]  +  [ 1/4 (2) ]  +  [ 1/4 (3) ]
    !
    !.........................................................................
    !
    ! --> hydro : matrix structure containing all blocs describing 
    !             saturation of the dangling bonds of an atom 
    !             in the super basis.
    !
    ! NOTE : This structure is allocated for each atom in the basis
    !        according to the number of its dangling bonds.
    !
    !.........................................................................
    !
    ! --> so_diag    : diagonal spin-orbit bloc matrix ( up-up or down-down )
    ! --> so_offdiag : non diagonal spin-orbit matrix ( up-down or down-up )
    ! --> tb_bloc    : Hamiltonian bloc matrix relevant to the coupling
    !                  between two atoms in the super basis.
    !
    ! NOTE : the dimension of these blocs is fixed to the maximum numbers
    !        of atomic orbitals which may be included on one atom.
    !        These are temporary data which contain the results
    !        of the Slater-Koster and spin-orbit routines.
    !
    ! The link between the maximum atomic orbitals basis and the one used
    ! in the simulation is achieved by an index vector :
    !
    ! --> a_ref : index of the row states in the reference list.
    ! --> b_ref : index of the column states in the reference list.
    !
    !=========================================================================

    INTEGER,       DIMENSION( 20) :: cpl_ind     ! Hard-coded maximum number
    TYPE (matrix), DIMENSION( 20) :: cpl, hydro  ! of couplings
    TYPE (matrix)                 :: h_diag      ! Matthias: to work on the onsite block

    TYPE(CSB)                     :: M, UM
    TYPE(CSR)                     :: TM

    COMPLEX ( dp ), DIMENSION(:,:), ALLOCATABLE :: so_diag
    COMPLEX ( dp ), DIMENSION(:,:), ALLOCATABLE :: so_off_diag
    COMPLEX ( dp ), DIMENSION(:,:), ALLOCATABLE :: U_bl
    REAL    ( dp ), DIMENSION(:,:), ALLOCATABLE :: tb_bloc
    REAL    ( dp ), DIMENSION(:,:), ALLOCATABLE :: intra_cpl, Q_bloc, quad_corr

    INTEGER,       DIMENSION(:), ALLOCATABLE :: a_ref
    INTEGER,       DIMENSION(:), ALLOCATABLE :: b_ref
    INTEGER,       DIMENSION(n_ref_states)   :: inv_ref
    !=========================================================================
    ! 
    ! Local variables from OUTP :
    !
    
    LOGICAL :: d_onsite_shift_flag
    LOGICAL :: potential_flag
    LOGICAL :: relat
    LOGICAL :: fuzzy, scale
    LOGICAL :: syst_rotated
    !LOGICAL :: alloy_random 
    LOGICAL :: ioutput_flag
    LOGICAL :: optmat
    TYPE (ion_basis), POINTER :: basis     
    TYPE (material_data), POINTER :: interface_data, current_mat
    TYPE (material_data), DIMENSION( : ), POINTER :: mat_data
    TYPE (material_data), DIMENSION( :,: ), POINTER :: int_data
    INTEGER, DIMENSION( : ), POINTER :: n_st
    INTEGER, DIMENSION( : ), POINTER :: n_dg
    REAL ( dp ), DIMENSION( : ), POINTER :: pot_data
    REAL ( dp ) :: d_H   ! scaling of interactions for empirical Hydrogens
    REAL ( dp ) :: E_H   ! scaling of on-site for empirical Hydrogens
    REAL ( dp ) :: dangling_shift
    INTEGER :: verbose   ! verobsity level
    INTEGER :: n_basis   !natoms 
    INTEGER :: n_spin    !spin variables
    INTEGER :: n_mat     !N. of materials
    INTEGER :: n_ham     !Hamiltonian size
    INTEGER :: poldir    !Polarization direction
    REAL ( dp ), DIMENSION( 3 )  :: veck    ! k-vector 
    REAL ( dp ), DIMENSION( 3 )  :: c_axis  ! c_axis for wurtzites
    REAL ( dp ), DIMENSION(3,3)  :: rec_latt
    !CHARACTER(5), DIMENSION(:), POINTER :: ref_couplings
    CHARACTER(1) :: sparse_format
    INTEGER(8) :: mem
    INTEGER :: ii

    !=========================================================================

    !COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE :: val_tmp
    !INTEGER,        DIMENSION(:), ALLOCATABLE :: ind_tmp

    !COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: values
    !INTEGER    , DIMENSION(:), ALLOCATABLE :: indeces
    !=========================================================================
    !
    ! ETB data :
    !
    ! --> onsite  : onsite energies for the current atom in the super basis.
    ! --> cpl_tmp : coupling energies for the current pair in the super basis.
    ! --> tb_pow  : coupling scaling exponents for the current pair.
    !
    ! NOTE : these data are temporary and correspond to the user input.
    !        They must be converted to the larger set of coupling elements
    !        corresponding to the reference list.
    !
    ! --> tb_cpl  : maximum coupling energies list for the current pair.
    !               It will be used in the Slater Koster routine.
    !
    !=========================================================================

    REAL ( dp ), DIMENSION(:), ALLOCATABLE   :: onsite, temp_onsite, ons_corr
    CHARACTER(5), DIMENSION(:), ALLOCATABLE  :: state
    REAL ( dp ), DIMENSION(:), ALLOCATABLE  :: tb_cpl, P_list, S_list, Q_list

    !=========================================================================
    ! Other local variables
    !=========================================================================
    TYPE(nearest_neighbours), POINTER :: current_near, near_near
    TYPE(parent_data), DIMENSION(:), POINTER :: p_parent
    TYPE(ion_orbit), DIMENSION(:), POINTER :: p_ion
    TYPE(pair_coupling), DIMENSION(:), POINTER :: p_pair
    INTEGER, DIMENSION( : ), POINTER :: ref_pair

    COMPLEX ( dp ) :: exp_fac, value

    REAL ( dp )    :: cutoff
    REAL ( dp ), DIMENSION( 3 ) :: D_bas, D_latt, D_tmp
    REAL ( dp ), DIMENSION( 3 ) :: cos_latt, u
    REAL ( dp ), DIMENSION( 2 ) :: so_cpl, so_corr
    REAL ( dp ), DIMENSION( 3 ) :: p_onsite ! p orbital onsite energies
    REAL ( dp ), DIMENSION( 5 ) :: d_onsite ! d orbital onsite energies
    REAL ( dp ), DIMENSION(:), ALLOCATABLE   :: dist_1nn_list, dist_2nn_list
    REAL ( dp ), DIMENSION(:), ALLOCATABLE   :: dd_rel_1nn_list, dd_rel_2nn_list
    REAL ( dp ), DIMENSION(:,:), ALLOCATABLE :: unit_disp_1nn_list, unit_disp_2nn_list

    REAL ( dp ) :: scaling_factor, scaled_factor_1, scaled_factor_2, phase

    INTEGER :: i_read, i_ind, n_ind, i_val, n_val, n_val_U, n_col, n_row, nnz
    INTEGER :: i_mat, i_a_mat, i_b_mat, i_ab_mat, i_ba_mat
    INTEGER :: i_a_ion, i_b_ion, n_pair, i_a_in_parent, i_ion_parent
    INTEGER :: i_pair, i_parent, j_parent, n_parent, n_neig, n_1nn
    INTEGER :: i_cpl, i_onsite, n_cpl, n_so, i_c, i_dg, i_dg_a
    INTEGER :: n_n, i_n, dec
    INTEGER :: i_a, i_b, i_b_old, n_st_a, n_st_b, sp_a, sp_b
    INTEGER :: i_spin_a, i_spin_b, i_state
    INTEGER :: a_st, b_st, ab_order
    INTEGER :: index_num, file_num, file_num2, ham_num, err
    INTEGER :: n_sparse, n_value, n_zeros_on, n_zeros_off, n_zeros_h
    INTEGER :: p_start, p_stop, offset
    INTEGER :: i_1nn_list, i_2nn_list, i_1n, i_2n, i_dump, n_orb_b, len_2nn_list
    REAL(sp) :: estimate_factor

    CHARACTER ( LEN = 200 ) :: file_name, file_name2
    CHARACTER ( LEN = 500 ) :: write_format
    CHARACTER ( LEN = 200 ) :: list 
    CHARACTER (3) :: str
    CHARACTER(ATOMLEN) :: name_a, name_b
    INTEGER, DIMENSION(:), ALLOCATABLE :: bl_arr, n_2nn_list

    ! the starting colum indices for each block
    INTEGER, DIMENSION(:), ALLOCATABLE :: col_start

    LOGICAL :: value_test, dg_bond, tmp_scale, ons_mix
    COMPLEX(dp) :: peierls_phase_factor ! Added for Peierls
    CHARACTER ( LEN = 10 ) :: scheme !the ETB scheme ('tan','jancu')
    CHARACTER ( LEN = 5 ), DIMENSION(n_intra) :: name_intra
    REAL(dp) :: dist_phase
    REAL(dp), DIMENSION(n_intra)   :: fac_C

    !=========================================================================
    ! Initializations
    !=========================================================================
    !write(*,*) '================================================='
    !write(*,*) 'Build Hamiltonian'
    !write(*,*) '================================================='

    verbose = upt%verbose

    cutoff = emach * 3.7D04

    !switch on/off onsite-mixing improvement for Jancu scheme
    ons_mix = .true.

    d_onsite_shift_flag = upt%d_onsite_shift_flag
    potential_flag = upt%potential_flag
    relat = upt%relat
    fuzzy = .true.
    scale = upt%scaling
    syst_rotated = upt%syst_rotated
    c_axis = upt%c_axis
    ioutput_flag = upt%ioutput_flag
    optmat = upt%optmat
    sparse_format = upt%sparse_format

    d_H = upt%d_H 
    E_H = upt%E_H 

    n_spin = upt%n_spin
    n_mat = upt%nr_mat
    mat_data => upt%materials
    int_data => upt%interfaces

    pot_data => upt%pot_data

    basis => upt%basis

    n_basis = upt%basis%n_basis

    !ref_couplings => upt%ref_couplings

    n_st => upt%basis%n_st
    n_dg => upt%basis%n_dg

    n_ham = n_spin * ( SUM( n_st(1:n_basis) ) )
    
    poldir = upt%poldir

    veck = upt%k_point
    rec_latt = upt%basis%rec_latt

    estimate_factor = 4.0*upt%estimate_factor


    i_val = 0
    n_val = 0        ! total number of sparse values    
    n_val_U = 0      ! total number of sparse values for U   
    n_row = 0        ! row counter
    n_col = 0        ! column counter
    n_zeros_on = 0   ! accidental zeors counter
    n_zeros_off = 0  ! accidental zeors counter
    n_zeros_h = 0    ! accidental zeors counter

    n_ind = 0        ! total number of sparse indexes    
    i_ind = 0        ! temporary counter of sparse indexes
    i_dg = 0       ! total dangling bond counter

    !=========================================================================
    ! ALLOCATE local arrays
    !=========================================================================
    ALLOCATE(so_diag(upt%n_ref_st, upt%n_ref_st), STAT = err)    
    ALLOCATE(so_off_diag(upt%n_ref_st, upt%n_ref_st), STAT = err) 
    ALLOCATE(tb_bloc(upt%n_ref_st, upt%n_ref_st), STAT = err)
    ALLOCATE(U_bl(upt%n_ref_st, upt%n_ref_st), STAT = err)     

    ALLOCATE(a_ref(upt%n_ref_st), STAT = err)    
    ALLOCATE(b_ref(upt%n_ref_st), STAT = err)    
    ALLOCATE(onsite(upt%n_ref_st), STAT = err)
    ALLOCATE(temp_onsite(upt%n_ref_st), STAT = err)
    ALLOCATE(ons_corr(upt%n_ref_st), STAT = err)
    ALLOCATE(tb_cpl(upt%n_ref_cpl), STAT = err)

    ALLOCATE(bl_arr(n_spin), STAT=err)

    so_diag = (0.d0,0.d0)
    so_off_diag = (0.d0,0.d0)
    U_bl = (0.d0,0.d0)

    do i_a = 1, upt%n_ref_st
       U_bl(i_a,i_a) = (1.d0, 0.d0)
    end do

    ALLOCATE( ref_pair(upt%n_ref_cpl), STAT = err )

    ! Initialize extra arrays for for Tan scheme
    ALLOCATE(P_list(upt%n_ref_cpl), STAT = err)
    ALLOCATE(S_list(upt%n_ref_cpl), STAT = err)
    ALLOCATE(Q_list(upt%n_ref_cpl), STAT = err)
    ALLOCATE(intra_cpl(upt%n_ref_st, upt%n_ref_st), STAT = err)
    ALLOCATE(Q_bloc(upt%n_ref_st, upt%n_ref_st), STAT = err)
    ALLOCATE(quad_corr(upt%n_ref_st, upt%n_ref_st), STAT = err)

    IF(err.ne.0) call alloc_error( 'TB_ham', 'TB_ham', 'local arrays' )
    
    ref_pair = 0

    if (id0 .and. verbose.gt.0) then
      select case(sparse_format)
      case('U')     
        write(*,*) '(TB ham) Build Upper Hamiltonian'
      case('L')     
        write(*,*) '(TB ham) Build Lower Hamiltonian'
      case('F')     
        write(*,*) '(TB ham) Build Full Hamiltonian'
      end select
    end if




    !**************************************************************************
    !==========================================================================
    !       Loop on TB Hamiltonian row part 1 : basis atoms a.
    !==========================================================================   


    !**************************************************************************
    ! Parallelization on atom loop:
    !**************************************************************************
    p_start = id*n_basis/num_procs + 1
    offset = n_spin * SUM( n_st( 1 : p_start - 1 ) )
    n_row = offset

    if (id.eq.num_procs-1) then
      p_stop = n_basis
    else    
      p_stop = (id+1)*n_basis/num_procs 
    endif


    !---------------------------------------------------------------------- 
    ! First calculate the starting column index of each block
    !----------------------------------------------------------------------
    ALLOCATE(col_start(n_spin*n_basis), STAT=err)
    col_start(1) = 1
    n_st_a = n_st(1)
    if (relat) then
      col_start(2) = 1 + n_st_a
    end if

    i_ind = n_st_a
    DO i_a = 2, n_basis
      i_n = (i_a-1)*n_spin + 1
      col_start(i_n) = col_start(i_n-1) + n_st_a
      n_st_a = n_st(i_a)
      if (i_ind .lt. n_st_a) then
        i_ind = n_st_a
      end if
      if (relat) then  
        col_start(i_n+1) = col_start(i_n) + n_st_a
      end if
    END DO
    
    ! in i_ind we have now temporarily the maximum number of states per atom

    ! the matrix dimension
    i_n = n_spin*SUM(n_st(1:n_basis))

  
    ! Estimate nnz (4 as typical number of neighbours +1) 
    !if (relat) then
    !   n_zeros_h =  2.0 * estimate_factor * (p_stop-p_start+1) * (n_spin * i_ind)**2
    !else
       n_zeros_h =  2 * estimate_factor * (p_stop-p_start+1) * (n_spin * i_ind)**2
    !endif

    if(id0 .and. verbose.gt.0) write(*,*) '(TB_ham) CPU '//id_s//' est non zero values= ',n_zeros_h


    call create_matrix(upt%ham, i_n, i_n, n_zeros_h) 
    upt%ham%sparse_fmt = sparse_format
    upt%ham%offset = offset
                         
    if (relat) then
       call create_matrix(upt%U, i_n, i_n, i_ind * (p_stop-p_start+1) * n_spin)
       upt%U%sparse_fmt = sparse_format
       upt%U%offset = offset
       if(id0 .and. verbose.gt.0) then
          write(*,*) '(TB_ham) CPU '//id_s//' est non zero values U= ',upt%U%nnz
       endif
    endif
    
 
    if(id0 .and. verbose.gt.0) then
      write(*,*) '(TB_ham) Relativistic:',upt%relat  
      write(*,*) '(TB_ham) Harrison Scaling:',upt%scaling
      write(*,*) '(TB_ham) d onsite shift:',upt%d_onsite_shift_flag
      write(*,*) '(TB_ham) hybrid passivation:',upt%hybrid_passivation
      write(*,*) '(TB_ham) External potential:',upt%potential_flag 
      write(*,'(a,3(f8.4))') ' (TB_ham) k-point:',upt%k_point
      write(*,*) '(TB_ham) d_H:',d_H
      write(*,*) '(TB_ham) E_H:',E_H
    endif



    basis_loop:DO i_a = p_start, p_stop 

       current_near => upt%nn_map(i_a)

       !----------------------------------------------------------------------
       ! Screen output of the atom a index
       !----------------------------------------------------------------------
       IF ( verbose.gt.0 .AND. MOD( i_a, 10000 ) .EQ. 0.0 ) THEN 
          WRITE ( *, * ) "(TB ham) atom a = ", i_a 
       END IF
       !----------------------------------------------------------------------
       !
       ! Get reference states index for the "a" atom :
       !
       ! This step is necessary to have a correct ordering of the states 
       ! during the calculation of matrix elements.
       !
       ! This is achieved through indentification of the "a" atom 
       ! in the ion data list "mat_data( : )%ion( : )"
       !
       ! It requires a call of the subroutine "check_ion" which returns 
       ! the material ( i_a_mat ) and ion ( i_a_ion ) index for atom "a" :
       !
       ! "a" atom <--> mat_data( i_a_mat )%ion( _ia_ion )
       !
       ! The ordered reference states index are then read from the field 
       ! "ind_ref" of the list.
       !       
       !----------------------------------------------------------------------
       ! Get the material index
       !----------------------------------------------------------------------
       name_a = basis%atomtypes(basis%type(i_a)) !LUAN
       i_a_mat = basis%mat(i_a)
       i_a_ion = basis%ion(i_a)
                 
       a_ref = 0
       p_ion => mat_data(i_a_mat)%ion
       scheme = TRIM(mat_data(i_a_mat)%scheme) !Tan
       n_st_a = n_st(i_a)

       IF ( relat ) THEN
         n_so = SIZE( p_ion(i_a_ion)%so_energy )
       END IF
       a_ref(1:n_st_a)= p_ion(i_a_ion)%ind_ref(1:n_st_a)
       inv_ref(1:n_ref_states) = p_ion(i_a_ion)%ind_ref_inv(1:n_ref_states)

       !----------------------------------------------------------------------
       ! Initialize n_neig (does not count Hydrogens,
       ! and each periodic copy once) 
       !----------------------------------------------------------------------
       n_neig = current_near%n_cpl

       n_1nn = 0 !number of actual 1NN of `a` (i.e. no dangling bonds)
       DO i_n = 1, SIZE(current_near%ind)
         IF ((current_near%order(i_n) .EQ. 1) .and. &
         .not.  is_dg_bond( i_n, basis, current_near ) ) THEN
           n_1nn = n_1nn + 1
         ENDIF
       ENDDO

       !----------------------------------------------------------------------
       ! Allocate the Hamiltonian block-matrices for atom a
       !----------------------------------------------------------------------
       bl_arr(:) = n_spin*n_neig
       call create_matrix(M, n_spin, n_basis, bl_arr)
       if (relat) call create_matrix(UM, n_spin, n_basis, bl_arr)

       DO i_cpl = 1, SIZE( current_near%ind )
          NULLIFY( cpl( i_cpl )%mat )
       END DO

       NULLIFY( h_diag%mat ) !Matthias
       ALLOCATE( h_diag%mat(n_st_a, n_st_a) )

       !write(*,*) i_a, current_near%n_cpl, SIZE( current_near%ind ) 

       !----------------------------------------------------------------------
       ! Initialization
       !----------------------------------------------------------------------
       i_b = 0         ! nearest neighbour index
       i_b_old = 0     ! previous nearest neighbour index
       i_cpl = 1       ! coupling block counter
       i_dg_a = 0      ! hydrogenated block counter
       onsite = 0.0d0  ! onsite
       h_diag%mat = 0.0! Matthias: diagonal block
       
       IF ( relat ) THEN
         so_cpl = 0.0D0  ! LUAN: SO
       END IF

       IF (scheme .EQ. 'tan') THEN
         P_list = 0.d0
         S_list = 0.d0
         Q_list = 0.d0
         Q_bloc = 0.d0
         quad_corr = 0.d0
         intra_cpl = 0.d0
       ENDIF
       


       !***********************************************************************
       !=======================================================================
       !        Loop on nearest neighbours atoms.
       !=======================================================================   
       !***********************************************************************

       ! For Tan scheme, perform preliminary loops over 1NN and 2NN of `a`
       ! to extract info that will be used later
       IF (scheme .EQ. 'tan') THEN

         ! Initiate 1NN_lists: these are about `b`, 1NN of `a`
         ALLOCATE( n_2nn_list(n_1nn), STAT = err) !list of numbers of 1NNs of `b`
         ALLOCATE( dist_1nn_list(n_1nn), STAT = err) !list of bondlengths a-b
         ALLOCATE( unit_disp_1nn_list(n_1nn,3), STAT = err) !list of unit displacement vectors from a->b
         ALLOCATE( dd_rel_1nn_list(n_1nn), STAT = err) !list of ratios (actual - average)/average bondlength
         !ALLOCATE( name_1nn_list(n_1nn), STAT = err) !list of atom names of `b`

         n_2nn_list = 0
         dist_1nn_list = 0.D0
         unit_disp_1nn_list = 0.D0
         dd_rel_1nn_list = 0.D0

         i_1nn_list = 0 !index reserved for scanning through 1NN_lists
         DO i_1n = 1, SIZE(current_near%ind) !scan through all neighbors of atom `a`
           IF (( current_near%order(i_1n) .EQ. 1) .and. &
                 .not.  is_dg_bond( i_1n, basis, current_near ) ) THEN
             
             i_1nn_list = i_1nn_list + 1
             dist_1nn_list(i_1nn_list) = current_near%dist(i_1n) !get a-b bond length
             D_latt = current_near%vec(:,i_1n) !lattice vector from a->b  
             i_b = current_near%ind(i_1n)
             unit_disp_1nn_list(i_1nn_list,1:3) = unit_displacement(basis, D_latt, i_a, i_b) !a->b displacement unit vec  
             ! name_1nn_list(i_1nn_list) = basis%atomtypes(basis%type(i_b)) !get atom type
  
             near_near => upt%nn_map(i_b) !Neighbor map of `b`

             DO i_2n = 1, SIZE(near_near%ind) !scan through all neighbors of `b`
               IF ( near_near%order(i_2n) .EQ. 1 .and. &
                    .not.  is_dg_bond( i_2n, basis, near_near ) ) THEN
                 n_2nn_list(i_1nn_list) = n_2nn_list(i_1nn_list) + 1
               ENDIF
             ENDDO

           ENDIF
    
         ENDDO
         ! Now n_2nn_list contains numbers of 1NN of `b`
         ! If a-b is a dangling bond, then the value is 0

         dd_rel_1nn_list(1:n_1nn) = dist_1nn_list(1:n_1nn)/(SUM(dist_1nn_list)/n_1nn) - 1.0D0

         len_2nn_list = SUM(n_2nn_list) !length of 2NN_lists (below)
         
         ! Initiate 2NN_lists: about 1NN of non-dangling `b`, similar to above 1NN_lists
         ! The order of the elements: (b1-c1, b1-c2,..., b2-d1, b2-d2,...)
         ! NOTE: one of 1NN of `b` coincides with `a`
         ALLOCATE(dist_2nn_list(len_2nn_list), STAT = err)
         ALLOCATE(unit_disp_2nn_list(len_2nn_list,3), STAT = err) !vector b->c;
         ALLOCATE(dd_rel_2nn_list(len_2nn_list), STAT = err)
         !ALLOCATE( name_2nn_list(len_2nn_list), STAT = err)

         dist_2nn_list = 0.0d0
         unit_disp_2nn_list = 0.0d0
         dd_rel_2nn_list = 0.d0

         i_2nn_list = 0 !reserved index to scan through the 2NN_lists
         ! Now scan once again
         DO i_1n = 1, SIZE(current_near%ind) !scan through NN of atom `a`
           IF ( current_near%order(i_1n) .EQ. 1 .AND. & !if 1NN ...
                .NOT. is_dg_bond( i_1n, basis, current_near ) ) THEN !and not dangling bond

             i_b = current_near%ind(i_1n)
             near_near => upt%nn_map(i_b) !Neighbor map of `b`
      
             DO i_2n = 1, SIZE(near_near%ind)
               
               IF ( near_near%order(i_2n) .EQ. 1 .and. &
                    .not.  is_dg_bond( i_2n, basis, near_near ) ) THEN

                 i_2nn_list = i_2nn_list + 1
                 i_c = near_near%ind(i_2n)
                 ! name_2nn_list(i_2nn_list) = basis%atomtypes(basis%type(i_c)) !get atom type
                 dist_2nn_list(i_2nn_list) = near_near%dist(i_2n) !get b-c bond length
                 D_latt = near_near%vec(:,i_2n) !lattice vector from b->c
                 unit_disp_2nn_list(i_2nn_list,1:3) = unit_displacement(basis, D_latt, i_b, i_c) !b->c displacement unit vec

               ENDIF
    
             ENDDO
  
           ENDIF
  
         ENDDO
  
         !IF (i_2nn_list > len_2nn_list) THEN !if more neighbors than expected, throw error then exit
         !  CALL throw_solve_exception(ERR_BOND_NUM)
         !ENDIF

         CALL get_dd_rel_2nn_list(n_2nn_list, dist_2nn_list, dd_rel_2nn_list)
         ! After the above quick scan, we know the necessary info about 
         ! all relative bond change and displacement direction upto 2NN shell

         ! Reset `i_1nn_list` and `i_2nn_list` to before the beginning of above lists
         i_1nn_list = 0
         i_2nn_list = 0

       ENDIF

   
       ! This is the actual loop over 1NN of `a` to build the coupling blocks
       nn_loop:DO i_n = 1, SIZE(current_near%ind) 
          
          !-------------------------------------------------------------------
          ! Get current neighbour's basis index
          !-------------------------------------------------------------------
                    
          i_b = current_near%ind( i_n )
          n_st_b = n_st(i_b)
          i_b_mat = basis%mat(i_b)
          i_b_ion = basis%ion(i_b)
          !-------------------------------------------------------------------
          ! If we jump to a new atom
          !-------------------------------------------------------------------
          IF ( i_b .NE. i_b_old ) THEN
             !----------------------------------------------------------------
             ! Increment coupling block counter if required
             !----------------------------------------------------------------
             IF ( ASSOCIATED( cpl( i_cpl )%mat ) ) i_cpl = i_cpl + 1
             !----------------------------------------------------------------
             ! Update previous element counter
             !----------------------------------------------------------------
             i_b_old = i_b
          END IF

          !-------------------------------------------------------------------
          !dangling bond test 
          !-------------------------------------------------------------------
          IF ( is_dg_bond( i_n, basis, current_near ) ) THEN  
             dg_bond = .true.
          ELSE
             dg_bond = .false.
          ENDIF

          !----------------------------------------------------------------
          ! Allocate work coupling array and update index
          !----------------------------------------------------------------
          IF ( .NOT. ASSOCIATED( cpl( i_cpl )%mat ) ) THEN
             
             ALLOCATE( cpl(i_cpl)%mat(n_st_a, n_st_b),  STAT = err )
             IF (err.NE.0) CALL alloc_error('TB_ham','sparse_ham','cpl%mat')
             
             cpl( i_cpl )%mat = (0.d0, 0.d0)
             
             cpl_ind( i_cpl ) = i_b
              
          END IF
          
          !----------------------------------------------------------------
          ! GET PAIR INFO AND COMPUTE PAIR-DEPENDENT PARAMETERS
          !----------------------------------------------------------------
          
          ! get the order of the current neighbours pair ( 0th, 1st, 2nd, ... )
          ab_order = current_near%order( i_n )
          
          IF (ab_order .NE. 0) THEN
             
             !----------------------------------------------------------------
             ! Get reference states index for atom b
             !----------------------------------------------------------------

             name_b = basis%atomtypes(basis%type(i_b))

             !-------------------------------------------------------------
             ! Calculate direction-cosinus'
             !-------------------------------------------------------------
             
             D_latt = current_near%vec(:,i_n)
             cos_latt = unit_displacement(basis, D_latt, i_a, i_b)
             
             !----------------------------------------------------------------
             ! Get pair-dependent ETB data for the a-b pair couplings
             !----------------------------------------------------------------
             ! IMPORTANT: we define here which set of ETB data are used.
             ! this is done by selecting a material index : i_ab_mat
             ! If the two atoms in the pair belong to the same material,
             ! we just have to read the relevant data.
             !----------------------------------------------------------------
             IF ( i_a_mat .EQ. i_b_mat ) THEN 
                
                !-------------------------------------------------------------
                ! Indentification of the "a" - "b" atom pair in the pair 
                ! data lists "mat_data( : )%parent%pair( : )"
                ! 
                ! pair "a" - "b" 
                ! <--> mat_data( i_ab_mat )%parent( i_parent )%pair( i_pair )
                !-------------------------------------------------------------
                interface_data => mat_data(i_a_mat)
                p_parent => interface_data%parent

                CALL check_pair(name_a, name_b, p_parent, ab_order, i_parent, i_pair)

                IF (i_parent .EQ. 0) THEN
                  print*,'atom',i_a,name_a,'in mat',i_a_mat
                  print*,'atom',i_b,name_b,'in mat',i_b_mat
                  call throw_solve_exception(ERR_HAM_UNPAIR)  
                ENDIF
                
             ELSE
                !-------------------------------------------------------------
                ! INTERFACE ATOMS. Works for nearest neighbours parameters
                ! at an interface like GaAs / AlAs with a common atom (As) 
                ! Simple case, since we can always use the parameters
                ! of one of the two materials ( as in a common atom alloy )
                ! It might be a problem when there is no common atom (Si/Ge)
                !
                ! POSSIBLE SOLUTION : 
                ! The couplings  between non-common atoms are taken from 
                ! a fictitious  parent of  SiGe, SiGe_cpl, so that all
                ! the relevant parameters are present for material "SiGe" : 
                ! cpl Si <-> Ge  in Si / SiGe interface.
                !   
                ! A better solution is to define specific interface models
                !-------------------------------------------------------------
                
                
                ! CALL check_interface_pair( basis, i_a, i_b, mat_data, ab_order, & 
                !                                   i_ab_mat, i_parent, i_pair )

                ! interface_data => mat_data(i_ab_mat)
                
                CALL check_interface_pair( basis, i_a, i_b, mat_data, ab_order, & 
                                                interface_data, i_parent, i_pair )
                
                IF (i_parent .EQ. 0) THEN 
                  CALL check_uncommon_interface_pair( basis, i_a, i_b, int_data, ab_order, &
                                                  interface_data, i_parent, i_pair )
                ENDIF
                

                IF (i_parent .EQ. 0) THEN
                   print*,'atom',i_a,name_a,' in material',i_a_mat
                   print*,'atom',i_b,name_b,' in material',i_b_mat
                   call throw_solve_exception(ERR_HAM_UNPAIR)
                ENDIF

                p_parent => interface_data%parent
                
             END IF
             
             ! Now we identified the pair and its corresponding material
             p_pair => p_parent(i_parent)%pair
             p_ion => p_parent(i_parent)%ion
             current_mat => interface_data

             IF ((scheme .EQ. 'jancu') .AND. (.NOT. ons_mix)) THEN !if ignore onsite-mixing then reset p_ion and current_mat
               p_ion => mat_data(i_a_mat)%ion
               current_mat => mat_data(i_a_mat)
             ENDIF
             
             IF (scheme .EQ. 'tan') THEN !calculate distance phase in Tan scheme
               dist_phase = current_near%dist(i_n) - p_pair(i_pair)%dist_ref &
                                                   + p_pair(i_pair)%delta_d
             ENDIF
             
             ! Calculate the weighted average onsite energies and SO energies
             ! and update intracouplings according to the 1NN neighbors
             IF ((ab_order .EQ. 1) .AND. (.NOT. dg_bond)) THEN ! only count 1NN and non-dangling
               
               ! Matthias: CHECK THIS PART --> This is needed because ion index for the same elements
               ! may be different in different materials
               ! Identify the ion-type index of `a` in the found parent material
               DO i_ion_parent = 1,SIZE(p_ion) 
                 IF (name_a .EQ. p_ion(i_ion_parent)%name) THEN 
                   i_a_in_parent = i_ion_parent
                   EXIT
                 ENDIF
                 !IF (name_b .EQ. p_ion(i_ion_parent)%name) THEN 
                 !  i_b_ion_in_parent = i_ion_parent
                 !ENDIF
               ENDDO
               
               !a_ref(1:n_st_a)= p_ion(i_a_in_parent)%ind_ref(1:n_st_a)
               !inv_ref(1:n_ref_states) = p_ion(i_a_in_parent)%ind_ref_inv(1:n_ref_states)
               !i_a_in_parent = i_a_ion

               temp_onsite(1:n_st_a) = p_ion(i_a_in_parent)%energy(1:n_st_a)
               
               ! Matthias: TODO check possible error in the following IF statement
               IF (scheme .EQ. 'jancu') THEN
                 IF (TRIM(p_ion( i_a_in_parent )%state( 1 )) .EQ. 's') THEN

                   ! onsite shift of s-orbitals to correct bandgap (Pecchia) 
                   temp_onsite(inv_ref(s)) = temp_onsite(inv_ref(s)) + &
                                             current_mat%gap_corr(1)
                 ENDIF
                 
                 ! Crystal-field splitting for strained zb
                 IF (TRIM(current_mat%cry) .EQ. 'zinc-blende' &  
                     .AND. d_onsite_shift_flag ) THEN
                    CALL crystal_field_zb(temp_onsite, p_ion, i_a_in_parent, c_axis, &
                         basis%strain(:,i_a), inv_ref)
                 ENDIF
                 
               ENDIF
               

               ! Accummulative sum for onsite
               ! NOTE: For Jancu scheme, this averaging in principle differs from the true onsite energies
               ! if there are dangling bonds in 1st neighbors of `a`, i.e. number of actual 1NN is less than the nominal one.
               ! This is unavoidable because we don't know the portions in which `a` and its 1NN contribute to its onsite energies
               onsite(1:n_st_a) = onsite(1:n_st_a) + &
                                  p_ion(i_a_in_parent)%offset/n_1nn + & !offset included Ev in .dat file already
                                  temp_onsite(1:n_st_a)/n_1nn

               ! Calculate the correction to onsite (including band offset) of `a`
               IF (scheme .EQ. 'tan') THEN
                 n_orb_b = SIZE(p_ion(i_a_in_parent)%deg)
                 ons_corr = 0.d0
                 CALL get_ons_corr(name_a, p_pair(i_pair), &
                                    p_ion(i_a_in_parent), dist_phase, &
                                    n_orb_b, ons_corr)

                 onsite(1:n_st_a) = onsite(1:n_st_a) + ons_corr(1:n_st_a)
               ENDIF
               
               ! Accummulative sum for SO
               IF ( relat ) THEN

                 so_cpl(1:n_so) = so_cpl(1:n_so) + &
                                  p_ion(i_a_in_parent)%so_energy(1:n_so)/n_1nn

                 IF (scheme .EQ. 'tan') THEN
                   so_corr = 0.d0 !Tan
                   CALL get_so_corr(name_a, p_pair(i_pair), n_so, so_corr)
                   so_cpl(1:n_so) = so_cpl(1:n_so) + so_corr(1:n_so)
                 ENDIF

               END IF
               
               ! Update intra-atomic couplings
               IF (scheme .EQ. 'tan') THEN
                 IF (name_a .NE. p_pair(i_pair)%name(1)) THEN !test whether `a` is cation
                   name_intra(1:n_intra) = p_pair(i_pair)%intra(1:n_intra)
                   fac_C(1:n_intra) = p_pair(i_pair)%fac_C(1:n_intra)
                 ELSE !or `a` is anion
                   name_intra(1:n_intra) = p_pair(i_pair)%intra(n_intra+1:n_intra+n_intra)
                   fac_C(1:n_intra) = p_pair(i_pair)%fac_C(n_intra+1:n_intra+n_intra)
                 ENDIF
                 
                 CALL update_intra_cpl(cos_latt, name_intra, fac_C, intra_cpl)
               ENDIF

             ENDIF
             

             !----------------------------------------------------------------
             ! Anion/Cation pairs are eventually swapped for Cation/Anion:
             ! e.g. Vpa,sc => Vsc,pa  
             !----------------------------------------------------------------
             CALL check_ref_pair(basis, i_a, i_b, &
                                 p_parent(i_parent)%pair(i_pair), &
                                 ref_pair )

             n_cpl = SIZE( p_parent( i_parent)%pair( i_pair )%cpl )

             !-------------------------------------------------------------
             ! Calculate exp( i . k . l )
             !-------------------------------------------------------------
             phase = pi*2.0d0*DOT_PRODUCT(veck,current_near%vec(:,i_n))
             
             IF (optmat) THEN
                D_tmp = MATMUL( basis%prim, current_near%vec(:,i_n) ) &
                        + D_basis( basis, i_b, i_a )
                ! Opt mat-element according to corrected TB-matrix
                ! PRB, 63, 201101(r) (2001)
                exp_fac = (0.d0,1.d0)*D_tmp(poldir)*EXP(j*phase)
                
             ELSE
                
                exp_fac = EXP(j*phase)
                
             END IF
             

             tb_cpl = 0.0D0

             IF (dg_bond) THEN
                
                i_dg_a = i_dg_a + 1
                i_dg = i_dg + 1

                ALLOCATE( hydro(i_dg_a)%mat(n_st_a, 1), STAT = err )
                IF (err.NE.0) CALL alloc_error('TB_ham','sparse_ham','hydro%mat' )                
                
                hydro(i_dg_a)%mat = 0.d0

                IF (scheme .EQ. 'tan') THEN
                  DO i_dump = 1, n_cpl
                    tb_cpl(ref_pair(i_dump)) = p_parent(i_parent)%pair(i_pair)%energy(i_dump)
                  ENDDO
                ELSEIF (scheme .EQ. 'jancu') THEN 
                  CALL harrison_scaling(.false., interface_data%type, &
                     p_parent, i_parent, i_pair, current_near%dist(i_n), &
                     interface_data%bowing,  interface_data%gap_corr, &
                     ref_pair, interface_data%alloy_random, n_cpl, tb_cpl)
                ENDIF
                
                tb_cpl = tb_cpl * d_H     

                CALL koster_slater(cos_latt, tb_cpl, tb_bloc)

                ! Calculate Peierls phase factor
                peierls_phase_factor = calculate_peierls_phase(basis%coord(:, i_a), basis%coord(:, i_b), &
                                                               magnetic_field_vector, gauge_choice)

                hydro( i_dg_a )%mat(:, 1) = exp_fac * peierls_phase_factor * tb_bloc(1:size(hydro( i_dg_a )%mat), 1)

             ELSE

                IF (scheme .EQ. 'tan') THEN
                  
                  i_1nn_list = i_1nn_list + 1

                  ! got bare intercouplings and the P, S, Q parameters
                  DO i_dump = 1, n_cpl
                    tb_cpl(ref_pair(i_dump)) = p_pair(i_pair)%energy(i_dump) * &
                           EXP( - p_pair(i_pair)%scaling(i_dump)*dist_phase)
                    P_list(ref_pair(i_dump)) = p_pair(i_pair)%fac_P(i_dump)
                    S_list(ref_pair(i_dump)) = p_pair(i_pair)%fac_S(i_dump)
                    Q_list(ref_pair(i_dump)) = p_pair(i_pair)%fac_Q(i_dump)
                  ENDDO
                  !IF (i_a .EQ. p_start) THEN
                  !  print*,'tb_cpl: ', tb_cpl
                  !ENDIF
                  ! add dipole correction for intercoupling parameters
                  CALL add_dip_corr(P_list, S_list, n_cpl, n_2nn_list, &
                                    unit_disp_1nn_list, dd_rel_1nn_list, &
                                    unit_disp_2nn_list, dd_rel_2nn_list, &
                                    i_1nn_list, i_2nn_list, tb_cpl)
                  !IF (i_a .EQ. p_start) THEN
                  !  print*,'tb_cpl: ', tb_cpl
                  !ENDIF
                  ! compute intercoupling matrix block
                  CALL koster_slater(cos_latt, tb_cpl, tb_bloc)

                  ! compute quadrupole correction ...
                  Q_bloc = 0.D0
                  CALL koster_slater(cos_latt, Q_list, Q_bloc)
                  
                  CALL get_quad_corr(Q_bloc, unit_disp_1nn_list, unit_disp_2nn_list, &
                                     n_2nn_list, i_1nn_list, i_2nn_list, quad_corr)
                  
                  ! ... then add it to intercoupling block 
                  tb_bloc = tb_bloc + quad_corr

                  ! done with this `b` neighbor,
                  ! move to before the next zone of the next `b` in 2nn-lists
                  i_2nn_list = i_2nn_list + n_2nn_list(i_1nn_list)
                  
                ELSEIF (scheme .EQ. 'jancu') THEN 

                  CALL harrison_scaling(scale, interface_data%type, &
                     p_parent, i_parent, i_pair, current_near%dist(i_n), &
                     interface_data%bowing,  interface_data%gap_corr, &
                     ref_pair, interface_data%alloy_random, n_cpl, tb_cpl)
                  
                  CALL koster_slater(cos_latt, tb_cpl, tb_bloc)

                ENDIF
                !IF (i_a .EQ. p_start) THEN
                !  print*,'tb_bloc: ', tb_bloc
                !ENDIF

                ! Calculate Peierls phase factor (after tb_bloc is computed)
                peierls_phase_factor = calculate_peierls_phase(basis%coord(:, i_a), basis%coord(:, i_b), &
                                                               magnetic_field_vector, gauge_choice)

                b_ref(1:n_st_b)= mat_data(i_b_mat)%ion(i_b_ion)%ind_ref(1:n_st_b)

                ! different atoms can have different orbital sets
                DO a_st = 1, n_st_a 
                  DO b_st = 1, n_st_b 
                    cpl(i_cpl)%mat(a_st,b_st) = cpl(i_cpl)%mat(a_st,b_st) &
                                                + exp_fac * peierls_phase_factor * tb_bloc(a_ref(a_st), b_ref(b_st))
                  END DO
                END DO
                
                !IF (i_a .EQ. p_start) THEN
                !  print*,"cpl%mat(1,:): ", cpl(i_cpl)%mat(1,:)
                !ENDIF
             END IF
             
          ELSE ! if b==a, remember current i_cpl as i_onsite
            i_onsite = i_cpl

          END IF  ! End of loop through neighbor map of `a`

       END DO nn_loop  
       
       !-------------------------------------------------------------
       ! Get final onsite terms then store them
       !------------------------------------------------------------
       
       IF (optmat) THEN
          ! intra atomic optical matrix elements
          !call intratomic_optics(poldir, tb_bloc)
          
          !DO a_st = 1, n_st( i_a )
          !   DO b_st = 1, n_st( i_a )
          
          !      cpl( i_cpl )%mat( a_st, b_st ) = &
          !           tb_bloc( a_ref( a_st ), a_ref( b_st ) )
          !   END DO
          !END DO
          
       ELSE    

          !----------------------------------------------------------------------
          ! Add potential energy
          !----------------------------------------------------------------------
          IF ( potential_flag ) THEN
             onsite(1:n_st_a) = onsite(1:n_st_a) - pot_data( i_a )
          ENDIF
         
          !----------------------------------------------------------------------
          ! Put numbers into diagonal block
          !----------------------------------------------------------------------
          DO a_st = 1, n_st_a !Matthias
             h_diag%mat(a_st,a_st) = onsite(a_st)*(1.d0,0.d0) !+ h_diag%mat(a_st,a_st)
          END DO

          !----------------------------------------------------------------------
          ! COMMENT ON THE ROTATION:
          ! In principle, the Slater-Koster framework is general no matter the
          ! choice of the Cartesian coordinates (x,y,z), provided that
          ! all the atomic-like orbitals of the same subshell, i.e. same (n,l),
          ! are associated with the same onsite parameter. For example, all three
          ! p-orbitals (similarly, all five d-orbitals) are equivalent, as how most
          ! of the ETB schemes treat zincblende crystals.
          ! However, for some other crystals, Jancu scheme and others do a
          ! trick on the onsite parameters of d-orbitals by dividing them into
          ! two subgroups, associated with two different onsite parameters.
          ! This is a phenomenological way to account for crystal-field splitting
          ! in first-nearest-neighbor schemes because of lacking 2NN info.
          ! (A similar trick should be done for p-orbitals.)
          ! This trick is not mathematically rigorous, but physically acceptable
          ! because, in spirit, it is equivalent to choosing a new basis set which
          ! adapts to the crystal symmetry instead of the full spherical symmetry
          ! of the isolated atom as for the atomic-like orbitals.
          ! Such a trick only makes sense in a reference Cartesian system where,
          ! e.g. z-axis is chosen along the growth direction in wurtzite case,
          ! or z-axis is chosen to be the normal to the material plane in 2D case.
          ! Thus, if the Cartesian system is not the reference Cartesian system,
          ! one need to perform a rotation of the onsite blocks from the reference
          ! Cartesian system to the actual Cartesian system used in the calculation.
          !
          ! In Tan scheme, the basis set are always the atomic-like orbitals,
          ! no specific choice of Cartesian system is assumed, so no rotation is needed.
          ! This is an advantage for the polytype transferability of the TB calculations
          ! because the same TB parameters can be used for different polytypes.
          !----------------------------------------------------------------------
          IF (scheme .EQ. 'jancu') THEN !Do rotation for Jancu scheme only
            ! d orbital rotation 
            IF (d_orbital_present(mat_data(i_a_mat),i_a_ion) ) THEN
               d_onsite = 0
               IF ( inv_ref(dxy) .gt. 0)   d_onsite(1) = onsite( inv_ref(dxy) )
               IF ( inv_ref(dyz) .gt. 0)   d_onsite(2) = onsite( inv_ref(dyz) )
               IF ( inv_ref(dzx) .gt. 0)   d_onsite(3) = onsite( inv_ref(dzx) )
               IF ( inv_ref(dx2y2) .gt. 0) d_onsite(4) = onsite( inv_ref(dx2y2) )
               IF ( inv_ref(dz2r2) .gt. 0) d_onsite(5) = onsite( inv_ref(dz2r2) )
               IF ((d_onsite(1) .ne. d_onsite(2)) .or. &
                   (d_onsite(1) .ne. d_onsite(3)) .or. &
                   (d_onsite(1) .ne. d_onsite(4)) .or. &
                   (d_onsite(1) .ne. d_onsite(5))) THEN
  
                  u = mat_data(i_a_mat)%ref_axis
                  call rotate_d(d_onsite, c_axis, u, inv_ref, h_diag%mat)
               END IF
            END IF
            
            ! p orbital rotation
            IF (p_orbital_present(mat_data(i_a_mat),i_a_ion)) THEN
              p_onsite = 0
              IF (inv_ref(px) .gt. 0) p_onsite(1) = onsite( inv_ref(px) )
              IF (inv_ref(py) .gt. 0) p_onsite(2) = onsite( inv_ref(py) )
              IF (inv_ref(pz) .gt. 0) p_onsite(3) = onsite( inv_ref(pz) )
              IF ((p_onsite(1) .ne. p_onsite(2)) .or. &
                  (p_onsite(1) .ne. p_onsite(3))) THEN
                  
                  u = mat_data(i_a_mat)%ref_axis
                  call rotate_p(p_onsite, c_axis, u, inv_ref, h_diag%mat)
              END IF
            END IF
          END IF

          ! Matthias: Now add the onsite block to cpl(i_onsite)
          cpl(i_onsite)%mat = cpl(i_onsite)%mat + h_diag%mat
          
          ! Matthias: Now add the intra-atomic couplings in case of Tan scheme
          IF (scheme .EQ. 'tan') THEN
            cpl(i_onsite)%mat = cpl(i_onsite)%mat + intra_cpl*(1.d0,0.d0)
          ENDIF
          
       END IF
       

       !!---------------------------------------------------------------
       ! PASSIVATION
       ! (Modified by Gabriele Penazzi, 11/09) Pecchia (11/2/2012)
       !!---------------------------------------------------------------
       IF (n_dg(i_a).gt.0) THEN
        IF (upt%hybrid_passivation) THEN 
          !call Klimeck_passiv(cpl(i_onsite)%mat, i_dg_a)
          DO i_dg_a = 1, n_dg(i_a)    
            ! Renormalize onsite taking E = 0
            DO b_st = 1, n_st_a
                DO a_st = 1, n_st_a
                  cpl(i_onsite)%mat(a_st,b_st) = cpl(i_onsite)%mat(a_st,b_st) - & 
                       hydro(i_dg_a)%mat(a_st,1)*conjg(hydro(i_dg_a)%mat(b_st,1))/E_H
                ENDDO
            ENDDO
          ENDDO
        ELSE
          DO i_cpl = 1, n_neig     
            IF (i_cpl.NE.i_onsite) THEN
                cpl(i_cpl)%mat = cpl(i_cpl)%mat * d_H
            ELSE    
                DO a_st = 1, n_st_a
                  cpl(i_onsite)%mat(a_st,a_st) = cpl(i_onsite)%mat(a_st,a_st) + E_H
                END DO
            ENDIF
          ENDDO        
        ENDIF
       ENDIF

       !----------------------------------------------------------------------
       ! Store TB Hamiltonian in CSB format 
       !----------------------------------------------------------------------
       spin_loop:DO i_spin_a = 1, n_spin
          
          sp_a = i_spin_a !+ n_spin * (i_a - 1) 
          !-------------------------------------------------------------------
          ! Calculate the relevant spin orbit matrixes
          !-------------------------------------------------------------------
          IF (relat) THEN
             IF (optmat) THEN            
                so_diag =(0.0D0,0D0); so_off_diag=(0.D0,0.D0)
             ELSE
                CALL spin_sparse( i_spin_a, so_cpl, so_diag, so_off_diag )
             END IF
          END IF


          col_loop:DO i_cpl = 1, n_neig
             i_b = cpl_ind( i_cpl )
      
             n_st_b = n_st(i_b) 

             ! up
             !    dn
             DO i_spin_b = 1, n_spin

                sp_b = n_spin * (i_cpl - 1) + i_spin_b  

                call create_matrix(M%Row(sp_a)%B(sp_b), n_st_a, n_st_b)

                IF (relat) THEN
                   call create_matrix(UM%Row(sp_a)%B(sp_b), n_st_a, n_st_b)
                   UM%Row(sp_a)%B(sp_b)%val = (0.d0, 0.d0)
                ENDIF

                IF ( i_b .EQ. i_a ) THEN  
                   IF (i_spin_a .EQ. i_spin_b) THEN
                      M%Row(sp_a)%B(sp_b)%val = cpl(i_cpl)%mat + so_diag(1:n_st_a,1:n_st_b) 
                   ELSE
                      M%Row(sp_a)%B(sp_b)%val = so_off_diag(1:n_st_a,1:n_st_b)
                      IF (i_spin_a.eq.1) THEN
                        IF (relat) UM%Row(sp_a)%B(sp_b)%val = -U_bl(1:n_st_a,1:n_st_b)
                      ELSE
                        IF (relat) UM%Row(sp_a)%B(sp_b)%val = U_bl(1:n_st_a,1:n_st_b) 
                      ENDIF
                   ENDIF  
                ELSE   
                   IF (i_spin_a .EQ. i_spin_b) THEN
                     M%Row(sp_a)%B(sp_b)%val = cpl(i_cpl)%mat
                   ELSE
                     M%Row(sp_a)%B(sp_b)%val = (0.d0, 0.d0)
                   ENDIF 
                ENDIF        

                ! store the column start position of the block
                M%Row(sp_a)%col(sp_b) = col_start(n_spin * (i_b-1) + i_spin_b)
                IF (relat) UM%Row(sp_a)%col(sp_b) = col_start(n_spin * (i_b-1) + i_spin_b)

                n_ind = n_ind + n_st_a*n_st_b

             END DO  !i_spin_b

          END DO col_loop

       END DO spin_loop  ! i_spin_a
       !----------------------------------------------------------------------
       ! Update the CSR Hamiltonian 
       !----------------------------------------------------------------------
       call csb_to_csr(M,TM,1.d-10)

       nnz = TM%nnz

       if (n_val+nnz.lt.n_zeros_h) then
          upt%ham%M(n_val+1: n_val+nnz) = TM%M(1:nnz)      
          upt%ham%Mj(n_val+1: n_val+nnz) = TM%Mj(1:nnz) 
          upt%ham%Mi(n_row+1:n_row+n_spin*n_st_a) = n_val + TM%Mi(1:n_spin*n_st_a) 
          n_val = n_val + nnz
          upt%ham%nnz = n_val
          upt%ham%Mi(n_row+n_spin*n_st_a+1) = n_val+1
       else    
          write(*,*) "STOP: need to enlarge the matrix"
          write(*,*) "n_val=",n_val
          write(*,*) "last block nnz=",nnz
          write(*,*) "M%srtrow M%endrow size: ", M%srtrow, M%endrow
          write(*,*) "size: ",size(M%Row(M%srtrow)%B(1)%val,1),size(M%Row(M%srtrow)%B(2)%val,1)
          stop
       endif

       call destroy_matrix(TM)

       if (relat) then
          call csb_to_csr(UM,TM,1.d-10)
          nnz = TM%nnz
          upt%U%M(n_val_U+1: n_val_U+nnz) = TM%M(1:nnz)      
          upt%U%Mj(n_val_U+1: n_val_U+nnz) = TM%Mj(1:nnz) 
          upt%U%Mi(n_row+1 : n_row+n_spin*n_st_a) = n_val_U + TM%Mi(1:n_spin*n_st_a) 
          n_val_U = n_val_U + nnz
          upt%U%nnz = n_val_U
          upt%U%Mi(n_row+n_spin*n_st_a+1) = n_val_U+1
          call destroy_matrix(TM)
       endif

       n_row = n_row + n_spin * n_st_a                

       !----------------------------------------------------------------------
       ! Deallocate work arrays
       !----------------------------------------------------------------------
       CALL destroy_matrix(M)
       if (relat) CALL destroy_matrix(UM)

       DO i_cpl = 1,  SIZE( current_near%ind )
          IF ( ASSOCIATED( cpl(i_cpl)%mat ) ) THEN
             DEALLOCATE( cpl( i_cpl )%mat, STAT = err )
             IF ( err .NE. 0 ) &
                  CALL dealloc_error( 'TB_ham', 'sparse_ham', 'cpl%mat' )
          END IF
       END DO

       DEALLOCATE(h_diag%mat, STAT = err) !Matthias
       IF ( err .NE. 0 ) &
        CALL dealloc_error( 'TB_ham', 'sparse_ham', 'h_diag%mat' )

       DO i_dg_a = 1, n_dg(i_a)

          DEALLOCATE( hydro( i_dg_a )%mat, STAT = err )
          IF ( err .NE. 0 ) &
               CALL dealloc_error( 'TB_ham', 'sparse_ham', 'hydro%mat' )

       END DO

       !----------------------------------------------------------------------
       ! Increment the nearest neighbours list
       !----------------------------------------------------------------------

    END DO basis_loop  ! end loop on i_a
    
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !  FINALIZE
    !----------------------------------------------------------------------
    IF ( verbose.gt.0 )  WRITE ( *, * ) "(TB ham) atom a = ", n_basis 
    
    n_value = n_val
    n_sparse = n_ind
    n_row = n_row - offset

    !write(*,*) '(TB_ham) CPU '//id_s//' Hamiltonian dim= ',n_ham
    !write(*,*) '(TB_ham) CPU '//id_s//' non zero values= ',n_val
    !write(*,*) '(TB_ham) CPU '//id_s//' non zeros of U = ',n_val_U


    ! deallocate any previously defined shift
    ! known problem: global shifts are not thread-safe. Move 'em in upt 
    call kill_shifts()
    
    call init_shifts(n_ham, offset)
    
    ! ADDITIONAL MATRIX INFO (Pecchia)
    mem = 16_8 * n_val 
    CALL memstr(mem,dec,str)
    if (id0 .and. verbose.gt.0) &
          write(*,'(a28,f10.2,a3)') ' (TB_ham) memory for values: ',real(mem)/real(dec),str
    
    mem = 4_8 * n_val 
    CALL memstr(mem,dec,str)    
    if (id0 .and. verbose.gt.0) &
          write(*,'(a28,f10.2,a3)') ' (TB_ham) memory for indeces:',real(mem)/real(dec),str

    if (id0 .and. verbose.gt.0) write(*,*) '(TB_ham) matrix done'
    if (id0 .and. verbose.gt.0) write(*,*)

    !if (id0 .and. verbose.gt.0) write(*,*) '(TB_ham) '//id_s//' Sort CSR ...'

    !call sort_csr(n_row, upt%ham%M, upt%ham%Mj, upt%ham%Mi, .true.)

    !if (relat) then
    !  call sort_csr(n_row, upt%U%M, upt%U%Mj, upt%U%Mi, .true.) 
    !endif
    
    DEALLOCATE(so_diag,so_off_diag,tb_bloc,onsite,temp_onsite,tb_cpl,ons_corr)
    DEALLOCATE(ref_pair,bl_arr,col_start)
    IF (scheme .EQ. 'tan') THEN
      DEALLOCATE(P_list,S_list,Q_list,Q_bloc,quad_corr,intra_cpl)
      DEALLOCATE(dist_1nn_list,unit_disp_1nn_list,dd_rel_1nn_list,n_2nn_list)
      DEALLOCATE(dist_2nn_list,unit_disp_2nn_list,dd_rel_2nn_list)
    ENDIF
 
  END SUBROUTINE sparse_ham


  !=================================================================================
  !**************************************************************
  ! TB parameters Harrison-scaling with distance
  ! Jancu scheme only
  ! If bond deviations are to be considered, scale parameters
  !**************************************************************
  ! Harrison 's scaling :
  !
  ! cpl_strained = cpl_unstrained *
  !               
  !  * ( dist_unstrained / dist_strained )^power
  !
  !************************************************************* 
  
  SUBROUTINE harrison_scaling(scale, type, parent, i_parent, i_pair, dist, &
                     bowing, gap_corr, ref_pair, alloy_random, n_cpl, tb_cpl)

    LOGICAL, INTENT(IN) :: scale
    CHARACTER ( LEN = 20 ), INTENT(IN) :: type
    TYPE(parent_data), DIMENSION(:), POINTER :: parent
    INTEGER, INTENT(IN) :: i_parent, i_pair
    REAL(dp), INTENT(IN) :: dist
    REAL(dp), DIMENSION(:) :: bowing
    REAL(dp), DIMENSION(2) :: gap_corr
    INTEGER, DIMENSION(:) :: ref_pair
    LOGICAL, INTENT(IN) :: alloy_random
    INTEGER, INTENT(IN) :: n_cpl

    REAL(dp), DIMENSION(:), INTENT(INOUT)  :: tb_cpl    
    !_______________________________________________________________________
    ! local variables
    REAL(dp),  DIMENSION( : ), POINTER :: energy1, energy2
    REAL(dp),  DIMENSION( : ), POINTER :: scaling1, scaling2
    REAL(dp) :: content1, content2
    REAL(dp) :: dist_ref1, dist_ref2, dist_ref
    REAL(dp) :: scaling_factor1, scaling_factor2, scaling_factor
    REAL(dp) :: scaled_factor_1, scaled_factor_2
    INTEGER :: i_c

  
    IF ( alloy_random .OR. TRIM(type).EQ.'simple' &
                      .OR. TRIM(type).EQ.'binary' ) THEN

       ! --> random alloy : do not average
       energy1 => parent(i_parent)%pair( i_pair )%energy
       scaling1 => parent(i_parent)%pair( i_pair )%scaling 
       dist_ref1 = parent(i_parent)%pair( i_pair )%dist_ref
        
       IF ( scale ) then             
          scaling_factor1 = dist_ref1 / dist
       ELSE
          scaling_factor1 = 1.d0
       END IF 

       DO i_c = 1, n_cpl          
          
          tb_cpl( ref_pair( i_c ) ) = energy1( i_c ) &
               * ( scaling_factor1 )**scaling1( i_c )
          
       ENDDO

       tb_cpl(sps) = tb_cpl(sps) * (1.d0+gap_corr(2))
       tb_cpl(pss) = tb_cpl(pss) * (1.d0+gap_corr(2))
       !print*,'tb_cpl(sss)',tb_cpl(sss)
    ELSE

       SELECT CASE( TRIM(type) )

       CASE('ternary')

          energy1 =>  parent(1)%pair( i_pair )%energy
          scaling1 => parent(1)%pair( i_pair )%scaling 
          dist_ref1 = parent(1)%pair( i_pair )%dist_ref
          content1  = parent(1)%content

          energy2 =>  parent(2)%pair( i_pair )%energy
          scaling2 => parent(2)%pair( i_pair )%scaling 
          dist_ref2 = parent(2)%pair( i_pair )%dist_ref
          content2  = parent(2)%content
       
          ! scale parent 1          
          IF (scale) THEN
             scaling_factor1 = dist_ref1 / dist
             scaling_factor2 = dist_ref2 / dist
          ELSE
             scaling_factor1 = 1.d0
             scaling_factor2 = 1.d0
          ENDIF

          DO i_c = 1, n_cpl

             !______________________________________________________________
             ! Jerome scaling: 
             ! V = x * V1 * (d1/d)^n1 + (1-x) * V2 * (d2/d)^n2 
             !           
             scaled_factor_1 = energy1(i_c) &
                  * (scaling_factor1)**scaling1(i_c) 
             
             scaled_factor_2 = energy2(i_c) &
                  * (scaling_factor2)**scaling2(i_c) 
             
             tb_cpl( ref_pair( i_c ) ) = &
                  content1 * scaled_factor_1 + content2 * scaled_factor_2 &
                  - bowing(1) * SQRT( content1 * content2 ) &
                  * ( scaled_factor_1 - scaled_factor_2 )
             !______________________________________________________________
             ! Aldo scaling:
             ! d0 =  x * d1 + (1-x) * d2
             !  n =  x * n1 + (1-x) * n2
             !  V = (x * V1 + (1-x) * V2) * (d0/d)^n
             !           
             !dist_ref = content1 * scaling_factor1 + content2 * scaling_factor2

             !scaling_factor = content1 * scaling1(i_c) + content2 * scaling2(i_c)             

             !tb_cpl( ref_pair( i_c) ) = (content1 * energy1(i_c) + content2 * energy2(i_c)) &
             !                    * dist_ref**scaling_factor
             !______________________________________________________________

          END DO

          tb_cpl(sps) = tb_cpl(sps) * (1.d0+gap_corr(2))
          tb_cpl(pss) = tb_cpl(pss) * (1.d0+gap_corr(2))


       CASE ( 'quaternary' )
          
          ! to do : implement alloying for quaternary pairs
          
       CASE DEFAULT
          
       END SELECT

    END IF
     
  END SUBROUTINE harrison_scaling


  !========================================================================
  ! Calculate the list of ratios delta_d/ave_d (Tan2016PRB)
  !========================================================================
  SUBROUTINE get_dd_rel_2nn_list(n_2nn, dist, dd)

    INTEGER, DIMENSION(:), INTENT(IN) :: n_2nn !`n_2nn_list`
    REAL(dp), DIMENSION(:), INTENT(IN) :: dist !`dist_2nn_list`
    REAL(dp), DIMENSION(:), INTENT(OUT) :: dd !`dd_rel_2nn_list`

    !local
    INTEGER :: i, i_start, i_end
    REAL(dp) :: ave_bond_b !averaged bondlength around `b`

    i_start = 1
    DO i = 1, SIZE(n_2nn) !scan through each `b`
      i_end = i_start+n_2nn(i) - 1 !end of the zone for this `b` in `dist_2nn_list`
      ave_bond_b = SUM(dist(i_start:i_end))/n_2nn(i)
      dd(i_start:i_end) = dist(i_start:i_end)/ave_bond_b - 1.0d0
      i_start = i_end + 1 !move to the next zone for the next `b`
    ENDDO

  END SUBROUTINE get_dd_rel_2nn_list


  !========================================================================
  ! Get corrections to onsite energies (including band offset)
  !========================================================================
  SUBROUTINE get_ons_corr(name,pair,ion,dist,n_orb_b,ons)

    CHARACTER(LEN = 10) :: scheme
    CHARACTER(ATOMLEN) :: name
    TYPE(pair_coupling) :: pair
    TYPE(ion_orbit) :: ion
    REAL(dp) :: dist 
    INTEGER :: n_orb_b
    REAL(dp), DIMENSION(:), INTENT(OUT) :: ons !onsite correction

    !local
    INTEGER :: i_state, ii, iii, i_offset
    

      IF (name .NE. pair%name(1)) THEN !test whether this atom is cation
        i_offset = 0
      ELSE !or is anion
        i_offset = n_orb_b
      ENDIF

      i_state = 0
      DO ii = 1, SIZE(ion%deg)
        DO iii = 1, ion%deg(ii)
          i_state = i_state + 1
          ons(i_state) = pair%fac_I(i_offset+ii) * &
                         EXP(-pair%l_I(i_offset+ii)*dist)
        ENDDO
      ENDDO

      ons(1:i_state) = ons(1:i_state) + &
                       pair%fac_O*EXP(-pair%l_O*dist)

  END SUBROUTINE get_ons_corr


  !========================================================================
  ! Get corrections to SO couplings
  !========================================================================
  SUBROUTINE get_so_corr(name,pair,n_so,so)

    CHARACTER(LEN = 10) :: scheme
    CHARACTER(ATOMLEN) :: name
    TYPE(pair_coupling) :: pair
    INTEGER :: n_so
    REAL(dp), DIMENSION(:), INTENT(OUT) :: so 
    
      IF (name .NE. pair%name(1)) THEN !test whether this atom is cation
        so(1:n_so) = pair%fac_so(1:n_so)
      ELSE !or is anion
        so(1:n_so) = pair%fac_so(n_so+1:n_so+n_so)
      ENDIF

  END SUBROUTINE get_so_corr


  !========================================================================
  ! Update intra-atomic coupling of `a` due to current pair `a`-`b`
  ! Tan scheme only
  !========================================================================
  SUBROUTINE update_intra_cpl(vec, name, C, intra)
    
    REAL(dp), DIMENSION(3) :: vec !unit displacement vector a->b
    CHARACTER(LEN = 5), DIMENSION(:) :: name !list of name of intracoupling
    REAL(dp), DIMENSION(:) :: C !list of C parameters
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: intra

    !local
    INTEGER :: ii
    CHARACTER(LEN = 5) :: n
    REAL(dp), DIMENSION(n_ref_states,n_ref_states) :: M1, M2 !M matrices
    
    CALL get_M_dip(vec, M1) !due to dipole component
    CALL get_M_quad(vec, M2) !due to quadrupole component
    
    DO ii = 1, SIZE(C)
      n = TRIM(name(ii))
      IF (n .EQ. 'spca' .OR. n .EQ. 'spac') &
         intra(s,px:pz) = intra(s,px:pz) + M1(s,px:pz)*C(ii)
      
      IF (n .EQ. 'sdca' .OR. n .EQ. 'sdac') &
         intra(s,dxy:dz2r2) = intra(s,dxy:dz2r2) + M2(s,dxy:dz2r2)*C(ii)

      IF (n .EQ. 'ppca' .OR. n .EQ. 'ppac') &
         intra(px:pz,px:pz) = intra(px:pz,px:pz) + M2(px:pz,px:pz)*C(ii)

      IF (n .EQ. 'ps*ca' .OR. n .EQ. 'ps*ac') &
         intra(px:pz,se) = intra(px:pz,se) + M1(px:pz,se)*C(ii)

      IF (n .EQ. 'pdca' .OR. n .EQ. 'pdac') &
         intra(px:pz,dxy:dz2r2) = intra(px:pz,dxy:dz2r2) + &
                                  M1(px:pz,dxy:dz2r2)*C(ii)

      IF (n .EQ. 's*dca' .OR. n .EQ. 's*dac') &
         intra(se,dxy:dz2r2) = intra(se,dxy:dz2r2) + &
                               M2(se,dxy:dz2r2)*C(ii)

      IF (n .EQ. 'ddca' .OR. n .EQ. 'ddac') &
         intra(dxy:dz2r2,dxy:dz2r2) = intra(dxy:dz2r2,dxy:dz2r2) + &
                                      M2(dxy:dz2r2,dxy:dz2r2)*C(ii)
    ENDDO
    
    intra(px:dz2r2,s:s) = transpose(intra(s:s,px:dz2r2))
    intra(se:dz2r2,px:pz) = transpose(intra(px:pz,se:dz2r2))
    intra(dxy:dz2r2,se:se) = transpose(intra(se:se,dxy:dz2r2))

  END SUBROUTINE update_intra_cpl


  !========================================================================
  ! Add corrections due to dipole component of potential to intercouplings
  ! Tan scheme only
  !========================================================================
  SUBROUTINE add_dip_corr(P, S, n_cpl, n_2nn, vec1, dd1, vec2, dd2, i1, i2, tb_cpl)

    REAL(dp), DIMENSION(:) :: P, S !lists of P, S parameters
    INTEGER :: n_cpl
    INTEGER, DIMENSION(:) :: n_2nn ! `n_2nn_list`
    INTEGER :: i1, i2 !i_1nn_list, i_2nn_list
    REAL(dp), DIMENSION(:,:) :: vec1, vec2 !`unit_disp_1nn_list`, `unit_disp_2nn_list`
    REAL(dp), DIMENSION(:) :: dd1, dd2 !`dd_rel_1nn_list` and `dd_rel_2nn_list`
    !NOTE: How averaged bondlength, and inherently, dd should be evaluated in cases of
    !      there are missing neighboring atoms is quite tricky. Here we assume ignoring
    !      missing atoms from the calculations is a good approximation.
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: tb_cpl
  
    !local
    REAL(dp), DIMENSION(3) :: vec_ik, vecdd_ik, vec_jk, vecdd_jk
    REAL(dp) :: p_sum, q_sum
    INTEGER :: i
    
    vec_ik = 0.d0
    vecdd_ik = 0.d0
    DO i = 1, SIZE(dd1) !scan through 1NN of `a`
      vec_ik = vec_ik + vec1(i,1:3)
      vecdd_ik = vecdd_ik + vec1(i,1:3)*dd1(i)
    ENDDO
    
    vec_jk = 0.d0
    vecdd_jk = 0.d0
    DO i = i2 + 1, i2 + n_2nn(i1)
      vec_jk = vec_jk + vec2(i,1:3)
      vecdd_jk = vecdd_jk + vec2(i,1:3)*dd2(i)
    ENDDO
    
    p_sum = + vec1(i1,1)*vec_ik(1) - vec1(i1,2)*vec_ik(2) + vec1(i1,3)*vec_ik(3) &!p_ij
            - vec1(i1,1)*vec_jk(1) + vec1(i1,2)*vec_jk(2) - vec1(i1,3)*vec_jk(3)  !p_ji
    !print*, 'p_sum: ', p_sum
    q_sum = + vec1(i1,1)*vecdd_ik(1) - vec1(i1,2)*vecdd_ik(2) + vec1(i1,3)*vecdd_ik(3) &!q_ij
            - vec1(i1,1)*vecdd_jk(1) + vec1(i1,2)*vecdd_jk(2) - vec1(i1,3)*vecdd_jk(3)  !q_ji
    !print*, 'q_sum: ', q_sum
    !NOTE: factor 4*pi/3 has been absorbed into p_sum and q_sum
    tb_cpl(1:n_cpl) = tb_cpl(1:n_cpl) + P(1:n_cpl)*p_sum + S(1:n_cpl)*q_sum

  END SUBROUTINE add_dip_corr


  !============================================================================
  ! Get corrections due to quadrupole component of potential to intercouplings
  ! Tan scheme only
  !============================================================================
  SUBROUTINE get_quad_corr(Q, vec1, vec2, n_2nn, i1, i2, quad)

    REAL(dp), DIMENSION(:,:) :: Q !Q_bloc
    REAL(dp), DIMENSION(:,:) :: vec1, vec2 !`unit_disp_1nn_list`, `unit_disp_2nn_list`
    INTEGER, DIMENSION(:) :: n_2nn !`n_2nn_list`
    INTEGER :: i1, i2 !current `i_1nn_list` and `i_2nn_list`
    REAL(dp), DIMENSION(n_ref_states,n_ref_states), INTENT(OUT) :: quad
  
    !local
    INTEGER :: i, j
    REAL(dp), DIMENSION(n_ref_states,n_ref_states) :: Mik, Mjk, Mik_sum, Mjk_sum

    Mik_sum = 0.d0
    Mjk_sum = 0.d0

    DO i = 1, SIZE(n_2nn) !scan through each `b` of `a`
      CALL get_M_quad(vec1(i,:), Mik)
      Mik_sum = Mik_sum + Mik
    ENDDO

    DO j = i2 + 1, i2 + n_2nn(i1) !scan through each `c` of current `b`
      CALL get_M_quad(vec2(j,:), Mjk)
      Mjk_sum = Mjk_sum + Mjk
    ENDDO

    quad = MATMUL(Mik_sum,Q) + MATMUL(Q,Mjk_sum)

  END SUBROUTINE get_quad_corr


  !============================================================================
  ! Calculate part of M matrix due to dipole component of potential 
  ! Tan scheme only
  !============================================================================
  SUBROUTINE get_M_dip(D, M)

    REAL(dp), DIMENSION(3), INTENT(IN) :: D !unit displacement vector
    REAL(dp), DIMENSION(n_ref_states,n_ref_states), INTENT(OUT) :: M

    !local
    REAL(dp) :: fac

    fac = 1.d0/4.d0/pi
    M = 0.d0

    !s-p, p-s
    M(s,px:pz) = D(1:3)*sqrt(3.d0)*fac
    M(px:pz,s:s) = transpose(M(s:s,px:pz))
  
    !p-se, se-p
    M(se,px:pz) = M(s,px:pz)
    M(px:pz,se) = M(px:pz,s)
  
    !p-d, d-p
    M(px,dxy  ) = 3*D(2)*fac/sqrt(5.d0)
    M(px,dyz  ) = 0.0D0
    M(px,dzx  ) = 3*D(3)*fac/sqrt(5.d0)
    M(px,dx2y2) = 3*D(1)*fac/sqrt(5.d0)
    M(px,dz2r2) = -sqrt(3.d0)*D(1)*fac/sqrt(5.d0)
    M(py,dxy  ) = 3*D(1)*fac/sqrt(5.d0)
    M(py,dyz  ) = 3*D(3)*fac/sqrt(5.d0)
    M(py,dzx  ) = 0.0D0
    M(py,dx2y2) = -3*D(2)*fac/sqrt(5.d0)
    M(py,dz2r2) = -sqrt(3.d0)*D(2)*fac/sqrt(5.d0)
    M(pz,dxy  ) = 0.0D0
    M(pz,dyz  ) = 3*D(2)*fac/sqrt(5.d0)
    M(pz,dzx  ) = 3*D(1)*fac/sqrt(5.d0)
    M(pz,dx2y2) = 0.0D0
    M(pz,dz2r2) = 2*sqrt(3.d0)*D(3)*fac/sqrt(5.d0)
  
    M(dxy:dz2r2,px:pz) = transpose(M(px:pz,dxy:dz2r2))

  END SUBROUTINE get_M_dip


  !============================================================================
  ! Calculate part of M matrix due to quadrupole component of potential 
  ! Tan scheme only
  !============================================================================
  SUBROUTINE get_M_quad(D, M)

    REAL(dp), DIMENSION(3), INTENT( IN ) :: D !unit displacement vector
    REAL(dp), DIMENSION(n_ref_states,n_ref_states), INTENT(OUT) :: M

    !local
    REAL(dp) :: fac

    fac = 1.d0/4.d0/pi
    M = 0.d0
    
    !s-d, d-s
    !M(s,dxy) = D(1)*D(2)*sqrt(15.d0)*fac
    !M(s,dyz) = D(2)*D(3)*sqrt(15.d0)*fac
    !M(s,dzx) = D(3)*D(1)*sqrt(15.d0)*fac
    !M(s,dx2y2) = 0.5d0*(D(1)*D(1) - D(2)*D(2))*sqrt(15.d0)*fac
    !M(s,dz2r2) = 0.5d0*(2*D(3)*D(3) - D(1)*D(1) - D(2)*D(2))*sqrt(5.d0)*fac
    !M(dxy:dz2r2,s:s) = transpose(M(s:s,dxy:dz2r2))
  
    !p-p
    M(px,px) = (2*D(1)**2 - D(2)**2 - D(3)**2)*fac
    M(px,py) = D(1)*D(2)*3*fac
    M(px,pz) = D(3)*D(1)*3*fac
    M(py,py) = (2*D(2)**2 - D(1)**2 - D(3)**2)*fac
    M(py,pz) = D(2)*D(3)*3*fac
    M(pz,pz) = (2*D(3)**2 - D(1)**2 - D(2)**2)*fac
    M(py,px) = M(px,py)
    M(pz,px) = M(px,pz)
    M(pz,py) = M(py,pz)
  
    !se-d, d-se, not sure
    !M(se,dxy:dz2r2) = M(s,dxy:dz2r2)
    !M(dxy:dz2r2,se) = M(dxy:dz2r2,s)
  
    !d-d
    M(dxy  ,dxy  ) = -(2*D(3)**2 - D(1)**2 - D(2)**2)*5*fac/7
    M(dxy  ,dyz  ) = D(1)*D(3)*15*fac/7
    M(dxy  ,dzx  ) = D(2)*D(3)*15*fac/7
    M(dxy  ,dx2y2) = 0.0D0
    M(dxy  ,dz2r2) = -D(1)*D(2)*sqrt(3.d0)*10*fac/7
    M(dyz  ,dxy  ) = M(dxy  ,dyz  )
    M(dzx  ,dxy  ) = M(dxy  ,dzx  )
    M(dx2y2,dxy  ) = M(dxy  ,dx2y2)
    M(dz2r2,dxy  ) = M(dxy  ,dz2r2)
  
    M(dyz  ,dyz  ) = -(2*D(1)**2 - D(2)**2 - D(3)**2)*5*fac/7
    M(dyz  ,dzx  ) = D(1)*D(2)*15*fac/7
    M(dyz  ,dx2y2) = -D(2)*D(3)*15*fac/7
    M(dyz  ,dz2r2) = D(2)*D(3)*sqrt(3.d0)*5*fac/7
    M(dzx  ,dyz  ) = M(dyz  ,dzx  )
    M(dx2y2,dyz  ) = M(dyz  ,dx2y2)
    M(dz2r2,dyz  ) = M(dyz  ,dz2r2)
  
    M(dzx  ,dzx  ) = -(2*D(2)**2 - D(1)**2 - D(3)**2)*5*fac/7
    M(dzx  ,dx2y2) = D(1)*D(3)*15*fac/7
    M(dzx  ,dz2r2) = D(1)*D(3)*sqrt(3.d0)*5*fac/7
    M(dx2y2,dzx  ) = M(dzx  ,dx2y2)
    M(dz2r2,dzx  ) = M(dzx  ,dz2r2)
  
    M(dx2y2,dx2y2) = -(2*D(3)**2 - D(1)**2 - D(2)**2)*5*fac/7
    M(dx2y2,dz2r2) = -(D(1)**2 - D(2)**2)/sqrt(3.d0)*15*fac/7
    M(dz2r2,dx2y2) = M(dx2y2,dz2r2)
  
    M(dz2r2,dz2r2) = (2*D(3)**2 - D(1)**2 - D(3)**2)*5*fac/7

  END SUBROUTINE get_M_quad




  !===========================================================================
  !
  ! Subroutine "koster_slater" :
  !
  ! calculate the Koster-Slater two center integrals for atomic orbitals. 
  !
  !===========================================================================
  !
  ! INPUT :
  !
  ! => cos_latt( 3 )  - REAL( dp ) array : direction cosine
  !                                              of the molecular axis.
  !
  ! => t - REAL( dp ) array : coupling matrix elements with respect
  !                                 to the molecular axis.
  !
  ! => no_d_flag - logical : true if no d-states are included.
  !
  !
  ! OUTPUT :
  !
  ! => v - REAL( dp ) array : matrix elements with respect
  !                                 to cartesian axes.
  !
  !===========================================================================
  !
  ! NOTE : order of states follows from state.dat
  !
  !    1  2   3   4   5    6    7    8     9     10
  !    s  px  py  pz  s*  dxy  dyz  dzx  dx2y2  dz2r2
  !
  !===========================================================================

  SUBROUTINE koster_slater( cos_latt, t, v)
 
    !=========================================================================
 
    ! input arguments :
 
    REAL ( dp ), DIMENSION( : ), INTENT( IN ) :: t
    REAL ( dp ), DIMENSION( 3 ), INTENT( IN ) :: cos_latt
 
    !_________________________________________________________________________
 
    ! output arguments :
 
    REAL ( dp ), DIMENSION(:,:) :: v
 
    !_________________________________________________________________________
 
    ! local variables :
 
    REAL ( dp ) :: al, am, an, al2, am2, an2, dsq3
 
    !=========================================================================
 
    al = cos_latt( 1 )
    am = cos_latt( 2 )
    an = cos_latt( 3 )
 
    al2 = al**2
    am2 = am**2
    an2 = an**2
 
    dsq3 = SQRT( 3.0d0 )
    !=========================================================================
    ! s-s
    v( s, s )   = t( sss )  

    ! s-p, p-s
    v( s, px )  =   al * t( sps )
    v( px, s )  = - al * t( pss )

    v( s, py )  =   am * t( sps )
    v( py, s )  = - am * t( pss )
!
    v( s, pz )  =   an * t( sps )
    v( pz, s )  = - an * t( pss )

    ! p-p
    v( px, px ) = al2 * t( pps ) + ( 1.0d0 - al2 ) * t( ppp )
    v( py, py ) = am2 * t( pps ) + ( 1.0d0 - am2 ) * t( ppp )
    v( pz, pz ) = an2 * t( pps ) + ( 1.0d0 - an2 ) * t( ppp )

    v( px, py ) = al * am * ( t( pps ) - t( ppp ) )
    v( py, px ) = al * am * ( t( pps ) - t( ppp ) )
!
    v( px, pz ) = al * an * ( t( pps ) - t( ppp ) )
    v( pz, px ) = al * an * ( t( pps ) - t( ppp ) )
!
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
!
    v( s, dyz )   = dsq3 * am * an * t( sds )
    v( dyz, s )   = dsq3 * am * an * t( dss )
!
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

    !========================================================================
 
  END SUBROUTINE koster_slater


  ! intratomic_optics: onsite matrix elements of p operator (or dipole moment) 
  !===========================================================================
  !
  ! NOTE : order of states follows from state.dat
  !
  !    1  2   3   4   5    6    7    8     9     10
  !    s  px  py  pz  s*  dxy  dyz  dzx  dx2y2  dz2r2
  !
  !===========================================================================

  SUBROUTINE intratomic_optics( poldir, v )
    ! output arguments :
    INTEGER :: poldir

    REAL ( dp ), DIMENSION(:,:), ALLOCATABLE :: v

    real(dp), parameter :: P0 = 10.d0

    v = 0.d0
 
    if (poldir.eq.1) THEN

       v( s, px )  = P0   
       v( px, s )  = -P0 
       v( se, px ) = 0
       v( px, se ) = 0

    else if (poldir.eq.2) THEN

       v( s, py )  = P0 
       v( py, s )  = -P0
       v( se, py ) = 0
       v( py, se ) = 0

    elseif (poldir.eq.3) THEN

       v( s, pz )  = P0 
       v( pz, s )  = -P0
       v( se, pz ) = 0
       v( pz, se ) = 0

    end if

  END SUBROUTINE intratomic_optics
  !==========================================================================


  SUBROUTINE spin_sparse( i_spin_a, so, diag, off_diag )

    !=========================================================================
    !
    !  build intra-atomic spin-orbit contribution to matrix <  L | H | L' > ,
    !
    !     L = ( alpha, l ), L' = ( alpha, l' ), for fixed atom alpha
    !
    !=========================================================================
    !
    ! INPUT :
    !
    ! => i_spin_a : integer - l
    ! => so - real array : l-s splitting parameter in eV
    !
    !=========================================================================
    !
    ! OUTPUT :
    !
    ! => diag, off_diag - real array :
    !
    !     complex spin-orbit Hamiltonian submatrix for the l-l' coupling
    !     with l = l' ( diag ) and l /= l' ( off_diag )
    !
    !=========================================================================

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------

    REAL ( dp ), DIMENSION( 2 ), INTENT( IN ) :: so

    INTEGER, INTENT( IN ) :: i_spin_a

    !-------------------------------------------------------------------------
    ! Output arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), DIMENSION(:,:),  INTENT( OUT ) :: diag, off_diag

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------

    COMPLEX ( dp ) :: ci, c1

    REAL ( dp ) :: so_p, so_d

    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------

    c1 = CMPLX( 1.0D0, 0.0D0 )
    ci = CMPLX( 0.0D0, 1.0D0 )

    so_p = so( 1 )
    !so_d = so( 2 )

    diag = CMPLX( 0.0D0, 0.0D0 )           
    off_diag = CMPLX( 0.0D0, 0.0D0 )

    !-------------------------------------------------------------------------
    ! Calculate the spin orbit matrix
    !-------------------------------------------------------------------------

    IF ( i_spin_a .EQ. 1 ) THEN

       !----------------------------------------------------------------------
       ! up - up
       !----------------------------------------------------------------------
 
       diag( px, py )  = -ci * so_p
       diag( py, px )  =  ci * so_p

       !diag( dzx, dyz )   = - ci * so_d
       !diag( dyz, dzx )   =   ci * so_d

       !diag( dx2y2, dxy ) = - 2.0d0 * ci * so_d
       !diag( dxy, dx2y2 ) =   2.0d0 * ci * so_d

       !-------------------------------------------------------------------
       ! up - down
       !-------------------------------------------------------------------

       off_diag( px, pz )  = -c1 * so_p
       off_diag( py, pz )  =  ci * so_p
       off_diag( pz, px )  =  c1 * so_p
       off_diag( pz, py )  = -ci * so_p

       !off_diag( dyz, dxy ) =   c1 * so_d
       !off_diag( dxy, dyz ) = - c1 * so_d
       !off_diag( dxy, dzx ) =   ci * so_d
       !off_diag( dzx, dxy ) = - ci * so_d

       !off_diag( dx2y2, dyz ) = - ci * so_d 
       !off_diag( dx2y2, dzx ) = - c1 * so_d
       !off_diag( dyz, dx2y2 ) =   ci * so_d
       !off_diag( dzx, dx2y2 ) =   c1 * so_d

       !off_diag( dz2r2, dyz ) = - ci * SQRT( 3.0d0 ) * so_d
       !off_diag( dz2r2, dzx ) =   c1 * SQRT( 3.0d0 ) * so_d
       !off_diag( dyz, dz2r2 ) =   ci * SQRT( 3.0d0 ) * so_d
       !off_diag( dzx, dz2r2 ) = - c1 * SQRT( 3.0d0 ) * so_d

    ELSE IF ( i_spin_a .EQ. 2 ) THEN

       !----------------------------------------------------------------------
       ! down - down
       !----------------------------------------------------------------------

       diag( px, py )  =   ci * so_p
       diag( py, px )  =  -ci * so_p

       !diag( dzx, dyz )   =   ci * so_d
       !diag( dyz, dzx )   = - ci * so_d

       !diag( dx2y2, dxy ) =   2.0d0 * ci * so_d
       !diag( dxy, dx2y2 ) = - 2.0d0 * ci * so_d

       !-------------------------------------------------------------------
       ! down - up (Careful! This is  a transpose matrix)
       !-------------------------------------------------------------------

       off_diag( px, pz )  =  c1 * so_p
       off_diag( py, pz )  =  ci * so_p
       off_diag( pz, px )  = -c1 * so_p
       off_diag( pz, py )  = -ci * so_p

       !off_diag( dyz, dxy ) = - c1 * so_d
       !off_diag( dxy, dyz ) =   c1 * so_d
       !off_diag( dxy, dzx ) =   ci * so_d
       !off_diag( dzx, dxy ) = - ci * so_d

       !off_diag( dx2y2, dyz ) = - ci * so_d
       !off_diag( dx2y2, dzx ) =   c1 * so_d
       !off_diag( dyz, dx2y2 ) =   ci * so_d
       !off_diag( dzx, dx2y2 ) = - c1 * so_d

       !off_diag( dz2r2, dyz ) = - ci * SQRT( 3.0d0 ) * so_d
       !off_diag( dz2r2, dzx ) = - c1 * SQRT( 3.0d0 ) * so_d
       !off_diag( dyz, dz2r2 ) =   ci * SQRT( 3.0d0 ) * so_d
       !off_diag( dzx, dz2r2 ) =   c1 * SQRT( 3.0d0 ) * so_d

    END IF

    !=========================================================================

  END SUBROUTINE spin_sparse


  SUBROUTINE rotate_p( p_onsite, c_axis, u, inv_ref, diag )

    !=========================================================================
    !
    !  Rotate p orbital onsite block to calculation coordinate system
    !
    !=========================================================================
    !
    ! INPUT :
    !
    ! => p_onsite : onsite p orbital energies in reference ordering
    ! => c_axis   : orientation of the reference axis in calculation system
    ! => u        : orientation of the reference axis for parametrisation
    ! => inv_ref  : map from reference state indices to active state indices
    !
    !=========================================================================
    !
    ! OUTPUT :
    !
    ! => diag  - complex array :
    !
    !     complex Hamiltonian submatrix
    !
    !=========================================================================

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------

    REAL ( dp ), DIMENSION( 3 ), INTENT( IN ) :: p_onsite
    REAL ( dp ) :: c_axis(3), u(3)
    INTEGER, DIMENSION(:), INTENT(IN) :: inv_ref
    !-------------------------------------------------------------------------
    ! Output arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), DIMENSION( :,: ), POINTER :: diag

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------
    COMPLEX ( dp ), DIMENSION( 3, 3 )  :: Rp, Hp, HR

    INTEGER :: k,l

    ! note the inversion of u and c_axis here w.r.t to d_rotation
    ! here, we use the rotation definition using the complex conjugate of the
    ! Wigner D matrix, with the rotation from u -> c_axis
    call p_rotation(u, c_axis, Rp)

    ! onsite diagonal hamiltonian block with order px, py, pz
    Hp(:,:) = CMPLX(0.d0, 0.d0)
    Hp(1,1) = CMPLX(p_onsite(1), 0.d0);
    Hp(2,2) = CMPLX(p_onsite(2), 0.d0);
    Hp(3,3) = CMPLX(p_onsite(3), 0.d0);
    
    call H_rotation(Hp,Rp,HR)
    !write(*,*) HR(1,1), " ", HR(2,2), " ", HR(3,3)

    if (inv_ref(px) .gt. 0) then
      diag( inv_ref(px), inv_ref(px) ) = HR(1, 1)
      if (inv_ref(py) .gt. 0) diag( inv_ref(px), inv_ref(py) ) = HR(1, 2)
      if (inv_ref(pz) .gt. 0) diag( inv_ref(px), inv_ref(pz) ) = HR(1, 3)
    end if

    if (inv_ref(py) .gt. 0) then
      if (inv_ref(px) .gt. 0) diag( inv_ref(py), inv_ref(px) ) = HR(2, 1)
      diag( inv_ref(py), inv_ref(py) ) = HR(2, 2)
      if (inv_ref(pz) .gt. 0) diag( inv_ref(py), inv_ref(pz) ) = HR(2, 3)
    end if

    if (inv_ref(pz) .gt. 0) then
      if (inv_ref(px) .gt. 0) diag( inv_ref(pz), inv_ref(px) ) = HR(3, 1)
      if (inv_ref(py) .gt. 0) diag( inv_ref(pz), inv_ref(py) ) = HR(3, 2)
      diag( inv_ref(pz), inv_ref(pz) ) = HR(3, 3)
    end if

   ! Remove almost 0 elements.....
    do k=1,size(diag,1)
       do l=1,size(diag,1)

          IF ( ABS(diag(k,l)).LT.1.0d-14 ) diag(k,l)=(0.d0,0.d0)

       end do
    end do

  END SUBROUTINE rotate_p



  SUBROUTINE rotate_d( d_onsite, c_axis, u, inv_ref, diag )

    !=========================================================================
    !
    !  build onsite contributions to the crystal field for the d-orbitals.
    !  The matrix < R_i L M | H | R_i L M' > , L = 2
    !
    !     for fixed atom R_i
    !
    !=========================================================================
    !
    ! INPUT :
    !
    ! => d_onsite : onsite d orbital energies in reference ordering
    ! => c_axis   : orientation of the reference axis in calculation system
    ! => u        : orientation of the reference axis for parametrisation
    ! => inv_ref  : map from reference state indices to active state indices
    !
    !=========================================================================
    !
    ! OUTPUT :
    !
    ! => diag  - complex array :
    !
    !     complex CF Hamiltonian submatrix 
    !
    !=========================================================================

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------

    REAL ( dp ), DIMENSION( 5 ), INTENT( IN ) :: d_onsite
    REAL ( dp ) :: c_axis(3), u(3)
    INTEGER, DIMENSION(:), INTENT(IN) :: inv_ref
    !-------------------------------------------------------------------------
    ! Output arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), DIMENSION( :,: ), POINTER :: diag

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------
    COMPLEX ( dp ), DIMENSION( 5, 5 )  :: Rd, Hd, HR

    COMPLEX ( dp ) :: j

    INTEGER :: k,l, perm(n_ref_states)

    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------

    j = CMPLX( 0.0D0, 1.0D0 )

    !-------------------------------------------------------------------------
    ! Calculate the crystal field matrix
    !-------------------------------------------------------------------------
    
    !  Natural d-state ordering:   |dxy>   |dyz>   |dz2>   |dxz>    |dx2y2>
 
    ! Get the represention for the rotation of the c-axis along <111>
    call d_rotation(c_axis, u, Rd)

    ! Obtain rotation into the tesseral states
    call tesseral_trans(Rd)
 

    ! Initialize the d-orbital diagonal matrix in the natural d-state tesseral
    ! basis set
    Hd(:,:) = CMPLX(0.d0, 0.d0)
    Hd(1,1) = CMPLX(d_onsite(1), 0.d0);
    Hd(2,2) = CMPLX(d_onsite(2), 0.d0);
    Hd(3,3) = CMPLX(d_onsite(5), 0.d0);
    Hd(4,4) = CMPLX(d_onsite(3), 0.d0);
    Hd(5,5) = CMPLX(d_onsite(4), 0.d0);
    ! 
    ! 
    call H_rotation(Hd,Rd,HR)

    !
    !
    !  Uptight ordering:           |dxy>   |dyz>   |dxz>   |dx2y2>  |dz2> 
    !
    !  ==>  we need a permutation array mapping the tesseral ordering
    !       into the UPTIGHT basis 
    !       
    !       example:
    !       
    !       dz2==5  =>  perm(5)=3  means that in the tesseral basis 
    !                              dz2 is the element 3  
    !
    perm(dxy)=1; perm(dyz)=2; perm(dz2r2)=3; perm(dzx)=4; perm(dx2y2)=5;

    ! TRANSFER TO MATRIX BLOCK:
    diag( inv_ref(dxy), inv_ref(dxy) )   = HR( perm(dxy) , perm(dxy) )
    diag( inv_ref(dxy), inv_ref(dyz) )   = HR( perm(dxy) , perm(dyz) ) 
    diag( inv_ref(dxy), inv_ref(dzx) )   = HR( perm(dxy) , perm(dzx) )
    diag( inv_ref(dxy), inv_ref(dx2y2) ) = HR( perm(dxy) , perm(dx2y2) )
    diag( inv_ref(dxy), inv_ref(dz2r2) ) = HR( perm(dxy) , perm(dz2r2) )
 
    diag( inv_ref(dyz), inv_ref(dxy) )   = diag( inv_ref(dxy), inv_ref(dyz) )
    diag( inv_ref(dyz), inv_ref(dyz) )   = HR( perm(dyz) , perm(dyz) ) 
    diag( inv_ref(dyz), inv_ref(dzx) )   = HR( perm(dyz) , perm(dzx) )
    diag( inv_ref(dyz), inv_ref(dx2y2) ) = HR( perm(dyz) , perm(dx2y2) )
    diag( inv_ref(dyz), inv_ref(dz2r2) ) = HR( perm(dyz) , perm(dz2r2) )

    diag( inv_ref(dzx), inv_ref(dxy) )   = diag( inv_ref(dxy), inv_ref(dzx) )
    diag( inv_ref(dzx), inv_ref(dyz) )   = diag( inv_ref(dyz), inv_ref(dzx) )
    diag( inv_ref(dzx), inv_ref(dzx) )   = HR( perm(dzx) , perm(dzx) )
    diag( inv_ref(dzx), inv_ref(dx2y2) ) = HR( perm(dzx) , perm(dx2y2) )
    diag( inv_ref(dzx), inv_ref(dz2r2) ) = HR( perm(dzx) , perm(dz2r2) )

    diag( inv_ref(dx2y2), inv_ref(dxy) )   = diag( inv_ref(dxy), inv_ref(dx2y2) )
    diag( inv_ref(dx2y2), inv_ref(dyz) )   = diag( inv_ref(dyz), inv_ref(dx2y2) )
    diag( inv_ref(dx2y2), inv_ref(dzx) )   = diag( inv_ref(dzx), inv_ref(dx2y2) )
    diag( inv_ref(dx2y2), inv_ref(dx2y2) ) = HR( perm(dx2y2) , perm(dx2y2) )
    diag( inv_ref(dx2y2), inv_ref(dz2r2) ) = HR( perm(dx2y2) , perm(dz2r2) )

    diag( inv_ref(dz2r2), inv_ref(dxy) )   = diag( inv_ref(dxy), inv_ref(dz2r2) )
    diag( inv_ref(dz2r2), inv_ref(dyz) )   = diag( inv_ref(dyz), inv_ref(dz2r2) )
    diag( inv_ref(dz2r2), inv_ref(dzx) )   = diag( inv_ref(dzx), inv_ref(dz2r2) )
    diag( inv_ref(dz2r2), inv_ref(dx2y2) ) = diag( inv_ref(dx2y2), inv_ref(dz2r2) )
    diag( inv_ref(dz2r2), inv_ref(dz2r2) ) = HR( perm(dz2r2) , perm(dz2r2) )

    ! Remove almost 0 elements.....
    do k=1,5
       do l=1,5
          
          IF ( ABS(diag(k,l)).LT.1.0d-14 ) diag(k,l)=(0.d0,0.d0)

       end do
    end do

  END SUBROUTINE rotate_d


  !--------------------------------------------------------------------------
  ! Add strain degeneracy lifting for zincblende
  ! Jancu scheme only
  !--------------------------------------------------------------------------
  SUBROUTINE crystal_field_zb(onsite, ion, i_ion, c_axis, strain, inv_ref )
    
    REAL(dp), DIMENSION(:) :: onsite
    TYPE(ion_orbit), DIMENSION(:), POINTER :: ion
    INTEGER :: i_ion
    REAL(dp), DIMENSION( 3 ) :: c_axis
    REAL(dp), DIMENSION( 3 ) :: strain
    INTEGER, DIMENSION(:), INTENT(IN) :: inv_ref

    !INTEGER :: i_state
    REAL(dp) :: d_e, b_d

    b_d = ion(i_ion)%b_d

    ! will now be rotated afterwards
    !if(c_axis(3).eq.1.d0) then
       d_e = strain(3) - (strain(1)+strain(2))*0.5d0
       onsite(inv_ref(dxy)) = onsite(inv_ref(dxy)) * (1 + 2 * b_d * d_e )
       onsite(inv_ref(dyz)) = onsite(inv_ref(dyz)) * (1 -  b_d * d_e )
       onsite(inv_ref(dzx)) = onsite(inv_ref(dzx)) * (1 -  b_d * d_e )

    !else if(c_axis(1).eq.1.d0) then
    !   d_e = strain(1) - (strain(2)+strain(3))*0.5d0
    !   onsite(dyz) = onsite(dyz) * (1 + 2 * b_d * d_e )
    !   onsite(dxy) = onsite(dxy) * (1 - b_d * d_e )
    !   onsite(dzx) = onsite(dzx) * (1 - b_d * d_e )

    !else if(c_axis(2).eq.1.d0) then

    !   d_e = strain(2) - (strain(1)+strain(3))*0.5d0
    !   onsite(dzx) = onsite(dzx) * (1 + 2 * b_d * d_e )
    !   onsite(dyz) = onsite(dyz) * (1 - b_d * d_e )
    !   onsite(dxy) = onsite(dxy) * (1 - b_d * d_e )     

    !endif


  END SUBROUTINE crystal_field_zb


  !--------------------------------------------------------------------------
  ! Add strain degeneracy lifting for zincblende (version of Zielinski2012)
  ! Jancu scheme only !NOT YET FINISHED, DON'T USE!
  !--------------------------------------------------------------------------
  SUBROUTINE crystal_field_zb2(onsite, ion, i_ion, v, v0, &
                               d, d0, n_bond, inv_ref )
    
    REAL(dp), DIMENSION(:) :: onsite
    TYPE(ion_orbit), DIMENSION(:), POINTER :: ion
    INTEGER :: i_ion, n_bond
    REAL(dp), DIMENSION(:,:) :: v !list of cos_latt
    REAL(dp), DIMENSION(:,:) :: v0 !list of normal cos_latt
    REAL(dp), DIMENSION(:) :: d !list of bondlength
    REAL(dp), DIMENSION(:) :: d0 !list of normal bondlength
    INTEGER, DIMENSION(:), INTENT(IN) :: inv_ref

    !local
    REAL(dp) :: split
    INTEGER :: i, N

    split = 0.d0
    N = 0
    DO i = 1, SIZE(d), n_bond+1
      split = split + (d(i)/d0(i)) * &
              ( v(i,3)/v0(i,3) - 0.5d0*(v(i,1)/v0(i,1) + v(i,2)/v0(i,2)) )
      N = N + 1
    ENDDO
    split = ion(i_ion)%b_d/N*split

    onsite(inv_ref(dxy)) = onsite(inv_ref(dxy)) + 2*split
    onsite(inv_ref(dyz)) = onsite(inv_ref(dyz)) - split
    onsite(inv_ref(dzx)) = onsite(inv_ref(dzx)) - split

  END SUBROUTINE crystal_field_zb2


  !===========================================================================
  !
  ! Subroutine "hydrogen_coupling"
  !
  !=========================================================================

  SUBROUTINE hydrogen_coupling( t )

    !=========================================================================

    ! input - output arguments :

    REAL ( dp ), DIMENSION( : ), INTENT( INOUT ) :: t

    !_________________________________________________________________________

    ! local variables :

    REAL ( dp ), DIMENSION( SIZE( t ) ) :: t_aux

    !=========================================================================

    t_aux = t

    t = 0.0D0

    t( sss )  = t_aux( sss )
    t( sps )  = t_aux( sps )
    t( pss )  = t_aux( pss )
    t( sses ) = t_aux( sses )
    t( sess ) = t_aux( sess )
    t( sds )  = t_aux( sds )
    t( dss )  = t_aux( dss )

    !=========================================================================

  END SUBROUTINE hydrogen_coupling

  !===========================================================================
  SUBROUTINE Klimeck_passiv(mat,i_dg_a)

    complex(dp), DIMENSION(:,:), POINTER :: mat
    integer :: i_dg_a

    real(dp), PARAMETER :: dangling_shift = 20.d0
    
    SELECT CASE(i_dg_a)
    CASE(1)
       mat(1:4,1:4) = mat(1:4,1:4) + dangling_shift/4.d0
    CASE(2)
    
       mat(s, s) = mat(s, s) + dangling_shift/2.d0
       !Putting this term on px, py or pz is the same
       mat(s, py) = mat(s, py) + dangling_shift/2.d0
       mat(py, s) = mat(py, s) + dangling_shift/2.d0
       mat(px, px) = mat(px, px) + dangling_shift/2.d0
       mat(px, pz) = mat(px, pz) + dangling_shift/2.d0
       mat(pz, px) = mat(pz, px) + dangling_shift/2.d0
       mat(py, py) = mat(py, py) + dangling_shift/2.d0
       mat(pz, pz) = mat(py, py) + dangling_shift/2.d0
       
    CASE(3)
              
       mat(s, s) = mat(s, s) + dangling_shift*3.d0/4.d0
       mat(s, px) = mat(s, px) + dangling_shift/4.d0
       mat(px, s) = mat(px, s) + dangling_shift/4.d0
       mat(s, py) = mat(s, py) + dangling_shift/4.d0
       mat(py, s) = mat(py, s) + dangling_shift/4.d0
       mat(s, pz) = mat(s, pz) + dangling_shift/(-4.d0)
       mat(pz, s) = mat(pz, s) + dangling_shift/(-4.d0)
       mat(px, px) = mat(px, px) + dangling_shift*3.d0/4.d0
       mat(px, py) = mat(px, py) + dangling_shift/(-4.d0)
       mat(py, px) = mat(py, px) + dangling_shift/(-4.d0)
       mat(px, pz) = mat(px, pz) + dangling_shift/4.d0
       mat(pz, px) = mat(pz, px) + dangling_shift/4.d0
       mat(py, py) = mat(py, py) + dangling_shift*3.d0/4.d0
       mat(py, pz) = mat(py, pz) + dangling_shift/4.d0
       mat(pz, py) = mat(pz, py) + dangling_shift/4.d0
       mat(pz, pz) = mat(py, py) + dangling_shift*3.d0/4.d0
       
    END SELECT

  END SUBROUTINE KLIMECK_PASSIV
  
!!$  SUBROUTINE TB_sparse_davidson( sparse_ind, sparse_mat, n_mat )
!!$    
!!$    !-------------------------------------------------------------------------
!!$    ! Input arguments
!!$    !----------------------------------------------------------------------
!!$
!!$    COMPLEX ( dp , DIMENSION( : ), INTENT( IN ) :: sparse_mat
!!$    INTEGER                  , DIMENSION( : ), INTENT( IN ) :: sparse_ind
!!$    INTEGER                                  , INTENT( IN ) :: n_mat
!!$
!!$    !-------------------------------------------------------------------------
!!$    ! Local variables
!!$    !-------------------------------------------------------------------------
!!$    
!!$    INTEGER :: i_ind, n_ind, i_val, i_row, i_col, i_targ
!!$    INTEGER :: N, dim
!!$
!!$    !-------------------------------------------------------------------------
!!$    ! Allocate new sparse matrixes
!!$    !-------------------------------------------------------------------------
!!$    
!!$    CALL setup_sparse_matrix( n_mat, SIZE( sparse_mat ) - n_mat )
!!$    
!!$    N = Mij( 1 ) - 2
!!$    dim = SIZE( Mij )
!!$    
!!$    !-------------------------------------------------------------------------
!!$    ! Initialization
!!$    !-------------------------------------------------------------------------
!!$    
!!$    i_val = 0
!!$    i_col = 0
!!$    i_targ = N + 2
!!$    
!!$    !-------------------------------------------------------------------------
!!$    ! Loop on sparse indexes
!!$    !-------------------------------------------------------------------------
!!$    
!!$    DO i_ind = 1, SIZE( sparse_ind )
!!$       
!!$       !----------------------------------------------------------------------
!!$       ! Test on index
!!$       !----------------------------------------------------------------------
!!$       
!!$       IF ( sparse_ind( i_ind ) .LT. 0 ) THEN
!!$          
!!$          !-------------------------------------------------------------------
!!$          ! Row index
!!$          !-------------------------------------------------------------------
!!$          
!!$          i_row = - sparse_ind( i_ind )
!!$          
!!$          !-------------------------------------------------------------------
!!$          ! Reset column counter
!!$          !-------------------------------------------------------------------
!!$          
!!$          Mij( i_row ) = i_targ
!!$          i_col = 0
!!$          
!!$          CYCLE
!!$          
!!$       ELSE
!!$          
!!$          !-------------------------------------------------------------------
!!$          ! Column index
!!$          !-------------------------------------------------------------------
!!$          
!!$          i_col = sparse_ind( i_ind )
!!$          
!!$          !-------------------------------------------------------------------
!!$          ! Update sparse element counter
!!$          !-------------------------------------------------------------------
!!$          
!!$          i_val = i_val + 1
!!$          
!!$          !-------------------------------------------------------------------
!!$          ! End test on sparse indexes
!!$          !-------------------------------------------------------------------
!!$          
!!$       END IF
!!$       
!!$       IF ( i_row .EQ. i_col ) THEN
!!$          
!!$          !------ diagonal element ----------
!!$          
!!$          M( i_row ) = sparse_mat( i_val )
!!$          
!!$       ELSE
!!$          
!!$          !------ offdiagonal element -------
!!$                       
!!$          Mij( i_targ ) = i_col
!!$          M( i_targ ) = sparse_mat( i_val )
!!$          i_targ = i_targ + 1
!!$                 
!!$       END IF
!!$              
!!$       !----------------------------------------------------------------------
!!$       ! End loop on sparse indexes
!!$       !----------------------------------------------------------------------
!!$
!!$    END DO
!!$    
!!$    Mij( N + 1 ) = i_targ
!!$
!!$    !=========================================================================
!!$
!!$  END SUBROUTINE TB_sparse_davidson


  !===========================================================================
  !
  ! Subroutine "store_ind" : store sparse index.
  !
  !=========================================================================
  
  SUBROUTINE store_ind_file( file_num, ind, i_ind, n_ind, ind_tmp )

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------

    INTEGER, INTENT( IN ) :: ind, file_num

    !-------------------------------------------------------------------------
    ! Input - output arguments
    !-------------------------------------------------------------------------

    INTEGER, DIMENSION( max_ind_D ), INTENT( INOUT ) :: ind_tmp
    
    INTEGER, INTENT( INOUT ) :: n_ind, i_ind

    !-------------------------------------------------------------------------
    ! Update counters
    !-------------------------------------------------------------------------

    i_ind = i_ind + 1
    ind_tmp( i_ind ) = ind
    
    !-------------------------------------------------------------------------
    ! If the indexes array is full, dump it in indexes file
    !-------------------------------------------------------------------------

    IF ( i_ind .EQ. max_ind_D ) THEN 
       WRITE ( file_num ) ind_tmp
       n_ind = n_ind + i_ind
       i_ind = 0
       
    END IF

    !=========================================================================

  END SUBROUTINE store_ind_file

  !=========================================================================
  
  SUBROUTINE store_ind_mem( ind, i_ind, n_ind, ind_tmp )

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------
    INTEGER, INTENT( IN ) :: ind

    INTEGER, INTENT( INOUT ) :: n_ind, i_ind

    !-------------------------------------------------------------------------
    ! Input - output arguments
    !-------------------------------------------------------------------------

    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_tmp
    

    !-------------------------------------------------------------------------
    ! Update counters
    !-------------------------------------------------------------------------

    i_ind = i_ind + 1
    ind_tmp( i_ind ) = ind
    
    !=========================================================================

  END SUBROUTINE store_ind_mem



  !===========================================================================
  !
  ! Subroutine "store_val" : store sparse index.
  !
  !=========================================================================
  
  SUBROUTINE store_val_file( file_num, value, i_val, n_val, val_tmp )

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), INTENT( IN ) :: value

    INTEGER,                 INTENT( IN ) :: file_num

    !-------------------------------------------------------------------------
    ! Input - output arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), DIMENSION( max_ind_D ), INTENT( INOUT ) :: &
         val_tmp
    
    INTEGER, INTENT( INOUT ) :: n_val, i_val

    !-------------------------------------------------------------------------
    ! Update counters
    !-------------------------------------------------------------------------

    i_val = i_val + 1
    val_tmp( i_val ) = value

    !-------------------------------------------------------------------------
    ! If the indexes array is full, dump it in indexes file
    !-------------------------------------------------------------------------

    IF ( i_val .EQ. max_ind_D ) THEN 
       WRITE ( file_num ) val_tmp
       n_val = n_val + i_val
       i_val = 0
    END IF

    !=========================================================================

  END SUBROUTINE store_val_file
  !=========================================================================
  
  SUBROUTINE store_val_mem( value, i_val, n_val, val_tmp )

    !-------------------------------------------------------------------------
    ! Input arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), INTENT( IN ) :: value
    INTEGER, INTENT( INOUT ) :: n_val, i_val

    !-------------------------------------------------------------------------
    ! Input - output arguments
    !-------------------------------------------------------------------------

    COMPLEX ( dp ), DIMENSION(:), ALLOCATABLE  ::  val_tmp
   

    !-------------------------------------------------------------------------
    ! Update counters
    !-------------------------------------------------------------------------

    i_val = i_val + 1
    val_tmp( i_val ) = value

  END SUBROUTINE store_val_mem


  !===========================================================================

  LOGICAL FUNCTION d_orbital_present(mat_data,i_ion)
    
    TYPE (material_data) :: mat_data
    integer :: i_ion, i_state
    character(STATELEN) :: state
    
    d_orbital_present = .false.
    
    
    DO i_state = 1, SIZE( mat_data%ion( i_ion )%energy )
     
       state = trim(mat_data%ion( i_ion )%state( i_state ))

       SELECT CASE ( state )
          
       CASE( 'dxy', 'dyz', 'dzx', 'dx2y2', 'dz2r2', 'd15', 'd' )
          
          d_orbital_present= .true.
          
       END SELECT
      
    END DO

  END FUNCTION d_orbital_present

  LOGICAL FUNCTION p_orbital_present(mat_data,i_ion)

    TYPE (material_data) :: mat_data
    integer :: i_ion, i_state
    character(STATELEN) :: state

    p_orbital_present = .false.


    DO i_state = 1, SIZE( mat_data%ion( i_ion )%energy )

       state = trim(mat_data%ion( i_ion )%state( i_state ))

       SELECT CASE ( state )

       CASE( 'px', 'py', 'pz', 'p' )

          p_orbital_present= .true.

       END SELECT

    END DO

  END FUNCTION p_orbital_present



!*****************************************************************************
  !=========================================================================
  ! Force to Hermitian by taking an average between Lower and Upper
  !=========================================================================     
  SUBROUTINE hermitianize_csr(ham)

    TYPE(CSR) :: ham

    ! Local variables:
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mi
    INTEGER,        DIMENSION(:), POINTER :: Mj
    CHARACTER(1)                  :: sparse_fmt

    INTEGER  :: n_ham
    INTEGER :: row, p, col
    COMPLEX( dp ) :: matel1, matel2, matav
     
    M => ham%M
    Mi => ham%Mi
    Mj => ham%Mj
    n_ham = ham%nrow

    SELECT CASE(ham%sparse_fmt)

    CASE('F','f')
    
       do row=1,n_ham
 
          do p = Mi(row), Mi(row+1) - 1 
             
             col = Mj(p)

             matel1 = sprs_element(M, Mj, Mi,row,col)
             matel2 = sprs_element(M, Mj, Mi,col,row)
            
             ! Either matel1 or matel2 == 0 
             if (abs(matel1).eq.0.d0 .neqv. abs(matel2).eq.0.d0) then

                !print*, row, col    
                !print*, matel1, matel2     
                matav=(0.d0,0.d0)     
                call set_sprs_element(M, Mj, Mi, row, col, matav)
                call set_sprs_element(M, Mj, Mi, col, row, matav)
                        
             ! matel1 .ne. matel2 
             elseif (abs(matel1-conjg(matel2)).gt.1d2*emach ) then
                     
                matav = (matel1 + matel2) / 2.d0

                !write(100,*) row, col, matel1,matel2    
                call set_sprs_element(M, Mj, Mi, row, col, matav)
                call set_sprs_element(M, Mj, Mi, col, row, matav)

             end if
             
          enddo
       enddo

    END SELECT


  END SUBROUTINE hermitianize_csr

  !=========================================================================
  ! Force to Hermitian by taking an average between Lower and Upper
  !=========================================================================     
  SUBROUTINE hermitianize_ex(ham)

    TYPE(CSR_ex) :: ham

    ! Local variables:
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mij
    CHARACTER(1)                  :: sparse_fmt

    INTEGER  :: n_ham
    INTEGER :: row, p, col
    COMPLEX( dp ) :: matel1, matel2, matav
     
    M => ham%M
    Mij => ham%Mij
    n_ham = ham%nrow

    SELECT CASE(ham%sparse_fmt)

    CASE('F','f')
    
       do row=1,n_ham
 
          do p = Mij(row), Mij(row+1) - 1 
             
             col = Mij(p)

             matel1 = sprs_element(M, Mij, row, col)
             matel2 = sprs_element(M, Mij, col, row)
             
             if( abs(matel1-conjg(matel2)).gt.1d2*emach ) then
                      
                matav = (matel1 + matel2) / 2.d0

                call set_sprs_element(M, Mij, row, col, matav)
                call set_sprs_element(M, Mij, col, row, matav)

             end if
             
          enddo
       enddo

    END SELECT


  END SUBROUTINE hermitianize_ex
  !=========================================================================
  ! Hermeticity Checks
  !=========================================================================     
  SUBROUTINE check_if_hermitian_csr(ham)

    TYPE(CSR) :: ham

    ! Local variables:
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mi
    INTEGER,        DIMENSION(:), POINTER :: Mj
    CHARACTER(1)                  :: sparse_fmt
    INTEGER  :: n_ham
 
    INTEGER :: row, p, col, file_num
    COMPLEX( dp ) :: matel1, matel2
    
    M => ham%M
    Mi => ham%Mi
    Mj => ham%Mj
    n_ham = ham%nrow
    
    CALL open_file( "herm_check.dat", file_num, "write", format_flag = .TRUE., &
         output_flag = .FALSE. ) 

    SELECT CASE(ham%sparse_fmt)
       !CONVERT MSR into full MSR matrix
    CASE('F','f')
    
       do row=1,n_ham
 
          do p = Mi(row), Mi(row+1) - 1 !row+1,n_ham
             
             col = Mj(p)

             matel1 = sprs_element(M, Mj, Mi, row, col)
             matel2 = sprs_element(M, Mj, Mi, col, row)
             
             if( abs(matel1-conjg(matel2)).gt.1d2*emach ) then
                      
                write(file_num,*) row,col,matel1                
                write(file_num,*) col,row,matel2
                write(file_num,*)


                !if(abs(matel1).lt.1.d-16) then
                !   call set_sprs_element(M, Mij,col,row,conjg(matel1))                   
                !elseif(abs(matel2).lt.1.d-16) then
                !   call set_sprs_element(M, Mij,row,col,conjg(matel2))                   
                !else       
                !   ! UPPER MATRIX IS COPIED ON LOWER MATRIX
                !   call set_sprs_element(M, Mij,col,row,conjg(matel1))
                !endif
                                
             end if
             
          enddo
       enddo

    END SELECT


    close(file_num)
   
    !WRITE ( * , * ) 'converting L/U sparse matrix to FULL sparse matrix'

    !call msr2msr_full
  END SUBROUTINE check_if_hermitian_csr

  SUBROUTINE check_if_hermitian_ex(ham)

    TYPE(CSR_ex) :: ham

    ! Local variables:
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mij
    CHARACTER(1)                  :: sparse_fmt
    INTEGER  :: n_ham
 
    INTEGER :: row, p, col, file_num
    COMPLEX( dp ) :: matel1, matel2
    
    M => ham%M
    Mij => ham%Mij
    n_ham = ham%nrow
    
    CALL open_file( "herm_check.dat", file_num, "write", format_flag = .TRUE., &
         output_flag = .FALSE. ) 

    SELECT CASE(ham%sparse_fmt)
       !CONVERT MSR into full MSR matrix
    CASE('F','f')
    
       do row=1,n_ham
 
          do p = Mij(row), Mij(row+1) - 1 !row+1,n_ham
             
             col = Mij(p)

             matel1 = sprs_element(M, Mij, row, col)
             matel2 = sprs_element(M, Mij, col, row)
             
             if( abs(matel1-conjg(matel2)).gt.1d2*emach ) then
                      
                write(file_num,*) row,col,matel1                
                write(file_num,*) col,row,matel2
                write(file_num,*)


                !if(abs(matel1).lt.1.d-16) then
                !   call set_sprs_element(M, Mij,col,row,conjg(matel1))                   
                !elseif(abs(matel2).lt.1.d-16) then
                !   call set_sprs_element(M, Mij,row,col,conjg(matel2))                   
                !else       
                !   ! UPPER MATRIX IS COPIED ON LOWER MATRIX
                !   call set_sprs_element(M, Mij,col,row,conjg(matel1))
                !endif
                                
             end if
             
          enddo
       enddo

    END SELECT


    close(file_num)
   
    !WRITE ( * , * ) 'converting L/U sparse matrix to FULL sparse matrix'

    !call msr2msr_full
  END SUBROUTINE check_if_hermitian_ex


!*****************************************************************************
  !=========================================================================
  ! Anti-Hermeticity Checks
  !=========================================================================     
  SUBROUTINE check_if_antihermitian(ham)

    TYPE(CSR) :: ham

    ! Local variables:
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER,        DIMENSION(:), POINTER :: Mj
    INTEGER,        DIMENSION(:), POINTER :: Mi
    CHARACTER(1)                  :: sparse_fmt
    INTEGER  :: n_ham
 
    INTEGER :: row, p, col
    COMPLEX( dp ) :: matel1, matel2
    
    M => ham%M
    Mi => ham%Mi
    Mj => ham%Mj
    n_ham = ham%nrow
    
    SELECT CASE(ham%sparse_fmt)
       !CONVERT MSR into full MSR matrix
    CASE('F','f')
    
       do row=1,n_ham
 
          do p = Mi(row), Mi(row+1) - 1 !row+1,n_ham
             
             col = Mj(p)

             matel1 = sprs_element(M, Mj, Mi, row, col)
             matel2 = sprs_element(M, Mj, Mi, col, row)
             
             if( abs(matel1+conjg(matel2)).gt.1d2*emach ) then
                      
                write(*,*) row,col,matel1                
                write(*,*) col,row,matel2
                write(*,*)


                !if(abs(matel1).lt.1.d-16) then
                !   call set_sprs_element(M, Mij,col,row,conjg(matel1))                   
                !elseif(abs(matel2).lt.1.d-16) then
                !   call set_sprs_element(M, Mij,row,col,conjg(matel2))                   
                !else       
                !   ! UPPER MATRIX IS COPIED ON LOWER MATRIX
                !   call set_sprs_element(M, Mij,col,row,conjg(matel1))
                !endif
                                
             end if
             
          enddo
       enddo

    END SELECT

   
    !WRITE ( * , * ) 'converting L/U sparse matrix to FULL sparse matrix'

    !call msr2msr_full
  END SUBROUTINE check_if_antihermitian
!*****************************************************************************
  !=========================================================================
  ! Anti-Hermeticity Checks
  !=========================================================================     
  SUBROUTINE check_if_real(ham)

    TYPE(CSR) :: ham

    ! Local variables:
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    INTEGER  :: n_ham
 
    INTEGER :: row, p, col
    COMPLEX( dp ) :: matel1, matel2
    
    M => ham%M
    n_ham = size(ham%M)
    
    DO p=1,n_ham
       IF (ABS(AIMAG(M(p))).gt.0.d0) THEN
          WRITE(*,*) 'IMAGINARY ELEMENT: ', p
       ENDIF
    END DO

  END SUBROUTINE check_if_real



  SUBROUTINE get_ion_orbitals(upt, i_a, orbital_ids)
    TYPE(OUPT), TARGET :: upt
    integer :: orbital_ids(*)


    TYPE (ion_basis), POINTER :: basis  
    TYPE (material_data), DIMENSION( : ), POINTER :: mat_data 
    integer :: i_a, i_a_mat, i_a_ion

    basis => upt%basis
    mat_data => upt%materials   

    CALL check_ion( basis, i_a, mat_data, i_a_mat, i_a_ion )
    
!write(*,*) size( orbital_ids ), size( mat_data( i_a_mat )%ion( i_a_ion )%ind_ref )

DO i_a = 1, size( mat_data( i_a_mat )%ion( i_a_ion )%ind_ref)
    orbital_ids(i_a) = mat_data( i_a_mat )%ion( i_a_ion )%ind_ref(i_a)
END DO


  END SUBROUTINE get_ion_orbitals
  ! ------------------------------------------------------------------------


  SUBROUTINE get_ion_block_size(upt, ion_block_size)
    TYPE(OUPT), TARGET :: upt
    integer :: ion_block_size(*)


    TYPE (ion_basis), POINTER :: basis  
    TYPE (material_data), DIMENSION( : ), POINTER :: mat_data 
    integer :: i_a, i_a_mat, i_a_ion, n_basis, i_n, n_spin

    basis => upt%basis
    mat_data => upt%materials   
    n_basis = upt%basis%n_basis
    n_spin = upt%n_spin

    DO i_a = 1, n_basis

       CALL check_ion( basis, i_a, mat_data, i_a_mat, i_a_ion )
       
       ion_block_size(i_a) = n_spin * SIZE(mat_data( i_a_mat )%ion( i_a_ion )%ind_ref)

    END DO


  END SUBROUTINE get_ion_block_size
  ! ------------------------------------------------------------------------
  
  SUBROUTINE get_local_strain(basis,i_a,near,parent,strain)

    TYPE(ion_basis), intent(in) :: basis
    INTEGER, intent(in) :: i_a
    TYPE(nearest_neighbours), intent(in) :: near
    TYPE(parent_data), DIMENSION(:), POINTER :: parent
    REAL(dp), DIMENSION(3), intent(out) :: strain

    ! locals.
    integer :: i_n, i_b, l, m    
    REAL(dp), DIMENSION(3) :: D_bas, D_latt, l0
    REAL(dp), DIMENSION(3,4) :: D_nn, D_tmp
    REAL(dp) :: dist_ref1, dist_ref2, dist_ref
    REAL(dp), DIMENSION(3,3) :: K, inv_K

    !dist_ref1 = parent(1)%pair(1)%dist_ref
    !dist_ref2 = parent(2)%pair(1)%dist_ref 
    !dist_ref = parent(1)%content*dist_ref1 + parent(2)%content*dist_ref2 
    
    dist_ref = parent(1)%pair(1)%dist_ref

    strain = 0.d0

    l = 0
    DO i_n = 1, SIZE( near%ind ) 

       IF(near%order(i_n) .ne. 0) THEN
          l = l+1
          i_b = near%ind(i_n)
          D_latt = MATMUL( basis%prim, near%vec( :, i_n ) )
          D_bas = D_basis( basis, i_b, i_a )
          D_nn(:,l) = D_latt + D_bas           
       ENDIF

    END DO

    D_tmp(:,1) = D_nn(:,1)-D_nn(:,2)
    D_tmp(:,2) = D_nn(:,3)-D_nn(:,4)
    D_tmp(:,3) = D_nn(:,3)-D_nn(:,2)    

    !print *, D_tmp(:,1)
    !print *, sqrt(dot_product(D_tmp(:,1),D_tmp(:,1)))
    !print *, D_tmp(:,2)
    !print *, sqrt(dot_product(D_tmp(:,2),D_tmp(:,2)))
    !print *, D_tmp(:,3)   
    !print *, sqrt(dot_product(D_tmp(:,3),D_tmp(:,3)))

    DO l=1,3
       DO m=1,3
          K(l,m) = D_tmp(m,l)*D_tmp(m,l)
       ENDDO
    ENDDO

    l0 = dist_ref*dist_ref*8.d0/3.d0
    print *,'ref distance=',sqrt(l0(1))


    print *, 'invert:'    
    print *, K(1,:)
    print *, K(2,:)    
    print *, K(3,:)

    call inverse(K, inv_K)

    print *, 'inverted:'    
    print *, inv_K(1,:)
    print *, inv_K(2,:)    
    print *, inv_K(3,:)

    strain = matmul(inv_K,l0)

    print*,'strain:'
    print*,strain    


  END SUBROUTINE get_local_strain
   
  subroutine inverse(K, inv_K)
    REAL(dp), DIMENSION(3,3) :: K, inv_K    
    REAL(dp) :: det

    det =+K(1,1)*(K(2,2)*K(3,3)-K(3,2)*K(2,3)) &
         -K(1,2)*(K(2,1)*K(3,3)-K(3,1)*K(2,3)) &
         +K(1,3)*(K(2,1)*K(3,2)-K(3,1)*K(2,2))

    print*,'det=',det

    inv_K(1,1) =  (K(2,2)*K(3,3)-K(3,2)*K(2,3))/det
    inv_K(1,2) = -(K(1,2)*K(3,3)-K(3,2)*K(1,3))/det
    inv_K(1,3) =  (K(1,2)*K(2,3)-K(2,2)*K(1,3))/det

    inv_K(2,1) = -(K(2,1)*K(3,3)-K(3,1)*K(2,3))/det
    inv_K(2,2) =  (K(1,1)*K(3,3)-K(3,1)*K(1,3))/det
    inv_K(2,3) = -(K(1,1)*K(2,3)-K(2,1)*K(1,3))/det
    inv_K(3,3) =  (K(1,1)*K(2,2)-K(1,2)*K(2,1))/det

  end subroutine inverse


END MODULE TB_ham




