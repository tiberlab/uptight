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
module upt_param 
  
  use precision, only : sp, dp
  use globals
  use type_defs
  use neighbours, only : nearest_neighbours
  use sparse_matrix

  implicit none
  private


  public :: OUPT, set_defaults

  !!* Parameters needed during UPT calculations  
  type OUPT

     INTEGER :: verbose
     INTEGER :: mpi_comm

     CHARACTER ( MST ) :: gen_filename
     CHARACTER ( MST ) :: gen_out
     CHARACTER ( MST ) :: state_file
     CHARACTER ( MST ) :: potential_file
     CHARACTER ( LST ) :: database_path
     CHARACTER ( LST ) :: work_path
     CHARACTER ( 1 )   :: sparse_format
     CHARACTER ( LST ) :: out_path
     CHARACTER ( LST ) :: load_path
     CHARACTER ( MST ) :: out_format

     type(TStructure) :: structure
     type(ion_basis)  :: basis
     integer :: nr_mat
     type(material_data), DIMENSION(:), POINTER :: materials
     integer :: nr_int
     type(material_data), DIMENSION(:,:), POINTER :: interfaces
     type(nearest_neighbours), DIMENSION(:), POINTER :: nn_map
     REAL ( dp ) :: d_H   ! Scaling of interactions for empirical Hydrogens
     REAL ( dp ) :: E_H   ! Scaling of on-site for empirical Hydrogens
     INTEGER :: n_spin    ! Spin variables 
     INTEGER :: poldir    ! Polarization direction for optical transitions
     INTEGER :: n_ref_st  ! number of reference states
     INTEGER :: n_ref_cpl ! number of reference couplings
     REAL ( dp ), DIMENSION( : ), POINTER :: pot_data

     REAL ( dp ) :: k_point(3)      ! k-point at which H is computed 
     REAL ( dp ) :: c_axis(3)       ! direction of c-axis for wurtzites

     LOGICAL :: d_onsite_shift_flag ! Onsite deformation for zincblend
     LOGICAL :: potential_flag      ! Adds atomic potential shifts
     LOGICAL :: relat               ! Spin-polarized Hamiltonian
     LOGICAL :: scaling             ! Harrison scaling of TB parameters
     LOGICAL :: syst_rotated        ! Tells if c-axis in wz is along (111)
     LOGICAL :: ioutput_flag        ! Performs i/o ?
     LOGICAL :: optmat              ! Computes optical matrix elements
     LOGICAL :: check_bondmap       ! Checks bondmap consistency
     LOGICAL :: hybrid_passivation  ! Set on the hybrid orbital passivation 

     REAL(dp) :: vol_fraction       ! cutoff isosurface    
     REAL(dp) :: grid_step          ! step for output grid

     CHARACTER(STATELEN),   DIMENSION( : ),  POINTER   :: ref_states  
     CHARACTER(CPLLEN),   DIMENSION( : ),  POINTER   :: ref_couplings  
     type(CSR_ex) :: ham2           ! Sparse Hamiltonian Matrix
     type(CSR)    :: ham            ! Sparse Hamiltonian Matrix
     type(CSR)    :: U              ! Time-reversal symmetry operator
 
     REAL ( dp ),   DIMENSION( : ),    POINTER     :: eigen_values
     COMPLEX ( dp ), DIMENSION( :,: ),    POINTER  :: eigen_vectors
     INTEGER, DIMENSION( : ), POINTER              :: particles

     ! Parameters for internal Lanczos Solver 
     INTEGER :: num_vb              ! number of valence eigenvalues 
     INTEGER :: num_cb              ! number of conduct eigenvalues 
     INTEGER :: start_vb            ! start valence eigenvalue (1) 
     INTEGER :: start_cb            ! start conduct eigenvalue (1)


     INTEGER :: min_iter            ! fast Lanczos steps 
     INTEGER :: long_iter           ! long (standard) L. steps
     INTEGER :: max_iter            ! max number of L. steps

     REAL(dp) :: lambda_cb          ! guess for conduction eigenv.
     REAL(dp) :: lambda_vb          ! guess for valence eigenv.
     REAL(dp) :: fast_tol           ! fast iteration tolerance
     REAL(dp) :: long_tol           ! long iteration tolerance
     REAL(dp) :: ort_tol            ! orthogonality tolerance
     
     INTEGER :: solver_flag         ! solver flavour (0: CPU, 1: GPU)
     LOGICAL :: dynamic             ! shift dyn for next eigenv 
     REAL(dp) :: bitoff             ! tiny offset in dynamic search
     LOGICAL :: seed_flag           ! random start or use old vector

     REAL(sp) :: estimate_factor    ! for estimating Hamil size

  end type OUPT

contains

 subroutine set_defaults(upt)
   type(OUPT) :: upt

   upt%verbose=0

   upt%gen_filename='str.upg'
   upt%gen_out='out.gen'
   upt%state_file='states.upt'
   upt%potential_file='pot_on_atoms.dat'
   upt%database_path='.' 
   upt%work_path='.'
   upt%sparse_format='F'
   upt%out_path='.'
   upt%load_path='.'
   upt%out_format='F'

   upt%d_H= 0.d1     ! Scaling of interactions for empirical Hydrogens
   upt%E_H= -200.0d0 ! Scaling of on-site for empirical Hydrogens
   upt%n_spin= 2     ! Spin variables 
   upt%poldir= 3     ! Polarization direction for optical transitions
   upt%n_ref_st= 10  ! number of reference states
   upt%n_ref_cpl=21  ! number of reference couplings
   upt%k_point(:) = (/ 0.d0 ,0.d0, 0.d0 /)     ! k-point at which H is computed 
   upt%c_axis(:) = (/ 0.d0 ,0.d0, 1.d0 /)      ! direction of c-axis for wurtzites

   upt%d_onsite_shift_flag = .false. ! Onsite deformation for zincblend
   upt%potential_flag = .false.      ! Adds atomic potential shifts
   upt%relat = .true.               ! Spin-polarized Hamiltonian
   upt%scaling = .true.            ! Harrison scaling of TB parameters
   upt%syst_rotated = .false.        ! Tells if c-axis in wz is along (111)
   upt%ioutput_flag = .true.       ! Performs i/o ?
   upt%optmat = .false.              ! Computes optical matrix elements
   upt%check_bondmap = .false.      ! Checks bondmap consistency
   upt%hybrid_passivation = .true.  ! Set on the hybrid orbital passivation 
                                    ! (PRB 69, 045316 2004)

   upt%vol_fraction = 0.60         ! cutoff isosurface    
   upt%grid_step = 0.5d0           ! step for output grid

   upt%num_vb = 2          ! number of valence eigenvalues 
   upt%num_cb = 2          ! number of conduct eigenvalues 
   upt%start_vb = 1        ! start valence eigenvalue (1) 
   upt%start_cb = 1        ! start conduct eigenvalue (1)

   upt%min_iter = 30       ! fast Lanczos steps 
   upt%long_iter = 32      ! long (standard) L. steps
   upt%max_iter = 10000    ! max number of L. steps

   upt%lambda_cb = 0.d0    ! guess for conduction eigenv.
   upt%lambda_vb = 0.d0    ! guess for valence eigenv.
   upt%fast_tol = 1.0d-1   ! fast iteration tolerance
   upt%long_tol = 1.0d-9   ! long iteration tolerance
   upt%ort_tol = 1.0d-6    ! orthogonality tolerance
                           
   upt%solver_flag =  0    ! perform lanczos on CPU
   upt%dynamic = .true.    ! shift dyn for next eigenv 
   upt%bitoff  = 0.1_dp    ! offset in dynamic search 
   upt%seed_flag = .true.  ! random start or use old vector

   upt%estimate_factor = 1.0

 end subroutine set_defaults

end module upt_param
