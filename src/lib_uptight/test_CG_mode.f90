! Test driver for coarse-graining (ICGN) against full diagonalization.
!
! Modes:
!   (1) Standard full diagonalization (no CG)
!   (2) Coupling CG  — partition_type = 0 (METIS/graph), needs num_blocks
!   (3) Nested CG    — partition_type = 2 (cube→row→plane), needs nested_axes/grids
!
! Flow (matches uptight.F90):
!   UPT_configure_coarse_graining → UPT_set_coarse_graining_partition
!     [+ UPT_set_coarse_graining_nested / nested_tiers if nested]
!   → upt_hamiltonian  (sparse_ham → icgn_prepare → swap_in)
!   → UPT_lapack / UPT_jd / UPT_lanczos  (solve + icgn_finalize_eigenvectors)
!
! AAD: modes 2/3 vs mode 1 reference (sorted ascending).
!   n_up:   n smallest positive evals;  n_down: n largest negative evals.
!
! Config file "input" (one value per line):
!   1:  structure file (.upg)
!   2:  relativistic (.true./.false.)
!   3:  Harrison scaling (.true./.false.)
!   4:  c-axis (3 floats)
!   5:  solver (LK / JD / LO)
!   6:  nVB
!   7:  nCB
!   8:  lambda_vb (eV)
!   9:  lambda_cb (eV)
!  10:  partition_type  (0 = coupling, 2 = nested; both modes are run if both usable)
!  11:  num_blocks      (coupling only; 0 = auto; ignored for nested)
!  12:  metis_imbalance
!  13:  core_emin (eV)
!  14:  core_emax (eV)
!  15:  top_buffer (eV)
!  16:  bottom_buffer (eV)
!  17:  outer_coupling_threshold (epsilon)
!  18:  add_core_outers  (0 or 1)
!  19:  neumann_order    (<0 = off)
!  20:  expansion_energy E0 (0.0 = auto = core midpoint)
!  21:  nested_axes  (3 ints in {1,2,3}, permutation; nested only)
!  22:  nested_grids (3 positive ints; nested only)
!  23:  layer_tolerance (Angstrom, <=0 = auto; nested only)
!  24:  n_up   (AAD)
!  25:  n_down (AAD)
program test_CG_mode

  USE precision
  USE globals,             only : MST, LST
  USE mpi_globals,         only : upt_mpi_init, upt_mpi_end
  USE upt_param,           only : OUPT
  USE struct_building,     only : init_structure, make_basis, init_basis, subs_dg_ions
  USE neighbours,          only : refine_neighbours_map, check_input_nn_list, &
                                  check_nn_map, write_neighbours_map
  USE input_data,          only : read_data
  USE type_defs,           only : write_basis, write_materials
  USE states_and_couplings, only : ref_states_and_couplings, init_n_st, &
                                   sort_states, set_max_order
  USE alloys,              only : init_mat_ion
  USE uptight,             only : UPT_configure_coarse_graining,  &
                                  UPT_set_coarse_graining_partition, &
                                  UPT_set_coarse_graining_nested, &
                                  UPT_set_coarse_graining_nested_tiers, &
                                  UPT_get_coarse_graining_info,    &
                                  upt_hamiltonian, upt_lapack, upt_jd, upt_lanczos
  USE sparse_matrix,       only : destroy_matrix
  USE clock,               only : set_clock, get_sclock

  IMPLICIT NONE

  TYPE(OUPT), TARGET  :: upt
  TYPE(OUPT), POINTER :: pupt

  INTEGER        :: i, err, n_ham, num_ev
  CHARACTER(LST) :: config_file
  CHARACTER(MST) :: solver_choice
  INTEGER        :: n_blocks, nVB, nCB, partition_type
  INTEGER        :: add_core_outers, neumann_order
  REAL(dp)       :: metis_imbalance
  REAL(dp)       :: core_emin, core_emax, top_buffer, bottom_buffer, epsilon, E0
  INTEGER        :: nested_axes(3), nested_grids(3)
  REAL(dp)       :: layer_tol
  INTEGER        :: n_up, n_down
  REAL(sp)       :: solve_time
  LOGICAL        :: ready
  INTEGER        :: orig_dim, red_dim, nb_out
  REAL(dp)       :: cut_frac
  REAL(dp), ALLOCATABLE :: ref_evals(:)
  ! Per-tier copies (nested uses same scalar on every tier in this test)
  REAL(dp) :: tier_emin(3), tier_emax(3), tier_top(3), tier_bot(3)
  REAL(dp) :: tier_eps(3), tier_E0(3)
  INTEGER  :: tier_outers(3), tier_orders(3)

  call upt_mpi_init(0)
  pupt => upt

  upt%database_path = './'
  upt%work_path     = './'
  upt%out_path      = './'
  upt%gen_out       = 'out.gen'
  upt%state_file    = 'states.data'
  upt%sparse_format = 'U'
  upt%verbose       = 10

  ! ---- read input -----------------------------------------------------------
  config_file = 'input'
  open(10, file=trim(config_file), status='old', action='read', iostat=err)
  if (err /= 0) then
     write(*,*) 'ERROR: cannot open input file: ', trim(config_file); stop 1
  end if
  read(10,*) upt%gen_filename
  read(10,*) upt%relat
  read(10,*) upt%scaling
  read(10,*) upt%c_axis(:)
  read(10,*) solver_choice
  read(10,*) nVB
  read(10,*) nCB
  read(10,*) upt%lambda_vb
  read(10,*) upt%lambda_cb
  read(10,*) partition_type
  read(10,*) n_blocks
  read(10,*) metis_imbalance
  read(10,*) core_emin
  read(10,*) core_emax
  read(10,*) top_buffer
  read(10,*) bottom_buffer
  read(10,*) epsilon
  read(10,*) add_core_outers
  read(10,*) neumann_order
  read(10,*) E0
  read(10,*) nested_axes(1), nested_axes(2), nested_axes(3)
  read(10,*) nested_grids(1), nested_grids(2), nested_grids(3)
  read(10,*) layer_tol
  read(10,*) n_up
  read(10,*) n_down
  close(10)

  if (abs(E0) < 1.0e-30_dp) E0 = 0.5_dp * (core_emin + core_emax)

  write(*,'(a)') '========================================'
  write(*,'(a,a)')     ' Structure:       ', trim(upt%gen_filename)
  write(*,'(a,l1)')    ' Relativistic:    ', upt%relat
  write(*,'(a,l1)')    ' Scaling:         ', upt%scaling
  write(*,'(a,a)')     ' Solver:          ', trim(solver_choice)
  write(*,'(a,i0)')    ' partition_type:  ', partition_type
  write(*,'(a,i0)')    ' num_blocks:      ', n_blocks
  write(*,'(a,2f8.3)') ' core window:     ', core_emin, core_emax
  write(*,'(a,2f8.3)') ' buffers top/bot: ', top_buffer, bottom_buffer
  write(*,'(a,es10.2)')' epsilon:         ', epsilon
  write(*,'(a,i0)')    ' add_core_outers: ', add_core_outers
  write(*,'(a,i0)')    ' neumann_order:   ', neumann_order
  write(*,'(a,f8.3)')  ' E0:              ', E0
  write(*,'(a,3i3)')   ' nested_axes:     ', nested_axes
  write(*,'(a,3i3)')   ' nested_grids:    ', nested_grids
  write(*,'(a,f8.4)')  ' layer_tol:       ', layer_tol
  write(*,'(a,i0,a,i0)') ' AAD bands +/ -: ', n_up, ' / ', n_down
  write(*,'(a)') '========================================'

  upt%structure%gen_filename = upt%gen_filename
  upt%d_onsite_shift_flag = .true.
  upt%potential_flag      = .false.
  upt%syst_rotated        = .false.
  upt%ioutput_flag        = .false.
  upt%optmat              = .false.
  upt%poldir  = 3
  upt%d_H     = 0.1d0
  upt%E_H     = -200.0d0
  upt%estimate_factor = 1.0
  upt%check_bondmap   = .false.
  upt%n_spin = merge(2, 1, upt%relat)

  upt%num_vb   = nVB;  upt%num_cb   = nCB
  upt%start_vb = 1;    upt%start_cb = 1
  upt%min_iter = 2;    upt%long_iter = 30;  upt%max_iter = 100000
  upt%fast_tol = 1.0d-1; upt%long_tol = 1.0d-10; upt%ort_tol = 1.0d-5
  upt%solver_flag = 0;  upt%dynamic = .true.
  upt%seed_flag   = .false.; upt%bitoff = 0.1_dp
  upt%k_point = (/ 0.0d0, 0.0d0, 0.0d0 /)
  upt%icgn_enabled = .false.

  ! ---- build structure once -------------------------------------------------
  write(*,'(a)') ' Building structure...'
  call set_machine_acc
  call init_structure(upt%verbose, upt%structure, upt%materials, upt%nr_mat, &
                      upt%interfaces, upt%nr_int)
  call check_input_nn_list(upt%structure)
  do i = 1, upt%nr_mat
     call read_data(upt%materials(i), upt%work_path, upt%database_path)
  end do
  call make_basis(upt%verbose, upt%structure, upt%basis)
  call ref_states_and_couplings(upt%ref_states, upt%n_ref_st, &
                                upt%ref_couplings, upt%n_ref_cpl)
  do i = 1, upt%nr_mat
     call sort_states(upt%materials(i), upt%ref_states, upt%ref_couplings)
  end do
  call set_max_order(upt%materials)
  do i = 1, upt%nr_mat; call init_mat_ion(upt%materials(i)); end do
  call write_materials(upt%nr_mat, upt%materials)
  call init_basis(upt%basis, upt%materials)
  call write_basis(upt%basis, 0)
  call refine_neighbours_map(upt%structure, upt%basis, upt%materials, upt%nn_map)
  call write_neighbours_map(upt%nn_map)
  call check_nn_map(upt%nn_map)
  call subs_dg_ions(upt%basis, upt%materials, upt%nn_map)
  call init_n_st(upt%basis, upt%materials)

  ! ---- MODE 1: standard (no CG) ---------------------------------------------
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 1: Standard full diagonalization'
  write(*,'(a)') '========================================'

  upt%icgn_enabled = .false.
  upt%verbose = 0
  call upt_hamiltonian(pupt)
  n_ham = upt%ham%nrow
  write(*,'(a,i0)') ' Full Hamiltonian dimension: ', n_ham

  num_ev = nVB + nCB
  allocate(upt%eigen_values(num_ev), upt%eigen_vectors(n_ham, num_ev), &
           upt%particles(num_ev), stat=err)
  upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0

  call set_clock()
  call run_solver(solver_choice, pupt)
  solve_time = get_sclock()

  write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
  write(*,'(a,f10.3)')  ' Total time:   ', solve_time
  write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)

  allocate(ref_evals(size(upt%eigen_values)))
  ref_evals = upt%eigen_values
  call sort_ascending(ref_evals)

  call write_eigenvalues('eigenvalues_standard.dat', upt%eigen_values, solve_time, &
       'STANDARD', n_ham, n_ham, 0.0_dp, ref_evals, 0, 0)
  deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
  call destroy_matrix(upt%ham)

  ! ---- MODE 2: coupling CG --------------------------------------------------
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 2: Coupling CG (partition_type = 0)'
  write(*,'(a)') '========================================'

  call UPT_configure_coarse_graining(upt, n_blocks, metis_imbalance, &
       core_emin, core_emax, top_buffer, bottom_buffer, epsilon, &
       add_core_outers, neumann_order, E0, .false., 1000, 1.0e-3_dp)
  call UPT_set_coarse_graining_partition(upt, 0, 0.0_dp)

  call set_clock()
  call upt_hamiltonian(pupt)
  call UPT_get_coarse_graining_info(upt, ready, orig_dim, red_dim, nb_out, cut_frac)
  if (.not. ready) then
     write(*,*) ' WARNING: coupling CG not ready, skip mode 2'
     call destroy_matrix(upt%ham)
  else
     write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
          '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(max(1,orig_dim),dp)), '% reduction)'
     write(*,'(a,i0)') ' Blocks used: ', nb_out

     num_ev = nVB + nCB
     allocate(upt%eigen_values(num_ev), upt%eigen_vectors(red_dim, num_ev), &
              upt%particles(num_ev), stat=err)
     upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0

     call run_solver(solver_choice, pupt)
     solve_time = get_sclock()

     write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
     write(*,'(a,f10.3)')  ' Total time:   ', solve_time
     write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
     call write_eigenvalues('eigenvalues_coupling.dat', upt%eigen_values, solve_time, &
          'COUPLING', orig_dim, red_dim, cut_frac, ref_evals, n_up, n_down)
     deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
     call destroy_matrix(upt%ham)
  end if

  ! ---- MODE 3: nested CG ----------------------------------------------------
  write(*,'(a)') '========================================'
  write(*,'(a)') ' MODE 3: Nested CG (partition_type = 2)'
  write(*,'(a)') '========================================'

  ! nblocks is ignored for nested; pass 0. Nested grids/axes drive partition.
  call UPT_configure_coarse_graining(upt, 0, metis_imbalance, &
       core_emin, core_emax, top_buffer, bottom_buffer, epsilon, &
       add_core_outers, neumann_order, E0, .false., 1000, 1.0e-3_dp)
  call UPT_set_coarse_graining_partition(upt, 2, layer_tol)
  call UPT_set_coarse_graining_nested(upt, nested_axes, nested_grids)

  tier_emin = core_emin; tier_emax = core_emax
  tier_top  = top_buffer; tier_bot = bottom_buffer
  tier_eps  = epsilon;    tier_E0  = E0
  tier_outers = add_core_outers; tier_orders = neumann_order
  call UPT_set_coarse_graining_nested_tiers(upt, tier_emin, tier_emax, &
       tier_top, tier_bot, tier_eps, tier_outers, tier_orders, tier_E0)

  call set_clock()
  call upt_hamiltonian(pupt)
  call UPT_get_coarse_graining_info(upt, ready, orig_dim, red_dim, nb_out, cut_frac)
  if (.not. ready) then
     write(*,*) ' WARNING: nested CG not ready, skip mode 3'
     call destroy_matrix(upt%ham)
  else
     write(*,'(a,i0,a,i0,a,f6.2,a)') ' Reduced: ', orig_dim, ' -> ', red_dim, &
          '  (', 100.0_dp*(1.0_dp - real(red_dim,dp)/real(max(1,orig_dim),dp)), '% reduction)'

     num_ev = nVB + nCB
     allocate(upt%eigen_values(num_ev), upt%eigen_vectors(red_dim, num_ev), &
              upt%particles(num_ev), stat=err)
     upt%eigen_values = 0.0d0; upt%eigen_vectors = (0.0d0,0.0d0); upt%particles = 0

     call run_solver(solver_choice, pupt)
     solve_time = get_sclock()

     write(*,'(a,i0)')     ' Bands found:  ', size(upt%eigen_values)
     write(*,'(a,f10.3)')  ' Total time:   ', solve_time
     write(*,'(a,2f10.4)') ' Energy range: ', minval(upt%eigen_values), maxval(upt%eigen_values)
     call write_eigenvalues('eigenvalues_nested.dat', upt%eigen_values, solve_time, &
          'NESTED', orig_dim, red_dim, cut_frac, ref_evals, n_up, n_down)
     deallocate(upt%eigen_values, upt%eigen_vectors, upt%particles)
     call destroy_matrix(upt%ham)
  end if

  if (allocated(ref_evals)) deallocate(ref_evals)
  call upt_mpi_end

contains

  subroutine run_solver(choice, p)
    character(*), intent(in) :: choice
    type(OUPT), pointer :: p
    select case (trim(choice))
    case ('LK'); call upt_lapack(p)
    case ('JD'); call upt_jd(p)
    case ('LO'); call upt_lanczos(p)
    case default
       write(*,*) 'Unknown solver: ', trim(choice); stop 1
    end select
  end subroutine run_solver

  subroutine sort_ascending(a)
    real(dp), intent(inout) :: a(:)
    integer  :: i, j
    real(dp) :: tmp
    do i = 2, size(a)
       tmp = a(i); j = i - 1
       do while (j >= 1 .and. a(j) > tmp)
          a(j+1) = a(j); j = j - 1
       end do
       a(j+1) = tmp
    end do
  end subroutine sort_ascending

  subroutine compute_aad(evals, ref, n_up, n_down, aad_up, aad_down)
    real(dp), intent(in)  :: evals(:), ref(:)
    integer,  intent(in)  :: n_up, n_down
    real(dp), intent(out) :: aad_up, aad_down
    integer  :: n_ev, n_ref, i, first_pos_ev, first_pos_ref, last_neg_ev, last_neg_ref

    n_ev  = size(evals); n_ref = size(ref)
    first_pos_ev  = n_ev  + 1
    first_pos_ref = n_ref + 1
    do i = 1, n_ev;  if (evals(i) >= 0.0_dp) then; first_pos_ev  = i; exit; end if; end do
    do i = 1, n_ref; if (ref(i)   >= 0.0_dp) then; first_pos_ref = i; exit; end if; end do
    last_neg_ev  = first_pos_ev  - 1
    last_neg_ref = first_pos_ref - 1

    aad_up = -1.0_dp
    if (n_up > 0) then
       if ((n_ev - first_pos_ev + 1 >= n_up) .and. (n_ref - first_pos_ref + 1 >= n_up)) then
          aad_up = 0.0_dp
          do i = 0, n_up - 1
             aad_up = aad_up + abs(evals(first_pos_ev + i) - ref(first_pos_ref + i))
          end do
          aad_up = aad_up / real(n_up, dp)
       end if
    end if

    aad_down = -1.0_dp
    if (n_down > 0) then
       if ((last_neg_ev >= n_down) .and. (last_neg_ref >= n_down)) then
          aad_down = 0.0_dp
          do i = 0, n_down - 1
             aad_down = aad_down + abs(evals(last_neg_ev - i) - ref(last_neg_ref - i))
          end do
          aad_down = aad_down / real(n_down, dp)
       end if
    end if
  end subroutine compute_aad

  subroutine write_eigenvalues(fname, evals, t, mode, ndim, nred, cut, ref, nu, nd)
    character(*), intent(in) :: fname, mode
    real(dp),     intent(in) :: evals(:), cut, ref(:)
    real(sp),     intent(in) :: t
    integer,      intent(in) :: ndim, nred, nu, nd

    integer  :: n, fu, ii
    real(dp), allocatable :: se(:)
    real(dp) :: aad_up, aad_down

    n = size(evals)
    allocate(se(n)); se = evals; call sort_ascending(se)

    aad_up = -1.0_dp; aad_down = -1.0_dp
    if ((nu > 0 .or. nd > 0) .and. size(ref) > 0) &
         call compute_aad(se, ref, nu, nd, aad_up, aad_down)

    open(newunit=fu, file=trim(fname), status='replace', action='write')
    write(fu,'(a,a)')     '# Mode: ', trim(mode)
    write(fu,'(a,i0)')    '# Full dimension:    ', ndim
    write(fu,'(a,i0)')    '# Reduced dimension: ', nred
    write(fu,'(a,f8.2,a)') '# Rank reduction: ', &
         100.0_dp*(1.0_dp - real(nred,dp)/real(max(1,ndim),dp)), ' %'
    write(fu,'(a,f10.4)') '# Cut fraction:      ', cut
    write(fu,'(a,f12.6)') '# Total time (s):    ', t
    write(fu,'(a,i0)')    '# Total bands:       ', n
    if (nu > 0) then
       if (aad_up >= 0.0_dp) then
          write(fu,'(a,i0,a,es14.6)') '# AAD ', nu, ' smallest positive evals (eV): ', aad_up
       else
          write(fu,'(a,i0,a)') '# AAD ', nu, ' smallest positive evals (eV): N/A'
       end if
    end if
    if (nd > 0) then
       if (aad_down >= 0.0_dp) then
          write(fu,'(a,i0,a,es14.6)') '# AAD ', nd, ' largest  negative evals (eV): ', aad_down
       else
          write(fu,'(a,i0,a)') '# AAD ', nd, ' largest  negative evals (eV): N/A'
       end if
    end if
    write(fu,'(a)') '#'
    write(fu,'(a)') '# Index    Energy(eV)'
    do ii = 1, n
       write(fu,'(i6,2x,f16.8)') ii, se(ii)
    end do
    close(fu)

    write(*,'(a,a)') ' Output: ', trim(fname)
    if (nu > 0 .and. aad_up   >= 0.0_dp) &
         write(*,'(a,i0,a,es12.4)') '  AAD ', nu, ' smallest pos: ', aad_up
    if (nd > 0 .and. aad_down >= 0.0_dp) &
         write(*,'(a,i0,a,es12.4)') '  AAD ', nd, ' largest  neg: ', aad_down
    deallocate(se)
  end subroutine write_eigenvalues

end program test_CG_mode
