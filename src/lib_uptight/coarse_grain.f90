! Coarse-grained tight-binding projection after Liu et al. (2022).
module coarse_grain
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use precision, only : dp
  use upt_param, only : OUPT, CGBlock
  use sparse_matrix, only : CSR, create_matrix, destroy_matrix
  use mpi_globals, only : num_procs
  implicit none
  private

  type CGPair
     integer :: a = 0, b = 0
     complex(dp), dimension(:,:), pointer :: v => null()
  end type CGPair

  public :: cg_configure, cg_prepare, cg_clear, cg_active, cg_lift
  public :: cg_get_info
  public :: icg_configure, icg_prepare, icg_clear, icg_active, icg_lift
  public :: icg_get_info
  public :: icgn_configure, icgn_prepare, icgn_clear, icgn_active, icgn_lift
  public :: icgn_get_info

  interface
     integer(c_int) function cg_metis_partition(nvtxs, xadj, adjncy, vwgt, &
          adjwgt, nparts, ufactor, seed, part) bind(C, name='upt_cg_metis_partition')
       import :: c_int
       integer(c_int), value :: nvtxs, nparts, ufactor, seed
       integer(c_int), intent(in) :: xadj(*), adjncy(*), vwgt(*), adjwgt(*)
       integer(c_int), intent(out) :: part(*)
     end function cg_metis_partition
  end interface

contains

  subroutine cg_configure(upt, enabled, nblocks, emin, emax, imbalance)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled
    integer, intent(in) :: nblocks
    real(dp), intent(in) :: emin, emax, imbalance
    call cg_clear(upt)
    upt%cg_enabled = enabled
    upt%cg_num_blocks = nblocks
    upt%cg_emin = emin
    upt%cg_emax = emax
    upt%cg_imbalance = imbalance
  end subroutine cg_configure

  logical function cg_active(upt)
    type(OUPT), intent(in) :: upt
    cg_active = upt%cg_enabled .and. upt%cg_ready
  end function cg_active

  subroutine cg_get_info(upt, ready, original_dim, reduced_dim, nblocks, cut_fraction)
    type(OUPT), intent(in) :: upt
    logical, intent(out) :: ready
    integer, intent(out) :: original_dim, reduced_dim, nblocks
    real(dp), intent(out) :: cut_fraction
    ready = upt%cg_ready
    original_dim = upt%cg_original_dim
    reduced_dim = upt%cg_reduced_dim
    nblocks = upt%cg_num_blocks
    cut_fraction = upt%cg_cut_fraction
  end subroutine cg_get_info

  subroutine cg_clear(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    if (associated(upt%cg_ham%M)) call destroy_matrix(upt%cg_ham)
    if (associated(upt%cg_U%M)) call destroy_matrix(upt%cg_U)
    if (associated(upt%cg_blocks)) then
       do i = 1, size(upt%cg_blocks)
          if (associated(upt%cg_blocks(i)%rows)) deallocate(upt%cg_blocks(i)%rows)
          if (associated(upt%cg_blocks(i)%eval)) deallocate(upt%cg_blocks(i)%eval)
          if (associated(upt%cg_blocks(i)%q)) deallocate(upt%cg_blocks(i)%q)
          if (associated(upt%cg_blocks(i)%evals_full)) deallocate(upt%cg_blocks(i)%evals_full)
          if (associated(upt%cg_blocks(i)%S_full)) deallocate(upt%cg_blocks(i)%S_full)
          if (associated(upt%cg_blocks(i)%retained_idx)) deallocate(upt%cg_blocks(i)%retained_idx)
       end do
       deallocate(upt%cg_blocks)
    end if
    upt%cg_ready = .false.
    upt%cg_original_dim = 0
    upt%cg_reduced_dim = 0
    upt%cg_cut_fraction = 0.0_dp
  end subroutine cg_clear

  subroutine cg_prepare(upt, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(out) :: ierr
    integer :: n, na, nb, i, j, k, p, q, status, nedge, maxedge
    integer :: r, c, br, bc, lr, lc, npair, total_ret, pos
    integer, allocatable :: atom_of(:), local_of(:), label(:), row_of(:)
    integer, allocatable :: counts(:), cursor(:), offsets(:), bsize(:)
    integer(c_int), allocatable :: xadj(:), adjncy(:), vwgt(:), adjwgt(:), part(:)
    real(dp), allocatable :: edge_weight(:)
    real(dp) :: max_weight, all_weight, cut_weight, workspace_mib
    type(CGPair), allocatable :: pairs(:)

    ierr = 0
    call cg_clear(upt)
    if (.not.upt%cg_enabled) return
    if (num_procs /= 1) then
       ierr = 1; write(*,*) '(coarse grain) MPI runs are not supported'; return
    end if
    na = upt%basis%n_basis
    n = upt%ham%nrow
    if (na < 1 .or. upt%cg_num_blocks < 1 .or. upt%cg_num_blocks > na) then
       ierr = 2; write(*,*) '(coarse grain) invalid number of blocks'; return
    end if
    if (upt%cg_emin >= upt%cg_emax .or. upt%cg_imbalance < 0.0_dp) then
       ierr = 3; write(*,*) '(coarse grain) invalid energy window or imbalance'; return
    end if
    if (.not.associated(upt%ham%M)) then
       ierr = 4; write(*,*) '(coarse grain) Hamiltonian is not initialized'; return
    end if
    if (upt%cg_num_blocks == 1 .and. upt%verbose > 0) then
       write(*,*) '(coarse grain) one block selected; energy window controls rank'
    end if

    allocate(atom_of(n), local_of(n), offsets(na+1), bsize(na))
    pos = 1
    do i = 1, na
       offsets(i) = pos
       bsize(i) = upt%n_spin * upt%basis%n_st(i)
       do j = 1, bsize(i)
          atom_of(pos) = i; local_of(pos) = j; pos = pos + 1
       end do
    end do
    offsets(na+1) = pos
    if (pos-1 /= n) then
       ierr = 5; write(*,*) '(coarse grain) atom/orbital mapping is inconsistent'; return
    end if

    ! Parallel graph edges are intentional: their summed weights are the
    ! Frobenius norm squared of an atom-to-atom Hamiltonian block.
    maxedge = max(1, upt%ham%nnz)
    allocate(edge_weight(maxedge), counts(na), label(na))
    nedge = 0; max_weight = 0.0_dp; counts = 0
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          nedge = nedge + 1
          if (nedge > maxedge) then
             ierr = 6; return
          end if
          edge_weight(nedge) = abs(upt%ham%M(k))**2
          counts(atom_of(r)) = counts(atom_of(r)) + 1
          counts(atom_of(c)) = counts(atom_of(c)) + 1
          max_weight = max(max_weight, edge_weight(nedge))
       end do
    end do
    if (nedge == 0 .or. max_weight == 0.0_dp) then
       ierr = 7; write(*,*) '(coarse grain) atom graph has no couplings'; return
    end if
    allocate(xadj(na+1), cursor(na), adjncy(2*nedge), adjwgt(2*nedge), vwgt(na), part(na))
    xadj(1) = 0_c_int
    do i = 1, na
       xadj(i+1) = xadj(i) + int(counts(i), c_int)
       cursor(i) = int(xadj(i)) + 1
       vwgt(i) = int(bsize(i), c_int)
    end do
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          br = atom_of(r); bc = atom_of(c)
          p = max(1, nint(abs(upt%ham%M(k))**2 / max_weight * 1000000.0_dp))
          adjncy(cursor(br)) = int(bc-1, c_int); adjwgt(cursor(br)) = int(p, c_int); cursor(br)=cursor(br)+1
          adjncy(cursor(bc)) = int(br-1, c_int); adjwgt(cursor(bc)) = int(p, c_int); cursor(bc)=cursor(bc)+1
       end do
    end do
    status = cg_metis_partition(int(na,c_int), xadj, adjncy, vwgt, adjwgt, &
         int(upt%cg_num_blocks,c_int), int(nint(1000.0_dp*upt%cg_imbalance),c_int), 42_c_int, part)
    if (status /= 0) then       ! METIS not available — fallback to connectivity-aware partition
       if (upt%verbose > 0) then
          write(*,*) '(coarse grain) METIS unavailable, using graph-BFS fallback partitioning'
       end if
       call cg_graph_partition(na, upt%cg_num_blocks, vwgt, xadj, adjncy, adjwgt, part)
    end if
    do i = 1, na
       label(i) = int(part(i)) + 1
    end do

    all_weight = 0.0_dp; cut_weight = 0.0_dp
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_weight = all_weight + abs(upt%ham%M(k))**2
          if (label(atom_of(r)) /= label(atom_of(c))) cut_weight = cut_weight + abs(upt%ham%M(k))**2
       end do
    end do
    upt%cg_cut_fraction = cut_weight / all_weight

    deallocate(counts)
    allocate(upt%cg_blocks(upt%cg_num_blocks), counts(upt%cg_num_blocks))
    counts = 0
    do i = 1, na
       counts(label(i)) = counts(label(i)) + bsize(i)
    end do
    do i = 1, upt%cg_num_blocks
       upt%cg_blocks(i)%nrow = counts(i)
       allocate(upt%cg_blocks(i)%rows(counts(i)))
    end do
    ! ZHEEVD needs the dense matrix plus work arrays.  This is deliberately
    ! only a prediction: the allocation itself remains inside diagonalize_block.
    workspace_mib = 16.0_dp * real(maxval(counts),dp)**2 / (1024.0_dp**2)
    if (upt%verbose > 0) write(*,'(a,i0,a,f10.2,a)') '(coarse grain) largest dense block ', &
         maxval(counts), ', matrix workspace at least ', workspace_mib, ' MiB'
    cursor = 0
    do i = 1, na
       br = label(i)
       do j = offsets(i), offsets(i+1)-1
          cursor(br) = cursor(br) + 1
          upt%cg_blocks(br)%rows(cursor(br)) = j
       end do
    end do

    allocate(row_of(n)); row_of = 0
    do i = 1, upt%cg_num_blocks
       do j = 1, upt%cg_blocks(i)%nrow
          row_of(upt%cg_blocks(i)%rows(j)) = j
       end do
       call diagonalize_block(upt, i, atom_of, row_of, ierr)
       if (ierr /= 0) return
       row_of(upt%cg_blocks(i)%rows) = 0
    end do
    do i = 1, upt%cg_num_blocks
       do j = 1, upt%cg_blocks(i)%nrow
          row_of(upt%cg_blocks(i)%rows(j)) = j
       end do
    end do
    total_ret = 0
    do i = 1, upt%cg_num_blocks
       total_ret = total_ret + upt%cg_blocks(i)%nret
    end do
    if (total_ret == 0) then
       ierr = 9; write(*,*) '(coarse grain) energy window retained no states'; return
    end if
    ! No check needed - we will compute ALL eigenvalues of reduced matrix
    upt%cg_original_dim = n; upt%cg_reduced_dim = total_ret

    call build_reduced_hamiltonian(upt, atom_of, label, row_of, pairs, npair, ierr)
    if (ierr /= 0) return
    call destroy_pairs(pairs)
    upt%cg_ready = .true.
    if (upt%verbose > 0) write(*,'(a,i0,a,i0,a,f8.4)') '(coarse grain) dimension ',n,' -> ',total_ret, &
         ', cut fraction ',upt%cg_cut_fraction
  end subroutine cg_prepare

  subroutine diagonalize_block(upt, ib, atom_of, local, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, atom_of(:), local(:)
    integer, intent(out) :: ierr
    integer :: n, i, k, r, c
    complex(dp), allocatable :: h(:,:)
    real(dp), allocatable :: w(:)
    ierr = 0; n = upt%cg_blocks(ib)%nrow
    if (n == 0) then
       upt%cg_blocks(ib)%nret = 0
       return
    end if
    allocate(h(n,n), w(n)); h = (0.0_dp,0.0_dp)
    do r = 1, upt%ham%nrow
       if (local(r) == 0) cycle
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (local(c) == 0) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          h(local(r),local(c)) = upt%ham%M(k)
          if (r /= c) h(local(c),local(r)) = conjg(upt%ham%M(k))
       end do
    end do
    call dense_eigh(h, w, ierr)
    if (ierr /= 0) return
    ! Store the FULL eigensystem — projection happens in build_reduced_hamiltonian
    ! using all nrow columns of S before the energy window truncation.
    allocate(upt%cg_blocks(ib)%evals_full(n), upt%cg_blocks(ib)%S_full(n,n))
    upt%cg_blocks(ib)%evals_full = w
    upt%cg_blocks(ib)%S_full     = h   ! columns are eigenvectors after zheevd
    ! Count how many states fall in the energy window (needed for reduced dim)
    upt%cg_blocks(ib)%nret = count(w >= upt%cg_emin .and. w <= upt%cg_emax)
    deallocate(h, w)
  end subroutine diagonalize_block

  subroutine dense_eigh(a, w, ierr)
    complex(dp), intent(inout) :: a(:,:)
    real(dp), intent(out) :: w(:)
    integer, intent(out) :: ierr
    integer :: n, info, lwork, lrwork, liwork
    complex(dp) :: workq(1)
    real(dp) :: rworkq(1)
    integer :: iworkq(1)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    n=size(w)
    if (n == 0) then; ierr=0; return; end if
    call zheevd('V','U',n,a,n,w,workq,-1,rworkq,-1,iworkq,-1,info)
    if (info /= 0) then; ierr=10; return; end if
    lwork=max(1,int(real(workq(1)))); lrwork=max(1,int(rworkq(1))); liwork=max(1,iworkq(1))
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    call zheevd('V','U',n,a,n,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    deallocate(work,rwork,iwork); ierr=info
  end subroutine dense_eigh

  logical function stored_entry(fmt, r, c)
    character(1), intent(in) :: fmt
    integer, intent(in) :: r, c
    select case (fmt)
    case ('U', 'u')
       stored_entry = (r <= c)
    case ('L', 'l')
       stored_entry = (r >= c)
    case default
       ! The Hamiltonian builder may emit a complete matrix even when the
       ! sparse format is labelled F.  Process one triangle here; the
       ! projected matrix is completed by its Hermitian conjugate below.
       stored_entry = (r <= c)
    end select
  end function stored_entry


  subroutine build_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(CGPair), allocatable, intent(out) :: pairs(:)
    integer, intent(out) :: npair, ierr
    integer :: i,j,k,r,c,a,b,ia,ib,nnz,pos,slot,nred,na,nb
    integer, allocatable :: roff(:), rowcount(:), next(:)
    integer, allocatable :: win_a(:), win_b(:)   ! indices of retained states within S_full
    complex(dp), allocatable :: g_full(:,:), g_ret(:,:)
    ierr=0; nred=upt%cg_reduced_dim

    ! --- Step 1: Build q/eval for each block from S_full by applying window ---
    ! (This is done here so S_full is available; S_full is freed at the end.)
    do i = 1, upt%cg_num_blocks
       na = upt%cg_blocks(i)%nrow
       if (na == 0 .or. upt%cg_blocks(i)%nret == 0) cycle
       ! collect indices of retained states
       allocate(win_a(upt%cg_blocks(i)%nret))
       k = 0
       do j = 1, na
          if (upt%cg_blocks(i)%evals_full(j) >= upt%cg_emin .and. &
              upt%cg_blocks(i)%evals_full(j) <= upt%cg_emax) then
             k = k + 1; win_a(k) = j
          end if
       end do
       allocate(upt%cg_blocks(i)%eval(upt%cg_blocks(i)%nret))
       allocate(upt%cg_blocks(i)%q(na, upt%cg_blocks(i)%nret))
       do j = 1, upt%cg_blocks(i)%nret
          upt%cg_blocks(i)%eval(j) = upt%cg_blocks(i)%evals_full(win_a(j))
          upt%cg_blocks(i)%q(:,j)  = upt%cg_blocks(i)%S_full(:, win_a(j))
       end do
       deallocate(win_a)
    end do

    allocate(roff(upt%cg_num_blocks+1)); roff(1)=1
    do i=1,upt%cg_num_blocks; roff(i+1)=roff(i)+upt%cg_blocks(i)%nret; end do
    allocate(pairs(max(1,upt%ham%nnz))); npair=0

    ! --- Step 2: Project inter-block couplings using FULL S matrices ---
    ! pairs(slot)%v has shape (nrow_a, nrow_b) — full projection first
    do r=1,upt%ham%nrow
       do k=upt%ham%Mi(r),upt%ham%Mi(r+1)-1
          c=upt%ham%Mj(k); a=label(atom_of(r)); b=label(atom_of(c))
          if(a==b) cycle
          if(upt%cg_blocks(a)%nrow==0 .or. upt%cg_blocks(b)%nrow==0) cycle
          if(.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          ia=min(a,b); ib=max(a,b)
          slot=pair_slot_full(pairs,npair,ia,ib,upt)
          if(slot==0) then; ierr=11; return; end if
          if(a < b) then
             ! S_a(local(r), :) is row local(r) of S_full_a
             call add_outer(pairs(slot)%v, upt%cg_blocks(a)%S_full(local(r),:), &
                  upt%cg_blocks(b)%S_full(local(c),:), upt%ham%M(k))
          else
             call add_outer(pairs(slot)%v, upt%cg_blocks(b)%S_full(local(c),:), &
                  upt%cg_blocks(a)%S_full(local(r),:), conjg(upt%ham%M(k)))
          end if
       end do
    end do

    ! --- Step 3: Cut window and store reduced pairs, then build CSR ---
    ! Convert pairs from full (nrow_a x nrow_b) to retained (nret_a x nret_b)
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       na = upt%cg_blocks(a)%nrow; nb = upt%cg_blocks(b)%nrow
       ! Build index arrays of retained states in a and b
       allocate(win_a(upt%cg_blocks(a)%nret), win_b(upt%cg_blocks(b)%nret))
       k = 0
       do j = 1, na
          if (upt%cg_blocks(a)%evals_full(j) >= upt%cg_emin .and. &
              upt%cg_blocks(a)%evals_full(j) <= upt%cg_emax) then
             k=k+1; win_a(k)=j
          end if
       end do
       k = 0
       do j = 1, nb
          if (upt%cg_blocks(b)%evals_full(j) >= upt%cg_emin .and. &
              upt%cg_blocks(b)%evals_full(j) <= upt%cg_emax) then
             k=k+1; win_b(k)=j
          end if
       end do
       ! Slice: g_ret(i,j) = g_full(win_a(i), win_b(j))
       g_full = pairs(i)%v   ! shape (nrow_a, nrow_b) — copy
       deallocate(pairs(i)%v)
       allocate(pairs(i)%v(upt%cg_blocks(a)%nret, upt%cg_blocks(b)%nret))
       do j = 1, upt%cg_blocks(b)%nret
          do k = 1, upt%cg_blocks(a)%nret
             pairs(i)%v(k,j) = g_full(win_a(k), win_b(j))
          end do
       end do
       deallocate(g_full, win_a, win_b)
    end do

    ! --- Step 4: Free S_full (no longer needed) ---
    do i = 1, upt%cg_num_blocks
       if (associated(upt%cg_blocks(i)%S_full)) then
          deallocate(upt%cg_blocks(i)%S_full)
          nullify(upt%cg_blocks(i)%S_full)
       end if
       if (associated(upt%cg_blocks(i)%evals_full)) then
          deallocate(upt%cg_blocks(i)%evals_full)
          nullify(upt%cg_blocks(i)%evals_full)
       end if
    end do

    ! --- Step 5: Build CSR reduced Hamiltonian (same structure as before) ---
    allocate(rowcount(nred),next(nred)); rowcount=1
    do i=1,npair
       a=pairs(i)%a; b=pairs(i)%b
       select case(upt%ham%sparse_fmt)
       case('F')
          rowcount(roff(a):roff(a+1)-1)=rowcount(roff(a):roff(a+1)-1)+upt%cg_blocks(b)%nret
          rowcount(roff(b):roff(b+1)-1)=rowcount(roff(b):roff(b+1)-1)+upt%cg_blocks(a)%nret
       case('L')
          rowcount(roff(b):roff(b+1)-1)=rowcount(roff(b):roff(b+1)-1)+upt%cg_blocks(a)%nret
       case default
          rowcount(roff(a):roff(a+1)-1)=rowcount(roff(a):roff(a+1)-1)+upt%cg_blocks(b)%nret
       end select
    end do
    nnz=sum(rowcount); call create_matrix(upt%cg_ham,nred,nred,nnz)
    upt%cg_ham%sparse_fmt=upt%ham%sparse_fmt; upt%cg_ham%Mi(1)=1
    do i=1,nred; upt%cg_ham%Mi(i+1)=upt%cg_ham%Mi(i)+rowcount(i); end do
    next=upt%cg_ham%Mi(1:nred)
    do a=1,upt%cg_num_blocks
       do i=1,upt%cg_blocks(a)%nret
          pos=next(roff(a)+i-1); upt%cg_ham%Mj(pos)=roff(a)+i-1; upt%cg_ham%M(pos)=upt%cg_blocks(a)%eval(i)
          next(roff(a)+i-1)=pos+1
       end do
    end do
    do i=1,npair
       a=pairs(i)%a; b=pairs(i)%b
       call emit_pair(upt%cg_ham,pairs(i),roff(a),roff(b),upt%ham%sparse_fmt,next)
    end do
    upt%cg_ham%nnz=nnz
    call create_matrix(upt%cg_U,nred,nred,nred)
    upt%cg_U%sparse_fmt='F'; upt%cg_U%Mi(1)=1
    do i=1,nred
       upt%cg_U%Mi(i)=i; upt%cg_U%Mj(i)=i; upt%cg_U%M(i)=(1.0_dp,0.0_dp)
    end do
    upt%cg_U%Mi(nred+1)=nred+1; upt%cg_U%nnz=nred
    deallocate(roff,rowcount,next)
  end subroutine build_reduced_hamiltonian

  integer function pair_slot(pairs,npair,a,b,upt)
    type(CGPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer,intent(in)::a,b
    type(OUPT),intent(in)::upt
    integer::i
    do i=1,npair
       if(pairs(i)%a==a .and. pairs(i)%b==b) then; pair_slot=i; return; end if
    end do
    npair=npair+1
    if(npair>size(pairs)) then; pair_slot=0; return; end if
    pairs(npair)%a=a; pairs(npair)%b=b
    allocate(pairs(npair)%v(upt%cg_blocks(a)%nret,upt%cg_blocks(b)%nret)); pairs(npair)%v=(0.0_dp,0.0_dp)
    pair_slot=npair
  end function pair_slot

  ! Variant that allocates v using the FULL block size (nrow x nrow),
  ! used during projection before the energy window is applied.
  integer function pair_slot_full(pairs,npair,a,b,upt)
    type(CGPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer,intent(in)::a,b
    type(OUPT),intent(in)::upt
    integer::i
    do i=1,npair
       if(pairs(i)%a==a .and. pairs(i)%b==b) then; pair_slot_full=i; return; end if
    end do
    npair=npair+1
    if(npair>size(pairs)) then; pair_slot_full=0; return; end if
    pairs(npair)%a=a; pairs(npair)%b=b
    allocate(pairs(npair)%v(upt%cg_blocks(a)%nrow, upt%cg_blocks(b)%nrow))
    pairs(npair)%v=(0.0_dp,0.0_dp)
    pair_slot_full=npair
  end function pair_slot_full

  subroutine add_outer(target, qr, qc, value)
    complex(dp), intent(inout) :: target(:,:)
    complex(dp), intent(in) :: qr(:), qc(:), value
    integer :: i,j
    do j=1,size(qc); do i=1,size(qr)
       target(i,j)=target(i,j)+conjg(qr(i))*value*qc(j)
    end do; end do
  end subroutine add_outer

  subroutine emit_pair(h,p,oa,ob,fmt,next)
    type(CSR),intent(inout)::h
    type(CGPair),intent(in)::p
    integer,intent(in)::oa,ob
    character(1),intent(in)::fmt
    integer,intent(inout)::next(:)
    integer::i,j,k
    if(fmt/='L') then
       do i=1,size(p%v,1); do j=1,size(p%v,2)
          k=next(oa+i-1); h%Mj(k)=ob+j-1; h%M(k)=p%v(i,j); next(oa+i-1)=k+1
       end do; end do
    end if
    if(fmt=='F' .or. fmt=='L') then
       do j=1,size(p%v,2); do i=1,size(p%v,1)
          k=next(ob+j-1); h%Mj(k)=oa+i-1; h%M(k)=conjg(p%v(i,j)); next(ob+j-1)=k+1
       end do; end do
    end if
  end subroutine emit_pair

  subroutine destroy_pairs(pairs)
    type(CGPair),allocatable,intent(inout)::pairs(:)
    integer::i
    if(.not.allocated(pairs)) return
    do i=1,size(pairs); if(associated(pairs(i)%v)) deallocate(pairs(i)%v); end do
    deallocate(pairs)
  end subroutine destroy_pairs

  subroutine cg_lift(upt, reduced, physical)
    type(OUPT),intent(in)::upt
    complex(dp),intent(in)::reduced(:,:)
    complex(dp),intent(out)::physical(:,:)
    integer::i,j,k,off
    physical=(0.0_dp,0.0_dp); off=1
    do i=1,size(upt%cg_blocks)
       if(upt%cg_blocks(i)%nret>0) then
          do j=1,size(upt%cg_blocks(i)%rows)
             physical(upt%cg_blocks(i)%rows(j),:)=matmul(upt%cg_blocks(i)%q(j,:),reduced(off:off+upt%cg_blocks(i)%nret-1,:))
          end do
       end if
       off=off+upt%cg_blocks(i)%nret
    end do
  end subroutine cg_lift

  ! ===========================================================================
  ! IMPROVED COARSE-GRAINING (core + buffer + level-1 acquaintance)
  ! ===========================================================================

  subroutine icg_configure(upt, enabled, nblocks, core_emin, core_emax, &
                            e_buffer, epsilon, imbalance)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled
    integer, intent(in) :: nblocks
    real(dp), intent(in) :: core_emin, core_emax, e_buffer, epsilon, imbalance
    call icg_clear(upt)
    upt%icg_enabled     = enabled
    upt%icg_num_blocks  = nblocks
    upt%icg_core_emin   = core_emin
    upt%icg_core_emax   = core_emax
    upt%icg_e_buffer    = e_buffer
    upt%icg_epsilon     = epsilon
    upt%icg_imbalance   = imbalance
  end subroutine icg_configure

  logical function icg_active(upt)
    type(OUPT), intent(in) :: upt
    icg_active = upt%icg_enabled .and. upt%icg_ready
  end function icg_active

  subroutine icg_get_info(upt, ready, original_dim, reduced_dim, nblocks, cut_fraction)
    type(OUPT), intent(in) :: upt
    logical, intent(out) :: ready
    integer, intent(out) :: original_dim, reduced_dim, nblocks
    real(dp), intent(out) :: cut_fraction
    ready         = upt%icg_ready
    original_dim  = upt%icg_original_dim
    reduced_dim   = upt%icg_reduced_dim
    nblocks       = upt%icg_num_blocks
    cut_fraction  = upt%icg_cut_fraction
  end subroutine icg_get_info

  subroutine icg_clear(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    if (associated(upt%icg_ham%M)) call destroy_matrix(upt%icg_ham)
    if (associated(upt%icg_U%M))   call destroy_matrix(upt%icg_U)
    if (associated(upt%icg_blocks)) then
       do i = 1, size(upt%icg_blocks)
          if (associated(upt%icg_blocks(i)%rows))       deallocate(upt%icg_blocks(i)%rows)
          if (associated(upt%icg_blocks(i)%eval))       deallocate(upt%icg_blocks(i)%eval)
          if (associated(upt%icg_blocks(i)%q))          deallocate(upt%icg_blocks(i)%q)
          if (associated(upt%icg_blocks(i)%evals_full)) deallocate(upt%icg_blocks(i)%evals_full)
          if (associated(upt%icg_blocks(i)%S_full))     deallocate(upt%icg_blocks(i)%S_full)
          if (associated(upt%icg_blocks(i)%retained_idx)) deallocate(upt%icg_blocks(i)%retained_idx)
       end do
       deallocate(upt%icg_blocks)
    end if
    upt%icg_ready        = .false.
    upt%icg_original_dim = 0
    upt%icg_reduced_dim  = 0
    upt%icg_cut_fraction = 0.0_dp
  end subroutine icg_clear

  ! ---------------------------------------------------------------------------
  ! icg_prepare: partition → diagonalize blocks (full S) → select states with
  !   improved criterion → project inter-block couplings → build icg_ham CSR
  ! ---------------------------------------------------------------------------
  subroutine icg_prepare(upt, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(out) :: ierr

    ! Reuse cg_prepare's partition + block-diag machinery by temporarily
    ! mapping icg_* params into cg_* fields, running cg_prepare, then
    ! replacing the keep mask with the improved criterion.
    ! We do NOT call cg_prepare directly to avoid coupling; instead we
    ! reproduce the minimal needed steps inline, calling the same helpers.

    integer :: n, na, nb_atoms, i, j, k, p, q_int, status, nedge, maxedge
    integer :: r, c, br, bc, npair, total_ret, pos
    integer, allocatable :: atom_of(:), local_of(:), label(:), row_of(:)
    integer, allocatable :: counts(:), cursor(:), offsets(:), bsize(:)
    integer(c_int), allocatable :: xadj(:), adjncy(:), vwgt(:), adjwgt(:), part(:)
    real(dp), allocatable :: edge_weight(:)
    real(dp) :: max_weight, all_weight, cut_weight, workspace_mib
    type(CGPair), allocatable :: pairs(:)

    ! keep_mask(b, i) : should state i of block b be retained?
    logical, allocatable :: is_core(:,:), keep_mask(:,:)

    integer :: ia, ib, slot, nred, nnz
    integer, allocatable :: roff(:), rowcount(:), next(:), win_a(:), win_b(:)
    complex(dp), allocatable :: g_full(:,:)
    real(dp) :: threshold

    ierr = 0
    call icg_clear(upt)
    if (.not. upt%icg_enabled) return
    if (num_procs /= 1) then
       ierr = 1; write(*,*) '(icg) MPI not supported'; return
    end if

    na = upt%basis%n_basis
    n  = upt%ham%nrow
    if (na < 1 .or. upt%icg_num_blocks < 1 .or. upt%icg_num_blocks > na) then
       ierr = 2; write(*,*) '(icg) invalid number of blocks'; return
    end if
    if (upt%icg_core_emin >= upt%icg_core_emax .or. upt%icg_e_buffer < 0.0_dp) then
       ierr = 3; write(*,*) '(icg) invalid core window or buffer'; return
    end if
    if (.not. associated(upt%ham%M)) then
       ierr = 4; write(*,*) '(icg) Hamiltonian not initialized'; return
    end if

    ! ---- Build atom→orbital mapping (same as cg_prepare) -------------------
    allocate(atom_of(n), local_of(n), offsets(na+1), bsize(na))
    pos = 1
    do i = 1, na
       offsets(i) = pos
       bsize(i)   = upt%n_spin * upt%basis%n_st(i)
       do j = 1, bsize(i)
          atom_of(pos) = i; local_of(pos) = j; pos = pos + 1
       end do
    end do
    offsets(na+1) = pos
    if (pos-1 /= n) then
       ierr = 5; write(*,*) '(icg) atom/orbital mapping inconsistent'; return
    end if

    ! ---- Build atom graph and METIS partition (same as cg_prepare) ----------
    maxedge = max(1, upt%ham%nnz)
    allocate(edge_weight(maxedge), counts(na), label(na))
    nedge = 0; max_weight = 0.0_dp; counts = 0
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          nedge = nedge + 1
          if (nedge > maxedge) then; ierr = 6; return; end if
          edge_weight(nedge) = abs(upt%ham%M(k))**2
          counts(atom_of(r)) = counts(atom_of(r)) + 1
          counts(atom_of(c)) = counts(atom_of(c)) + 1
          max_weight = max(max_weight, edge_weight(nedge))
       end do
    end do
    if (nedge == 0 .or. max_weight == 0.0_dp) then
       ierr = 7; write(*,*) '(icg) atom graph has no couplings'; return
    end if
    allocate(xadj(na+1), cursor(na), adjncy(2*nedge), adjwgt(2*nedge), vwgt(na), part(na))
    xadj(1) = 0_c_int
    do i = 1, na
       xadj(i+1) = xadj(i) + int(counts(i), c_int)
       cursor(i)  = int(xadj(i)) + 1
       vwgt(i)    = int(bsize(i), c_int)
    end do
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          br = atom_of(r); bc = atom_of(c)
          p  = max(1, nint(abs(upt%ham%M(k))**2 / max_weight * 1000000.0_dp))
          adjncy(cursor(br)) = int(bc-1,c_int); adjwgt(cursor(br)) = int(p,c_int); cursor(br)=cursor(br)+1
          adjncy(cursor(bc)) = int(br-1,c_int); adjwgt(cursor(bc)) = int(p,c_int); cursor(bc)=cursor(bc)+1
       end do
    end do
    status = cg_metis_partition(int(na,c_int), xadj, adjncy, vwgt, adjwgt, &
         int(upt%icg_num_blocks,c_int), int(nint(1000.0_dp*upt%icg_imbalance),c_int), 42_c_int, part)
    if (status /= 0) then
       if (upt%verbose > 0) write(*,*) '(icg) METIS unavailable, using graph-BFS fallback partitioning'
       call cg_graph_partition(na, upt%icg_num_blocks, vwgt, xadj, adjncy, adjwgt, part)
    end if
    do i = 1, na; label(i) = int(part(i)) + 1; end do

    ! ---- cut fraction -------------------------------------------------------
    all_weight = 0.0_dp; cut_weight = 0.0_dp
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_weight = all_weight + abs(upt%ham%M(k))**2
          if (label(atom_of(r)) /= label(atom_of(c))) &
               cut_weight = cut_weight + abs(upt%ham%M(k))**2
       end do
    end do
    upt%icg_cut_fraction = cut_weight / all_weight

    ! ---- Allocate blocks and fill rows arrays --------------------------------
    deallocate(counts)
    allocate(upt%icg_blocks(upt%icg_num_blocks), counts(upt%icg_num_blocks))
    counts = 0
    do i = 1, na; counts(label(i)) = counts(label(i)) + bsize(i); end do
    do i = 1, upt%icg_num_blocks
       upt%icg_blocks(i)%nrow = counts(i)
       allocate(upt%icg_blocks(i)%rows(counts(i)))
    end do
    cursor = 0
    do i = 1, na
       br = label(i)
       do j = offsets(i), offsets(i+1)-1
          cursor(br) = cursor(br) + 1
          upt%icg_blocks(br)%rows(cursor(br)) = j
       end do
    end do

    ! ---- Diagonalize each block fully (store S_full) -------------------------
    allocate(row_of(n)); row_of = 0
    do i = 1, upt%icg_num_blocks
       do j = 1, upt%icg_blocks(i)%nrow
          row_of(upt%icg_blocks(i)%rows(j)) = j
       end do
       ! Reuse diagonalize_block but operating on icg_blocks:
       ! We replicate inline for icg_blocks (can't pass icg vs cg distinction).
       call icg_diagonalize_block(upt, i, row_of, ierr)
       if (ierr /= 0) return
       row_of(upt%icg_blocks(i)%rows) = 0
    end do

    ! ---- Improved keep mask: core + buffer + acquaintance -------------------
    allocate(is_core(upt%icg_num_blocks, maxval(counts)), &
             keep_mask(upt%icg_num_blocks, maxval(counts)))
    is_core   = .false.
    keep_mask = .false.

    ! Step 1 & 2: core and buffer
    do i = 1, upt%icg_num_blocks
       do j = 1, upt%icg_blocks(i)%nrow
          if (.not. associated(upt%icg_blocks(i)%evals_full)) cycle
          associate(e => upt%icg_blocks(i)%evals_full(j))
            if (e >= upt%icg_core_emin .and. e <= upt%icg_core_emax) then
               is_core(i,j)   = .true.
               keep_mask(i,j) = .true.
            else if (e >= upt%icg_core_emin - upt%icg_e_buffer .and. &
                     e <= upt%icg_core_emax + upt%icg_e_buffer) then
               keep_mask(i,j) = .true.
            end if
          end associate
       end do
    end do

    ! Step 3: level-1 acquaintance via inter-block coupling in eigenbasis
    ! We need the transformed couplings g_ab = S_a^T * V_ab * S_b.
    ! Build them on the fly from S_full and the physical Hamiltonian.
    threshold = upt%icg_epsilon * upt%icg_e_buffer

    ! Rebuild row_of for all blocks simultaneously (needed for coupling loop)
    row_of = 0
    do i = 1, upt%icg_num_blocks
       do j = 1, upt%icg_blocks(i)%nrow
          row_of(upt%icg_blocks(i)%rows(j)) = j
       end do
    end do

    ! For each inter-block pair (a,b) compute g = S_a^† V_ab S_b in blocks,
    ! then check coupling from core states of a to unselected states of b and vice versa.
    ! We accumulate |g(i,j)|^2 per (block_a state i, block_b state j).
    ! To avoid allocating a full dense g for every pair, we use a scalar
    ! accumulation: for each physical CSR entry (r→c) crossing block boundary,
    ! add contribution conj(S_a(r_local,i)) * H(r,c) * S_b(c_local,j) to g(i,j).
    ! We only need to check if |g(i,j)|^2 > threshold for core i.
    ! Strategy: for each inter-block CSR entry (r,c), iterate over all core
    ! states i of block_a and all unselected states j of block_b and accumulate.
    ! This is O(n_core * n_unselected) per entry — potentially expensive for
    ! large blocks; acceptable for the problem sizes targeted.
    ! We use a dense (na_core × nb_all) temporary matrix per pair.

    ! Simpler: collect full g_ab dense then check. Allocate per pair.
    ! We track pairs already processed with a simple visited array.
    allocate(pairs(max(1, upt%ham%nnz))); npair = 0
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          if (upt%icg_blocks(ia)%nrow == 0 .or. upt%icg_blocks(ib)%nrow == 0) cycle
          slot = pair_slot_icg(pairs, npair, min(ia,ib), max(ia,ib), upt)
          if (slot == 0) then; ierr = 11; return; end if
          if (ia < ib) then
             call add_outer(pairs(slot)%v, &
                  upt%icg_blocks(ia)%S_full(row_of(r),:), &
                  upt%icg_blocks(ib)%S_full(row_of(c),:), upt%ham%M(k))
          else
             call add_outer(pairs(slot)%v, &
                  upt%icg_blocks(ib)%S_full(row_of(c),:), &
                  upt%icg_blocks(ia)%S_full(row_of(r),:), conjg(upt%ham%M(k)))
          end if
       end do
    end do

    ! Now scan each pair: for core states of block a, find unselected states of b
    ! with |g|^2 > threshold, and mark them as acquaintances.
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       ! core in a → unselected in b
       do j = 1, upt%icg_blocks(ia)%nrow
          if (.not. is_core(ia, j)) cycle
          do k = 1, upt%icg_blocks(ib)%nrow
             if (keep_mask(ib, k)) cycle
             if (abs(pairs(i)%v(j,k))**2 > threshold) keep_mask(ib, k) = .true.
          end do
       end do
       ! core in b → unselected in a (g_ba = g_ab^†)
       do k = 1, upt%icg_blocks(ib)%nrow
          if (.not. is_core(ib, k)) cycle
          do j = 1, upt%icg_blocks(ia)%nrow
             if (keep_mask(ia, j)) cycle
             if (abs(pairs(i)%v(j,k))**2 > threshold) keep_mask(ia, j) = .true.
          end do
       end do
    end do
    call destroy_pairs(pairs)

    ! ---- Apply keep_mask to set nret and allocate q/eval --------------------
    total_ret = 0
    do i = 1, upt%icg_num_blocks
       upt%icg_blocks(i)%nret = count(keep_mask(i, 1:upt%icg_blocks(i)%nrow))
       total_ret = total_ret + upt%icg_blocks(i)%nret
    end do
    if (total_ret == 0) then
       ierr = 9; write(*,*) '(icg) no states retained'; return
    end if
    upt%icg_original_dim = n; upt%icg_reduced_dim = total_ret

    ! Build q/eval for each block using the keep_mask
    do i = 1, upt%icg_num_blocks
       if (upt%icg_blocks(i)%nret == 0) cycle
       allocate(upt%icg_blocks(i)%eval(upt%icg_blocks(i)%nret))
       allocate(upt%icg_blocks(i)%q(upt%icg_blocks(i)%nrow, upt%icg_blocks(i)%nret))
       allocate(upt%icg_blocks(i)%retained_idx(upt%icg_blocks(i)%nret))
       j = 0
       do k = 1, upt%icg_blocks(i)%nrow
          if (.not. keep_mask(i, k)) cycle
          j = j + 1
          upt%icg_blocks(i)%eval(j)         = upt%icg_blocks(i)%evals_full(k)
          upt%icg_blocks(i)%q(:, j)         = upt%icg_blocks(i)%S_full(:, k)
          upt%icg_blocks(i)%retained_idx(j) = k
       end do
    end do
    deallocate(is_core, keep_mask)

    ! ---- Rebuild row_of for all blocks (needed by build_icg_reduced_ham) ----
    row_of = 0
    do i = 1, upt%icg_num_blocks
       do j = 1, upt%icg_blocks(i)%nrow
          row_of(upt%icg_blocks(i)%rows(j)) = j
       end do
    end do

    ! ---- Build reduced Hamiltonian using full S projection then cut ---------
    call build_icg_reduced_hamiltonian(upt, atom_of, label, row_of, pairs, npair, ierr)
    if (ierr /= 0) return
    call destroy_pairs(pairs)

    ! ---- Free S_full after projection ---------------------------------------
    do i = 1, upt%icg_num_blocks
       if (associated(upt%icg_blocks(i)%S_full))     deallocate(upt%icg_blocks(i)%S_full)
       if (associated(upt%icg_blocks(i)%evals_full)) deallocate(upt%icg_blocks(i)%evals_full)
       nullify(upt%icg_blocks(i)%S_full, upt%icg_blocks(i)%evals_full)
    end do

    upt%icg_ready = .true.
    if (upt%verbose > 0) write(*,'(a,i0,a,i0,a,f8.4)') &
         '(icg) dimension ', n, ' -> ', total_ret, ', cut fraction ', upt%icg_cut_fraction

  end subroutine icg_prepare

  ! Diagonalize block ib of icg_blocks (identical logic to diagonalize_block but for icg).
  subroutine icg_diagonalize_block(upt, ib, local, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, local(:)
    integer, intent(out) :: ierr
    integer :: n, k, r, c
    complex(dp), allocatable :: h(:,:)
    real(dp), allocatable :: w(:)
    ierr = 0; n = upt%icg_blocks(ib)%nrow
    if (n == 0) then; upt%icg_blocks(ib)%nret = 0; return; end if
    allocate(h(n,n), w(n)); h = (0.0_dp, 0.0_dp)
    do r = 1, upt%ham%nrow
       if (local(r) == 0) cycle
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (local(c) == 0) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          h(local(r), local(c)) = upt%ham%M(k)
          if (r /= c) h(local(c), local(r)) = conjg(upt%ham%M(k))
       end do
    end do
    call dense_eigh(h, w, ierr)
    if (ierr /= 0) return
    allocate(upt%icg_blocks(ib)%evals_full(n), upt%icg_blocks(ib)%S_full(n,n))
    upt%icg_blocks(ib)%evals_full = w
    upt%icg_blocks(ib)%S_full     = h
    deallocate(h, w)
  end subroutine icg_diagonalize_block

  ! pair_slot variant for icg: allocates v(nrow_a, nrow_b).
  integer function pair_slot_icg(pairs, npair, a, b, upt)
    type(CGPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(in) :: a, b
    type(OUPT), intent(in) :: upt
    integer :: i
    do i = 1, npair
       if (pairs(i)%a == a .and. pairs(i)%b == b) then; pair_slot_icg = i; return; end if
    end do
    npair = npair + 1
    if (npair > size(pairs)) then; pair_slot_icg = 0; return; end if
    pairs(npair)%a = a; pairs(npair)%b = b
    allocate(pairs(npair)%v(upt%icg_blocks(a)%nrow, upt%icg_blocks(b)%nrow))
    pairs(npair)%v = (0.0_dp, 0.0_dp)
    pair_slot_icg = npair
  end function pair_slot_icg

  ! Build reduced Hamiltonian for ICG — same algorithm as build_reduced_hamiltonian
  ! but operating on icg_blocks and icg_ham.
  subroutine build_icg_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(CGPair), allocatable, intent(out) :: pairs(:)
    integer, intent(out) :: npair, ierr
    integer :: i, j, k, r, c, a, b, ia, ib, nnz, pos, slot, nred
    integer, allocatable :: roff(:), rowcount(:), next(:)
    complex(dp), allocatable :: g_full(:,:)
    ierr = 0; nred = upt%icg_reduced_dim
    allocate(roff(upt%icg_num_blocks+1)); roff(1) = 1
    do i = 1, upt%icg_num_blocks; roff(i+1) = roff(i) + upt%icg_blocks(i)%nret; end do
    allocate(pairs(max(1, upt%ham%nnz))); npair = 0

    ! Project with full S
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
          if (a == b) cycle
          if (upt%icg_blocks(a)%nrow == 0 .or. upt%icg_blocks(b)%nrow == 0) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          ia = min(a,b); ib = max(a,b)
          slot = pair_slot_icg(pairs, npair, ia, ib, upt)
          if (slot == 0) then; ierr = 11; return; end if
          if (a < b) then
             call add_outer(pairs(slot)%v, &
                  upt%icg_blocks(a)%S_full(local(r),:), &
                  upt%icg_blocks(b)%S_full(local(c),:), upt%ham%M(k))
          else
             call add_outer(pairs(slot)%v, &
                  upt%icg_blocks(b)%S_full(local(c),:), &
                  upt%icg_blocks(a)%S_full(local(r),:), conjg(upt%ham%M(k)))
          end if
       end do
    end do

    ! Slice to retained states using retained_idx
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       g_full = pairs(i)%v
       deallocate(pairs(i)%v)
       allocate(pairs(i)%v(upt%icg_blocks(a)%nret, upt%icg_blocks(b)%nret))
       do j = 1, upt%icg_blocks(b)%nret
          do k = 1, upt%icg_blocks(a)%nret
             pairs(i)%v(k,j) = g_full(upt%icg_blocks(a)%retained_idx(k), &
                                       upt%icg_blocks(b)%retained_idx(j))
          end do
       end do
       deallocate(g_full)
    end do

    ! Build CSR
    allocate(rowcount(nred), next(nred)); rowcount = 1
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       select case(upt%ham%sparse_fmt)
       case('F')
          rowcount(roff(a):roff(a+1)-1) = rowcount(roff(a):roff(a+1)-1) + upt%icg_blocks(b)%nret
          rowcount(roff(b):roff(b+1)-1) = rowcount(roff(b):roff(b+1)-1) + upt%icg_blocks(a)%nret
       case('L')
          rowcount(roff(b):roff(b+1)-1) = rowcount(roff(b):roff(b+1)-1) + upt%icg_blocks(a)%nret
       case default
          rowcount(roff(a):roff(a+1)-1) = rowcount(roff(a):roff(a+1)-1) + upt%icg_blocks(b)%nret
       end select
    end do
    nnz = sum(rowcount); call create_matrix(upt%icg_ham, nred, nred, nnz)
    upt%icg_ham%sparse_fmt = upt%ham%sparse_fmt; upt%icg_ham%Mi(1) = 1
    do i = 1, nred; upt%icg_ham%Mi(i+1) = upt%icg_ham%Mi(i) + rowcount(i); end do
    next = upt%icg_ham%Mi(1:nred)
    do a = 1, upt%icg_num_blocks
       do i = 1, upt%icg_blocks(a)%nret
          pos = next(roff(a)+i-1)
          upt%icg_ham%Mj(pos) = roff(a)+i-1
          upt%icg_ham%M(pos)  = upt%icg_blocks(a)%eval(i)
          next(roff(a)+i-1) = pos + 1
       end do
    end do
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       call emit_pair(upt%icg_ham, pairs(i), roff(a), roff(b), upt%ham%sparse_fmt, next)
    end do
    upt%icg_ham%nnz = nnz
    call create_matrix(upt%icg_U, nred, nred, nred)
    upt%icg_U%sparse_fmt = 'F'; upt%icg_U%Mi(1) = 1
    do i = 1, nred
       upt%icg_U%Mi(i) = i; upt%icg_U%Mj(i) = i; upt%icg_U%M(i) = (1.0_dp, 0.0_dp)
    end do
    upt%icg_U%Mi(nred+1) = nred + 1; upt%icg_U%nnz = nred
    deallocate(roff, rowcount, next)
  end subroutine build_icg_reduced_hamiltonian

  ! Lift ICG eigenvectors from reduced basis back to physical space.
  subroutine icg_lift(upt, reduced, physical)
    type(OUPT), intent(in) :: upt
    complex(dp), intent(in) :: reduced(:,:)
    complex(dp), intent(out) :: physical(:,:)
    integer :: i, j, off
    physical = (0.0_dp, 0.0_dp); off = 1
    do i = 1, size(upt%icg_blocks)
       if (upt%icg_blocks(i)%nret > 0) then
          do j = 1, size(upt%icg_blocks(i)%rows)
             physical(upt%icg_blocks(i)%rows(j),:) = &
                  matmul(upt%icg_blocks(i)%q(j,:), reduced(off:off+upt%icg_blocks(i)%nret-1,:))
          end do
       end if
       off = off + upt%icg_blocks(i)%nret
    end do
  end subroutine icg_lift


  ! ==========================================================================
  ! ICGN (improved CG + Neumann self-energy correction) routines
  !
  ! Strategy: icgn_prepare is a near-copy of icg_prepare operating on
  ! icgn_* OUPT fields.  After building the reduced Hamiltonian (identical
  ! P-space selection), we extract the H_PQ coupling matrices from the
  ! already-built 'pairs' array (which holds the sliced block-eigenbasis
  ! couplings between retained and discarded states) and add the Neumann
  ! self-energy  Sigma = H_PQ * (E0-D)^-1 * [W*(E0-D)^-1]^order * H_QP
  ! directly into the dense form of icgn_ham before converting back to CSR.
  ! ==========================================================================

  subroutine icgn_configure(upt, enabled, nblocks, core_emin, core_emax, &
       e_buffer, epsilon, selfenergy_order, E0, imbalance, &
       check_convergence, pi_maxiter, pi_tol)
    type(OUPT), intent(inout) :: upt
    logical, intent(in) :: enabled, check_convergence
    integer, intent(in) :: nblocks, selfenergy_order, pi_maxiter
    real(dp), intent(in) :: core_emin, core_emax, e_buffer, epsilon, E0, imbalance, pi_tol
    upt%icgn_enabled = enabled
    upt%icgn_num_blocks = nblocks
    upt%icgn_core_emin = core_emin
    upt%icgn_core_emax = core_emax
    upt%icgn_e_buffer = e_buffer
    upt%icgn_epsilon = epsilon
    upt%icgn_selfenergy_order = selfenergy_order
    upt%icgn_E0 = E0
    upt%icgn_imbalance = imbalance
    upt%icgn_check_convergence = check_convergence
    upt%icgn_pi_maxiter = pi_maxiter
    upt%icgn_pi_tol     = pi_tol
    upt%icgn_sigma_T2   = -1.0_dp
    upt%icgn_pi_converged = .true.
    upt%icgn_ready = .false.
  end subroutine icgn_configure

  logical function icgn_active(upt)
    type(OUPT), intent(in) :: upt
    icgn_active = upt%icgn_enabled .and. upt%icgn_ready
  end function icgn_active

  subroutine icgn_get_info(upt, ready, orig_dim, red_dim, nblocks, cut_frac, &
       sigma_T2, pi_converged)
    type(OUPT), intent(in) :: upt
    logical, intent(out) :: ready, pi_converged
    integer, intent(out) :: orig_dim, red_dim, nblocks
    real(dp), intent(out) :: cut_frac, sigma_T2
    ready        = upt%icgn_ready
    orig_dim     = upt%icgn_original_dim
    red_dim      = upt%icgn_reduced_dim
    nblocks      = upt%icgn_num_blocks
    cut_frac     = upt%icgn_cut_fraction
    sigma_T2     = upt%icgn_sigma_T2
    pi_converged = upt%icgn_pi_converged
  end subroutine icgn_get_info

  subroutine icgn_clear(upt)
    type(OUPT), intent(inout) :: upt
    integer :: i
    if (associated(upt%icgn_ham%M)) call destroy_matrix(upt%icgn_ham)
    if (associated(upt%icgn_U%M))   call destroy_matrix(upt%icgn_U)
    if (associated(upt%icgn_blocks)) then
       do i = 1, size(upt%icgn_blocks)
          if (associated(upt%icgn_blocks(i)%rows))         deallocate(upt%icgn_blocks(i)%rows)
          if (associated(upt%icgn_blocks(i)%eval))         deallocate(upt%icgn_blocks(i)%eval)
          if (associated(upt%icgn_blocks(i)%q))            deallocate(upt%icgn_blocks(i)%q)
          if (associated(upt%icgn_blocks(i)%evals_full))   deallocate(upt%icgn_blocks(i)%evals_full)
          if (associated(upt%icgn_blocks(i)%S_full))       deallocate(upt%icgn_blocks(i)%S_full)
          if (associated(upt%icgn_blocks(i)%retained_idx)) deallocate(upt%icgn_blocks(i)%retained_idx)
       end do
       deallocate(upt%icgn_blocks)
    end if
    upt%icgn_ready        = .false.
    upt%icgn_original_dim = 0
    upt%icgn_reduced_dim  = 0
    upt%icgn_cut_fraction = 0.0_dp
  end subroutine icgn_clear

  ! ---------------------------------------------------------------------------
  ! icgn_prepare: same P-space selection as icg_prepare, but after building the
  ! reduced Hamiltonian we add a Neumann-series self-energy correction for the
  ! discarded Q states. 'pairs' holds the full (nrow_a x nrow_b) coupling
  ! matrices in the block eigenbasis; we reuse them to extract H_PQ and H_QQ.
  ! ---------------------------------------------------------------------------
  subroutine icgn_prepare(upt, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(out) :: ierr

    integer :: n, na, nb_atoms, i, j, k, p, status, nedge, maxedge
    integer :: r, c, br, bc, npair, total_ret, pos
    integer, allocatable :: atom_of(:), local_of(:), label(:), row_of(:)
    integer, allocatable :: counts(:), cursor(:), offsets(:), bsize(:)
    integer(c_int), allocatable :: xadj(:), adjncy(:), vwgt(:), adjwgt(:), part(:)
    real(dp), allocatable :: edge_weight(:)
    real(dp) :: max_weight, all_weight, cut_weight
    type(CGPair), allocatable :: pairs(:)

    logical, allocatable :: is_core(:,:), keep_mask(:,:)

    integer :: ia, ib, slot, nred, nnz
    integer, allocatable :: roff(:), rowcount(:), next(:), win_a(:), win_b(:)
    complex(dp), allocatable :: g_full(:,:)
    real(dp) :: threshold

    ! For Neumann self-energy
    integer :: nb, ord, ip, iq, gp_row, gp_col, q_idx, qb_idx
    integer, allocatable :: block_of_state(:), local_of_state(:), ret_offset(:)
    real(dp), allocatable :: evals_q(:)
    complex(dp), allocatable :: H_dense(:,:), Sigma(:,:)
    complex(dp), allocatable :: amp_vec(:), next_vec(:)
    real(dp) :: E0_used, resolvent_val
    complex(dp) :: amp_val, contrib
    ! pq_list: list of (global_p_in_red, global_q_flat, coupling_value)
    integer, allocatable :: pq_p(:), pq_q(:)
    complex(dp), allocatable :: pq_v(:)
    ! qq_list: Q-Q coupling edges in block eigenbasis
    integer, allocatable :: qq_i(:), qq_j(:)
    complex(dp), allocatable :: qq_v(:)
    integer :: npq, nqq, cnt_pq, cnt_qq
    integer :: nstates_total
    ! block→reduced-row offset
    integer, allocatable :: roff_blk(:)
    ! local state index → retained index (0 if discarded)
    integer, allocatable :: local_ret_idx(:,:)

    ierr = 0
    call icgn_clear(upt)
    if (.not. upt%icgn_enabled) return
    if (num_procs /= 1) then
       ierr = 1; write(*,*) '(icgn) MPI not supported'; return
    end if

    na = upt%basis%n_basis
    n  = upt%ham%nrow
    if (na < 1 .or. upt%icgn_num_blocks < 1 .or. upt%icgn_num_blocks > na) then
       ierr = 2; write(*,*) '(icgn) invalid number of blocks'; return
    end if
    if (upt%icgn_core_emin >= upt%icgn_core_emax .or. upt%icgn_e_buffer < 0.0_dp) then
       ierr = 3; write(*,*) '(icgn) invalid core window or buffer'; return
    end if
    if (.not. associated(upt%ham%M)) then
       ierr = 4; write(*,*) '(icgn) Hamiltonian not initialized'; return
    end if

    ! ---- Build atom→orbital mapping (identical to icg_prepare) -------------
    allocate(atom_of(n), local_of(n), offsets(na+1), bsize(na))
    pos = 1
    do i = 1, na
       offsets(i) = pos
       bsize(i)   = upt%n_spin * upt%basis%n_st(i)
       do j = 1, bsize(i)
          atom_of(pos) = i; local_of(pos) = j; pos = pos + 1
       end do
    end do
    offsets(na+1) = pos
    if (pos-1 /= n) then
       ierr = 5; write(*,*) '(icgn) atom/orbital mapping inconsistent'; return
    end if

    ! ---- METIS partition (identical to icg_prepare) -------------------------
    maxedge = max(1, upt%ham%nnz)
    allocate(edge_weight(maxedge), counts(na), label(na))
    nedge = 0; max_weight = 0.0_dp; counts = 0
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          nedge = nedge + 1
          if (nedge > maxedge) then; ierr = 6; return; end if
          edge_weight(nedge) = abs(upt%ham%M(k))**2
          counts(atom_of(r)) = counts(atom_of(r)) + 1
          counts(atom_of(c)) = counts(atom_of(c)) + 1
          max_weight = max(max_weight, edge_weight(nedge))
       end do
    end do
    if (nedge == 0 .or. max_weight == 0.0_dp) then
       ierr = 7; write(*,*) '(icgn) atom graph has no couplings'; return
    end if
    allocate(xadj(na+1), cursor(na), adjncy(2*nedge), adjwgt(2*nedge), vwgt(na), part(na))
    xadj(1) = 0_c_int
    do i = 1, na
       xadj(i+1) = xadj(i) + int(counts(i), c_int)
       cursor(i)  = int(xadj(i)) + 1
       vwgt(i)    = int(bsize(i), c_int)
    end do
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          br = atom_of(r); bc = atom_of(c)
          p  = max(1, nint(abs(upt%ham%M(k))**2 / max_weight * 1000000.0_dp))
          adjncy(cursor(br)) = int(bc-1,c_int); adjwgt(cursor(br)) = int(p,c_int); cursor(br)=cursor(br)+1
          adjncy(cursor(bc)) = int(br-1,c_int); adjwgt(cursor(bc)) = int(p,c_int); cursor(bc)=cursor(bc)+1
       end do
    end do
    status = cg_metis_partition(int(na,c_int), xadj, adjncy, vwgt, adjwgt, &
         int(upt%icgn_num_blocks,c_int), int(nint(1000.0_dp*upt%icgn_imbalance),c_int), 42_c_int, part)
    if (status /= 0) then
       if (upt%verbose > 0) write(*,*) '(icgn) METIS unavailable, using graph-BFS fallback partitioning'
       call cg_graph_partition(na, upt%icgn_num_blocks, vwgt, xadj, adjncy, adjwgt, part)
    end if
    do i = 1, na; label(i) = int(part(i)) + 1; end do

    all_weight = 0.0_dp; cut_weight = 0.0_dp
    do r = 1, n
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (atom_of(r) == atom_of(c)) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          all_weight = all_weight + abs(upt%ham%M(k))**2
          if (label(atom_of(r)) /= label(atom_of(c))) &
               cut_weight = cut_weight + abs(upt%ham%M(k))**2
       end do
    end do
    upt%icgn_cut_fraction = cut_weight / all_weight

    ! ---- Allocate icgn_blocks and fill rows arrays --------------------------
    deallocate(counts)
    nb = upt%icgn_num_blocks
    allocate(upt%icgn_blocks(nb), counts(nb))
    counts = 0
    do i = 1, na; counts(label(i)) = counts(label(i)) + bsize(i); end do
    do i = 1, nb
       upt%icgn_blocks(i)%nrow = counts(i)
       allocate(upt%icgn_blocks(i)%rows(counts(i)))
    end do
    cursor = 0
    do i = 1, na
       br = label(i)
       do j = offsets(i), offsets(i+1)-1
          cursor(br) = cursor(br) + 1
          upt%icgn_blocks(br)%rows(cursor(br)) = j
       end do
    end do

    ! ---- Diagonalize each block fully (store S_full, evals_full) ------------
    allocate(row_of(n)); row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
       call icgn_diagonalize_block(upt, i, row_of, ierr)
       if (ierr /= 0) return
       row_of(upt%icgn_blocks(i)%rows) = 0
    end do

    ! ---- Improved keep-mask: core + buffer + acquaintance -------------------
    allocate(is_core(nb, maxval(counts)), keep_mask(nb, maxval(counts)))
    is_core = .false.; keep_mask = .false.

    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          if (.not. associated(upt%icgn_blocks(i)%evals_full)) cycle
          associate(e => upt%icgn_blocks(i)%evals_full(j))
            if (e >= upt%icgn_core_emin .and. e <= upt%icgn_core_emax) then
               is_core(i,j)   = .true.
               keep_mask(i,j) = .true.
            else if (e >= upt%icgn_core_emin - upt%icgn_e_buffer .and. &
                     e <= upt%icgn_core_emax + upt%icgn_e_buffer) then
               keep_mask(i,j) = .true.
            end if
          end associate
       end do
    end do

    threshold = upt%icgn_epsilon * upt%icgn_e_buffer
    row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
    end do

    ! Level-1 acquaintance: same as icg_prepare
    allocate(pairs(max(1, upt%ham%nnz))); npair = 0
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          if (upt%icgn_blocks(ia)%nrow == 0 .or. upt%icgn_blocks(ib)%nrow == 0) cycle
          slot = pair_slot_icgn(pairs, npair, min(ia,ib), max(ia,ib), upt)
          if (slot == 0) then; ierr = 11; return; end if
          if (ia < ib) then
             call add_outer(pairs(slot)%v, &
                  upt%icgn_blocks(ia)%S_full(row_of(r),:), &
                  upt%icgn_blocks(ib)%S_full(row_of(c),:), upt%ham%M(k))
          else
             call add_outer(pairs(slot)%v, &
                  upt%icgn_blocks(ib)%S_full(row_of(c),:), &
                  upt%icgn_blocks(ia)%S_full(row_of(r),:), conjg(upt%ham%M(k)))
          end if
       end do
    end do

    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       do j = 1, upt%icgn_blocks(ib)%nrow
          if (is_core(ia, j)) then
             do k = 1, upt%icgn_blocks(ib)%nrow
                if (.not. keep_mask(ib, k)) then
                   if (abs(pairs(i)%v(j,k))**2 > threshold) keep_mask(ib, k) = .true.
                end if
             end do
          end if
       end do
       do k = 1, upt%icgn_blocks(ib)%nrow
          if (is_core(ib, k)) then
             do j = 1, upt%icgn_blocks(ia)%nrow
                if (.not. keep_mask(ia, j)) then
                   if (abs(pairs(i)%v(j,k))**2 > threshold) keep_mask(ia, j) = .true.
                end if
             end do
          end if
       end do
    end do
    call destroy_pairs(pairs)

    ! ---- Apply keep_mask to nret, q, eval, retained_idx --------------------
    total_ret = 0
    do i = 1, nb
       upt%icgn_blocks(i)%nret = count(keep_mask(i, 1:upt%icgn_blocks(i)%nrow))
       total_ret = total_ret + upt%icgn_blocks(i)%nret
    end do
    if (total_ret == 0) then
       ierr = 9; write(*,*) '(icgn) no states retained'; return
    end if
    upt%icgn_original_dim = n
    upt%icgn_reduced_dim  = total_ret

    do i = 1, nb
       if (upt%icgn_blocks(i)%nret == 0) cycle
       allocate(upt%icgn_blocks(i)%eval(upt%icgn_blocks(i)%nret))
       allocate(upt%icgn_blocks(i)%q(upt%icgn_blocks(i)%nrow, upt%icgn_blocks(i)%nret))
       allocate(upt%icgn_blocks(i)%retained_idx(upt%icgn_blocks(i)%nret))
       j = 0
       do k = 1, upt%icgn_blocks(i)%nrow
          if (.not. keep_mask(i, k)) cycle
          j = j + 1
          upt%icgn_blocks(i)%eval(j)         = upt%icgn_blocks(i)%evals_full(k)
          upt%icgn_blocks(i)%q(:, j)         = upt%icgn_blocks(i)%S_full(:, k)
          upt%icgn_blocks(i)%retained_idx(j) = k
       end do
    end do
    deallocate(is_core, keep_mask)

    ! ---- Rebuild row_of for reduced-ham build -------------------------------
    row_of = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nrow
          row_of(upt%icgn_blocks(i)%rows(j)) = j
       end do
    end do

    ! ---- Build reduced Hamiltonian (same as ICG), keeping pairs for Sigma --
    call build_icgn_reduced_hamiltonian(upt, atom_of, label, row_of, pairs, npair, ierr)
    if (ierr /= 0) return

    ! =====================================================================
    ! Self-energy correction  Sigma_{p,p'} via Neumann series
    ! pairs(s)%v is now (nret_a x nret_b): H_PP sliced couplings.
    ! We need H_PQ: from the *unsliced* g_full before retained_idx slicing.
    ! Strategy: re-project directly from S_full to extract P-Q and Q-Q.
    ! =====================================================================

    nred = total_ret
    allocate(H_dense(nred, nred))
    H_dense = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Fill H_dense from icgn_ham CSR (may be triangular — symmetrize immediately)
    do r = 1, upt%icgn_ham%nrow
       do k = upt%icgn_ham%Mi(r), upt%icgn_ham%Mi(r+1)-1
          c = upt%icgn_ham%Mj(k)
          H_dense(r, c) = upt%icgn_ham%M(k)
          H_dense(c, r) = conjg(upt%icgn_ham%M(k))
       end do
    end do

    ! Build roff: global reduced-row offset per block
    allocate(roff_blk(nb+1)); roff_blk(1) = 1
    do i = 1, nb; roff_blk(i+1) = roff_blk(i) + upt%icgn_blocks(i)%nret; end do

    ! Build local_ret_idx(block, local_state) → retained index in block (0=discarded)
    allocate(local_ret_idx(nb, maxval([(upt%icgn_blocks(i)%nrow, i=1,nb)])))
    local_ret_idx = 0
    do i = 1, nb
       do j = 1, upt%icgn_blocks(i)%nret
          local_ret_idx(i, upt%icgn_blocks(i)%retained_idx(j)) = j
       end do
    end do

    ! Count total block states for flat Q-space indexing
    nstates_total = sum([(upt%icgn_blocks(i)%nrow, i=1,nb)])
    allocate(block_of_state(nstates_total), local_of_state(nstates_total), &
         evals_q(nstates_total), ret_offset(nb))
    pos = 0
    do i = 1, nb
       ret_offset(i) = pos
       do j = 1, upt%icgn_blocks(i)%nrow
          pos = pos + 1
          block_of_state(pos) = i
          local_of_state(pos) = j
          evals_q(pos) = upt%icgn_blocks(i)%evals_full(j)
       end do
    end do

    ! Build P-Q edge list: (global_red_row_p, flat_q_idx, coupling g(p,q))
    ! These come from pairs before slicing. Re-project with S_full.
    ! For each inter-block physical CSR entry (r,c) where r in P, c in Q:
    cnt_pq = 0
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          ia = label(atom_of(r)); ib = label(atom_of(c))
          if (ia == ib) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          ! We count all P-Q edges across the full g matrix
          cnt_pq = cnt_pq + upt%icgn_blocks(ia)%nret * &
               (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
          cnt_pq = cnt_pq + upt%icgn_blocks(ib)%nret * &
               (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret)
       end do
    end do
    ! Upper bound — just pre-allocate generously using npair info
    ! Exact P-Q: for each pair, nret_a * (nrow_b - nret_b) + nret_b * (nrow_a - nret_a)
    cnt_pq = 0
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       cnt_pq = cnt_pq + upt%icgn_blocks(ia)%nret * &
            (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
       cnt_pq = cnt_pq + upt%icgn_blocks(ib)%nret * &
            (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret)
    end do
    ! pairs%v now holds the *sliced* P-P block. We need the g_full from
    ! the projection step. Since we already freed g_full inside
    ! build_icgn_reduced_hamiltonian, we need to re-project.
    ! Re-project inter-block couplings to get the full g matrices.
    allocate(pq_p(max(1,cnt_pq)), pq_q(max(1,cnt_pq)))
    allocate(pq_v(max(1,cnt_pq)))
    npq = 0

    ! Re-project: for each coupling pair (a,b), build g_full(nrow_a, nrow_b) fresh.
    ! Then extract P-Q entries.
    do i = 1, npair
       ia = pairs(i)%a; ib = pairs(i)%b
       allocate(g_full(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
       g_full = cmplx(0.0_dp, 0.0_dp, kind=dp)
       ! Accumulate S_a^H * V_ab * S_b from physical ham entries
       do r = 1, upt%ham%nrow
          if (label(atom_of(r)) /= ia) cycle
          do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
             c = upt%ham%Mj(k)
             if (label(atom_of(c)) /= ib) cycle
             if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
             call add_outer(g_full, &
                  upt%icgn_blocks(ia)%S_full(row_of(r),:), &
                  upt%icgn_blocks(ib)%S_full(row_of(c),:), upt%ham%M(k))
          end do
       end do
       ! Extract P(ia)-Q(ib) entries
       do j = 1, upt%icgn_blocks(ia)%nret
          do k = 1, upt%icgn_blocks(ib)%nrow
             if (local_ret_idx(ib, k) /= 0) cycle ! k is in P
             if (abs(g_full(upt%icgn_blocks(ia)%retained_idx(j), k)) < 1.0e-14_dp) cycle
             npq = npq + 1
             if (npq > size(pq_p)) then
                call grow_int_array(pq_p, 2*size(pq_p))
                call grow_int_array(pq_q, 2*size(pq_q))
                call grow_cx_array(pq_v, 2*size(pq_v))
             end if
             pq_p(npq) = roff_blk(ia) + j - 1
             pq_q(npq) = ret_offset(ib) + k
             pq_v(npq) = g_full(upt%icgn_blocks(ia)%retained_idx(j), k)
          end do
       end do
       ! Extract P(ib)-Q(ia) entries (g_ba = g_ab^†)
       do j = 1, upt%icgn_blocks(ib)%nret
          do k = 1, upt%icgn_blocks(ia)%nrow
             if (local_ret_idx(ia, k) /= 0) cycle ! k is in P
             if (abs(g_full(k, upt%icgn_blocks(ib)%retained_idx(j))) < 1.0e-14_dp) cycle
             npq = npq + 1
             if (npq > size(pq_p)) then
                call grow_int_array(pq_p, 2*size(pq_p))
                call grow_int_array(pq_q, 2*size(pq_q))
                call grow_cx_array(pq_v, 2*size(pq_v))
             end if
             pq_p(npq) = roff_blk(ib) + j - 1
             pq_q(npq) = ret_offset(ia) + k
             pq_v(npq) = conjg(g_full(k, upt%icgn_blocks(ib)%retained_idx(j)))
          end do
       end do
       deallocate(g_full)
    end do

    ! Build Q-Q edge list if order >= 1 or convergence check is requested
    if (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence) then
       ! Both directions: cnt_qq * 2 (upper bound)
       cnt_qq = 0
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          cnt_qq = cnt_qq + 2 * (upt%icgn_blocks(ia)%nrow - upt%icgn_blocks(ia)%nret) * &
               (upt%icgn_blocks(ib)%nrow - upt%icgn_blocks(ib)%nret)
       end do
       allocate(qq_i(max(1,cnt_qq)), qq_j(max(1,cnt_qq)))
       allocate(qq_v(max(1,cnt_qq)))
       nqq = 0
       do i = 1, npair
          ia = pairs(i)%a; ib = pairs(i)%b
          allocate(g_full(upt%icgn_blocks(ia)%nrow, upt%icgn_blocks(ib)%nrow))
          g_full = cmplx(0.0_dp, 0.0_dp, kind=dp)
          do r = 1, upt%ham%nrow
             if (label(atom_of(r)) /= ia) cycle
             do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
                c = upt%ham%Mj(k)
                if (label(atom_of(c)) /= ib) cycle
                if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
                call add_outer(g_full, &
                     upt%icgn_blocks(ia)%S_full(row_of(r),:), &
                     upt%icgn_blocks(ib)%S_full(row_of(c),:), upt%ham%M(k))
             end do
          end do
          ! Store BOTH directions (j->k and k->j with conj(v))
          do j = 1, upt%icgn_blocks(ia)%nrow
             if (local_ret_idx(ia, j) /= 0) cycle ! j in P, skip
             do k = 1, upt%icgn_blocks(ib)%nrow
                if (local_ret_idx(ib, k) /= 0) cycle ! k in P, skip
                if (abs(g_full(j, k)) < 1.0e-14_dp) cycle
                nqq = nqq + 1
                if (nqq > size(qq_i)) then
                   call grow_int_array(qq_i, 2*size(qq_i))
                   call grow_int_array(qq_j, 2*size(qq_j))
                   call grow_cx_array(qq_v, 2*size(qq_v))
                end if
                qq_i(nqq) = ret_offset(ia) + j
                qq_j(nqq) = ret_offset(ib) + k
                qq_v(nqq) = g_full(j, k)
                ! Reverse direction: k -> j, conj(v)
                nqq = nqq + 1
                if (nqq > size(qq_i)) then
                   call grow_int_array(qq_i, 2*size(qq_i))
                   call grow_int_array(qq_j, 2*size(qq_j))
                   call grow_cx_array(qq_v, 2*size(qq_v))
                end if
                qq_i(nqq) = ret_offset(ib) + k
                qq_j(nqq) = ret_offset(ia) + j
                qq_v(nqq) = conjg(g_full(j, k))
             end do
          end do
          deallocate(g_full)
       end do
    else
       nqq = 0
    end if

    ! Free S_full now that projections are done
    do i = 1, nb
       if (associated(upt%icgn_blocks(i)%S_full))     deallocate(upt%icgn_blocks(i)%S_full)
       if (associated(upt%icgn_blocks(i)%evals_full)) deallocate(upt%icgn_blocks(i)%evals_full)
       nullify(upt%icgn_blocks(i)%S_full, upt%icgn_blocks(i)%evals_full)
    end do
    call destroy_pairs(pairs)

    ! ---- Convergence check: estimate ||T||_2 = ||(E0-D)^-1 W||_2 -----------
    ! Uses qq_i/qq_j/qq_v already built. Only meaningful when order >= 1.
    upt%icgn_sigma_T2     = -1.0_dp
    upt%icgn_pi_converged = .true.
    E0_used = upt%icgn_E0
    if (E0_used == 0.0_dp) E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp
    if (upt%icgn_check_convergence .and. nqq > 0) then
       call icgn_power_iteration(evals_q, nstates_total, &
            qq_i, qq_j, qq_v, nqq, E0_used, &
            upt%icgn_pi_maxiter, upt%icgn_pi_tol, &
            upt%icgn_sigma_T2, upt%icgn_pi_converged)
       if (.not. upt%icgn_pi_converged) then
          write(*,'(a,i0,a)') &
               '  (icgn) WARNING: power iteration did not converge in ', &
               upt%icgn_pi_maxiter, ' sweeps — ||T||_2 estimate may be unreliable.'
       end if
       if (upt%icgn_sigma_T2 >= 1.0_dp) then
          write(*,'(a,f10.6,a)') &
               '  (icgn) WARNING: ||T||_2 = ', upt%icgn_sigma_T2, &
               ' >= 1 — Neumann series not guaranteed to converge at this E0.'
          if (upt%icgn_selfenergy_order == 0) write(*,'(a)') &
               '  (icgn) NOTE: order=0 does not use the Q-Q coupling W, so this' // &
               ' warning is informational only for this run.'
       else if (upt%verbose > 0) then
          write(*,'(a,f10.6,a,l1)') &
               '  (icgn) ||T||_2 = ', upt%icgn_sigma_T2, &
               ', power-iter converged: ', upt%icgn_pi_converged
          if (upt%icgn_selfenergy_order == 0) write(*,'(a)') &
               '  (icgn) NOTE: order=0 ignores W — norm above is for diagnostics only.'
       end if
    end if

    ! ---- Compute Sigma and add to H_dense ----------------------------------
    allocate(Sigma(nred, nred))
    Sigma = cmplx(0.0_dp, 0.0_dp, kind=dp)
    allocate(amp_vec(nstates_total), next_vec(nstates_total))

    E0_used = upt%icgn_E0
    if (E0_used == 0.0_dp) E0_used = (upt%icgn_core_emin + upt%icgn_core_emax) / 2.0_dp

    ! For each unique P-state, propagate into Q-space and close back
    do ip = 1, nred
       amp_vec = cmplx(0.0_dp, 0.0_dp, kind=dp)
       ! Order-0 propagation: amp_vec(q) = g_{ip,q} * resolvent(q)
       do i = 1, npq
          if (pq_p(i) /= ip) cycle
          q_idx = pq_q(i)
          resolvent_val = 1.0_dp / (E0_used - evals_q(q_idx))
          amp_vec(q_idx) = amp_vec(q_idx) + pq_v(i) * resolvent_val
       end do

       do ord = 0, upt%icgn_selfenergy_order
          if (ord > 0) then
             next_vec = cmplx(0.0_dp, 0.0_dp, kind=dp)
             do i = 1, nqq
                q_idx  = qq_i(i)
                amp_val = amp_vec(q_idx)
                if (abs(amp_val) < 1.0e-14_dp) cycle
                qb_idx = qq_j(i)
                resolvent_val = 1.0_dp / (E0_used - evals_q(qb_idx))
                next_vec(qb_idx) = next_vec(qb_idx) + amp_val * qq_v(i) * resolvent_val
             end do
             amp_vec = next_vec
             if (maxval(abs(amp_vec)) < 1.0e-14_dp) exit
          end if

          ! Close chain: for each Q state with nonzero amp, sum over P' connected to Q
          do i = 1, npq
             q_idx = pq_q(i)
             amp_val = amp_vec(q_idx)
             if (abs(amp_val) < 1.0e-14_dp) cycle
             gp_col = pq_p(i)
             contrib = amp_val * conjg(pq_v(i))
             Sigma(ip, gp_col) = Sigma(ip, gp_col) + contrib
          end do
       end do
    end do

    deallocate(amp_vec, next_vec)
    deallocate(pq_p, pq_q, pq_v)
    if (upt%icgn_selfenergy_order >= 1 .or. upt%icgn_check_convergence) then
       if (allocated(qq_i)) deallocate(qq_i, qq_j, qq_v)
    end if
    deallocate(block_of_state, local_of_state, evals_q, ret_offset)
    deallocate(roff_blk, local_ret_idx)

    ! Add Sigma into H_dense and re-Hermitize
    H_dense = H_dense + Sigma
    do i = 1, nred
       do j = i+1, nred
          H_dense(i,j) = (H_dense(i,j) + conjg(H_dense(j,i))) / 2.0_dp
          H_dense(j,i) = conjg(H_dense(i,j))
       end do
    end do
    deallocate(Sigma)

    ! Rebuild icgn_ham from H_dense
    call destroy_matrix(upt%icgn_ham)
    nnz = count(abs(H_dense) > 1.0e-14_dp)
    call create_matrix(upt%icgn_ham, nred, nred, nnz)
    upt%icgn_ham%sparse_fmt = 'F'  ! always full — H_dense is complete after Sigma+Hermitianize
    upt%icgn_ham%Mi(1) = 1
    k = 0
    do i = 1, nred
       do j = 1, nred
          if (abs(H_dense(i,j)) > 1.0e-14_dp) then
             k = k + 1
             upt%icgn_ham%Mj(k) = j
             upt%icgn_ham%M(k)  = H_dense(i,j)
          end if
       end do
       upt%icgn_ham%Mi(i+1) = k + 1
    end do
    upt%icgn_ham%nnz = k
    deallocate(H_dense)

    ! Build icgn_U (identity)
    call create_matrix(upt%icgn_U, nred, nred, nred)
    upt%icgn_U%sparse_fmt = 'F'; upt%icgn_U%Mi(1) = 1
    do i = 1, nred
       upt%icgn_U%Mj(i) = i; upt%icgn_U%M(i) = (1.0_dp, 0.0_dp)
       upt%icgn_U%Mi(i+1) = i + 1
    end do
    upt%icgn_U%nnz = nred

    upt%icgn_ready = .true.
    if (upt%verbose > 0) write(*,'(a,i0,a,i0,a,f8.4,a,i0)') &
         '(icgn) dimension ', n, ' -> ', total_ret, &
         ', cut fraction ', upt%icgn_cut_fraction, ', order ', upt%icgn_selfenergy_order

    ! Cleanup
    deallocate(atom_of, local_of, offsets, bsize, label, row_of)
    deallocate(counts, edge_weight, xadj, adjncy, adjwgt, vwgt, part, cursor)

  end subroutine icgn_prepare

  ! Diagonalize block ib of icgn_blocks (same logic as icg_diagonalize_block).
  subroutine icgn_diagonalize_block(upt, ib, local, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: ib, local(:)
    integer, intent(out) :: ierr
    integer :: nn, k, r, c
    complex(dp), allocatable :: h(:,:)
    real(dp), allocatable :: w(:)
    ierr = 0; nn = upt%icgn_blocks(ib)%nrow
    if (nn == 0) then; upt%icgn_blocks(ib)%nret = 0; return; end if
    allocate(h(nn,nn), w(nn)); h = (0.0_dp, 0.0_dp)
    do r = 1, upt%ham%nrow
       if (local(r) == 0) cycle
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k)
          if (local(c) == 0) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          h(local(r), local(c)) = upt%ham%M(k)
          if (r /= c) h(local(c), local(r)) = conjg(upt%ham%M(k))
       end do
    end do
    call dense_eigh(h, w, ierr)
    if (ierr /= 0) return
    allocate(upt%icgn_blocks(ib)%evals_full(nn), upt%icgn_blocks(ib)%S_full(nn,nn))
    upt%icgn_blocks(ib)%evals_full = w
    upt%icgn_blocks(ib)%S_full     = h
    deallocate(h, w)
  end subroutine icgn_diagonalize_block

  ! pair_slot variant for icgn: allocates v(nrow_a, nrow_b) using icgn_blocks.
  integer function pair_slot_icgn(pairs, npair, a, b, upt)
    type(CGPair), intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(in) :: a, b
    type(OUPT), intent(in) :: upt
    integer :: i
    do i = 1, npair
       if (pairs(i)%a == a .and. pairs(i)%b == b) then; pair_slot_icgn = i; return; end if
    end do
    npair = npair + 1
    if (npair > size(pairs)) then; pair_slot_icgn = 0; return; end if
    pairs(npair)%a = a; pairs(npair)%b = b
    allocate(pairs(npair)%v(upt%icgn_blocks(a)%nrow, upt%icgn_blocks(b)%nrow))
    pairs(npair)%v = (0.0_dp, 0.0_dp)
    pair_slot_icgn = npair
  end function pair_slot_icgn

  ! Build reduced Hamiltonian for ICGN (same algorithm as build_icg_reduced_hamiltonian).
  subroutine build_icgn_reduced_hamiltonian(upt, atom_of, label, local, pairs, npair, ierr)
    type(OUPT), intent(inout) :: upt
    integer, intent(in) :: atom_of(:), label(:)
    integer, intent(inout) :: local(:)
    type(CGPair), allocatable, intent(inout) :: pairs(:)
    integer, intent(inout) :: npair
    integer, intent(out) :: ierr
    integer :: i, j, k, r, c, a, b, ia, ib, nnz, pos, slot, nred
    integer, allocatable :: roff(:), rowcount(:), next(:)
    complex(dp), allocatable :: g_full(:,:)
    ierr = 0; nred = upt%icgn_reduced_dim
    allocate(roff(upt%icgn_num_blocks+1)); roff(1) = 1
    do i = 1, upt%icgn_num_blocks; roff(i+1) = roff(i) + upt%icgn_blocks(i)%nret; end do
    if (allocated(pairs)) call destroy_pairs(pairs)
    allocate(pairs(max(1, upt%ham%nnz))); npair = 0

    ! Project with full S
    do r = 1, upt%ham%nrow
       do k = upt%ham%Mi(r), upt%ham%Mi(r+1)-1
          c = upt%ham%Mj(k); a = label(atom_of(r)); b = label(atom_of(c))
          if (a == b) cycle
          if (upt%icgn_blocks(a)%nrow == 0 .or. upt%icgn_blocks(b)%nrow == 0) cycle
          if (.not. stored_entry(upt%ham%sparse_fmt, r, c)) cycle
          ia = min(a,b); ib = max(a,b)
          slot = pair_slot_icgn(pairs, npair, ia, ib, upt)
          if (slot == 0) then; ierr = 11; return; end if
          if (a < b) then
             call add_outer(pairs(slot)%v, &
                  upt%icgn_blocks(a)%S_full(local(r),:), &
                  upt%icgn_blocks(b)%S_full(local(c),:), upt%ham%M(k))
          else
             call add_outer(pairs(slot)%v, &
                  upt%icgn_blocks(b)%S_full(local(c),:), &
                  upt%icgn_blocks(a)%S_full(local(r),:), conjg(upt%ham%M(k)))
          end if
       end do
    end do

    ! Slice to retained states using retained_idx
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       g_full = pairs(i)%v
       deallocate(pairs(i)%v)
       allocate(pairs(i)%v(upt%icgn_blocks(a)%nret, upt%icgn_blocks(b)%nret))
       do j = 1, upt%icgn_blocks(b)%nret
          do k = 1, upt%icgn_blocks(a)%nret
             pairs(i)%v(k,j) = g_full(upt%icgn_blocks(a)%retained_idx(k), &
                                       upt%icgn_blocks(b)%retained_idx(j))
          end do
       end do
       deallocate(g_full)
    end do

    ! Build CSR
    allocate(rowcount(nred), next(nred)); rowcount = 1
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       select case(upt%ham%sparse_fmt)
       case('F')
          rowcount(roff(a):roff(a+1)-1) = rowcount(roff(a):roff(a+1)-1) + upt%icgn_blocks(b)%nret
          rowcount(roff(b):roff(b+1)-1) = rowcount(roff(b):roff(b+1)-1) + upt%icgn_blocks(a)%nret
       case('L')
          rowcount(roff(b):roff(b+1)-1) = rowcount(roff(b):roff(b+1)-1) + upt%icgn_blocks(a)%nret
       case default
          rowcount(roff(a):roff(a+1)-1) = rowcount(roff(a):roff(a+1)-1) + upt%icgn_blocks(b)%nret
       end select
    end do
    nnz = sum(rowcount); call create_matrix(upt%icgn_ham, nred, nred, nnz)
    upt%icgn_ham%sparse_fmt = upt%ham%sparse_fmt; upt%icgn_ham%Mi(1) = 1
    do i = 1, nred; upt%icgn_ham%Mi(i+1) = upt%icgn_ham%Mi(i) + rowcount(i); end do
    next = upt%icgn_ham%Mi(1:nred)
    do a = 1, upt%icgn_num_blocks
       do i = 1, upt%icgn_blocks(a)%nret
          pos = next(roff(a)+i-1)
          upt%icgn_ham%Mj(pos) = roff(a)+i-1
          upt%icgn_ham%M(pos)  = upt%icgn_blocks(a)%eval(i)
          next(roff(a)+i-1) = pos + 1
       end do
    end do
    do i = 1, npair
       a = pairs(i)%a; b = pairs(i)%b
       call emit_pair(upt%icgn_ham, pairs(i), roff(a), roff(b), upt%ham%sparse_fmt, next)
    end do
    upt%icgn_ham%nnz = nnz
    deallocate(roff, rowcount, next)
  end subroutine build_icgn_reduced_hamiltonian

  ! Lift ICGN eigenvectors from reduced basis back to physical space.
  subroutine icgn_lift(upt, reduced, physical)
    type(OUPT), intent(in) :: upt
    complex(dp), intent(in) :: reduced(:,:)
    complex(dp), intent(out) :: physical(:,:)
    integer :: i, j, off
    physical = (0.0_dp, 0.0_dp); off = 1
    do i = 1, size(upt%icgn_blocks)
       if (upt%icgn_blocks(i)%nret > 0) then
          do j = 1, size(upt%icgn_blocks(i)%rows)
             physical(upt%icgn_blocks(i)%rows(j),:) = &
                  matmul(upt%icgn_blocks(i)%q(j,:), reduced(off:off+upt%icgn_blocks(i)%nret-1,:))
          end do
       end if
       off = off + upt%icgn_blocks(i)%nret
    end do
  end subroutine icgn_lift

  ! ---------------------------------------------------------------------------
  ! Power iteration to estimate ||T||_2, T = (E0-D)^{-1} W
  ! where D = diag(evals_q restricted to Q) and W = Q-Q coupling matrix
  ! stored as the edge list (qq_i, qq_j, qq_v, nqq), BOTH directions.
  !
  ! Algorithm: standard power method on T^dagger T (same as Julia prototype).
  !   v_0   = random unit vector on Q states with at least one Q-Q edge
  !   u     = T * v_k         (apply_T)
  !   w     = T^dagger * u    (apply_T_adjoint)
  !   v_{k+1} = w / ||w||
  !   sigma_est = ||u|| = ||T v_k||   (Rayleigh quotient)
  !   stop when ||v_{k+1} - v_k|| < tol or k == maxiter
  !
  ! Uses dense arrays indexed by flat state index (faster than Dict).
  ! Only states appearing in the Q-Q graph need nonzero entries.
  subroutine icgn_power_iteration(evals_q, nstates_total, &
       qq_i_arr, qq_j_arr, qq_v_arr, nqq, E0, maxiter, tol, &
       sigma_out, converged_out)
    real(dp),     intent(in)  :: evals_q(:)
    integer,      intent(in)  :: nstates_total
    integer,      intent(in)  :: qq_i_arr(:), qq_j_arr(:)
    complex(dp),  intent(in)  :: qq_v_arr(:)
    integer,      intent(in)  :: nqq, maxiter
    real(dp),     intent(in)  :: E0, tol
    real(dp),     intent(out) :: sigma_out
    logical,      intent(out) :: converged_out

    integer  :: k, e, qi, qj
    real(dp) :: nrm, sigma_est, diff2, res_qi
    complex(dp) :: amp
    complex(dp), allocatable :: v(:), u(:), w(:), v_next(:)

    sigma_out     = 0.0_dp
    converged_out = .true.

    if (nqq == 0) return

    allocate(v(nstates_total), u(nstates_total), w(nstates_total), v_next(nstates_total))

    ! ---- initialise v as random unit vector on Q-Q support -----------------
    call random_number_cx(v, nstates_total, qq_i_arr, nqq)
    nrm = sqrt(real(dot_product(v, v), kind=dp))
    if (nrm < 1.0e-14_dp) then; sigma_out = 0.0_dp; converged_out = .true.; return; end if
    v = v / nrm

    converged_out = .false.
    sigma_est     = 0.0_dp

    do k = 1, maxiter

       ! u = T * v = (E0 - D)^{-1} W v
       ! Step 1: w_tmp = W * v  (sparse matvec over qq edges, both directions stored)
       u = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do e = 1, nqq
          qi = qq_i_arr(e); qj = qq_j_arr(e)
          u(qj) = u(qj) + qq_v_arr(e) * v(qi)
       end do
       ! Step 2: scale by resolvent  u(q) = u(q) / (E0 - evals_q(q))
       do qi = 1, nstates_total
          if (abs(u(qi)) < 1.0e-300_dp) cycle
          res_qi = 1.0_dp / (E0 - evals_q(qi))
          u(qi) = u(qi) * res_qi
       end do

       ! sigma_est = ||u|| = ||T v||
       sigma_est = sqrt(real(dot_product(u, u), kind=dp))

       ! w = T^dagger * u = W^dagger * (E0-D)^{-1} * u
       ! Step 1: scale by resolvent (E0-D is real diagonal, self-adjoint)
       w = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do qi = 1, nstates_total
          if (abs(u(qi)) < 1.0e-300_dp) cycle
          w(qi) = u(qi) / (E0 - evals_q(qi))
       end do
       ! Step 2: apply W^dagger via conjugated edge walk
       v_next = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do e = 1, nqq
          qi = qq_i_arr(e); qj = qq_j_arr(e)
          v_next(qi) = v_next(qi) + conjg(qq_v_arr(e)) * w(qj)
       end do

       ! normalize v_{k+1}
       nrm = sqrt(real(dot_product(v_next, v_next), kind=dp))
       if (nrm < 1.0e-14_dp) then
          sigma_out = 0.0_dp; converged_out = .true.; return
       end if
       v_next = v_next / nrm

       ! convergence: ||v_{k+1} - v_k||
       diff2 = real(dot_product(v_next - v, v_next - v), kind=dp)
       v = v_next

       if (sqrt(diff2) < tol) then
          converged_out = .true.
          exit
       end if
    end do

    sigma_out = sigma_est
    deallocate(v, u, w, v_next)
  end subroutine icgn_power_iteration

  ! Fill v with random complex values on the support of qq edges, zeros elsewhere
  subroutine random_number_cx(v, n, qi_arr, nqq)
    complex(dp), intent(out) :: v(:)
    integer,     intent(in)  :: n, nqq, qi_arr(:)
    real(dp) :: rr, ri
    integer  :: e
    v = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do e = 1, nqq
       call random_number(rr); call random_number(ri)
       v(qi_arr(e)) = cmplx(rr - 0.5_dp, ri - 0.5_dp, kind=dp)
    end do
  end subroutine random_number_cx

  ! Helper: grow integer array
  subroutine grow_int_array(arr, new_size)
    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    integer, allocatable :: tmp(:)
    allocate(tmp(new_size))
    tmp(1:size(arr)) = arr
    call move_alloc(tmp, arr)
  end subroutine grow_int_array

  ! Helper: grow complex array
  subroutine grow_cx_array(arr, new_size)
    complex(dp), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    complex(dp), allocatable :: tmp(:)
    allocate(tmp(new_size))
    tmp(1:size(arr)) = arr
    call move_alloc(tmp, arr)
  end subroutine grow_cx_array

  ! ============================================================================
  ! cg_graph_partition: connectivity-aware fallback partition when METIS is
  ! unavailable.
  !
  ! Algorithm: weighted greedy graph growing (BFS-seeded, priority-queue-free).
  ! Guarantees:
  !   (1) Every atom in a block is reachable from the block seed via bonds —
  !       so no block is a disjoint set of atoms.
  !   (2) All orbitals of an atom go to the same block (works at the atom level,
  !       so the orbital → block mapping is done after by the caller).
  !   (3) Approximate balance: blocks grow until their total orbital count
  !       reaches the target weight (n_orbitals / nblocks).  Remaining atoms
  !       are appended to the last unfilled block.
  !
  ! Inputs:
  !   na         — number of atoms
  !   nblocks    — desired number of blocks
  !   vwgt(na)   — vertex weight of each atom (= number of orbitals, bsize(i))
  !   xadj(na+1) — CSR row pointers of atom adjacency graph (0-based)
  !   adjncy(*)  — CSR column indices (0-based atom indices)
  !   adjwgt(*)  — CSR edge weights (integer, Hamiltonian coupling strength)
  !
  ! Output:
  !   part(na)   — block index (0-based, in [0, nblocks-1]) for each atom
  ! ============================================================================
  subroutine cg_graph_partition(na, nblocks, vwgt, xadj, adjncy, adjwgt, part)
    use, intrinsic :: iso_c_binding, only : c_int
    integer,            intent(in)  :: na, nblocks
    integer(c_int),     intent(in)  :: vwgt(na)
    integer(c_int),     intent(in)  :: xadj(na+1), adjncy(*), adjwgt(*)
    integer(c_int),     intent(out) :: part(na)

    ! --- local ---
    integer :: i, j, atom, nb_atom, blk, seed
    integer :: target_wt, cur_wt, total_wt
    integer :: qhead, qtail, qsize
    integer, allocatable :: queue(:)     ! BFS queue (atom indices, 1-based)
    logical, allocatable :: visited(:)  ! atom already assigned?
    integer :: best_deg, deg, adj_start, adj_end
    integer :: max_adj_wt, wval

    ! ----- compute target weight per block -----
    total_wt = 0
    do i = 1, na
       total_wt = total_wt + int(vwgt(i))
    end do
    target_wt = (total_wt + nblocks - 1) / nblocks   ! ceiling

    allocate(visited(na), queue(na))
    visited = .false.
    part    = int(nblocks - 1, c_int)   ! default: last block (catches unvisited atoms)

    do blk = 0, nblocks - 2    ! assign blocks 0 .. nblocks-2; last gets remainder

       ! ---- pick seed: unvisited atom with highest degree ----
       seed     = -1
       best_deg = -1
       do i = 1, na
          if (visited(i)) cycle
          deg = int(xadj(i+1) - xadj(i))   ! number of inter-atom bonds
          if (deg > best_deg) then
             best_deg = deg
             seed     = i
          end if
       end do
       if (seed == -1) exit    ! all atoms assigned already

       ! ---- BFS grow from seed until block weight reaches target ----
       part(seed) = int(blk, c_int)
       visited(seed) = .true.
       queue(1) = seed
       qhead = 1; qtail = 1; qsize = 1
       cur_wt = int(vwgt(seed))

       do while (qhead <= qtail .and. cur_wt < target_wt)
          atom = queue(qhead); qhead = qhead + 1

          ! Add all unvisited neighbours of this atom to the block, in
          ! decreasing edge-weight order (greedy: strongest bond first).
          ! Cost: O(degree^2) per atom — degree is small for tight-binding.
          adj_start = int(xadj(atom)) + 1    ! convert 0-based xadj to 1-based
          adj_end   = int(xadj(atom+1))

          do while (cur_wt < target_wt)
             ! pick the unvisited neighbour with maximum edge weight
             max_adj_wt = -1
             nb_atom    = -1
             do j = adj_start, adj_end
                i = int(adjncy(j)) + 1    ! 0-based → 1-based
                if (.not. visited(i)) then
                   wval = int(adjwgt(j))
                   if (wval > max_adj_wt) then
                      max_adj_wt = wval
                      nb_atom    = i
                   end if
                end if
             end do
             if (nb_atom == -1) exit     ! no more unvisited neighbours of this atom

             visited(nb_atom) = .true.
             part(nb_atom)    = int(blk, c_int)
             cur_wt           = cur_wt + int(vwgt(nb_atom))
             qtail            = qtail + 1
             queue(qtail)     = nb_atom
          end do
          if (cur_wt >= target_wt) exit
       end do

    end do  ! blk loop

    ! ---- any unvisited atom (disconnected components or remainder) → last block ----
    ! (part already initialised to nblocks-1)

    deallocate(visited, queue)
  end subroutine cg_graph_partition

end module coarse_grain
