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
    if (status /= 0) then       ! METIS not available - fallback to simple sequential partitioning
       if (upt%verbose > 0) then
          write(*,*) '(coarse grain) METIS unavailable, using sequential partitioning'
       end if
       ! Simple partition: divide atoms sequentially into blocks
       do i = 1, na
          part(i) = int((i-1) * upt%cg_num_blocks / na, c_int)
       end do
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

end module coarse_grain
