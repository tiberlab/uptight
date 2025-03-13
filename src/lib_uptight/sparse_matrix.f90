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
MODULE sparse_matrix
  USE globals, only : LST
  USE mpi_globals
  USE precision
  USE exceptions
  implicit none
  private

  public :: CSR_ex, CSR, COO, ELL, DNS, CSB, CSR_real, CSR_real_sp_mx_prec
  public :: create_matrix, destroy_matrix, sprs_to_csr, sprs_to_coo
  public :: coo_to_csr, csb_to_csr
  public :: write_sprs_to_coo, write_csr_to_coo
  public :: read_coo 
  public :: get_nrow, get_nnz, write_coo, write_csr, analyze
  public :: drop, zcheck_nnz, drop_sp_mx_prec, zcheck_nnz_sp_mx_prec
  public :: split_matrix, split_matrix_sp_mx_prec
  !===========================================================================
  !
  ! Matrix stored in row-indexed sparse storage mode.
  ! Numerical recipes p.71
  ! Note: complex version
  !
  ! [ D1 V1 0  V2 ]
  ! [ V3 D2 0  0  ]
  ! [ V4 0  D3 0  ]
  ! [ 0  V5 V6 D4 ]
  ! 
  !     [  1  2  3  4  5  6  7  8  9 10 11] 
  ! M = [ D1 D2 D3 D4 ** V1 V2 V3 V4 V5 V6]
  ! Mij=[  6  8  9 10 12  2  4  1  1  2  3]
  !
  ! nrow = Mij(1) - 2
  ! nnz = Mij(nrow + 1) - 2
  !
  !___________________________________________________________________________
  ! max_num_off_diag  ...  Maximum number offdiagonal elements
  !
  !___________________________________________________________________________
  !
  ! sa, ija   (1..N+1+number_of_off_diagonal_elements)
  !
  !===========================================================================
  ! Additional infos
  TYPE CSR_ex

     INTEGER :: nrow, ncol, offset
     CHARACTER(1) :: sparse_fmt     ! L(ower),U(pper),F(ull)
     !
     COMPLEX ( dp ), DIMENSION( : ), POINTER :: M
     INTEGER,        DIMENSION( : ), POINTER :: Mij
     
  END TYPE CSR_ex
  !===========================================================================
  !
  ! Matrix stored in csr row-indexed sparse storage mode.
  !
  ! [ D1 V1 0  V2 ]
  ! [ V3 D2 0  0  ]
  ! [ V4 0  D3 0  ]
  ! [ 0  V5 V6 D4 ]
  ! 
  !     [  1  2  3  4  5  6  7  8  9 10 ] 
  ! M = [ D1 V1 V2 V3 D2 V4 D3 V5 V6 D4 ]  
  ! Mj= [  1  2  4  1  2  1  3  2  3  4 ]
  ! Mi= [  1  4  6  8 11 ] 
  !
  ! nrow = SIZE(Mi) - 1
  ! nnz = Mi(nrow + 1) - 1
  !
  TYPE CSR
     INTEGER :: nrow, ncol, offset
     CHARACTER(1) :: sparse_fmt     ! L(ower),U(pper),F(ull)
     INTEGER :: nnz
     !
     COMPLEX ( dp ), DIMENSION( : ), POINTER :: M
     INTEGER,        DIMENSION( : ), POINTER :: Mi  !ROW POINTER
     INTEGER,        DIMENSION( : ), POINTER :: Mj  !COL INDEX
  END TYPE CSR

  TYPE CSR_real
     INTEGER :: nrow, ncol, offset
     CHARACTER(1) :: sparse_fmt     ! L(ower),U(pper),F(ull)
     INTEGER :: nnz
     !
     REAL ( dp ),    DIMENSION( : ), POINTER :: M
     INTEGER,        DIMENSION( : ), POINTER :: Mi  !ROW POINTER
     INTEGER,        DIMENSION( : ), POINTER :: Mj  !COL INDEX
  END TYPE CSR_real

  TYPE CSR_real_sp_mx_prec
     INTEGER :: nrow, ncol, offset
     CHARACTER(1) :: sparse_fmt     ! L(ower),U(pper),F(ull)
     INTEGER :: nnz
     !
     REAL (sp),      DIMENSION( : ), POINTER :: M
     INTEGER,        DIMENSION( : ), POINTER :: Mi  !ROW POINTER
     INTEGER,        DIMENSION( : ), POINTER :: Mj  !COL INDEX
  END TYPE CSR_real_sp_mx_prec

  !===========================================================================

  TYPE COO
     INTEGER :: nrow, ncol
     !
     COMPLEX ( dp ), DIMENSION( : ), POINTER :: M
     INTEGER,        DIMENSION( : ), POINTER :: Mi  ! ROW
     INTEGER,        DIMENSION( : ), POINTER :: Mj  ! COL
  END TYPE COO

  TYPE ELL
     INTEGER :: nrow, ncol, offset
     !
     COMPLEX (dp), DIMENSION( :, : ), POINTER :: M
     INTEGER,      DIMENSION( :, : ), POINTER :: Mj  ! COL
     INTEGER,      DIMENSION( : ), POINTER :: rowsize  ! ROW
  END TYPE ELL

  TYPE DNS
      COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: val
  END TYPE DNS   

  !===========================================================================
  TYPE RowBlocks
      TYPE(DNS), DIMENSION(:), ALLOCATABLE :: B
      INTEGER, DIMENSION(:), ALLOCATABLE :: col
  END TYPE RowBlocks

  TYPE CSB
      INTEGER :: nrow, ncol, srtrow, endrow    
      INTEGER, DIMENSION(:), ALLOCATABLE :: nbl
      TYPE(RowBlocks), DIMENSION(:), ALLOCATABLE :: Row      
  END TYPE CSB
  !===========================================================================

  INTERFACE create_matrix
    module procedure create_matrix_csr
    module procedure create_matrix_csr_real
    module procedure create_matrix_csr_real_sp_mx_prec
    module procedure create_matrix_coo
    module procedure create_matrix_csb
    module procedure create_matrix_dns
  END INTERFACE

  INTERFACE destroy_matrix
    module procedure destroy_matrix_csrex
    module procedure destroy_matrix_csr
    module procedure destroy_matrix_csr_real
    module procedure destroy_matrix_csr_real_sp_mx_prec
    module procedure destroy_matrix_coo
    module procedure destroy_matrix_csb
  END INTERFACE

  INTERFACE get_nrow
    module procedure get_nrow_csrex
    module procedure get_nrow_csr
  END INTERFACE

  INTERFACE get_nnz
    module procedure get_nnz_csrex
    module procedure get_nnz_csr
  END INTERFACE

  INTERFACE read_coo
    module procedure read_zcoomat
  END INTERFACE
          

  contains

  subroutine destroy_matrix_csrex(mat)

    TYPE(CSR_ex) :: mat
    
    if(associated(mat%M)) deallocate(mat%M)
    if(associated(mat%Mij)) deallocate(mat%Mij)

  end subroutine destroy_matrix_csrex

  
  subroutine destroy_matrix_csr(mat)

    TYPE(CSR) :: mat
    
    if(associated(mat%M)) deallocate(mat%M)
    if(associated(mat%Mj)) deallocate(mat%Mj)
    if(associated(mat%Mi)) deallocate(mat%Mi)    

  end subroutine destroy_matrix_csr
  
  subroutine destroy_matrix_csr_real(mat)

    TYPE(CSR_real) :: mat
    
    if(associated(mat%M)) deallocate(mat%M)
    if(associated(mat%Mj)) deallocate(mat%Mj)
    if(associated(mat%Mi)) deallocate(mat%Mi)    

  end subroutine destroy_matrix_csr_real

  subroutine destroy_matrix_csr_real_sp_mx_prec(mat)

    TYPE(CSR_real_sp_mx_prec) :: mat
    
    if(associated(mat%M)) deallocate(mat%M)
    if(associated(mat%Mj)) deallocate(mat%Mj)
    if(associated(mat%Mi)) deallocate(mat%Mi)    

  end subroutine destroy_matrix_csr_real_sp_mx_prec


  subroutine create_matrix_csr(mat,nrow,ncol,nnz)

    TYPE(CSR) :: mat
    integer :: nrow,ncol,nnz, ierr

    allocate(mat%M(nnz),stat=ierr)
    allocate(mat%Mj(nnz),stat=ierr)
    allocate(mat%Mi(nrow+1),stat=ierr)
    if (ierr.ne.0) then
        write(*,*) 'create_matrix ALLOCATION ERROR'
        call throw_solve_exception(ERR_GENERAL) 
    endif
    Mat%nrow = nrow    
    Mat%ncol = ncol    
    Mat%nnz  = nnz
    Mat%offset = 0     

  end subroutine create_matrix_csr

  subroutine create_matrix_csr_real(mat,nrow,ncol,nnz)
    TYPE(CSR_real) :: mat
    integer :: nrow,ncol,nnz, ierr

    allocate(mat%M(nnz),stat=ierr)
    allocate(mat%Mj(nnz),stat=ierr)
    allocate(mat%Mi(nrow+1),stat=ierr)
    if (ierr.ne.0) then
        write(*,*) 'create_matrix ALLOCATION ERROR'
        call throw_solve_exception(ERR_GENERAL) 
    endif
    Mat%nrow = nrow    
    Mat%ncol = ncol    
    Mat%nnz  = nnz
    Mat%offset = 0     

  end subroutine create_matrix_csr_real

  subroutine create_matrix_csr_real_sp_mx_prec(mat,nrow,ncol,nnz)

    TYPE(CSR_real_sp_mx_prec) :: mat
    integer :: nrow,ncol,nnz, ierr

    allocate(mat%M(nnz),stat=ierr)
    allocate(mat%Mj(nnz),stat=ierr)
    allocate(mat%Mi(nrow+1),stat=ierr)
    if (ierr.ne.0) then
        write(*,*) 'create_matrix ALLOCATION ERROR'
        print*,'create_matrix error function'
        call throw_solve_exception(ERR_GENERAL) 
    endif
    Mat%nrow = nrow    
    Mat%ncol = ncol    
    Mat%nnz  = nnz
    Mat%offset = 0     

  end subroutine create_matrix_csr_real_sp_mx_prec


  subroutine create_matrix_coo(mat,nrow,ncol,nnz)

    TYPE(COO) :: mat
    integer :: nrow,ncol,nnz,ierr

    allocate(mat%M(nnz),stat=ierr)
    allocate(mat%Mi(nnz),stat=ierr)
    allocate(mat%Mj(nnz),stat=ierr)
    if (ierr.ne.0) then
        write(*,*) 'create_matrix ALLOCATION ERROR'
        call throw_solve_exception(ERR_GENERAL) 
    endif
    Mat%nrow = nrow    
    Mat%ncol = ncol    
    !Mat%nnz  = nnz     

  end subroutine create_matrix_coo
  
  subroutine destroy_matrix_coo(mat)

    TYPE(COO) :: mat
    
    if(associated(mat%M)) deallocate(mat%M)
    if(associated(mat%Mj)) deallocate(mat%Mj)
    if(associated(mat%Mi)) deallocate(mat%Mi)    

  end subroutine destroy_matrix_coo

  
  subroutine create_matrix_dns(M,nrow,ncol)
     TYPE(DNS) :: M
     integer :: nrow, ncol

     integer :: ierr

     ierr=0
     IF (.not.allocated(M%val)) allocate(M%val(nrow,ncol), stat=ierr)
     if (ierr.ne.0) then
        write(*,*) 'create_matrix (DNS) ALLOCATION ERROR'
        call throw_solve_exception(ERR_GENERAL) 
     endif

  end subroutine create_matrix_dns

  subroutine create_matrix_csb(M,nrow,ncol,nbl,szbl)
     TYPE(CSB) :: M
     INTEGER :: nrow, ncol
     INTEGER, DIMENSION(:), OPTIONAL :: nbl
     INTEGER, OPTIONAL :: szbl

     INTEGER :: i, j, ierr

     allocate(M%Row(nrow),stat=ierr)
     allocate(M%nbl(nrow),stat=ierr)
     if (ierr.ne.0) then
        write(*,*) 'create_matrix (Row) ALLOCATION ERROR'
        call throw_solve_exception(ERR_GENERAL) 
     endif

     M%nrow = nrow
     M%ncol = ncol
     M%srtrow = 1
     M%endrow = nrow

     IF (present(nbl)) THEN
       M%nbl = nbl

       do i=1,nrow
         allocate(M%Row(i)%B(nbl(i)),stat=ierr)
         allocate(M%Row(i)%col(nbl(i)),stat=ierr)
         if (ierr.ne.0) then
            write(*,*) 'create_matrix (Row',i,') ALLOCATION ERROR'
            call throw_solve_exception(ERR_GENERAL) 
         endif
         if (present(szbl)) then
            do j=1,nbl(i)
               call create_matrix(M%Row(i)%B(j),szbl,szbl)
            enddo
         endif   
       end do

     END IF

  end subroutine create_matrix_csb

  subroutine destroy_matrix_csb(M)
     TYPE(CSB) :: M

     INTEGER :: i, j

     do i = 1, M%nrow
       do j = 1, M%nbl(i)
          if(allocated(M%Row(i)%B(j)%val)) deallocate(M%Row(i)%B(j)%val)
       enddo
       if(allocated(M%Row(i)%B)) deallocate(M%Row(i)%B)
       if(allocated(M%Row(i)%col)) deallocate(M%Row(i)%col)
     enddo

     if(allocated(M%Row)) deallocate(M%Row)
     if(allocated(M%nbl)) deallocate(M%nbl)

  end subroutine destroy_matrix_csb
  !----------------------------------------------------------------------
  !Convert a sparse matrix from internal row indexed storage
  !to CSR storage (Fortran indexing) 
  !----------------------------------------------------------------------
  SUBROUTINE sprs_to_csr(M, Mij, A, JA, IA)

    !IN data
    INTEGER, DIMENSION(:), POINTER :: Mij
    COMPLEX(dp), DIMENSION(:), POINTER :: M
    COMPLEX(dp) :: A(*)
    INTEGER :: JA(*), IA(*)

    !Work data
    INTEGER nnz, nrow, j, i, l

    !     [  1  2  3  4  5  6  7  8  9 10 11] 
    ! M = [ D1 D2 D3 D4 ** V1 V2 V3 V4 V5 V6]
    ! Mij=[  r o w p n t  | c o l i n d     ]

    nrow = Mij(1) - 2
    nnz = Mij(nrow + 1) - 2

    j = 1

    DO i = 1, nrow

       IA(i) = j
       A(j) = M(i)
       JA(j) = i

       j = j + 1

       DO l = Mij(i), Mij(i + 1) - 1

          A(j) = M(l)
          JA(j) = Mij(l)

          j = j + 1

       ENDDO

    ENDDO

    IA(nrow + 1) = nnz + 1

  END SUBROUTINE sprs_to_csr

  !----------------------------------------------------------------------
  !Convert a COO sparse matrix to CSR storage (Fortran indexing) 
  !----------------------------------------------------------------------
  SUBROUTINE coo_to_csr(M_in,M_out)
    type(COO) :: M_in
    type(CSR) :: M_out

    integer :: i,j,k, k0, nrow, nnz, iad
    integer, dimension(:), pointer :: ir, jc, cind, rpnt
    complex(dp), dimension(:), pointer :: a, ao    
    complex(dp) :: x

    ir => M_in%Mi
    jc => M_in%Mj
    a => M_in%M

    nrow = M_in%nrow
    nnz = size(M_in%M)

    ! WE ASSUME CSR HAS BEEN ALLOCATED !
    
    ao => M_out%M
    cind => M_out%Mj
    rpnt => M_out%Mi


    ! NOW WE CAN DO THE JOB

    do k=1,nrow+1
       rpnt(k) = 0
    end do

    ! determine row-lengths.
    do k= 1, nnz
       rpnt(ir(k)) = rpnt(ir(k))+1
    end do

    ! starting position of each row..
    k = 1
    do j=1,nrow+1
       k0 = rpnt(j)
       rpnt(j) = k
       k = k+k0
    end do

    !  go through the structure  once more. Fill in output matrix.
    do k=1, nnz
       i = ir(k)
       j = jc(k)
       x = a(k)
       iad = rpnt(i)
       ao(iad) =  x
       cind(iad) = j
       rpnt(i) = iad+1
    end do
    ! shift back iao
    do j=nrow,1,-1
       rpnt(j+1) = rpnt(j)
    end do
    rpnt(1) = 1

  END SUBROUTINE coo_to_csr

  !----------------------------------------------------------------------
  !Convert a CSR sparse matrix to COO storage (Fortran indexing) 
  !----------------------------------------------------------------------
  SUBROUTINE csr_to_coo(M_in,M_out)
    type(CSR) :: M_in
    type(COO) :: M_out
    
    
    integer :: i,k, k1, k2, nnz, nrow

    nnz = M_in%nnz
    nrow = M_in%nrow


    do k = 1, nnz
       M_out%M(k) = M_in%M(k)
    end do

    do k=1,nnz
       M_out%Mj(k) = M_in%Mj(k)
    end do
    
    do i = nrow, 1, -1
       k1 = M_in%Mi(i+1)-1
       k2 = M_in%Mi(i)
       do k = k1, k2, -1
          M_out%Mi(k) = i
       enddo
    enddo
    
  end subroutine csr_to_coo
  
  !---------------------------------------------------------------------------------
  !Write the Uptight hamiltonian to ASCII file directly in a coo format
  !---------------------------------------------------------------------------------
  subroutine write_sprs_to_coo(M, Mij, output)

  !IN data
  INTEGER, DIMENSION(:), POINTER :: Mij
  COMPLEX(dp), DIMENSION(:), POINTER :: M
  CHARACTER(LST) :: output

  COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: A
  INTEGER, DIMENSION(:), ALLOCATABLE :: JA, IA
  INTEGER nnz

  nnz = Mij(Mij(1) - 1) - 2
  
  ALLOCATE(A(nnz)); ALLOCATE(JA(nnz)); ALLOCATE(IA(nnz));

  call sprs_to_coo(M, Mij, A, JA, IA)

  call write_coo(output, A, JA, IA, (Mij(1)-2), nnz)

  deallocate(A); deallocate(JA); deallocate(IA);


  end subroutine write_sprs_to_coo


  !---------------------------------------------------------------------------------
  !Write the Uptight hamiltonian to ASCII file directly in a coo format
  !---------------------------------------------------------------------------------
  subroutine write_csr_to_coo(M_csr,workpath)
    type(CSR) :: M_csr
    character(LST) :: workpath

    INTEGER :: nrow, nnz

    type(COO) :: M_coo

    nrow = M_csr%nrow
    nnz = M_csr%nnz

    call create_matrix(M_COO,nrow,nrow,nnz)
    
    call csr_to_coo(M_csr, M_COO)
    
    call write_coo(workpath, M_COO%M, M_COO%Mj, M_COO%Mi, nrow, nnz)
    
    call destroy_matrix(M_COO)

  end subroutine write_csr_to_coo

  
  !---------------------------------------------------------------------------------
  !Convert a sparse matrix from internal row indexed storage
  !to COO storage (Fortran indexing)
  !---------------------------------------------------------------------------------
  SUBROUTINE sprs_to_coo(M, Mij, A, JA, IA)

    !IN data
    INTEGER, DIMENSION(:), POINTER :: Mij
    COMPLEX(dp), DIMENSION(:), POINTER :: M
    COMPLEX(dp) :: A(*)
    INTEGER :: JA(*), IA(*)


    INTEGER nnz, nrow, j, i, l, file_id

    nrow = Mij(1) - 2

    nnz = Mij(Mij(1) - 1) - 2

    !ALLOCATE(A(nnz)); ALLOCATE(JA(nnz)); ALLOCATE(IA(nnz))

    j = 1

    DO i = 1, nrow

       IA(j) = i
       A(j) = M(i)
       JA(j) = i

       j = j + 1

       DO l = Mij(i), Mij(i + 1) - 1

          IA(j) = i
          A(j) = M(l)
          JA(j) = Mij(l)

          j = j + 1

       ENDDO

    ENDDO

  END SUBROUTINE sprs_to_coo
  !*********************************************************************

  SUBROUTINE write_coo(outpath,A, JA, IA, nrow, nnz)
    character(LST) :: outpath
    COMPLEX(dp) :: A(*)
    INTEGER :: JA(*), IA(*)
    INTEGER :: nrow, nnz

    INTEGER i, file_id
    !IF(.not.allocated(A)) then
    !   write(*,*) 'ERROR: trying to print an unallocated matrix'
    !   return
    !END IF
    !Start writing on file
    file_id = 90
    OPEN(unit=file_id, file = trim(outpath)//'Hr.m', form = 'FORMATTED')
    write(file_id,*) '% Size =', nrow, nrow 
    write(file_id,*) '% Nonzeros =', nnz
    write(file_id,*) 'zzz = zeros(',nnz,',3)'
    write(file_id,*) 'zzz = ['
    !NOTE: nnz and nrow are inverted in csr and coo conversion routines 
    !(they served different wrappers....)
    DO i = 1,nnz
       WRITE(file_id,*) IA(i), JA(i), REAL(A(i))
    ENDDO
    write(file_id,*) '];'
    write(file_id,*) 'Mat_0 = spconvert(zzz);'
    CLOSE(file_id)

    OPEN(unit=file_id, file = trim(outpath)//'Hi.m', form = 'FORMATTED')
    write(file_id,*) '% Size =', nrow, nrow 
    write(file_id,*) '% Nonzeros =', nnz
    write(file_id,*) 'zzz = zeros(',nnz,',3)'
    write(file_id,*) 'zzz = ['
    !NOTE: nnz and nrow are inverted in csr and coo conversion routines 
    !(they served different wrappers....)
    DO i = 1,nnz
       WRITE(file_id,*) IA(i), JA(i), aimag(A(i))
    ENDDO
    write(file_id,*) '];'
    write(file_id,*) 'Mat_2 = spconvert(zzz);'
    CLOSE(file_id)

  END SUBROUTINE write_coo

  !*********************************************************************
  SUBROUTINE write_csr(outpath, A, JA, IA)
    character(LST) :: outpath
    COMPLEX(dp), DIMENSION(:) :: A
    INTEGER, DIMENSION(:) :: JA, IA

    INTEGER :: nrow, nnz
    INTEGER i, file_id

    nrow = size(IA)-1
    nnz = size(JA)
    !IF(.not.allocated(A)) then
    !   write(*,*) 'ERROR: trying to print an unallocated matrix'
    !   return
    !END IF
    print*,'printing H on ',trim(outpath)//'H.dat'
    !Start writing on file
    file_id = 90
    OPEN(unit=file_id, file = trim(outpath)//'H.dat', access = 'STREAM')
    !NOTE: nnz and nrow are inverted in csr and coo conversion routines 
    !(they served different wrappers....)
    WRITE(file_id) nrow, nnz
    WRITE(file_id) A
    WRITE(file_id) JA 
    WRITE(file_id) IA 

    CLOSE(file_id)

  END SUBROUTINE write_csr
  !*********************************************************************
 
  function zcheck_nnz(A_csr, i1, i2, j1, j2, sorted) result(nnz)
    
    !*********************************************************************************
    !                                                                                |
    !Input:                                                                          |
    !A_csr: sparse CSR matrix                                                        |
    !i1: starting row                                                                |
    !i2: ending row                                                                  |
    !j1: starting column                                                             |
    !j2: ending column                                                               |
    !sorted: if working with sorted matrix put sorted=1 (faster)                     |
    !                                                                                |
    !Output:                                                                         |
    !nzval: non zero values found in submatrix specified by i1,i2,j1,j2              |
    !                                                                                |
    !*********************************************************************************

    type(CSR_real) :: A_csr
    integer :: i1,i2,j1,j2,sorted,nnz
    integer :: i,j

    !Check i1, i2, j1, j2 validity
    if ((i1.lt.1).or.(i2.gt.A_csr%nrow).or.(j2.lt.j1) &
        & .or.(i2.lt.i1).or.(j1.lt.1).or.(j2.gt.A_csr%ncol)) then
       STOP 'Error in check_nzval: wrong row or column specification'
    endif

    nnz=0

    IF (A_csr%nnz.NE.0) THEN
       do i=1, A_csr%nnz
          
          if (abs(A_csr%M(i)).gt.0) then
             nnz=nnz+1
          ENDIF
          
       end do
       
    ENDIF
    
  end function zcheck_nnz


  function zcheck_nnz_sp_mx_prec(A_csr, i1, i2, j1, j2, sorted) result(nnz)

    !*********************************************************************************
    !                                                                                |
    !Input:                                                                          |
    !A_csr: sparse CSR matrix                                                        |
    !i1: starting row                                                                |
    !i2: ending row                                                                  |
    !j1: starting column                                                             |
    !j2: ending column                                                               |
    !sorted: if working with sorted matrix put sorted=1 (faster)                     |
    !                                                                                |
    !Output:                                                                         |
    !nzval: non zero values found in submatrix specified by i1,i2,j1,j2              |
    !                                                                                |
    !*********************************************************************************

    type(CSR_real_sp_mx_prec) :: A_csr
    integer :: i1,i2,j1,j2,sorted,nnz
    integer :: i,j

    nnz=0

    IF (A_csr%nnz.NE.0) THEN
        do i=1, A_csr%nnz
   
            if (abs(A_csr%M(i)).gt.0) then
            nnz=nnz+1
            ENDIF

        end do

   ENDIF

  end function zcheck_nnz_sp_mx_prec
        



  ! ******************************************************************

  ! NOTE: WHEN merging to || CUDA take from Walter bug fixes (sp)
  !       
  !       merge check_nnz as well 
  ! ******************************************************************  

  ! ******************************************************************
  subroutine drop(A)
     type(CSR_real) :: A
     type(CSR_real) :: C
     type(CSR_real) :: B
     integer :: nnz, i, k, r, j, p
     
     nnz = zcheck_nnz(A,1,A%nrow,1,A%ncol,1)

     !  print*,'zcheck_nnz done'
  
     call create_matrix(B,A%nrow,A%ncol,nnz)

     !    print*,'call create matrix done'

     k = 1
     r = 1
     
     !    print*,'i m here'
     do i=1, A%nrow

       B%Mi(r) = k
       do j = A%Mi(i), (A%Mi(i+1)-1)

           
     ! if( (abs(A%M(i)).gt.0)) then

           !  B%M(k) = A%M(i)
           !  B%Mj(k) = A%Mj(i)

          if( (abs(A%M(j)).gt.0)) then

             B%M(k) = A%M(j)
             B%Mj(k) = A%Mj(j)
             k = k + 1
	   endif
        enddo

         r = r + 1

     enddo
     !    print*,'B%nrow B%nrow+1', B%nrow, B%Mi(B%nrow-2), B%Mi(B%nrow-1), B%Mi(B%nrow), B%Mi(B%nrow+1)
     !     print*,'r k', r, k
     B%Mi(r)=k

     do p=r, A%nrow+1
        B%Mi(p) = k
     end do

     !	do i=1, 50
     !	print*,B%Mi(i)
     !	end do
     
     !        print*,'out of 2 do'
     
     call destroy_matrix(A)
     
     !    print*,'destroy_matrix done'
     
     call create_matrix(A,B%nrow,B%ncol,B%nnz)
     
     !    print*,'create matrix2 done'
     
     do k=1, B%nnz
        A%M(k) = B%M(k)
        !A%Mi = B%Mi
        A%Mj(k) = B%Mj(k)
     end do
     
     do k=1, B%nrow+1
        A%Mi(k) = B%Mi(k)
     end do
     !     print*,'here i m inside sparse_matrix'
     call destroy_matrix(B)
     
   end subroutine drop


   subroutine drop_sp_mx_prec(A)
     !INTEGER, DIMENSION( : ), POINTER :: row_offset
     type(CSR_real_sp_mx_prec) :: A
     type(CSR_real_sp_mx_prec) :: C
     type(CSR_real_sp_mx_prec) :: B
     integer :: nnz, i, k, r, j, p

     !write(*,*) "calling check nnz"
     
    ! if(id .ne. num_procs-1) then
     nnz = zcheck_nnz_sp_mx_prec(A,1,A%nrow,1,A%ncol,1)
     !else if( id .eq. num_procs-1) then
     !nnz = zcheck_nnz(A,1, A%nrow, 1, A%ncol,1)
     !endif

      !write(*,*) "calling check nnz done"

     !call mpi_barrier(upt_comm, ierr)


     call create_matrix(B,A%nrow,A%ncol,nnz)

     !write(*,*) "create matrix done"

     !call mpi_barrier(upt_comm, ierr)
    

     k = 1
     r = 1
     
     !write(*,*) "print id", id 
     !write(*,*) "shift_init, shift_end", shift_init, shift_end, id  
     if(id .ne. num_procs-1) then
     do i=shift_init, shift_end
     !do i=row_offset(id+1)+1, row_offset(id+2)

        !write(*,*) "inside first do", i, id

       B%Mi(r) = k
       do j = A%Mi(i), (A%Mi(i+1)-1)
       
      ! print*, 'j', j, id
 
         !write(*,*) "inside first do", j, id

     ! if( (abs(A%M(i)).gt.0)) then
           !  B%M(k) = A%M(i)
           !  B%Mj(k) = A%Mj(i)

          if( (abs(A%M(j)).gt.0)) then

             B%M(k) = A%M(j)
             B%Mj(k) = A%Mj(j)
             k = k + 1
	   endif
        enddo

         r = r + 1

     enddo
     end if
 
     
     
     if(id .eq. num_procs-1) then
     do i=shift_init, shift_end
     !do i=row_offset(id+1)+1, A%nrow

       B%Mi(r) = k
       do j = A%Mi(i), (A%Mi(i+1)-1)
       
       !print*, 'j', j, id
 
         

     ! if( (abs(A%M(i)).gt.0)) then
           !  B%M(k) = A%M(i)
           !  B%Mj(k) = A%Mj(i)

          if( (abs(A%M(j)).gt.0)) then

             B%M(k) = A%M(j)
             B%Mj(k) = A%Mj(j)
             k = k + 1
	   endif
        enddo

         r = r + 1

     enddo
     end if

     !print*,'B%nrow B%nrow+1', B%nrow, B%Mi(B%nrow-2), B%Mi(B%nrow-1), B%Mi(B%nrow), B%Mi(B%nrow+1), id
    
     B%Mi(r)=k

     do p = r, A%nrow+1
        B%Mi(p)=k
     end do

     call destroy_matrix(A)

    !print*,'destroy_matrix done', id

     call create_matrix(A,B%nrow,B%ncol,B%nnz)

    !print*,'create matrix2 done', id

     do k=1, B%nnz
     A%M(k) = B%M(k)
     !A%Mi = B%Mi
     A%Mj(k) = B%Mj(k)
     end do

     do k=1, B%nrow+1
     A%Mi(k) = B%Mi(k)
     end do

     !print*,'here i m inside sparse_matrix'
     call destroy_matrix(B)

   end subroutine drop_sp_mx_prec

  !*********************************************************************

  SUBROUTINE split_matrix(H,H_real,H_imag)
    TYPE(CSR) :: H
    TYPE(CSR_real) :: H_real, H_imag
    INTEGER :: numberofvalues, n_ham, i
    
    numberofvalues = size(H%M)
    n_ham = size(H%Mi)-1

    call create_matrix(H_real, n_ham, n_ham, numberofvalues)
    call create_matrix(H_imag, n_ham, n_ham, numberofvalues)

    do i=1, numberofvalues
       H_real%M(i) = real(H%M(i))
       H_real%Mj(i) = H%Mj(i)
       H_imag%M(i) = AIMAG(H%M(i))
       H_imag%Mj(i) = H%Mj(i)
    end do
    do i=1, n_ham+1
       H_real%Mi(i) = H%Mi(i)
       H_imag%Mi(i) = H%Mi(i)
    end do
    call drop(H_imag)
    
  END SUBROUTINE split_matrix

  SUBROUTINE split_matrix_sp_mx_prec(H,H_real,H_imag)
    TYPE(CSR) :: H
    TYPE(CSR_real_sp_mx_prec) :: H_real, H_imag
    INTEGER :: numberofvalues, n_ham, i
    !INTEGER, DIMENSION( : ), POINTER :: row_offset
    
    numberofvalues = size(H%M)
    n_ham = size(H%Mi)-1

    call create_matrix(H_real, n_ham, n_ham, numberofvalues)
    call create_matrix(H_imag, n_ham, n_ham, numberofvalues)
  

    do i=1, numberofvalues
       H_real%M(i) = real(H%M(i))
       H_real%Mj(i) = H%Mj(i)
       H_imag%M(i) = AIMAG(H%M(i))
       H_imag%Mj(i) = H%Mj(i)
    end do
    do i=1, n_ham+1
       H_real%Mi(i) = H%Mi(i)
       H_imag%Mi(i) = H%Mi(i)
    end do

    call drop_sp_mx_prec(H_imag)
    
  END SUBROUTINE split_matrix_sp_mx_prec

  !*********************************************************************
  SUBROUTINE get_nrow_csr(mat, dim)
    TYPE(CSR) :: mat    
    INTEGER :: dim
    dim = mat%nrow
  END SUBROUTINE get_nrow_csr

  SUBROUTINE get_nrow_csrex(mat, dim)
    TYPE(CSR_ex) :: mat    
    INTEGER :: dim
    dim = mat%nrow
  END SUBROUTINE get_nrow_csrex


  !*********************************************************************
  SUBROUTINE get_nnz_csr(ham, dim)
    TYPE(CSR) :: ham   
    INTEGER :: nrow,dim
    
    nrow=ham%nrow
    if (associated(ham%M)) then
       dim = ham%Mi(nrow+1)-1
    else
       dim = 0
    endif
  END SUBROUTINE get_nnz_csr

  SUBROUTINE get_nnz_csrex(ham, dim)
    TYPE(CSR_ex) :: ham    
    INTEGER :: dim
    if (associated(ham%M)) then
       dim = ham%Mij(ham%Mij(1) - 1) - 2
    else
       dim = 0
    endif
  END SUBROUTINE get_nnz_csrex


  !*********************************************************************

  SUBROUTINE analyze(ham)
    TYPE(CSR) :: ham  

    INTEGER :: nrow, nnz, i, lrow, minrow, maxrow
    INTEGER, DIMENSION(:), ALLOCATABLE :: stat

    call get_nnz_csr(ham, nnz)

    write(*,*) 'nnz: ',nnz

    nrow = ham%nrow

    minrow = 1000000000
    maxrow = 0
    do i = 1, nrow
       lrow = ham%Mi(i+1) - ham%Mi(i)
       if (lrow.lt.minrow) minrow = lrow
       if (lrow.gt.maxrow) maxrow = lrow
    enddo
  
    write(*,*) 'min row length:',minrow
    write(*,*) 'max row length:',maxrow

    allocate(stat(maxrow))
    stat = 0

    do i=1, nrow
       lrow = ham%Mi(i+1) - ham%Mi(i)
       stat(lrow) = stat(lrow) + 1
    enddo

    do i=minrow,maxrow
      write(*,*) i,'nrows:',stat(i)
    enddo

    deallocate(stat)
  END SUBROUTINE analyze

  SUBROUTINE csr2ell(A, B)
    TYPE(csr) :: A
    TYPE(ell) :: B 

    INTEGER :: nrow, nnz, i, k, str, stp 
    INTEGER :: lrow, minrow, maxrow

    nrow = A%nrow
    minrow = 1000000000
    maxrow = 0
    do i = 1, nrow
       lrow = A%Mi(i+1) - A%Mi(i)
       if (lrow.lt.minrow) minrow = lrow
       if (lrow.gt.maxrow) maxrow = lrow
    enddo

    B%nrow = A%nrow
    B%ncol = A%ncol
    B%offset = A%offset

    allocate(B%M(maxrow,nrow))
    allocate(B%Mj(maxrow,nrow))
    allocate(B%rowsize(nrow)) 

    B%M=(0.d0,0.d0)

    do i = 1, nrow
       str = A%Mi(i)
       stp = A%Mi(i+1) - 1
       B%rowsize(i) = stp-str+1
       do k = str, stp
          B%M(k,i) = A%M(k)
          B%Mj(k,i) = A%Mj(k)
       enddo
    enddo

   END SUBROUTINE csr2ell

   SUBROUTINE csb_to_csr(M, A, cutoff)
      TYPE(CSB) :: M
      TYPE(CSR) :: A
      real(dp) :: cutoff

      INTEGER :: nnz, nrows, i, j, szc, szr, cum_szc, k, r, c, col
      COMPLEX(dp) :: val

      ! Compute nnz, nrows, cols
      nrows = 0 
      nnz = 0
      DO i = M%srtrow, M%endrow
        szr=size(M%Row(i)%B(1)%val,1)
        DO j = 1, M%nbl(i)        
           szc=size(M%Row(i)%B(j)%val,2)
           col = M%Row(i)%col(j)
           !cols(col) = szc
           DO c = 1, szc
              DO r = 1, szr
                 val = M%Row(i)%B(j)%val(r,c) 
 
                 if (abs(val).gt.cutoff) nnz = nnz + 1 
                
              END DO
           END DO
        END DO   
        nrows = nrows + szr
      ENDDO
      
      call create_matrix(A,nrows,M%ncol,nnz)

!      k=1
!      nrows = 1
!      A%Mi(1) = 1
!      ! Loop on the row of blocks
!      DO i = M%srtrow, M%endrow
!        szr = size(M%Row(i)%B(1)%val,1)
!        ! Loop on the rows of each block
!        DO r = 1, szr
!           ! Loop on the col of blocks
!           DO j = 1, M%nbl(i)        
!              col = M%Row(i)%col(j)
!              ! Loop on the cols of each block
!              DO c = 1, szr
!                 val = M%Row(i)%B(j)%val(r,c) 
!                 ! Store only large values and 
!                 ! Compute the stride of each blk assuming 
!                 ! all blocks have the same size (sz).
!                 if (abs(val).gt.cutoff) then
!                    A%M(k) = val
!                    A%Mj(k) = szr*(col-1) + c 
!                    k = k + 1 
!                 endif    
!              END DO
!           END DO
!           nrows = nrows + 1
!           A%Mi(nrows) = k
!        ENDDO
!      ENDDO   
!
!
!
!

      k=1
      nrows = 1
      A%Mi(1) = 1
      ! Loop on the row of blocks
      DO i = M%srtrow, M%endrow
        szr = size(M%Row(i)%B(1)%val,1)
        ! Loop on the rows of each block
        DO r = 1, szr
           ! Loop on the col of blocks
           DO j = 1, M%nbl(i)        
              ! The first row index of the block
              col = M%Row(i)%col(j)
              ! Loop on the cols of each block
              szc = size(M%Row(i)%B(j)%val,2)
              DO c = 1, szc
                 val = M%Row(i)%B(j)%val(r,c) 
                 ! Store only large values and 
                 ! Compute the stride of each blk assuming 
                 ! all blocks have the same size (sz).
                 if (abs(val).gt.cutoff) then
                    A%M(k) = val
                    A%Mj(k) = col + c - 1
!                    A%Mj(k) = szc*(col-1) + c 
                    k = k + 1 
                 endif    
              END DO
           END DO
           nrows = nrows + 1
           A%Mi(nrows) = k
        ENDDO
      ENDDO   


   END SUBROUTINE csb_to_csr
    
   subroutine read_zcoomat(iounit,A)
      integer :: iounit
      real(dp), dimension(:,:) :: A
      
      integer :: i1,i2
      complex(dp) :: mat_el

      do
         read (iounit,end=100) i1,i2,mat_el
         A(i1,i2)=mat_el 
      enddo
100   continue
     
    end subroutine read_zcoomat


END MODULE sparse_matrix
