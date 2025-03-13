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
module lapack_diag

  use precision
  use errors
  use exceptions

  implicit none
  private

  public :: assemble_full_Hamiltonian, diagonalize_ham_i
  public :: diagonalize_ham

contains

   subroutine Diagonalize_ham(HAM, N, eigval)
     !N - dimension of the Hamiltonian so that ham is NxN matrix
     complex(dp), DIMENSION(:,:) :: ham     
     integer  :: N
     ! outputs
     REAL( dp ), DIMENSION(:) :: eigval

     ! LOCALS
     complex(dp),dimension(:),allocatable   :: WORK
     real(dp), dimension(:),allocatable     :: RWORK
     integer  :: LWORK, INFO, err

     character (len  = 1 )       :: storage

     storage = 'L'

     LWORK = 2*N
     allocate(WORK(LWORK),STAT=err)
     allocate(RWORK(3*N-2),STAT=err)
     IF(err .NE. 0) CALL alloc_error('lapack', 'diagonalize_ham', 'work')

     write(*,*) SIZE(ham)
     write(*,*) SIZE(eigval)

     write(*,*) 'call zheev'

     call ZHEEV('V',storage, N, ham, N, eigval, WORK,LWORK,RWORK,INFO)


     if ((INFO.ne.0)) call throw_solve_exception(ERR_LAPACK_DIAG)

     deallocate(WORK,RWORK)     

   end subroutine Diagonalize_ham


   subroutine Diagonalize_ham_i(HAM, N, Emin, Emax, eigval, eigvec, num_conv)
     !N - dimension of the Hamiltonian so that ham is NxN matrix
     !num_conv - number of eigenvalues found between Emin and Emax
     !Emin and Emax define range for the eigenvalues
     complex(dp), dimension(:,:), allocatable :: ham     
     integer  :: N
     real(dp)  :: Emin,Emax
     ! outputs
     COMPLEX( dp ), DIMENSION( :, : ), POINTER :: eigvec
     REAL   ( dp ), DIMENSION( : ),    POINTER :: eigval
     INTEGER :: num_conv

     ! LOCALS
     complex(dp),dimension(:),allocatable   :: WORK
     real(dp), dimension(:),allocatable     :: RWORK
     real(dp), dimension(:),allocatable     :: W
     integer,dimension(:),allocatable      :: IWORK
     integer, dimension(:),allocatable     :: IFAIL
     integer                               :: INFO,i1,i2,M
     integer                     :: Nmin,Nmax
     complex(dp)                  :: au1,au2,au3
     real(dp)                     :: DLAMCH,ACC
     character (len  = 1 )       :: storage

     storage = 'L'
     
     allocate(WORK(2*N))
     allocate(RWORK(7*N))
     !allocate(W(N))
     allocate(IWORK(5*N))
     allocate(IFAIL(N))

     storage='L'
     ACC=1d-8

     call ZHEEVX('V','V',storage, N, ham, N, Emin, Emax, Nmin, Nmax,&
                     ACC,  num_conv, eigval, eigvec, N, WORK, 2*N,&
                     RWORK, IWORK, IFAIL, INFO)

     if ((INFO.ne.0))  call throw_solve_exception(ERR_LAPACK_DIAG)
  
     deallocate(WORK,RWORK,IWORK,IFAIL)

!!$        do i2=1,num_conv
!!$           au3=(0.d0,0.d0)
!!$           au1=(0.d0,0.d0)
!!$           do i1=1,N
!!$              au1=au1+(abs(eigvec(i1,i2)))**2.d0
!!$              au3=au3+eigvec(i1,i2)
!!$           enddo
!!$
!!$           au3=au3/abs(au3)
!!$           au1=sqrt(au1)*au3
!!$           au2=0_dp
!!$
!!$           do i1=1,N
!!$              eigvec(i1,i2)=eigvec(i1,i2)/au1
!!$           enddo
!!$
!!$        enddo

    
  
   end subroutine Diagonalize_ham_i

   subroutine assemble_full_Hamiltonian(M,Mij,nrows,ham)
     !subrouutine creates full Hamiltonian matrix 
     !from a sparse matrix representation "a" and "ab"
     !the code is based on  sparse_ham routine algorithm
     COMPLEX ( dp ), DIMENSION(:), POINTER :: M
     INTEGER,        DIMENSION(:), POINTER :: Mij
     INTEGER                               :: nrows

     ! OUTPUT
     complex(dp), dimension(:,:), allocatable :: ham

     ! LOCALS
     integer :: i_val, i_a, i_b, row, col, err
     
     ALLOCATE( ham(nrows,nrows), STAT= err )
     IF(err .NE. 0) CALL alloc_error('lapack', 'assemble', 'ham')


     ham = (0.d0, 0.d0)
     i_val = 0

     DO i_a = 1, Mij(1)-2
        ham(i_a, i_a) = M(i_a)

        IF(Mij(i_a) .NE. Mij(i_a+1)) THEN
           DO i_b = Mij(i_a), Mij(i_a+1)-1
              ham(i_a, Mij(i_b)) = M(i_b )
              ham(Mij(i_b), i_a) = conjg(M(i_b))
           END DO     
        END IF

     END DO
 
     open (99,file='H.dat')
     DO i_a=1,Mij(1)-2

        DO i_b=1,i_a
           if(abs(ham(i_a,i_b)).gt.1e-10) then
              WRITE (99,*) i_a, i_b, ham(i_a, i_b)
           endif
        END DO

     ENDDO
     close(99)

     !WRITE (*,*) Mij(1)-2, Mij(Mij(1)-1)

   end subroutine assemble_full_Hamiltonian

 end module lapack_diag
