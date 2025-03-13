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
MODULE lapack_driver

  USE precision
  USE upt_param
  USE input_output
  use errors
  USE exceptions

  implicit none
  private
  
  public :: lapack

  private :: assemble_dense_H, diagonalize_ham_i
  private :: diagonalize_ham

contains

    subroutine lapack(upt)

      type(oupt) :: upt

      ! LOCALS 
      integer :: num_ev, n_ham,  i_k, i, err
      complex(dp), DIMENSION(:,:), ALLOCATABLE :: ham, hamh
      real(dp), DIMENSION(:), POINTER :: tmp_eigvals
      complex(dp), DIMENSION(:,:), POINTER :: p_eigvectors
      integer :: num_conv, file_num
      CHARACTER(200) :: file_name

      n_ham = upt%ham%nrow
      num_ev = upt%num_vb + upt%num_cb

      !write(*,*) '(lapack) make dense H', n_ham

      ALLOCATE( ham(n_ham,n_ham), STAT= err )
      !ALLOCATE( hamh(n_ham,n_ham), STAT= err )
      IF(err .NE. 0) CALL alloc_error('lapack', 'assemble', 'ham')

      allocate(tmp_eigvals(n_ham), STAT = err)
      IF (err.NE.0) CALL alloc_error('lapack driver','main','eigen_vectors')  

      tmp_eigvals = 0.D0

      !write(*,*) '(lapack) diagonalize',upt%ham%Mij(1)-2
     
      call assemble_dense_H(upt%ham%M, upt%ham%Mj, upt%ham%Mi, &
                                             upt%ham%sparse_fmt, n_ham, ham)

      !hamh = conjg(transpose(ham))
      !if (any( ham .ne. hamh)) stop 'not hermitian!' 
      !open(1000, file='HAM.dat')
      !call printmat_fc(1000,'F12.5',ham,n_ham,n_ham)
      !close(1000)

      !write(*,*) '(lapack) diagonalize'
 
      call diagonalize_ham(HAM, n_ham, tmp_eigvals)
      
      !do i_k = 1,n_ham
      !  write (*,*) i_k, " : ", tmp_eigvals(i_k)
      !end do

      i_k = 1
      DO WHILE (tmp_eigvals(i_k).lt.upt%lambda_vb)
        i_k = i_k + 1
        IF (i_k .eq. size(tmp_eigvals)) THEN
          exit
        END IF
      END DO

      i_k = i_k - 1

      
      IF (i_k .lt. upt%num_vb) THEN
        write(*,*) "Number of found valence band eigenvalues is less than requested: ", i_k, " < ", upt%num_vb
        write(*,*) "Check the atomistic structure, or try to increase value of guess_valence. Current value: ", upt%lambda_vb
        call throw_solve_exception(ERR_LAPACK_DIAG)
      END IF
      write(*,*) "Number of bands lower than guess_valence: ", i_k
      

      !do i = i_k, i_k-upt%num_vb, -1
      !   print*, real(dot_product(ham(:,i),ham(:,i_k)))
      !enddo

      !write(*,*) '(lapack) store valence', i_k, i_k-upt%num_vb+1
      if (.not. associated(upt%eigen_values)) then
         allocate(upt%eigen_values(num_ev), STAT = err)
      end if

      if (.not. associated(upt%eigen_vectors)) then
         allocate(upt%eigen_vectors(n_ham,num_ev), STAT = err)
      end if

      if (.not. associated(upt%particles)) then
         allocate(upt%particles(num_ev), STAT = err)
      end if

      upt%particles(1:upt%num_vb) = -1
      upt%particles(upt%num_vb+1:num_ev) = 1

      DO i = 1, upt%num_vb

         upt%eigen_values(i) = tmp_eigvals(i_k)
         upt%eigen_vectors(:,i) = ham(:,i_k)
         i_k = i_k - 1

         
      END DO

      i_k = i_k + upt%num_vb + 1

      !write(*,*) '(lapack) store conduction', i_k, i_k+upt%num_cb-1

      DO i = upt%num_vb+1, upt%num_vb+upt%num_cb

        if (i_k <= n_ham) then
         upt%eigen_values(i) = tmp_eigvals(i_k)
         upt%eigen_vectors(:,i) = ham(:,i_k)
         i_k = i_k + 1
        end if
         
      END DO
     

      !call diagonalize_ham_i(HAM, n_ham, upt%lambda_vb, upt%lambda_cb, & 
      !                        upt%eigen_values, upt%eigen_vectors, num_conv)


      !write(*,*) '(lapack) found',num_conv,'eigenvalues'

      deallocate(ham,tmp_eigvals)
 
    end subroutine lapack

      
   subroutine Diagonalize_ham(HAM, N, eigval)
     !N - dimension of the Hamiltonian so that ham is NxN matrix
     complex(dp), DIMENSION(:,:) :: ham     
     integer  :: N
     ! outputs
     REAL( dp ), DIMENSION(:) :: eigval

     ! LOCALS
     complex(dp),dimension(:),allocatable   :: WORK
     real(dp), dimension(:),allocatable     :: RWORK
     integer, dimension(:), allocatable     :: IWORK
     integer  :: LWORK, LRWORK, LIWORK, INFO, err

     character (len  = 1 )       :: storage

     storage = 'U'

     LWORK = 2*N + N*N
     LRWORK = 1 + 5*N + 2*N*N 
     LIWORK = 3 + 5*N

     allocate(WORK(LWORK),STAT=err)
     allocate(RWORK(LRWORK),STAT=err)
     allocate(IWORK(LIWORK),STAT=err)
     IF(err .NE. 0) CALL alloc_error('lapack', 'diagonalize_ham', 'work')

     !write(*,*) '(lapack) zheev', ham(N,N)
      
     !call ZHEEV('V',storage, N, ham, N, eigval, WORK,LWORK,RWORK,INFO)
     ! Divide and conquer 
     call ZHEEVD('V',storage, N, ham, N, eigval, WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)


     if ((INFO.ne.0))  call throw_solve_exception(ERR_LAPACK_DIAG)

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

     if ((INFO.ne.0)) call throw_solve_exception(ERR_LAPACK_DIAG)

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

   subroutine assemble_dense_H(M,Mj,Mi,sp_fmt,nrows,ham)
     !subrouutine creates full Hamiltonian matrix 
     !from a sparse matrix representation "a" and "ab"
     !the code is based on  sparse_ham routine algorithm
     COMPLEX ( dp ), DIMENSION(:), POINTER :: M
     INTEGER,        DIMENSION(:), POINTER :: Mj
     INTEGER,        DIMENSION(:), POINTER :: Mi
     CHARACTER(1)                          :: sp_fmt
     INTEGER                               :: nrows

     ! OUTPUT
     complex(dp), dimension(:,:), allocatable :: ham

     ! LOCALS
     integer :: i_val, i_a, i_b, row, col, err

     ham = (0.d0, 0.d0)

     i_val = 0

     SELECT CASE( TRIM(sp_fmt) )

     CASE( 'U', 'L', 'u', 'l' )

        DO i_a = 1, SIZE(Mi)-1
           
           !ham(i_a, i_a) = M(i_a)
           
           !IF(Mij(i_a) .NE. Mij(i_a+1)) THEN
              DO i_b = Mi(i_a), Mi(i_a+1)-1
                 ham(i_a, Mj(i_b)) = M(i_b )
                 ham(Mj(i_b), i_a) = conjg(M(i_b))
              END DO
           !END IF
           
        END DO
     
     CASE ( 'F', 'f' )

        DO i_a = 1, SIZE(Mi)-1
           
           !ham(i_a, i_a) = M(i_a)
           
           !IF(Mij(i_a) .NE. Mij(i_a+1)) THEN
              DO i_b = Mi(i_a), Mi(i_a+1)-1
                 ham(i_a, Mj(i_b)) = M(i_b )
              END DO
           !END IF
           
        END DO
        
     END SELECT


     !open (99,file='H.dat')
     !DO i_a=1, Mij(1)-2

     !   DO i_b=1, Mij(1)-2
     !      if(abs(ham(i_a,i_b)).gt.1e-10) then
     !         WRITE (99,*) i_a, i_b, ham(i_a, i_b)
     !      endif
     !   END DO

     !ENDDO
     !close(99)

     !WRITE (*,*) Mij(1)-2, Mij(Mij(1)-1)

   end subroutine assemble_dense_H

   subroutine printmat_fc(iu,FRM,gs,n,k)

      integer :: n,k,n1,k1
      complex(dp) :: gs(n,n)
      integer :: i1,i2,iu,length
      character(6) :: FRM
      !character(26) :: STR
      
      !STR = "'(a2,"//FRM//",a1,"//FRM//",a2)'"

      !print*, STR

      do i1 = 1, n 
         do i2 = 1, k

            write(iu,'(a2,F8.5,a1,F8.5,a2)', advance='NO') &
                 ' (',real(gs(i1,i2)),',',aimag(gs(i1,i2)),') '

         enddo
         write (iu,*) " "

      enddo

    end subroutine printmat_fc


END MODULE lapack_driver
