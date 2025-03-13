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
MODULE feast_driver

USE precision
USE upt_param
USE input_output
USE errors
USE sparse_matrix
USE clock
USE exceptions

IMPLICIT NONE
PRIVATE

PUBLIC :: feast

CONTAINS

SUBROUTINE feast(upt)

type(oupt) :: upt
!COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: A
!INTEGER, DIMENSION(:), ALLOCATABLE :: JA, IA
INTEGER :: nnz, nrow
!INTEGER, DIMENSION(:), POINTER :: pMij
INTEGER, DIMENSION(64) :: feastparam
character(len=1) :: UPLO
REAL(dp) :: epsout, Emin, Emax
INTEGER :: loop, M0, info, M, i, err, num_ev
REAL(dp), DIMENSION(:), ALLOCATABLE :: E, res
complex(dp),dimension(:,:),allocatable :: X ! eigenvectors

UPLO=upt%sparse_format
! pMij => upt%ham%mij
! nnz = pMij(pMij(1) - 1) - 1
! nrow = pMij(1) - 2
! ALLOCATE(A(nnz), stat=err); 
! ALLOCATE(JA(nnz), stat=err); 
! ALLOCATE(IA( nrow + 1), stat=err);
! if(err.ne.0) call throw_init_exception(ERR_ALLOC_ERR) 
 
!WRITE(*,*) 'Converting to csr'

!  call sprs_to_csr(upt%ham%m, pMij, A, JA, IA)

WRITE(*,*) 'Init feast'

call feastinit(feastparam)

 feastparam(1)=1
 
 Emin = upt%lambda_vb 
 Emax= upt%lambda_cb 
 M0= upt%num_cb
 allocate(e(1:M0), stat=err)         ! Eigenvalue
 allocate(X(1:nrow,1:M0), stat=err)  ! Eigenvectors
 allocate(res(1:M0), stat=err)       ! Residual (if needed)
 if(err.ne.0) call throw_init_exception(ERR_ALLOC_ERR) 


call set_clock()

CALL zfeast_sst(UPLO, nrow,  upt%ham%M, upt%ham%Mi, upt%ham%Mj, &
                  feastparam, epsout,loop,Emin,Emax,M0,E,X,M,res,info)

print *,'FEAST OUTPUT INFO',info
if (info==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     print *,'SIMULATION TIME'
     call write_clock()
     print *,'# Search interval [Emin,Emax]',Emin,Emax
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),res(i)
     enddo
else
   call throw_solve_exception(ERR_FEAST_DIAG)
endif

upt%num_vb=0
upt%num_cb=M
num_ev = M

 if (associated(upt%eigen_values)) then
   if( size(upt%eigen_values).lt.num_ev ) then
       deallocate(upt%eigen_values)
       allocate(upt%eigen_values(num_ev), STAT = err)
   end if
 else
  allocate(upt%eigen_values(num_ev), STAT = err)
 end if

 if (associated(upt%eigen_vectors)) then
   if( size(upt%eigen_vectors,2).lt.num_ev ) then
       deallocate(upt%eigen_vectors)
       allocate(upt%eigen_vectors(nrow,num_ev), STAT = err)
   end if
 else
   allocate(upt%eigen_vectors(nrow,num_ev), STAT = err)
 end if

upt%eigen_values(1:m)=E(1:m)
upt%eigen_vectors(:,1:m)=X(:,1:m)

deallocate(e)     ! Eigenvalue
deallocate(X)     ! Eigenvectors
deallocate(res)   ! Residual (if needed)
!DEALLOCATE(A); DEALLOCATE(JA); DEALLOCATE(IA);
!Allocating Hamiltonian in CSR format

END SUBROUTINE feast

 
END MODULE feast_driver
