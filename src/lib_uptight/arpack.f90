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
MODULE arpack

USE precision
USE upt_param
USE errors
USE sparse_numrec

IMPLICIT NONE
PRIVATE

PUBLIC :: arpack_dr, arpack_zh

CONTAINS

  SUBROUTINE arpack_dr(upt)

    type(oupt) :: upt
    
    ! ---------------------------------------------------------------------
    ! local variables
    ! ---------------------------------------------------------------------
    integer :: num_ev, num_vb, num_cb, n_ham, err, i_k, file_num, i
    REAL ( dp ),   DIMENSION( : ),  POINTER     :: p_eigen_values
    COMPLEX ( dp ),  DIMENSION( :, : ), POINTER :: p_eigen_vectors
    
    !----------------------------------------------------------------------

    num_cb = upt%num_cb
    num_vb = upt%num_vb
    num_ev = num_vb + num_cb
    n_ham = upt%ham%nrow
    i_k = 1

    allocate( upt%eigen_values(num_ev), STAT = err)
    allocate( upt%eigen_vectors(n_ham,num_ev), STAT = err)

    IF ( upt%num_cb.gt.0 ) THEN

       CALL sprs_shift(upt%ham%M, upt%ham%Mij, upt%lambda_cb)

       p_eigen_values => upt%eigen_values(num_vb+1:num_ev)
       p_eigen_vectors => upt%eigen_vectors(:,num_vb+1:num_ev)

       CALL arpack_zh(n_ham, num_cb, upt%long_iter, upt%max_iter, upt%long_tol, & 
            p_eigen_values, p_eigen_vectors, upt%ham%M, upt%ham%Mij, upt%ham%sparse_fmt)

       CALL sprs_shift(upt%ham%M, upt%ham%Mij, -upt%lambda_cb)

       upt%eigen_values(num_vb+1:num_ev) = &
                    sqrt(upt%eigen_values(num_vb+1:num_ev)) + upt%lambda_cb

    END IF

    IF ( upt%num_vb.gt.0 ) THEN

       CALL sprs_shift(upt%ham%M, upt%ham%Mij, upt%lambda_vb)

       p_eigen_values => upt%eigen_values(1:num_vb)
       p_eigen_vectors => upt%eigen_vectors(:,1:num_vb)

       CALL arpack_zh(n_ham, num_vb, upt%long_iter, upt%max_iter, upt%long_tol, &
            p_eigen_values, p_eigen_vectors, upt%ham%M, upt%ham%Mij, upt%ham%sparse_fmt)

       CALL sprs_shift(upt%ham%M, upt%ham%Mij, -upt%lambda_vb)

       upt%eigen_values(1:num_vb) = &
                    -sqrt(upt%eigen_values(1:num_vb)) + upt%lambda_vb

    END IF
    
  END SUBROUTINE arpack_dr
  
  !====================================================================================

  SUBROUTINE arpack_zh(maxn,maxnev,maxncv,max_arn,tol,eigen_values,eigen_vectors, &
                                                                          M,Mij,sp_fmt)
    !     %----------------------------------------------------------%
    !     | Define maximum dimensions                                |
    !     | for all arrays.                                          |
    !     | MAXN:   dimension of the matrix A                        |
    !     | MAXNEV: Maximum N of eigenvalues computed (0<NEV<N-1)    |
    !     | MAXNCV: Maximum N of Column Vectors per iteration        |
    !     | MAX_ARN: Maximum N of Arnoldi Iterations per Vector      |
    !     |                                                          |
    !     | NCV - NEV >= 2                                           |
    !     %----------------------------------------------------------%
    INTEGER maxn, maxnev, maxncv, max_arn
    REAL ( dp ) :: tol
    REAL ( dp ), DIMENSION(:), POINTER :: eigen_values ! d
    COMPLEX ( dp ), DIMENSION(:,:), POINTER :: eigen_vectors ! v
    INTEGER, DIMENSION(:), POINTER :: Mij
    COMPLEX ( dp ), DIMENSION(:), POINTER :: M
    CHARACTER(1)  :: sp_fmt
    COMPLEX ( dp ), DIMENSION(:), POINTER :: p_eig_vectors
    !   parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
    !
    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
    !
    INTEGER   ldv
    INTEGER, DIMENSION(11) :: iparam
    INTEGER, DIMENSION(14) :: ipntr
    LOGICAL, DIMENSION(maxncv) :: select
    COMPLEX ( dp ), DIMENSION(:,:), POINTER :: eig_vectors
    COMPLEX ( dp ), DIMENSION(:), POINTER :: temp
    COMPLEX ( dp ), DIMENSION(:), POINTER :: eig_values, workd! d
    COMPLEX ( dp ), DIMENSION(:), POINTER :: p_ax, workev, resid, workl
    COMPLEX ( dp ), DIMENSION(:), POINTER :: p_workd
    REAL ( dp ), DIMENSION(:), POINTER :: rwork
    REAL ( dp ), DIMENSION(:,:), POINTER :: rd
    !
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    !
    CHARACTER  bmat*1, which*2 
    INTEGER   ido, nev, ncv, lworkl, info, j, ierr, nconv
    INTEGER   maxitr, ishfts, mode, err, n
    COMPLEX ( dp ) sigma
    LOGICAL  rvec,  startflag
    
    REAL( dp) :: dznrm2, dlapy2 
    EXTERNAL dznrm2, zaxpy, dlapy2
    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    ! 
    
    !    ALLOCATE WORK SPACE 
    
    ALLOCATE( eig_values(maxnev), STAT = err )
    ALLOCATE( eig_vectors(maxn,maxncv), STAT = err )
    ALLOCATE( temp(maxn), STAT = err )
    ALLOCATE( p_ax(maxn), STAT = err )
    ALLOCATE( workd(3*maxn), STAT = err )
    ALLOCATE( workev(3*maxncv), STAT = err )
    ALLOCATE( resid(maxn), STAT = err )
    ALLOCATE( workl(3*maxncv*maxncv+5*maxncv), STAT = err )
    ALLOCATE( rwork(maxncv), STAT = err )
    ALLOCATE( rd(maxncv,3), STAT = err )

    IF (err.ne.0) call alloc_error( 'arpack', 'arpack', 'work' ) 

    nev   = maxnev
    ncv   = maxncv
    ldv = maxn
    n = maxn
    bmat  = 'I'
    which = 'SM'
    !
    !     %---------------------------------------------------%
    !     | The work array WORKL is used in ZNAUPD as         | 
    !     | workspace.  Its dimension LWORKL is set as        |
    !     | illustrated below.  The parameter TOL determines  |
    !     | the stopping criterion. If TOL<=0, machine        |
    !     | precision is used.  The variable IDO is used for  |
    !     | reverse communication, and is initially set to 0. |
    !     | Setting INFO=0 indicates that a random vector is  |
    !     | generated to start the ARNOLDI iteration.         | 
    !     %---------------------------------------------------%
    !
    lworkl  = 3*ncv**2+5*ncv 
    ido    = 0
    info   = 0
    !
    !     %---------------------------------------------------%
    !     | This program uses exact shift with respect to     |
    !     | the current Hessenberg matrix (IPARAM(1) = 1).    |
    !     | IPARAM(3) specifies the maximum number of Arnoldi |
    !     | iterations allowed.  Mode 1 of ZNAUPD is used     |
    !     | (IPARAM(7) = 1). All these options can be changed |
    !     | by the user. For details see the documentation in |
    !     | ZNAUPD.                                           |
    !     %---------------------------------------------------%
    !
    ishfts = 1
    maxitr = max_arn
    mode   = 1
    
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(4) = 1
    iparam(7) = mode 
    !
    !     %-------------------------------------------%
    !     | M A I N   L O O P (Reverse communication) | 
    !     %-------------------------------------------%
    !
    startflag = .TRUE.
    DO WHILE (ido .eq. -1 .or. ido .eq. 1 .or. startflag)
       startflag = .FALSE.
       !
       !        %---------------------------------------------%
       !        | Repeatedly call the routine ZNAUPD and take |
       !        | actions indicated by parameter IDO until    |
       !        | either convergence is indicated or maxitr   |
       !        | has been exceeded.                          |
       !        %---------------------------------------------%
       !
       
       CALL znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
            eig_vectors, ldv, iparam, &
            ipntr, workd, workl, lworkl, rwork, info )
       
       
       IF (ido .eq. -1 .or. ido .eq. 1) then
          !
          !           %-------------------------------------------%
          !           | Perform matrix vector multiplication      |
          !           |                y <--- OP*x                |
          !           | The user should supply his/her own        |
          !           | matrix vector multiplication routine here |
          !           | that takes workd(ipntr(1)) as the input   |
          !           | vector, and return the matrix vector      |
          !           | product to workd(ipntr(2)).               | 
          !           %-------------------------------------------%
          !
          
          p_workd => workd(ipntr(1):ipntr(1)+maxn-1)

          call sprs_ax( M, Mij, sp_fmt, p_workd, temp ) 
         
          p_workd => workd(ipntr(2):ipntr(2)+maxn-1)

          call sprs_ax( M, Mij, sp_fmt, temp, p_workd ) 
          
       END IF
       
       !
       !           %-----------------------------------------%
       !           | L O O P   B A C K to call ZNAUPD again. |
       !           %-----------------------------------------%
       !

    END DO

    WRITE (*,*) 'END OF ARPACK LOOP'
    !     %----------------------------------------%
    !     | Either we have convergence or there is |
    !     | an error.                              |
    !     %----------------------------------------%
    !
    IF ( info .lt. 0 ) THEN
       !
       !        %--------------------------%
       !        | Error message, check the |
       !        | documentation in ZNAUPD  |
       !        %--------------------------%
       !
       WRITE(*,*) ' '
       WRITE(*,*) ' Error with _naupd, info = ', info
       WRITE(*,*) ' Check the documentation of _naupd'
       WRITE(*,*) ' '
       
    ELSE
       !
       !        %-------------------------------------------%
       !        | No fatal errors occurred.                 |
       !        | Post-Process using ZNEUPD.                |
       !        |                                           |
       !        | Computed eigenvalues may be extracted.    |
       !        |                                           |
       !        | Eigenvectors may also be computed now if  |
       !        | desired.  (indicated by rvec = .true.)    |
       !        %-------------------------------------------%
       !
       rvec = .true.
       
       call zneupd (rvec, 'A', select, eig_values, eig_vectors, ldv, &
                                   sigma, workev, bmat, maxn, which, &
            nev, tol, resid, ncv, eig_vectors, ldv, iparam, ipntr, &
                                   workd, workl, lworkl, rwork, ierr)
       
       !        %----------------------------------------------%
       !        | Eigenvalues are returned in the one          |
       !        | dimensional array D.  The corresponding      |
       !        | eigenvectors are returned in the first NCONV |
       !        | (=IPARAM(5)) columns of the two dimensional  | 
       !        | array V if requested.  Otherwise, an         |
       !        | orthogonal basis for the invariant subspace  |
       !        | corresponding to the eigenvalues in D is     |
       !        | returned in V.                               |
       !        %----------------------------------------------%
       !
       IF ( ierr .ne. 0) THEN
          ! 
          !           %------------------------------------%
          !           | Error condition:                   |
          !           | Check the documentation of ZNEUPD. |
          !           %------------------------------------%
          !
          WRITE(*,*) ' '
          WRITE(*,*) ' Error with _neupd, info = ', ierr
          WRITE(*,*) ' Check the documentation of _neupd. '
          WRITE(*,*) ' '
          STOP

       ELSE
          
          nconv = iparam(5)
          DO j=1, nconv
             !
             !               %---------------------------%
             !               | Compute the residual norm |
             !               |                           |
             !               |   ||  A*x - lambda*x ||   |
             !               |                           |
             !               | for the NCONV accurately  |
             !               | computed eigenvalues and  |
             !               | eigenvectors.  (iparam(5) |
             !               | indicates how many are    |
             !               | accurate to the requested |
             !               | tolerance)                |
             !               %---------------------------%
             !
             p_eig_vectors => eig_vectors(:,j)

             call sprs_ax( M, Mij, sp_fmt, p_eig_vectors, temp ) 

             call sprs_ax(M, Mij, sp_fmt, temp, p_ax)

             call zaxpy(maxn, -eig_values(j), p_eig_vectors, 1, p_ax, 1)

             rd(j,1) = real(eig_values(j))
             rd(j,2) = aimag(eig_values(j))
             rd(j,3) = dznrm2(maxn, p_ax, 1)
             rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))

          END DO
          eigen_values = rd(:,1)
          !
          !            %-----------------------------%
          !            | Display computed residuals. |
          !            %-----------------------------%
          !
          call dmout(6, nconv, 3, rd, maxncv, -6, &
                            'Ritz values (Real, Imag) and relative residuals')
          
       END IF
       !             %-------------------------------------------%
       !             | Print additional convergence information. |
       !             %-------------------------------------------%
       If ( iparam(3) .gt. maxitr ) THEN
          WRITE (*,*) ' '
          WRITE (*,*) ' Maximum number of iterations reached.'
          WRITE (*,*) ' '
       ELSE IF ( info .eq. 3) THEN
          WRITE (*,*) ' ' 
          WRITE (*,*) ' No shifts could be applied during implicit Arnoldi &
                                                  update, try increasing NCV.'
          WRITE (*,*) ' '
       END IF
       
       WRITE (*,*) ' '
       WRITE (*,*) 'ARPACK SUMMARY'
       WRITE (*,*) '=============='
       WRITE (*,*) ' '
       WRITE (*,*) ' Size of the matrix is ', maxn
       WRITE (*,*) ' The number of Ritz values requested is ', nev
       WRITE (*,*) ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
       WRITE (*,*) ' What portion of the spectrum: ', which
       WRITE (*,*) ' The number of converged Ritz values is ', nconv 
       WRITE (*,*) ' The number of Implicit Arnoldi update', ' iterations taken is ', iparam(3)
       WRITE (*,*) ' The number of OP*x is ', iparam(9)
       WRITE (*,*) ' The convergence criterion is ', tol
       WRITE (*,*) ' '
       
    END IF
    
    !DEALLOCATE WORK SPACE 
    DEALLOCATE( eig_values, STAT = err )
    DEALLOCATE( eig_vectors, STAT = err )
    DEALLOCATE( temp, STAT = err )
    DEALLOCATE( p_ax, STAT = err )
    DEALLOCATE( workd, STAT = err )
    DEALLOCATE( workev, STAT = err )
    DEALLOCATE( resid, STAT = err )
    DEALLOCATE( workl, STAT = err )
    DEALLOCATE( rwork, STAT = err )
    DEALLOCATE( rd, STAT = err )
    
  END SUBROUTINE arpack_zh

END MODULE arpack
