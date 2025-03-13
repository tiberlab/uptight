  !=========================================================================================
  !Copyright (c) 2009, The Regents of the University of Massachusetts, Amherst.
  !Developed by E. Polizzi
  !All rights reserved.
  !
  !Redistribution and use in source and binary forms, with or without modification, 
  !are permitted provided that the following conditions are met:
  !
  !1. Redistributions of source code must retain the above copyright notice, this list of conditions 
  !   and the following disclaimer.
  !2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
  !   and the following disclaimer in the documentation and/or other materials provided with the distribution.
  !3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
  !    products derived from this software without specific prior written permission.
  !
  !THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
  !BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  !ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
  !EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
  !SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
  !LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  !IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !==========================================================================================
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST REVERSE COMMUNICATION INTERFACES !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! List of routines:
  !-------------------
  ! feastinit
  ! checkfeastparam
  ! scheck_rci_input
  ! dcheck_rci_input
  ! dfeast_rci
  ! zfeast_rci
  ! sfeast_rci
  ! cfeast_rci
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  include 'point_gauss_legendre.f90' 




  subroutine feastinit(feastparam)
    !  Purpose
    !  =======
    !
    !  Define the default values for the input FEAST parameters.
    !
    !  Arguments
    !  =========
    !
    !  feastparam  (output) INTEGER(64): FEAST parameters
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    integer,dimension(64) :: feastparam
    feastparam=0
    feastparam(1)=0 !com
    feastparam(2)=8 !nbe
    feastparam(3)=12!tol double precision
    feastparam(4)=20!maxloop
    feastparam(5)=0 !IS
    feastparam(6)=1 !resid
    feastparam(7)=5 !tol single precision
    feastparam(10)=1 ! omp (number of threads)
    feastparam(11)=2 ! algo (type of algorithm 1, 2 for zfeast_rci)
    feastparam(12)=0 ! customize eigenvalue solver for rci interface (undocumented)
  end subroutine feastinit



  subroutine checkfeastparam(fp,info)
    !  Purpose 
    !  =======
    !  Error handling for input FEAST parameters.
    !  Check the values for the input FEAST parameters, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  fp   (input) INTEGER(64) : FEAST parameters
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    integer,dimension(64) :: fp
    integer :: info
    integer:: i
    logical :: test
    integer,parameter :: max=13
    integer, dimension(max):: tnbe=(/3,4,5,6,8,10,12,16,20,24,32,40,48/)


    if ((fp(1)/=0).and.(fp(1)/=1)) info=101
    test=.true.
    do i=1,max
       if (fp(2)==tnbe(i)) test=.false. 
    enddo
    if (test) info=102
    if (fp(3)<0) info=103
    if (fp(4)<0) info=104
    if ((fp(5)/=0).and.(fp(5)/=1)) info=105
    if ((fp(6)/=0).and.(fp(6)/=1)) info=106
    if ((fp(11)/=1).and.(fp(11)/=2)) info=111
    if ((fp(12)/=0).and.(fp(12)/=1)) info=112

  end subroutine checkfeastparam




  subroutine dcheck_rci_input(Emin,Emax,M0,N,info)
    !  Purpose 
    !  =======
    !  Error handling for the FEAST RCI double precision interfaces input parameters.
    !  Check the values of Emin, Emax, M0,N, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  Emin,Emax   (input) REAL DOUBLE PRECISION: search interval
    !  M0          (input) INTEGER: Size subspace
    !  N           (input) INTEGER: Size system
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    double precision :: Emin,Emax
    integer :: N,M0,info

    if (Emin>=Emax) info=200 ! problem with Emin, Emax
    if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
    if (N<=0) info=202 ! problem with N

  end subroutine dcheck_rci_input





  subroutine scheck_rci_input(Emin,Emax,M0,N,info)
    !  Purpose 
    !  =======
    !  Error handling for the FEAST RCI single precision interfaces input parameters.
    !  Check the values of Emin, Emax, M0,N, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  Emin,Emax   (input) REAL SINGLE PRECISION: search interval
    !  M0          (input) INTEGER: Size subspace
    !  N           (input) INTEGER: Size system
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    real :: Emin,Emax
    integer :: N,M0,info

    if (Emin>=Emax) info=200 ! problem with Emin, Emax
    if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
    if (N<=0) info=202 ! problem with N

  end subroutine scheck_rci_input





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine dfeast_rci(ijob,N,Ze,work,workc,Aq,Sq,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Ax=eBx and standard Ax=eX eigenvalue problems
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE  
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) REAL DOUBLE PRECISION (N,M0):  Workspace 
    !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
    !  Aq,Sq      (input/output) REAL DOUBLE PRECISION (M0,M0) : Worspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  feastparam (input/output) INTEGER(64) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if feastparam(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option feastparam(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    double precision, dimension(N,*) ::work
    complex(kind=(kind(1.0d0))), dimension(N,*):: workc
    double precision, dimension(M0,*) ::Aq,Sq
    integer,dimension(64) :: feastparam
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    double precision,dimension(*)  :: lambda
    double precision,dimension(N,*):: q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    !! parameters
    double precision, Parameter :: pi=3.1415926535897932d0
    double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
    double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
    integer*8,parameter :: fout =6
    !! variable for FEAST
    integer :: i,m_min,m_max,e,Mf
    integer,dimension(4) :: iseed
    double precision :: theta,r,Emid
    complex(kind=(kind(1.0d0))) :: jac
    double precision ::xe,we ! Gauss-Legendre
    double precision, dimension(:,:),pointer :: Sqo
    logical :: testconv
    double precision :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    double precision, dimension(:),pointer :: work_loc
    integer,dimension(:),pointer :: ipiv
    integer :: lwork,lwork_loc,info_lap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ijob==-1) then
       info=0
       call checkfeastparam(feastparam,info)
       call dcheck_rci_input(Emin,Emax,M0,N,info)
       if (info/=0) then
          ijob=0
          return
       endif
       M0=min(M0,N)
       feastparam(23)=M0 
       feastparam(21)=0
       if (feastparam(1)==1) then
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '*********** FEAST- BEGIN **********************', 1)
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2)
          call wwrite(fout, 'Size subspace', 1)  
          call wwrite(fout, '\t', -2) 
          call wwrite(fout,M0,2)
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '#Loop | #Eig  |       Trace           |     Error-trace', 1)  
          call wwrite(fout, '\n', -2)
       endif
    end if

    if (feastparam(21)==0) then
       loop=0
       feastparam(21)=1 ! prepare reentry
       if (feastparam(5)==0) then !!! random vectors
          iseed=(/56,890,3456,2333/)
          call DLARNV(3,iseed,N*M0,work)
       elseif (feastparam(5)==1) then !!!!!! q is the initial guess
          ijob=40
          return
       end if
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)<=3) then

       if (feastparam(21)==-1) then  !!! no contour is used 
          call DLACPY( 'F', N, M0,work , N, q, N )
       else

          if (feastparam(21)==1) then !! we start a new contour integration
             q(1:N,1:M0)=DZERO
             feastparam(20)=1
             feastparam(21)=2
             ijob=20 ! just initialization 
          end if


          do e=feastparam(20),feastparam(2) !!!! loop over the contour

             if ((feastparam(21)==2).and.(ijob==20)) then !!Factorize the linear system (complex) (zS-A)
                call dset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wdcos(theta)+r*(DZERO,DONE)*wdsin(theta)
                ijob=10 ! for fact
                feastparam(20)=e
                return
             endif

             if ((feastparam(21)==2).and.(ijob==10)) then !!Solve the linear system (complex) (zS-A)q=v
                call ZLACP2( 'F', N, M0,work , N, workc, N )
                feastparam(21)=3 ! preparing reentry
                ijob=20 ! for solve
                feastparam(20)=e
                return
             endif

!!!!!! Add contribution to the integral integral
             if (feastparam(21)==3) then              
                call dset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0 
                jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                q(1:N,1:M0)=q(1:N,1:M0)+ba*(-DONE/pi)*we*dble(jac*workc(1:N,1:M0))
                feastparam(21)=2
             endif
          end do
       endif

       feastparam(21)=4 
    end if ! feastparam(21)<=3




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq =eq Sq xq
!!!!!!! with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!for Aq=> Aq=Q^T A Q 
    if (feastparam(21)==4) then
       feastparam(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif
    if (feastparam(21)==5) then ! Aq=Q^T A Q       
       call DGEMM('T','N',M0,M0,N,DONE,q(1,1),N,work(1,1),N,DZERO,Aq,M0) ! create new leading dimension for Aq
       feastparam(21)=6
    endif

!!!!!!!!!for  Sq=> Sq=Q^T S Q
    if (feastparam(21)==6) then
       feastparam(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec B*q => work
    end if
    if (feastparam(21)==7) then ! Bq=Q^T B Q
       call DGEMM('T','N',M0,M0,N,DONE,q(1,1),N,work(1,1),N,DZERO,Sq,M0) ! create new leading dimension for Sq
       feastparam(21)=8
    endif



    if (feastparam(21)==8) then
       mode=M0 ! save value of M0
       if (feastparam(12)==1) then ! customize eigenvalue solver
          feastparam(21)=9 ! preparing reentry - could return new value of M0 if reduced subspace is needed
          ijob=50
          return
       endif
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem using LAPACK!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)==8) then
       JOBZ='V'
       UPLO='L'
       info_lap=1
       i=1
       LWORK_LOC=3*M0-1 !! for lapack eig reduced system
       call wallocate_1d(WORK_LOC,LWORK_LOC,info)
       call wallocate_2d(Sqo,M0,M0,info)
       Mf=M0
       do while (info_lap/=0)
          i=i+1
          if (i>10) then
             if (feastparam(1)==1) then
                call wwrite(fout, 'problem reduced system', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             info=-3
             ijob=0
             return
          end if

          call DLACPY( 'F', Mf, Mf,Sq , M0, Sqo, M0 )
          call dsygv(1,JOBZ,UPLO,Mf,Aq,M0,Sqo,M0,lambda,work_loc,Lwork_loc,info_lap)

          if ((info_lap<=Mf).and.(info_lap/=0)) then
             info=-3
             ijob=0
             return
          end if


          if (info_lap/=0) then
             Mf=info_lap-Mf-1            
             if (feastparam(1)==1) then 
                call wwrite(fout, 'Resize subspace', 1)  
                call wwrite(fout, '\t', -2) 
                call wwrite(fout,Mf,2)
                call wwrite(fout, '\n', -2)
             end if
          end if
       end do

       call wdeallocate_1d(WORK_LOC)
       call wdeallocate_2d(Sqo)
       feastparam(21)=9
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Ritz values/vectors !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    testconv=.true. ! by default


    if (feastparam(21)==9) then
       if (feastparam(12)==1) Mf=M0   ! redefine value
       M0=mode ! reload value 
       testconv=.false.
       mode=0
       do i=1,Mf
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) mode=mode+1
       enddo

       if (mode==0) then ! no mode detected in the interval
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: No eigenvalues have been found in the proposed search interval', 1)
             call wwrite(fout, '\n', -2)
          endif
          ! return here with info error
          info=1
          ijob=0
          return
       endif

       if ((mode==feastparam(23)).and.(mode/=N)) then
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: Size subspace M0 too small', 1)  
             call wwrite(fout, '\n', -2)
          end if
          ! return here with info error
          info=3
          ijob=0
          return
       endif


       m_min=1
       do i=1,Mf
          if (lambda(i)<Emin) m_min=i+1
       enddo
       m_max=m_min+mode-1
       trace=sum(lambda(m_min:m_max))

       if (loop>0) then
          epsout=(abs(trace-epsout)/abs(epsout))
          if (epsout/=DZERO) then
             if (log10(epsout)<(-feastparam(3))) testconv=.true.
          else
             testconv=.true.
          end if
       end if


       if (feastparam(1)==1) then 
          call wwrite(fout,loop,2)
          call wwrite(fout, '\t', -2)
          call wwrite(fout,mode,2)
          call wwrite(fout, '\t', -2)
          call wwrite(fout,trace,4)
          if (loop>0) then
             call wwrite(fout, '\t', -2) 
             call wwrite(fout,epsout,4)
          endif
          call wwrite(fout, '\n', -2) 

          if (testconv) then
             call wwrite(fout, '==>FEAST has converged (reaches desired tolerance)', 1)  
             call wwrite(fout, '\n', -2) 
          end if
       end if


       if (.not.testconv) then
          epsout=trace
          if (loop==feastparam(4)) then
             if (feastparam(1)==1) then 
                call wwrite(fout, '==>FEAST did not converge (#loop reaches maximum)', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             ! return here with info error
             info=2
             ijob=0
             return
          end if
       endif

    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (testconv) then  !!! final eigenvectors/eigenvalues

       if (feastparam(21)==9) then
!!! shift lambda
          if (m_min/=1) then
             do i=1,mode
                lambda(i)=lambda(m_min+i-1)
             enddo
          end if

!!!! what are the vectors
          call DLACPY( 'F', N, Mf,q , N, work, N )
          !! option - shifted
          call DGEMM('N','N',N,Mf-m_min+1,Mf,DONE,work(1,1),N,Aq(1,m_min),M0,DZERO,q(1,1),N) 
          if (m_min>1) call DGEMM('N','N',N,m_min-1,Mf,DONE,work(1,1),N,Aq(1,1),M0,DZERO,q(1,Mf-m_min+2),N)           
          M0=Mf 
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(6)==1) then !!! compute residual

          if (feastparam(21)==9) then
             M0=mode !! temporary value (m_min and Mf are saved)
             mode=Mf !! save value of Mf
             feastparam(21)=10 ! preparing reentry
             ijob=30 
             return  ! mat-vec A*q => work
          endif

          if (feastparam(21)==10) then
             call ZLACP2( 'F', N, M0,work , N, workc, N )
             feastparam(21)=11 ! preparing reentry
             ijob=40 
             return  ! mat-vec S*q => work
          endif

          if (feastparam(21)==11) then
             Mf=mode !! recover value
             mode=M0 !! recover vlaue
             do i=1,mode
                res(i)=sum(abs(dble(workc(1:N,i))-lambda(i)*work(1:N,i)))/sum(abs(dble(workc(1:N,i))))
             end do
             M0=Mf ! new subspace
             ijob=0
          endif

       else
          ijob=0
       end if

!!!!!!!!!!!!!!!!!!!!

       if (ijob==0) then
          if (feastparam(1)==1) then
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '*********** FEAST- END*************************', 1)
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '\n', -2)
          endif
       end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else !!! need refinement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(21)==9) then

          call DLACPY( 'F', N, Mf,q , N, work, N )
          !! option - non shifted
          call DGEMM('N','N',N,Mf,Mf,DONE,work(1,1),N,Aq(1,1),M0,DZERO,q,N) ! eigenvectors with no refinement
          M0=Mf ! new search dimension
!!!!!!!!!!! here q are the eigenvectors, work is the result on the first integration 
          feastparam(21)=1   ! prepare reentry- reloop (with contour integration)
          !feastparam(21)=-1 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=40  ! mat-vec=> S*q => work
          return  
       end if
    end if

  end subroutine dfeast_rci




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine zfeast_rci(ijob,N,Ze,work,workc,zAq,zSq,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Ax=eBx and standard Ax=eX eigenvalue problems
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE  
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
    !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
    !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Worspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  feastparam (input/output) INTEGER(64) : FEAST parameters
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if feastparam(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option feastparam(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(64) :: feastparam
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    double precision,dimension(*)  :: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    !! parameters
    double precision, Parameter :: pi=3.1415926535897932d0
    double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
    double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
    integer*8,parameter :: fout =6
    !! variable for FEAST
    integer :: i,m_min,m_max,e,Mf,nbe2,algo
    integer,dimension(4) :: iseed
    double precision :: theta,r,Emid
    complex(kind=(kind(1.0d0))) :: jac,aux
    double precision ::xe,we ! Gauss-Legendre
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: zSqo
    logical :: testconv
    double precision :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    double precision, dimension(:),pointer :: work_loc
    complex(kind=(kind(1.0d0))), dimension(:),pointer :: zwork_loc
    integer,dimension(:),pointer :: ipiv
    integer :: lwork,lwork_loc,info_lap
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
    algo=feastparam(11)
    if (ijob==-1) then
       info=0
       call checkfeastparam(feastparam,info)
       call dcheck_rci_input(Emin,Emax,M0,N,info)
       if (info/=0) then
          ijob=0
          return
       endif
       M0=min(M0,N)
       feastparam(23)=M0 
       feastparam(21)=0
       if (feastparam(1)==1) then
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '*********** FEAST- BEGIN **********************', 1)
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, 'Size subspace', 1)  
          call wwrite(fout, '\t', -2) 
          call wwrite(fout,M0,2)
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '#Loop | #Eig  |       Trace           |     Error-Trace', 1)  
          call wwrite(fout, '\n', -2)
       endif
    end if


    if (feastparam(21)==0) then
       loop=0
       if (algo==1) then ! half contour
          feastparam(22)=feastparam(2)
       else ! full contour (nbe is only the number of points in the half-contour)
          feastparam(22)=feastparam(2)*2
       end if
       feastparam(21)=1 ! prepare reentry
       if (feastparam(5)==0) then !!! random vectors
          iseed=(/56,890,3456,2333/)
          call ZLARNV(3,iseed,N*M0,work)
       elseif (feastparam(5)==1) then !!!!!! q is the initial guess
          ijob=40
          return
       end if
    endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)<=3) then

       if (feastparam(21)==-1) then  !!! no contour is used 
          call ZLACPY( 'F', N, M0,work , N, q, N )
       else

          if (feastparam(21)==1) then !! we start a new contour integration 
             q(1:N,1:M0)=ZEROC
             feastparam(20)=1
             feastparam(21)=2
             ijob=-2 ! just initialization 
          end if

          do e=feastparam(20),feastparam(22) !!!! loop over the contour 

             if ((feastparam(21)==2).and.(ijob==-2)) then !!Factorize the linear system (complex) (zS-A)

                if (algo==1) then !! half-contour
                   call dset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                else ! full contour
                   if (e<=feastparam(2)) then
                      i=e
                   else
                      i=e-feastparam(2)
                   end if
                   call dset_point_gauss_legendre(feastparam(2),i,xe,we) !! Gauss-points  
                end if
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wdcos(theta)+r*(DZERO,DONE)*wdsin(theta)
                if (algo==2) then
                   if (e>feastparam(2)) Ze=conjg(Ze)
                end if
                ijob=10 ! for fact
                feastparam(20)=e
                return
             endif

             if ((feastparam(21)==2).and.(ijob==10)) then !!Solve the linear system (complex) (zS-A)q=v --> v real part
                call ZLACPY( 'F', N, M0,work , N, workc, N )
                if (algo==1) then ! half-contour 
                   feastparam(21)=2 ! preparing reentry
                else ! full contour
                   feastparam(21)=3 ! preparing reentry
                end if
                ijob=20 ! for solve
                return
             endif


             if (algo==1) then ! half-contour 
                if ((feastparam(21)==2).and.(ijob==20)) then!!!! Solve the linear system (complex) (zS-A)^Tq=v
                   call dset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                   theta=ba*xe+ab
                   r=(Emax-Emin)/2.0d0 
                   jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                   aux=-ba*(ONEC/(2.0d0*pi))*we*jac
                   call ZAXPY(N*M0,aux,workc,1,q,1)
                   call ZLACPY( 'F', N, M0,work , N, workc, N )
                   ijob=21 ! for solve with transpose
                   feastparam(21)=3 ! preparing reentry
                   return
                end if
                if (feastparam(21)==3) then
                   call dset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                   theta=ba*xe+ab
                   r=(Emax-Emin)/2.0d0 
                   jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                   aux=-ba*(ONEC/(2.0d0*pi))*we*conjg(jac)
                   call ZAXPY(N*M0,aux,workc,1,q,1)
                   feastparam(21)=2
                   ijob=-2 ! just for identification
                end if

             else !! full contour ! algo=2
                if (feastparam(21)==3) then
                   if (e<=feastparam(2)) then
                      i=e
                   else
                      i=e-feastparam(2)
                   end if
                   call dset_point_gauss_legendre(feastparam(2),i,xe,we) !! Gauss-points  
                   theta=ba*xe+ab
                   r=(Emax-Emin)/2.0d0 
                   jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                   if (e>feastparam(2)) jac=conjg(jac)
                   aux=-ba*(ONEC/(2.0d0*pi))*we*jac
                   call ZAXPY(N*M0,aux,workc,1,q,1)
                   feastparam(21)=2
                   ijob=-2 ! just for identification
                end if
             end if
          end do
       endif

       feastparam(21)=4 
    end if ! feastparam(21)<=3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq =eq Sq xq
!!!!!!! with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!for Aq=> Aq=Q^T A Q 
    if (feastparam(21)==4) then
       feastparam(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif
    if (feastparam(21)==5) then 
       call ZGEMM('C','N',M0,M0,N,ONEC,q(1,1),N,work(1,1),N,ZEROC,zAq,M0) ! new leading diemsnion for zAq
       feastparam(21)=6
    endif

!!!!!!!!!for  Sq=> Sq=Q^T S Q
    if (feastparam(21)==6) then 
       feastparam(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec B*q => work
    end if
    if (feastparam(21)==7) then 
       call ZGEMM('C','N',M0,M0,N,ONEC,q(1,1),N,work(1,1),N,ZEROC,zSq,M0) ! new leading dimension for zSq
       feastparam(21)=8
    endif


    if (feastparam(21)==8) then
       mode=M0 ! save value of M0
       if (feastparam(12)==1) then ! customize eigenvalue solver
          feastparam(21)=9 ! preparing reentry - could return new value of M0 if reduced subspace is needed
          ijob=50
          return
       endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)==8) then
       JOBZ='V'
       UPLO='L'
       info_lap=1
       i=1
       LWORK_LOC=2*M0-1 !! for lapack eig reduced system
       call wallocate_1z(zWORK_LOC,LWORK_LOC,info)
       call wallocate_1d(WORK_LOC,3*M0-2,info)
       call wallocate_2z(zSqo,M0,M0,info)
       Mf=M0
       do while (info_lap/=0)
          i=i+1
          if (i>10) then
             if (feastparam(1)==1) then 
                call wwrite(fout, 'problem reduced system', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             info=-3
             ijob=0
             return
          end if

          call ZLACPY( 'F', Mf, Mf,zSq , M0, zSqo, M0 )
          call ZHEGV(1, JOBZ, UPLO, Mf, zAq, M0, zSqo, M0, lambda, zWORK_loc,Lwork_loc, WORK_loc, INFO_lap)

          if ((info_lap<=Mf).and.(info_lap/=0)) then
             info=-3
             ijob=0
             return
          end if

          if (info_lap/=0) then
             Mf=info_lap-Mf-1
             if (feastparam(1)==1) then
                call wwrite(fout, 'Resize subspace', 1)  
                call wwrite(fout, '\t', -2) 
                call wwrite(fout,Mf,2)
                call wwrite(fout, '\n', -2)
             end if
          end if
       end do
       M0=Mf
       call wdeallocate_2z(zSqo)
       call wdeallocate_1z(zwork_loc)
       call wdeallocate_1d(work_loc)
       feastparam(21)=9
    end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Ritz values/vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    testconv=.true. ! by default

    if (feastparam(21)==9) then
       if (feastparam(12)==1) Mf=M0 ! recover value
       M0=mode  ! recover value
       testconv=.false.
       mode=0
       do i=1,Mf
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) mode=mode+1
       enddo

       if (mode==0) then ! no mode detected in the interval
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: No eigenvalues have been found in the proposed search interval', 1)
             call wwrite(fout, '\n', -2)
          endif
          ! return here with info error
          info=1
          ijob=0
          return
       end if

       if ((mode==feastparam(23)).and.(mode/=N)) then
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: Size subspace M0 too small', 1)  
             call wwrite(fout, '\n', -2)
          end if
          ! return here with info error
          info=3
          ijob=0
          return
       endif


       m_min=1
       do i=1,Mf
          if (lambda(i)<Emin) m_min=i+1
       enddo
       m_max=m_min+mode-1
       trace=sum(lambda(m_min:m_max))

       if (loop>0) then
          epsout=(abs(trace-epsout)/abs(epsout))
          if (epsout/=DZERO) then
             if (log10(epsout)<(-feastparam(3))) testconv=.true.
          else
             testconv=.true.
          end if
       end if


       if (feastparam(1)==1) then
          call wwrite(fout,loop,2)
          call wwrite(fout, '\t', -2) 
          call wwrite(fout,mode,2)
          call wwrite(fout, '\t', -2)
          call wwrite(fout,trace,4)
          if (loop>0) then
             call wwrite(fout, '\t', -2) 
             call wwrite(fout,epsout,4)
          endif
          call wwrite(fout, '\n', -2) 

          if (testconv) then
             call wwrite(fout, '==>FEAST has converged (reaches desired tolerance)', 1)  
             call wwrite(fout, '\n', -2) 
          end if
       end if


       if (.not.testconv) then
          epsout=trace
          if (loop==feastparam(4)) then
             if (feastparam(1)==1) then
                call wwrite(fout, '==>FEAST did not converge (#loop reaches maximum)', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             ! return here with info error
             info=2
             ijob=0
             return
          end if
       endif

    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (testconv) then  !!! final eigenvectors/eigenvalues


       if (feastparam(21)==9) then
!!! shift lambda
          if (m_min/=1) then
             do i=1,mode
                lambda(i)=lambda(m_min+i-1)
             enddo
          end if
!!!! what are the vectors
          call ZLACPY( 'F', N, Mf,q , N, work, N )
          !! option - shifted
          call ZGEMM('N','N',N,Mf-m_min+1,Mf,ONEC,work(1,1),N,zAq(1,m_min),M0,ZEROC,q(1,1),N) 
          if (m_min>1) call ZGEMM('N','N',N,m_min-1,Mf,ONEC,work(1,1),N,zAq(1,1),M0,ZEROC,q(1,Mf-m_min+2),N)           
          M0=Mf 
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(6)==1) then !!! compute residual

          if (feastparam(21)==9) then
             M0=mode !! temporary value (m_min and Mf are saved)
             mode=Mf !! save value of Mf
             feastparam(21)=10 ! preparing reentry
             ijob=30 
             return  ! mat-vec A*q => work
          endif

          if (feastparam(21)==10) then
             call ZLACPY( 'F', N, M0,work , N, workc, N )
             feastparam(21)=11 ! preparing reentry
             ijob=40 
             return  ! mat-vec S*q => work
          endif

          if (feastparam(21)==11) then
             Mf=mode !! recover value
             mode=M0 !! recover vlaue
             do i=1,mode
                res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs((workc(1:N,i))))
             end do
             M0=Mf ! new subspace
             ijob=0
          endif

       else
          ijob=0
       end if

!!!!!!!!!!!!!!!

       if (ijob==0) then
          if (feastparam(1)==1) then
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '*********** FEAST- END*************************', 1)
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '\n', -2)
          endif
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else !!! need refinement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       if (feastparam(21)==9) then
          work(1:N,1:Mf)=q(1:N,1:Mf)
          !! option - non shifted
          call ZGEMM('N','N',N,Mf,Mf,ONEC,work(1,1),N,zAq(1,1),M0,ZEROC,q,N) ! eigenvectors with no refinement
          M0=Mf ! new search dimension
!!!!!!!!!!! here q are the eigenvectors, work is the result on the first integration 
          feastparam(21)=1   ! prepare reentry- reloop (with contour integration)
          !feastparam(21)=-1 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=40  ! mat-vec=> S*q => work
          return  
       end if

    end if

  end subroutine zfeast_rci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine sfeast_rci(ijob,N,Ze,work,workc,Aq,Sq,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Ax=eBx and standard Ax=eX eigenvalue problems
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE  
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) REAL SINGLE PRECISION (N,M0):  Workspace 
    !  workc      (input/output) COMPLEX SINGLE PRECISION (N,M0):  Workspace 
    !  Aq,Sq      (input/output) REAL SINGLE PRECISION (M0,M0) : Worspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  feastparam (input/output) INTEGER(64) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if feastparam(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option feastparam(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex  :: Ze
    real, dimension(N,*) ::work
    complex , dimension(N,*):: workc
    real, dimension(M0,*) ::Aq,Sq
    integer,dimension(64) :: feastparam
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    real,dimension(*)  :: lambda
    real,dimension(N,*):: q
    integer :: mode
    real,dimension(*) :: res
    integer :: info
    !! parameters
    real, Parameter :: pi=3.1415926535897932e0
    real, Parameter :: SONE=1.0E0, SZERO=0.0E0
    complex ,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
    real, parameter :: ba=-pi/2.0e0, ab=pi/2.0e0
    integer*8,parameter :: fout =6
    !! variable for FEAST
    integer :: i,m_min,m_max,e,Mf
    integer,dimension(4) :: iseed
    real :: theta,r,Emid
    complex  :: jac
    real ::xe,we ! Gauss-Legendre
    real, dimension(:,:),pointer :: Sqo
    logical :: testconv
    real :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    real, dimension(:),pointer :: work_loc
    integer,dimension(:),pointer :: ipiv
    integer :: lwork,lwork_loc,info_lap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ijob==-1) then
       info=0
       call checkfeastparam(feastparam,info)
       call scheck_rci_input(Emin,Emax,M0,N,info)
       if (info/=0) then
          ijob=0
          return
       endif
       M0=min(M0,N)
       feastparam(23)=M0 
       feastparam(21)=0
       if (feastparam(1)==1) then
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '*********** FEAST- BEGIN **********************', 1)
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2)
          call wwrite(fout, 'Size subspace', 1)  
          call wwrite(fout, '\t', -2) 
          call wwrite(fout,M0,2)
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '#Loop | #Eig |  Trace    |   Error-trace', 1)  
          call wwrite(fout, '\n', -2)
       endif
    end if

    if (feastparam(21)==0) then
       loop=0
       feastparam(21)=1 ! prepare reentry
       if (feastparam(5)==0) then !!! random vectors
          iseed=(/56,890,3456,2333/)
          call SLARNV(3,iseed,N*M0,work)
       elseif (feastparam(5)==1) then !!!!!! q is the initial guess
          ijob=40
          return
       end if
    endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)<=3) then

       if (feastparam(21)==-1) then  !!! no contour is used  
          call SLACPY( 'F', N, M0,work , N, q, N )
       else
          if (feastparam(21)==1) then !! we start a new contour integration 
             q(1:N,1:M0)=SZERO
             feastparam(20)=1
             feastparam(21)=2
             ijob=20 ! just initialization 
          end if


          do e=feastparam(20),feastparam(2) !!!! loop over the contour

             if ((feastparam(21)==2).and.(ijob==20)) then !!Factorize the linear system (complex) (zS-A)
                call sset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0E0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wscos(theta)+r*(SZERO,SONE)*wssin(theta)
                ijob=10 ! for fact
                feastparam(20)=e
                return
             endif

             if ((feastparam(21)==2).and.(ijob==10)) then !!Solve the linear system (complex) (zS-A)q=v
                call CLACP2( 'F', N, M0,work , N, workc, N )
                feastparam(21)=3 ! preparing reentry
                ijob=20 ! for solve
                feastparam(20)=e
                return
             endif

!!!!!! Add contribution to the integral integral
             if (feastparam(21)==3) then              
                call sset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0E0 
                jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                q(1:N,1:M0)=q(1:N,1:M0)+ba*(-SONE/pi)*we*real(jac*workc(1:N,1:M0))
                feastparam(21)=2
             endif
          end do
       endif

       feastparam(21)=4 
    end if ! feastparam(21)<=3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq =eq Sq xq
!!!!!!! with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!for Aq=> Aq=Q^T A Q 
    if (feastparam(21)==4) then
       feastparam(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif
    if (feastparam(21)==5) then ! Aq=Q^T A Q       
       call SGEMM('T','N',M0,M0,N,SONE,q(1,1),N,work(1,1),N,SZERO,Aq,M0) ! create new leading dimension for Aq
       feastparam(21)=6
    endif

!!!!!!!!!for  Sq=> Sq=Q^T S Q
    if (feastparam(21)==6) then
       feastparam(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec B*q => work
    end if
    if (feastparam(21)==7) then ! Bq=Q^T B Q
       call SGEMM('T','N',M0,M0,N,SONE,q(1,1),N,work(1,1),N,SZERO,Sq,M0) ! create new leading dimension for Sq
       feastparam(21)=8
    endif


    if (feastparam(21)==8) then
       mode=M0 ! save value of M0
       if (feastparam(12)==1) then ! customize eigenvalue solver
          feastparam(21)=9 ! preparing reentry - could return new value of M0 if reduced subspace is needed
          ijob=50
          return
       endif
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem using LAPACK!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)==8) then
       JOBZ='V'
       UPLO='L'
       info_lap=1
       i=1
       LWORK_LOC=3*M0-1 !! for lapack eig reduced system
       call wallocate_1s(WORK_LOC,LWORK_LOC,info)
       call wallocate_2s(Sqo,M0,M0,info)
       Mf=M0
       do while (info_lap/=0)
          i=i+1          
          if (i>10) then
                if (feastparam(1)==1) then 
                call wwrite(fout, 'problem reduced system', 1)  
                call wwrite(fout, '\n', -2) 
             end if
              info=-3
              ijob=0
              return
          end if

          call SLACPY( 'F', Mf, Mf,Sq , M0, Sqo, M0 )
          call ssygv(1,JOBZ,UPLO,Mf,Aq,M0,Sqo,M0,lambda,work_loc,Lwork_loc,info_lap)

          if ((info_lap<=Mf).and.(info_lap/=0)) then
             info=-3
             ijob=0
             return
          end if

          if (info_lap/=0) then
             Mf=info_lap-Mf-1            
             if (feastparam(1)==1) then 
                call wwrite(fout, 'Resize subspace', 1)  
                call wwrite(fout, '\t', -2) 
                call wwrite(fout,Mf,2)
                call wwrite(fout, '\n', -2)
             end if
          end if
       end do

       call wdeallocate_1s(WORK_LOC)
       call wdeallocate_2s(Sqo)
       feastparam(21)=9
    end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Ritz values/vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    testconv=.true. ! by default
    if (feastparam(21)==9) then
       if (feastparam(12)==1) Mf=M0   ! redefine value
       M0=mode ! reload value (important in case ijob=50)
       testconv=.false.
       mode=0
       do i=1,Mf
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) mode=mode+1
       enddo

       if (mode==0) then ! no mode detected in the interval
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: No eigenvalues have been found in the proposed search interval', 1)
             call wwrite(fout, '\n', -2)
          end if
          ! return here with info error
          info=1
          ijob=0
          return
       endif

       if ((mode==feastparam(23)).and.(mode/=N)) then
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: Size subspace M0 too small', 1)  
             call wwrite(fout, '\n', -2)
          end if
          ! return here with info error
          info=3
          ijob=0
          return
       endif


       m_min=1
       do i=1,Mf
          if (lambda(i)<Emin) m_min=i+1
       enddo
       m_max=m_min+mode-1
       trace=sum(lambda(m_min:m_max))

       if (loop>0) then
          epsout=(abs(trace-epsout)/abs(epsout))
          if (epsout/=SZERO) then
             if (log10(epsout)<(-feastparam(7))) testconv=.true.
          else
             testconv=.true.
          end if
       end if


       if (feastparam(1)==1) then 
          call wwrite(fout,loop,2)
          call wwrite(fout, '\t', -2)
          call wwrite(fout,mode,2)
          call wwrite(fout, '\t', -2)
          call wwrite(fout,trace,3)
          if (loop>0) then
             call wwrite(fout, '\t', -2) 
             call wwrite(fout,epsout,3)
          endif
          call wwrite(fout, '\n', -2) 

          if (testconv) then
             call wwrite(fout, '==>FEAST has converged (reaches desired tolerance)', 1)  
             call wwrite(fout, '\n', -2) 
          end if
       end if


       if (.not.testconv) then
          epsout=trace
          if (loop==feastparam(4)) then
             if (feastparam(1)==1) then 
                call wwrite(fout, '==>FEAST did not converge (#loop reaches maximum)', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             ! return here with info error
             info=2
             ijob=0
             return
          end if
       endif

    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (testconv) then  !!! final eigenvectors/eigenvalues

       if (feastparam(21)==9) then
!!! shift lambda
          if (m_min/=1) then
             do i=1,mode
                lambda(i)=lambda(m_min+i-1)
             enddo
          end if
!!!! what are the vectors
          call SLACPY( 'F', N, Mf,q , N, work, N )
          !! option - shifted
          call SGEMM('N','N',N,Mf-m_min+1,Mf,SONE,work(1,1),N,Aq(1,m_min),M0,SZERO,q(1,1),N) 
          if (m_min>1) call SGEMM('N','N',N,m_min-1,Mf,SONE,work(1,1),N,Aq(1,1),M0,SZERO,q(1,Mf-m_min+2),N)           
          M0=Mf 
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(6)==1) then !!! compute residual
          if (feastparam(21)==9) then
             M0=mode !! temporary value (m_min and Mf are saved)
             mode=Mf !! save value of Mf
             feastparam(21)=10 ! preparing reentry
             ijob=30 
             return  ! mat-vec A*q => work
          endif

          if (feastparam(21)==10) then
             call CLACP2( 'F', N, M0,work , N, workc, N )
             feastparam(21)=11 ! preparing reentry
             ijob=40 
             return  ! mat-vec S*q => work
          endif

          if (feastparam(21)==11) then
             Mf=mode !! recover value
             mode=M0 !! recover vlaue
             do i=1,mode
                res(i)=sum(abs(real(workc(1:N,i))-lambda(i)*work(1:N,i)))/sum(abs(real(workc(1:N,i))))
             end do
             M0=Mf ! new subspace
             ijob=0
          endif

       else
          ijob=0
       end if

!!!!!!!!!!!!!!!!!!!

       if (ijob==0) then
          if (feastparam(1)==1) then
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '*********** FEAST- END*************************', 1)
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '\n', -2)
          endif
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else !!! need refinement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(21)==9) then
          call SLACPY( 'F', N, Mf,q , N, work, N )
          !! option - non shifted
          call SGEMM('N','N',N,Mf,Mf,SONE,work(1,1),N,Aq(1,1),M0,SZERO,q,N) ! eigenvectors with no refinement
          M0=Mf ! new search dimension
!!!!!!!!!!! here q are the eigenvectors, work is the result on the first integration 
          feastparam(21)=1   ! prepare reentry- reloop (with contour integration)
          !feastparam(21)=-1 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=40  ! mat-vec=> S*q => work
          return  
       end if
    end if

  end subroutine sfeast_rci




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







  subroutine cfeast_rci(ijob,N,Ze,work,workc,zAq,zSq,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Ax=eBx and standard Ax=eX eigenvalue problems
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE  
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) COMPLEX SINGLE PRECISION (N,M0):  Workspace 
    !  workc      (input/output) COMPLEX SINGLE PRECISION (N,M0):  Workspace 
    !  zAq,zSq    (input/output) COMPLEX SINGLE PRECISION (M0,M0) : Worspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  feastparam (input/output) INTEGER(64) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if feastparam(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option feastparam(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009
    ! ====================================================================
    implicit none
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex  :: Ze
    complex ,dimension(N,*):: work,workc
    complex ,dimension(M0,*):: zAq,zSq
    integer,dimension(64) :: feastparam
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    real,dimension(*)  :: lambda
    complex ,dimension(N,*):: q
    integer :: mode
    real,dimension(*) :: res
    integer :: info
    !! parameters
    real, Parameter :: pi=3.1415926535897932e0
    real, Parameter :: SONE=1.0E0, SZERO=0.0E0
    complex ,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
    real, parameter :: ba=-pi/2.0e0, ab=pi/2.0e0
    integer*8,parameter :: fout =6
    !! variable for FEAST
    integer :: i,m_min,m_max,e,Mf,nbe2,algo
    integer,dimension(4) :: iseed
    real :: theta,r,Emid
    complex  :: jac,aux
    real ::xe,we ! Gauss-Legendre
    complex , dimension(:,:),pointer :: zSqo
    logical :: testconv
    real :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    real, dimension(:),pointer :: work_loc
    complex , dimension(:),pointer :: zwork_loc
    integer,dimension(:),pointer :: ipiv
    integer :: lwork,lwork_loc,info_lap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    algo=feastparam(11)
    if (ijob==-1) then
       info=0
       call checkfeastparam(feastparam,info)
       call scheck_rci_input(Emin,Emax,M0,N,info)
       if (info/=0) then
          ijob=0
          return
       endif
       M0=min(M0,N)
       feastparam(23)=M0 
       feastparam(21)=0
       if (feastparam(1)==1) then
          call wwrite(fout, '\n', -2)
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '*********** FEAST- BEGIN **********************', 1)
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, '***********************************************', 1)  
          call wwrite(fout, '\n', -2) 
          call wwrite(fout, 'Size subspace', 1)  
          call wwrite(fout, '\t', -2) 
          call wwrite(fout,M0,2)
          call wwrite(fout, '\n', -2)  
          call wwrite(fout, '#Loop | #Eig |  Trace    |   Error-trace', 1)  
          call wwrite(fout, '\n', -2)
       endif
    end if


    if (feastparam(21)==0) then
       loop=0
       if (algo==1) then ! half contour
          feastparam(22)=feastparam(2)
       else ! full contour (nbe is only the number of points in the half-contour)
          feastparam(22)=feastparam(2)*2
       end if
       feastparam(21)=1 ! prepare reentry
       if (feastparam(5)==0) then !!! random vectors
          iseed=(/56,890,3456,2333/)
          call CLARNV(3,iseed,N*M0,work)
       elseif (feastparam(5)==1) then !!!!!! q is the initial guess
          ijob=40
          return
       end if
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)<=3) then

       if (feastparam(21)==-1) then  !!! no contour is used
          call CLACPY( 'F', N, M0,work , N, q, N )
       else

          if (feastparam(21)==1) then !! we start a new contour integration 
             q(1:N,1:M0)=ZEROC
             feastparam(20)=1
             feastparam(21)=2
             ijob=-2 ! just initialization 
          end if

          do e=feastparam(20),feastparam(22) !!!! loop over the contour 

             if ((feastparam(21)==2).and.(ijob==-2)) then !!Factorize the linear system (complex) (zS-A)

                if (algo==1) then !! half-contour
                   call sset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                else ! full contour
                   if (e<=feastparam(2)) then
                      i=e
                   else
                      i=e-feastparam(2)
                   end if
                   call sset_point_gauss_legendre(feastparam(2),i,xe,we) !! Gauss-points  
                end if
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0E0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wscos(theta)+r*(SZERO,SONE)*wssin(theta)
                if (algo==2) then
                   if (e>feastparam(2)) Ze=conjg(Ze)
                end if
                ijob=10 ! for fact
                feastparam(20)=e
                return
             endif

             if ((feastparam(21)==2).and.(ijob==10)) then !!Solve the linear system (complex) (zS-A)q=v --> v real part
                call CLACPY( 'F', N, M0,work , N, workc, N )
                if (algo==1) then ! half-contour 
                   feastparam(21)=2 ! preparing reentry
                else ! full contour
                   feastparam(21)=3 ! preparing reentry
                end if
                ijob=20 ! for solve
                return
             endif


             if (algo==1) then ! half-contour 
                if ((feastparam(21)==2).and.(ijob==20)) then!!!! Solve the linear system (complex) (zS-A)^Tq=v
                   call sset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                   theta=ba*xe+ab
                   r=(Emax-Emin)/2.0E0 
                   jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                   aux=-ba*(ONEC/(2.0E0*pi))*we*jac
                   call CAXPY(N*M0,aux,workc,1,q,1)
                   call CLACPY( 'F', N, M0,work , N, workc, N )
                   ijob=21 ! for solve with transpose
                   feastparam(21)=3 ! preparing reentry
                   return
                end if
                if (feastparam(21)==3) then
                   call sset_point_gauss_legendre(feastparam(2),e,xe,we) !! Gauss-points 
                   theta=ba*xe+ab
                   r=(Emax-Emin)/2.0E0 
                   jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                   aux=-ba*(ONEC/(2.0E0*pi))*we*conjg(jac)
                   call CAXPY(N*M0,aux,workc,1,q,1)
                   feastparam(21)=2
                   ijob=-2 ! just for identification
                end if

             else !! full contour ! algo=2
                if (feastparam(21)==3) then
                   if (e<=feastparam(2)) then
                      i=e
                   else
                      i=e-feastparam(2)
                   end if
                   call sset_point_gauss_legendre(feastparam(2),i,xe,we) !! Gauss-points  
                   theta=ba*xe+ab
                   r=(Emax-Emin)/2.0E0 
                   jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                   if (e>feastparam(2)) jac=conjg(jac)
                   aux=-ba*(ONEC/(2.0E0*pi))*we*jac
                   call CAXPY(N*M0,aux,workc,1,q,1)
                   feastparam(21)=2
                   ijob=-2 ! just for identification
                end if
             end if
          end do
       endif

       feastparam(21)=4 
    end if ! feastparam(21)<=3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq =eq Sq xq
!!!!!!! with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!for Aq=> Aq=Q^T A Q 
    if (feastparam(21)==4) then
       feastparam(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif
    if (feastparam(21)==5) then
       call CGEMM('C','N',M0,M0,N,ONEC,q(1,1),N,work(1,1),N,ZEROC,zAq,M0) ! new leading diemsnion for zAq
       feastparam(21)=6
    endif

!!!!!!!!!for  Sq=> Sq=Q^T S Q
    if (feastparam(21)==6) then 
       feastparam(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec B*q => work
    end if
    if (feastparam(21)==7) then 
       call CGEMM('C','N',M0,M0,N,ONEC,q(1,1),N,work(1,1),N,ZEROC,zSq,M0) ! new leading dimension for zSq
       feastparam(21)=8
    endif


    if (feastparam(21)==8) then
       mode=M0 ! save value of M0
       if (feastparam(12)==1) then ! customize eigenvalue solver
          feastparam(21)=9 ! preparing reentry - could return new value of M0 if reduced subspace is needed
          ijob=50
          return
       endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (feastparam(21)==8) then
       JOBZ='V'
       UPLO='L'
       info_lap=1
       i=1
       LWORK_LOC=2*M0-1 !! for lapack eig reduced system
       call wallocate_1c(zWORK_LOC,LWORK_LOC,info)
       call wallocate_1s(WORK_LOC,3*M0-2,info)
       call wallocate_2c(zSqo,M0,M0,info)
       Mf=M0
       do while (info_lap/=0)
          i=i+1
          if (i>10) then
             if (feastparam(1)==1) then 
                call wwrite(fout, 'problem reduced system', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             info=-3
             ijob=0
             return
          end if

          call CLACPY( 'F', Mf, Mf,zSq , M0, zSqo, M0 )
          call CHEGV(1, JOBZ, UPLO, Mf, zAq, M0, zSqo, M0, lambda, zWORK_loc,Lwork_loc, WORK_loc, INFO_lap)

          if ((info_lap<=Mf).and.(info_lap/=0)) then
             info=-3
             ijob=0
             return
          end if

          if (info_lap/=0) then
             Mf=info_lap-Mf-1
             if (feastparam(1)==1) then
                call wwrite(fout, 'Resize subspace', 1)  
                call wwrite(fout, '\t', -2) 
                call wwrite(fout,Mf,2)
                call wwrite(fout, '\n', -2)
             end if

          end if
       end do
       M0=Mf
       call wdeallocate_2c(zSqo)
       call wdeallocate_1c(zwork_loc)
       call wdeallocate_1s(work_loc)
       feastparam(21)=9
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Ritz values/vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    testconv=.true. ! by default

    if (feastparam(21)==9) then
       if (feastparam(12)==1) Mf=M0 ! recover value
       M0=mode  ! recover value
       testconv=.false.
       mode=0
       do i=1,Mf
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) mode=mode+1
       enddo

       if (mode==0) then ! no mode detected in the interval
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: No eigenvalues have been found in the proposed search interval', 1)
             call wwrite(fout, '\n', -2)
          endif
          ! return here with info error
          info=1
          ijob=0
          return
       end if

       if ((mode==feastparam(23)).and.(mode/=N)) then
          if (feastparam(1)==1) then
             call wwrite(fout, '==>WARNING: Size subspace M0 too small', 1)  
             call wwrite(fout, '\n', -2)
          end if
          ! return here with info error
          info=3
          ijob=0
          return
       endif


       m_min=1
       do i=1,Mf
          if (lambda(i)<Emin) m_min=i+1
       enddo
       m_max=m_min+mode-1
       trace=sum(lambda(m_min:m_max))

       if (loop>0) then
          epsout=(abs(trace-epsout)/abs(epsout))
          if (epsout/=SZERO) then
             if (log10(epsout)<(-feastparam(7))) testconv=.true.
          else
             testconv=.true.
          end if
       end if


       if (feastparam(1)==1) then
          call wwrite(fout,loop,2)
          call wwrite(fout, '\t', -2) 
          call wwrite(fout,mode,2)
          call wwrite(fout, '\t', -2)
          call wwrite(fout,trace,3)
          if (loop>0) then
             call wwrite(fout, '\t', -2) 
             call wwrite(fout,epsout,3)
          endif
          call wwrite(fout, '\n', -2) 

          if (testconv) then
             call wwrite(fout, '==>FEAST has converged (reaches desired tolerance)', 1)  
             call wwrite(fout, '\n', -2) 
          end if
       end if

       if (.not.testconv) then
          epsout=trace
          if (loop==feastparam(4)) then
             if (feastparam(1)==1) then
                call wwrite(fout, '==>FEAST did not converge (#loop reaches maximum)', 1)  
                call wwrite(fout, '\n', -2) 
             end if
             ! return here with info error
             info=2
             ijob=0
             return
          end if
       endif
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (testconv) then  !!! final eigenvectors/eigenvalues


       if (feastparam(21)==9) then
!!! shift lambda
          if (m_min/=1) then
             do i=1,mode
                lambda(i)=lambda(m_min+i-1)
             enddo
          end if
!!!! what are the vectors
          call CLACPY( 'F', N, Mf,q , N, work, N )
          !! option - shifted
          call CGEMM('N','N',N,Mf-m_min+1,Mf,ONEC,work(1,1),N,zAq(1,m_min),M0,ZEROC,q(1,1),N) 
          if (m_min>1) call CGEMM('N','N',N,m_min-1,Mf,ONEC,work(1,1),N,zAq(1,1),M0,ZEROC,q(1,Mf-m_min+2),N)           
          M0=Mf 
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(6)==1) then !!! compute residual

          if (feastparam(21)==9) then
             M0=mode !! temporary value (m_min and Mf are saved)
             mode=Mf !! save value of Mf
             feastparam(21)=10 ! preparing reentry
             ijob=30 
             return  ! mat-vec A*q => work
          endif

          if (feastparam(21)==10) then
             call CLACPY( 'F', N, M0,work , N, workc, N )
             feastparam(21)=11 ! preparing reentry
             ijob=40 
             return  ! mat-vec S*q => work
          endif

          if (feastparam(21)==11) then
             Mf=mode !! recover value
             mode=M0 !! recover vlaue
             do i=1,mode
                res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs((workc(1:N,i))))
             end do
             M0=Mf! new subspace
             ijob=0
          endif

       else
          ijob=0
       end if
!!!!!!!!!!!!!!!!!


       if (ijob==0) then
          if (feastparam(1)==1) then
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '*********** FEAST- END*************************', 1)
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '***********************************************', 1)  
             call wwrite(fout, '\n', -2) 
             call wwrite(fout, '\n', -2)
          endif
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else !!! need refinement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (feastparam(21)==9) then
          work(1:N,1:Mf)=q(1:N,1:Mf)
          !! option - non shifted
          call CGEMM('N','N',N,Mf,Mf,ONEC,work(1,1),N,zAq(1,1),M0,ZEROC,q,N) ! eigenvectors with no refinement
          M0=Mf ! new search dimension
!!!!!!!!!!! here q are the eigenvectors, work is the result on the first integration 
          feastparam(21)=1   ! prepare reentry- reloop (with contour integration)
          !feastparam(21)=-1 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=40  ! mat-vec=> S*q => work
          return  
       end if
    end if

  end subroutine cfeast_rci
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
