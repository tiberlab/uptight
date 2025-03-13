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
!!!!!!!!!!!!!!!!! FEAST PREDEFINED SPARSE INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! List of routines:
  !-------------------
  ! dfeast_sge
  ! dfeast_sst
  ! zfeast_sge
  ! zfeast_sst
  ! sfeast_sge
  ! sfeast_sst
  ! cfeast_sge
  ! cfeast_sst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  include 'lsprim.f90' !! Sparse primitives





  subroutine dfeast_sge(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE CSR FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(64) :: feastparam
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    double precision,dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::workc,caux
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k
    integer :: opt,nnza,nnzb,nnz


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SGE', -INFO+100 )
       RETURN
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    infoloc=0

    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if


       opt=1
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 

       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1d(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb =>  sb(1:nnzb)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnzb)




    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1d(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call dcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)


    endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
    opt=1
    call daddcsr(N,opt,DONE,ssa,sisa,sjsa,DONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1z(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=2
    call daddcsr(N,opt,DONE,ssa,sisa,sjsa,DONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    symbfact=.true.


    ijob=-1 ! initialization
    do while (ijob/=0)


       call dfeast_rci(ijob,N,Ze,work,workc,Aq,Sq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)

       select case(ijob)

       case(10) !! Factorize (ZeB-A)

          opt=3
          call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!!!!!!!!!!!!!!!


          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=6 ! complex and symmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc


          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,1:M0) result in work(1:N,1:M0)

          call dcsrmm('U',N,N,M0,DONE,ssa,sisa,sjsa,X,DZERO,work)

       case(40) !! perform multiplication B*x(1:N,1:M0) result in work(1:N,1:M0)

          call dcsrmm('U',N,N,M0,DONE,ssb,sisb,sjsb,X,DZERO,work)


       end select
    end do

    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2d(work)
    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1d(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif



  end subroutine dfeast_sge



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dfeast_sst(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: SPARSE CSR FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    double precision,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(64) :: feastparam
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    double precision,dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::workc,caux
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k
    integer :: opt,nnza,nnzb,nnz



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SST', -INFO+100 )
       RETURN
    END IF

    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1d(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=DONE
    enddo
    sisb(n+1)=nnzb+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=1
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 

       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if


       opt=2
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)


    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr



       nnza=isa(n+1)-1
       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call dcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)


    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=1
    call daddcsr(N,opt,DONE,ssa,sisa,sjsa,DONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1z(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=2
    call daddcsr(N,opt,DONE,ssa,sisa,sjsa,DONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    symbfact=.true.


    ijob=-1 ! initialization


    do while (ijob/=0)


       call dfeast_rci(ijob,N,Ze,work,workc,Aq,Sq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! Factorize (ZeB-A)


          opt=3
          call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!!!!!!!!!!!!!!!



          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=6 ! complex and symmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,1:M0) result in work(1:N,1:M0)

          call dcsrmm('U',N,N,M0,DONE,ssa,sisa,sjsa,X,DZERO,work)

       case(40) !! perform multiplication B*x(1:N,1:M0) result in work(1:N,1:M0)

          call DLACPY( 'F', N, M0,X , N, work, N )

       end select
    end do

    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
    end if

    call wdeallocate_1d(ssb)
    call wdeallocate_1i(sisb)
    call wdeallocate_1i(sjsb)



  end subroutine dfeast_sst


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine zfeast_sge(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE CSR FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex(kind=kind(1.0d0)),dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(64) :: feastparam
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    complex(kind=kind(1.0d0)),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! csr-upper format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_SGE', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb => sb(1:nnza)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnza)

    else !! upper-csr or lower-csr to full csr


       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1z(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1z(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call zhcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)


    endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work1,N,M0,infoloc)
    call wallocate_2z(work2,N,M0,infoloc)

    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    opt=1
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1z(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=2
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    symbfact=.true.

    !feastparam(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
    feastparam(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


    ijob=-1 ! initialization




    do while (ijob/=0)


       call zfeast_rci(ijob,N,Ze,work1,work2,zAq,zSq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! Factorize (ZeB-A)


          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz

          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=3 ! complex and unsymmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2

          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if



       case(30) !! perform multiplication A*x(1:N,1:M0) result in work1(1:N,1:M0)

          call zhcsrmm('F',N,N,M0,ONEC,ssa,sisa,sjsa,X,ZEROC,work1)

       case(40) !! perform multiplication B*x(1:N,1:M0) result in work1(1:N,1:M0)

          call zhcsrmm('F',N,N,M0,ONEC,ssb,sisb,sjsb,X,ZEROC,work1)

       end select
    end do

    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work1)
    call wdeallocate_2z(work2)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif



  end subroutine zfeast_sge



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine zfeast_sst(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE CSR FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex(kind=kind(1.0d0)),dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(64) :: feastparam
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    integer :: M0
    double precision,dimension(*)  :: E
    complex(kind=kind(1.0d0)),dimension(N,*):: X
    integer :: mode
    double precision,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    double precision :: ddum1
!!!!! csr-upper format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_SST', -INFO+100 )
       RETURN
    END IF

    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1z(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=ONEC
    enddo
    sisb(n+1)=nnzb+1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) !may be overestimated
       call wallocate_1z(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)


    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work1,N,M0,infoloc)
    call wallocate_2z(work2,N,M0,infoloc)

    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=1
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1z(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
    opt=2
    call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    symbfact=.true.

    !feastparam(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
    feastparam(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


    ijob=-1 ! initialization


    do while (ijob/=0)


       call zfeast_rci(ijob,N,Ze,work1,work2,zAq,zSq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)



       select case(ijob)

       case(10) !! Factorize (ZeB-A)


          opt=3
          call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz


          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=3 ! complex and unsymmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2


          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,1:M0) result in work1(1:N,1:M0)

          call zhcsrmm('F',N,N,M0,ONEC,ssa,sisa,sjsa,X,ZEROC,work1)


       case(40) !! perform multiplication B*x(1:N,1:M0) result in work1(1:N,1:M0)


          call ZLACPY( 'F', N, M0,X , N, work1, N )

       end select
    end do


    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work1)
    call wdeallocate_2z(work2)
    call wdeallocate_2z(caux)

    call wdeallocate_1z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
    end if

    call wdeallocate_1z(ssb)
    call wdeallocate_1i(sisb)
    call wdeallocate_1i(sjsb)




  end subroutine zfeast_sst



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine sfeast_sge(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE CSR FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(64) :: feastparam
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex  :: Ze
    complex ,dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex ,dimension(:,:),pointer ::workc,caux
    real, dimension(:,:),pointer ::work,Aq,Sq
    real :: norm,nzero
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex ,parameter :: ONEC=(SONE,SZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! csr-upper format
    real,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k
    integer :: opt,nnza,nnzb,nnz


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SGE', -INFO+100 )
       RETURN
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    infoloc=0

    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if


       opt=1
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call scsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 

       call wallocate_1s(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1s(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call scsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb =>  sb(1:nnzb)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnzb)




    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1s(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1s(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call scsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call scsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)


    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2s(Aq,M0,M0,infoloc)
    call wallocate_2s(Sq,M0,M0,infoloc)
    call wallocate_2s(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
    opt=1
    call saddcsr(N,opt,SONE,ssa,sisa,sjsa,SONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1c(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=2
    call saddcsr(N,opt,SONE,ssa,sisa,sjsa,SONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if



    symbfact=.true.

    !feastparam(12)=1


    ijob=-1 ! initialization
    do while (ijob/=0)


       call sfeast_rci(ijob,N,Ze,work,workc,Aq,Sq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! Factorize (ZeB-A)

          opt=3
          call csaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!!!!!!!!!!!!!!!

          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=6 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
!!!!!!!! single precision pardiso (MKL only)
             IPARM(28)=1 ! pardiso single precision
             !          iparm(11)=1 !turn on MPS scaling
             !          iparm(13)=1 !weighted matching
             !          iparm(21)=0 !1x1 diagonal pivoting
             !          iparm(8)=-10 !turn on iterative refinement with double precision
             !          IPARM(10)=6
!!!!!!!!!!!!!!!

             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if
             symbfact=.false.
          end if





       case(20) !! solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc


          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
!!!!!!!! in case of NAN returned with single precision (and infoloc==0) ==> propose clean exit for FEAST
          if (infoloc==0) then
             infoloc=-2
             if (((real(workc(1,1))>0.0).or.(real(workc(1,1))<=0.0))) infoloc=0
          endif
!!!!!!!!!!
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(30) !! perform multiplication A*x(1:N,1:M0) result in work(1:N,1:M0)


          call scsrmm('U',N,N,M0,SONE,ssa,sisa,sjsa,X,SZERO,work)

       case(40) !! perform multiplication B*x(1:N,1:M0) result in work(1:N,1:M0)

          call scsrmm('U',N,N,M0,SONE,ssb,sisb,sjsb,X,SZERO,work)

       end select
    end do

    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2s(work)
    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1s(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif



  end subroutine sfeast_sge



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sfeast_sst(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: SPARSE CSR FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(64) :: feastparam
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex  :: Ze
    complex ,dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex ,dimension(:,:),pointer ::workc,caux
    real, dimension(:,:),pointer ::work,Aq,Sq
    real :: norm,nzero
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex ,parameter :: ONEC=(SONE,SZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! csr-upper format
    real,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k
    integer :: opt,nnza,nnzb,nnz


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SST', -INFO+100 )
       RETURN
    END IF

    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1s(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=SONE
    enddo
    sisb(n+1)=nnzb+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=1
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 

       call wallocate_1s(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)


!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       call wallocate_1s(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call scsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)


    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call wallocate_2s(Aq,M0,M0,infoloc)
    call wallocate_2s(Sq,M0,M0,infoloc)
    call wallocate_2s(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=1
    call saddcsr(N,opt,SONE,ssa,sisa,sjsa,SONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1c(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=2
    call saddcsr(N,opt,SONE,ssa,sisa,sjsa,SONE,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    symbfact=.true.


    ijob=-1 ! initialization


    do while (ijob/=0)


       call sfeast_rci(ijob,N,Ze,work,workc,Aq,Sq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! Factorize (ZeB-A)


          opt=3
          call csaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!!!!!!!!!!!!!!!



          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=6 ! complex and symmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             IPARM(28)=1 ! pardiso single precision (MKL)
             !IPARM(10)=6

             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc


          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc(:,1:M0),caux(:,1:M0),infoloc)
!!!!!!!! in case of NAN returned with single precision (and infoloc==0) ==> propose clean exit for FEAST
          if (infoloc==0) then
             infoloc=-2
             if (((real(workc(1,1))>0.0).or.(real(workc(1,1))<=0.0))) infoloc=0
          endif
!!!!!!!!!!


          if (infoloc/=0) then
             info=-2
             return
          end if


       case(30) !! perform multiplication A*x(1:N,1:M0) result in work(1:N,1:M0)

          call scsrmm('U',N,N,M0,SONE,ssa,sisa,sjsa,X,SZERO,work)

       case(40) !! perform multiplication B*x(1:N,1:M0) result in work(1:N,1:M0)

          call SLACPY( 'F', N, M0,X , N, work, N )


       end select
    end do

    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2s(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
    end if

    call wdeallocate_1s(ssb)
    call wdeallocate_1i(sisb)
    call wdeallocate_1i(sjsb)



  end subroutine sfeast_sst



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cfeast_sge(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE CSR FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B   
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*),target:: sa,sb
    integer,dimension(*),target:: isa,jsa,isb,jsb
    integer,dimension(64) :: feastparam
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex  :: Ze
    complex ,dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex ,dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    real :: norm,nzero
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex ,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! csr-upper format
    complex,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_SGE', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)

       nnzb=isb(n+1)-1
       ssb => sb(1:nnza)
       sisb => isb(1:n+1)
       sjsb => jsb(1:nnza)



    else !! upper-csr or lower-csr to full csr


       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1c(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1c(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call chcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call chcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)


    endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work1,N,M0,infoloc)
    call wallocate_2c(work2,N,M0,infoloc)

    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    opt=1
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1c(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=2
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    symbfact=.true.

    !feastparam(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
    feastparam(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


    ijob=-1 ! initialization

    do while (ijob/=0)


       call cfeast_rci(ijob,N,Ze,work1,work2,zAq,zSq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! Factorize (ZeB-A)


          opt=3
          call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!!!!!!!!!!!!!!!



          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=3 ! complex and unsymmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             IPARM(28)=1 ! pardiso single precision (MKL)
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2

          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
!!!!!!!! in case of NAN returned with single precision (and infoloc==0) ==> propose clean exit for FEAST
          if (infoloc==0) then
             infoloc=-2
             if (((real(work2(1,1))>0.0).or.(real(work2(1,1))<=0.0))) infoloc=0
          endif
!!!!!!!!!!


          if (infoloc/=0) then
             info=-2
             return
          end if


       case(30) !! perform multiplication A*x(1:N,1:M0) result in work1(1:N,1:M0)

          call chcsrmm('F',N,N,M0,ONEC,ssa,sisa,sjsa,X,ZEROC,work1)

       case(40) !! perform multiplication B*x(1:N,1:M0) result in work1(1:N,1:M0)

          call chcsrmm('F',N,N,M0,ONEC,ssb,sisb,sjsb,X,ZEROC,work1)

       end select
    end do

    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work1)
    call wdeallocate_2c(work2)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1c(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif



  end subroutine cfeast_sge



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cfeast_sst(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: SPARSE CSR FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A   
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*),target:: sa
    integer,dimension(*),target:: isa,jsa
    integer,dimension(64) :: feastparam
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc
    complex  :: Ze
    complex ,dimension(:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    complex ,dimension(:,:),pointer ::work1,work2,caux,zAq,zSq
    real :: norm,nzero
    real,parameter :: SONE=1.0e0, SZERO=0.0e0
    complex ,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
!!!!!for pardiso
    integer*8,dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum,nrhs
    real :: ddum1
!!!!! csr-upper format
    complex,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    logical :: symbfact
    integer :: i,k 
    integer :: opt,nnza,nnzb,nnz


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_SST', -INFO+100 )
       RETURN
    END IF

    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! identity B matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzb=n
    call wallocate_1i(sisb,n+1,infoloc)
    call wallocate_1c(ssb,nnzb,infoloc)
    call wallocate_1i(sjsb,nnzb,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    do i=1,n
       sisb(i)=i
       sjsb(i)=i
       ssb(i)=ONEC
    enddo
    sisb(n+1)=nnzb+1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       nnza=isa(n+1)-1
       ssa => sa(1:nnza)
       sisa => isa(1:n+1)
       sjsa => jsa(1:nnza)



    else !! upper-csr or lower-csr to full csr


       nnza=2*(isa(n+1)-1) !may be overestimated
       call wallocate_1c(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call chcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)


    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! FEAST 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work1,N,M0,infoloc)
    call wallocate_2c(work2,N,M0,infoloc)

    call wallocate_1i(isaz,n+1,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    opt=1
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    call wallocate_1c(saz,nnz,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
    opt=2
    call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

    symbfact=.true.

    !feastparam(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
    feastparam(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


    ijob=-1 ! initialization


    do while (ijob/=0)


       call cfeast_rci(ijob,N,Ze,work1,work2,zAq,zSq,feastparam,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)



       select case(ijob)

       case(10) !! Factorize (ZeB-A)


          opt=3
          call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz,isaz,jsaz) !! get saz
!!!!!!!!!!!!!!!

          if (symbfact) then
!!!!!!!!!! analysis step (symbolic factorizarion for pardiso)
             MAXFCT=1
             MNUM=1
             MTYPE=3 ! complex and unsymmetric 
             call pardisoinit(PT,MTYPE,IPARM)
             IPARM(3)=feastparam(10) !omp_num_threads !! openmp -number of threads
             IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
             MSGLVL=0 !0- no output, 1- output
             PHASE=11 
             IPARM(28)=1 ! pardiso single precision (MKL)
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

             symbfact=.false.
          end if



       case(20) !! solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2


          PHASE=23
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2(:,1:M0),caux(:,1:M0),infoloc)
!!!!!!!! in case of NAN returned with single precision (and infoloc==0) ==> propose clean exit for FEAST
          if (infoloc==0) then
             infoloc=-2
             if (((real(work2(1,1))>0.0).or.(real(work2(1,1))<=0.0))) infoloc=0
          endif
!!!!!!!!!!

          if (infoloc/=0) then
             info=-2
             return
          end if



       case(30) !! perform multiplication A*x(1:N,1:M0) result in work1(1:N,1:M0)


          call chcsrmm('F',N,N,M0,ONEC,ssa,sisa,sjsa,X,ZEROC,work1)


       case(40) !! perform multiplication B*x(1:N,1:M0) result in work1(1:N,1:M0)

          call CLACPY( 'F', N, M0,X , N, work1, N )


       end select
    end do



    PHASE=-1 
    call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz,isaz,jsaz,idum,M0,IPARM,MSGLVL,work2,caux,infoloc)
    if (infoloc/=0) then
       info=-2
       return
    end if

    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work1)
    call wdeallocate_2c(work2)
    call wdeallocate_2c(caux)

    call wdeallocate_1c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)
    end if

    call wdeallocate_1c(ssb)
    call wdeallocate_1i(sisb)
    call wdeallocate_1i(sjsb)




  end subroutine cfeast_sst



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
