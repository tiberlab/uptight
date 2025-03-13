/*
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
*/




extern void dfeast_sge_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_sst_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_sge_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_sst_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

extern void sfeast_sge_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void sfeast_sst_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_sge_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_sst_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);



void dfeast_sge(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void dfeast_sst(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void zfeast_sge(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void zfeast_sst(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}



void sfeast_sge(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void sfeast_sst(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void cfeast_sge(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void cfeast_sst(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}




void DFEAST_SGE(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SST(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void ZFEAST_SGE(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void ZFEAST_SST(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam[64],double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}



void SFEAST_SGE(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void SFEAST_SST(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void CFEAST_SGE(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_sge_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}

void CFEAST_SST(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam[64],float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_sst_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}



