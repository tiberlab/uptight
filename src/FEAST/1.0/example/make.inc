############################################################################
#  PATHs to Compilers and libraries
############################################################################ 

########## C compiler and flags (for C examples)
CC     = gcc
ICFLAGS = -O3

########## F90 compiler and flags (for Fortran examples)
FORT  =  ifort
FFLAGS = -O3

######### PATH to FEAST (main directory)- can be defined as shell variable
FEASTROOT=/home/gpenazzi/Develop/FEAST/1.0

########  PATH to INCLUDE FEAST files (do not change)
#INCL=-I/$(FEASTROOT)/include


######### PATH to libraries: FEAST (path to compiled FEAST library for your architecture)
#########                    PARDISO (mandatory if sparse examples are compiled)
#########                    LAPACK (mandatory) 
#########                    BLAS (mandatory)
######### Remark: Intel MKL includes Pardiso, Lapack and Blas

LOCLIBS = -L/$(FEASTROOT)/lib/x86_64 -L/usr/pack/intel_mkl-10.2-gp/lib/em64t -Wl,-rpath,/usr/pack/intel_mkl-10.2-gp/lib/em64t

######## Libraries PARDISO (if sparse examples are compiled), LAPACK, BLAS, (and lm for C examples)

#LIBS=-lmkl_solver -lmkl_intel_lp64  -Wl,--start-group -lmkl_intel_thread -lmkl_core  -Wl,--end-group -liomp5 -lpthread -lm
#LIBS=-lmkl_solver -Wl,--start-group -lmkl_intel -lmkl_intel_thread -lmkl_core  -Wl,--end-group -liomp5 -lpthread -lm
LIBS= -lmkl_solver_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

