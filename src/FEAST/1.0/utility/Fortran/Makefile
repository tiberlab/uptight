help:
	@echo "Usage: make {all clean} "
	@echo "!!!!Please correct accordingly compiler and libraries paths, change compiler options,etc.  " 
	@echo " in the "User Specifications region" in file Makefile !!!!"
	@echo

###########################################
# USER SPECIFICATIONS #####################
###########################################

#==============================================
# Specify Compiler
#==============================================
FC = ifort

#==============================================
# Specify OPTIMIZATION flags
#==============================================
FLAGS =-O3  

#==============================================
# Specify PATH to FEAST (main directory)- can be defined as shell variable
#==============================================
FEASTROOT=/home/polizzi/FEAST/1.0

#==============================================
# Specify libraries path (FEAST+PARDISO+LAPACK+BLAS) 
#==============================================
LOCLIBS = -L/$(FEASTROOT)/lib/x64 -L/opt/intel/mkl/10.2.1.017/lib/em64t/  

#==============================================================
# Specify all the libraries (feast - pardiso- lapack- blas) mkl==pardiso - lapack -blas  
#==============================================================
LIB=$(LOCLIBS) -lfeast_sparse -lmkl_solver -lmkl_intel_lp64  -Wl,--start-group -lmkl_intel_thread -lmkl_core  -Wl,--end-group -liomp5 -lpthread 

###########################################
# END USER SPECIFICATIONS #################
###########################################

#==============================================================
# Include file path (unchanged)
#==============================================================
INCL = -I/$(FEASTROOT)/include 


#============================================================
# COMPILE and LINK
#============================================================
EXAMPLES = driver_feast_sparse 


all: examples 


examples: 
	@for file in $(EXAMPLES); \
	do \
		echo $(FC)  $(FLAGS) $(INCL) -c $$file.f90 ;\
		$(FC)  $(FLAGS) $(INCL) -c $$file.f90 ;\
		echo $(FC)   -o $$file $$file.o $(LIB) ;\
		$(FC)   -o $$file $$file.o $(LIB) ;\
	done


#=========================================================================
# Generate for each Sub-Program its Object File: do NOT link them yet
#==========================================================
.SUFFIXES: .f90 .mod
%.o: %.f90
	$(FC) $(FLAGS) -c $<


#==========================================================
# Clean up directory: delete object files and modules
#==========================================================
clean: 
	-@rm  $(EXAMPLES) *.o
