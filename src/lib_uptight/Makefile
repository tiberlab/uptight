SHELL = /bin/bash

.SUFFIXES:
.SUFFIXES: .f .f90 .F90 .o

include ../make.common

HERE = $(shell pwd)
FEASTSRC = ../FEAST/1.0/src
LFEAST = #-L ../FEAST/1.0/lib/$(ARCH1) -lfeast -lfeast_sparse
LIBARPACK = #../../lib_arpack/libarpack_debian_ifort_mkl.a 


FPPFLAGS += -D__REVISION='"$(shell git describe --always --tags --dirty )"'
FPPFLAGS += -D__MODIFIED='"$(shell git log -n 1 | grep Date | awk '{print $$3, $$4, $$5, $$6}' )"'
FPPFLAGS += -D__COMPDATE='"$(shell date +'%F')"'

#LIBMPI=-L$(MPILIB)

# Link against LAPACK and BLAS libraries
LIBMKLIB += -llapack -lblas

#LIBS =	$(LIBMKL) $(LIBMPI) $(LIBCUDA) 

ifdef UPT_MPI
  FPPFLAGS += -DUPT_MPI 
  FFLAGS += $(INCLUDE)
  CFLAGS += $(INCLUDE)
  #export OMPI_FC=$(FC)
  #export OMPI_FCFLAGS=$(FCFLAGS) $(OMPFLAGS) $(MPI_INCL)
  #export OMPI_LDFLAGS=$(OMPFLAGS)
  #export OMPI_LIBS=$(LIBS)
  #FC = $(MPIDIR)/bin/$(MPIFC) 
endif	



SOURCES = precision.f90 \
 	  constants.f90 \
          clock.f90 \
 	  globals.f90 \
 	  rcm_module.f90 \
 	  allocation.f90 \
 	  vectors.f90 \
 	  supercell.f90 \
 	  listtypes.f90 \
 	  type_defs.f90 \
 	  input_output.f90 \
		errors.f90 \
 	  input_data.f90 \
          checks.f90 \
 	  states_and_couplings.f90 \
 	  lists.f90 \
 	  sort.f90 \
 	  neighbours.f90 \
 	  struct_building.f90 \
          sparse_matrix.f90 \
 	  upt_param.f90 \
          alloys.f90 \
 	  latt_fun.f90 \
 	  rota.f90 \
 	  sparse_matrix.f90 \
 	  TB_ham.f90 \
 	  lanczos_driver.f90 \
 	  JD_driver.f90 \
 	  lapack_driver.f90 \
 	  savemofile.f90 \
 	  f_exceptions.f90 #\
 	  #feast_driver.f90 

SOURCES2 = mpi_globals.F90 \
	   exceptions.F90 \
	   lanczos.F90 \
	   JD.F90 \
           sparse_numrec.F90 \
	   uptight.F90
						 

FEASTCSOURCES = f90_noruntime.c \
                feast_banded_wc.c \
                feast_dense_wc.c \
                feast_rci_wc.c \
                feast_sparse_wc.c 

FEASTFSOURCES = feast_banded.f90 \
                feast_dense.f90 \
                feast_rci.f90 \
                feast_sparse.f90 

OBJECTS = $(SOURCES:.f90=.o) $(SOURCES2:.F90=.o)

ifdef CUDA
  OBJECTS += lanczos_cuda.o JD_N_PC_Mxprec_cuda.o
  ifdef UPT_MPI
    OBJECTS += JD_N_PC_Mxprec_ll_cuda.o
  endif
endif
ifdef UPT_MPI
  OBJECTS += JD_N_PC_Mxprec_ll_cpu.o
endif

#$(FEASTCSOURCES:.c=.o) $(FEASTFSOURCES:.f90=.o)  
# UPTIGHT_linux_mod_sige: $(OBJECTS)
# UPTIGHT_linux_sige_LAPACK: $(OBJECTS)
# UPTIGHT_linux_sige_QW_LAPACK: $(OBJECTS)


all: PARAMETERIZER

# UPT_TEST6: test6.o f_exceptions.o $(OBJECTS)  
#		$(FC) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFEAST) $(LFLAGS) $(LIBCUDA) $(LIBMKLIB)

# UPT_TEST7: test7.o f_exceptions.o $(OBJECTS) 
#		$(FC) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFEAST) $(LFLAGS) $(LIBCUDA) $(LIBMKLIB)

PARAMETERIZER: parameterizer.o f_exceptions.o $(OBJECTS) 
		$(FC) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFEAST) $(LFLAGS) $(LIBCUDA) $(LIBMKLIB)

# pUPT: test8.o f_exceptions.o  $(OBJECTS) 
#		@$(FC) $(SHOW) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFEAST) $(LIBMPI) $(LIBCUDA) $(LIBMKLIB)
#		$(FC) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFEAST) $(LIBMPI) $(LIBCUDA) $(LIBMKLIB)


getdist: $(OBJECTS) f_exceptions.o getdist.o 
		$(FC) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFEAST) $(LFLAGS) 

mo: $(OBJECTS) f_exceptions.o mo.o 
		$(FC) $(LDFLAGS) -o $@ $^ $(LIBARPACK) $(LFLAGS) $(LIBCUDA) $(LIBMKLIB) 




FEAST:
	(cd $(FEASTSRC); $(MAKE) ARCH=$(ARCH1) libfeast.a)

	
uptight.a: $(OBJECTS) 
	ar cq $@ $^


ifdef CUDA

CUDAARCHS = -arch=sm_$(word 1,$(CUDA_VER))
CUDAARCHS += $(foreach a,$(CUDA_VER),-gencode=arch=compute_$(a),code=sm_$(a))

lanczos_cuda.o: lanczos_cuda.cu
	$(NVCC) $(CUDAARCHS) \
	       	lanczos_cuda.cu -c -w \
	 	-I$(CUDADIR)/include -I$(CUDADIR)/src $(NVFLAGS) \
	 	$(MKLINC)  

ifdef UPT_MPI
JD_N_PC_Mxprec_cuda.o: JD_N_PC_Mxprec_cuda.cu
	$(NVCC) $(CUDAARCHS) \
		JD_N_PC_Mxprec_cuda.cu -c -w \
	 	-I$(CUDADIR)/include -I$(CUDADIR)/src $(NVFLAGS) \
	 	$(MKLINC) 

JD_N_PC_Mxprec_ll_cuda.o: JD_N_PC_Mxprec_ll_cuda.cu
	$(NVCC) $(CUDAARCHS) \
		JD_N_PC_Mxprec_ll_cuda.cu -c -w \
	 	-I$(CUDADIR)/include -I$(CUDADIR)/src $(NVFLAGS) \
	 	$(MKLINC) -I$(MPIDIR)/include 

endif
endif

ifdef UPT_MPI
JD_N_PC_Mxprec_ll_cpu.o: JD_N_PC_Mxprec_ll_cpu.c
	$(CC) $(CFLAGS) -c JD_N_PC_Mxprec_ll_cpu.c $(COPENMP) $(MKLINC) -o $@
endif

%.o : %.F90
	$(CPP) $(FPPFLAGS) $< > $*.f90
	$(FC) -c $(FFLAGS) -o $@ $*.f90 
	rm -f $*.f90

.f90.o:
	$(FC) -c $(FFLAGS) -o $@ $< 


#%.o: %.F90
#	gcc -E $(FPP) $< > $*.f90 
#	$(FC) -nocpp -c $(FFLAGS) -o $@ $*.f90
#%.o: %.f90
#	$(FC) -c $(FFLAGS) $(OPENMP) -o $@ $< 
#%.o: %.f
#	$(FC) -c $(FFLAGS) $(OPENMP) -o $@ $<


$(SOURCES):

.PHONY: clean check
clean:
	rm -f *.a *.o; rm -f *.mod UPT_TEST6 UPT_TEST7 PARAMETERIZER pUPT

check:
	@echo $(have_cuda)

tar_source:
	tar cvf uptight.tar ${SOURCES} Makefile


include ./Makefile.dep
