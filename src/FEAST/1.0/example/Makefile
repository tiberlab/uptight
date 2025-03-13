include ./make.inc

help:
	@echo "Usage: make {all allC allF dense banded sparse rall rallC rallF rdense rbanded rsparse"
	@echo "                             clean cleanC cleanF clean_dense clean_banded clean_sparse}"
	@echo "************* Compilation ************************ "
	@echo "make all    -- compile all C and F90 examples "
	@echo "make allC   -- compile all C examples "
	@echo "make allF   -- compile all F90 examples "
	@echo "make dense  -- compile all C and F90 examples for feast_dense interface "
	@echo "make banded -- compile all C and F90 examples for feast_banded interface "
	@echo "make sparse -- compile all C and F90 examples for feast_sparse interface "
	@echo "************* Both Compilation and Run ********** "
	@echo "make rall    -- compile and run all C and F90 examples "
	@echo "make rallC   -- compile and run all C examples "
	@echo "make rallF   -- compile and run all F90 examples "
	@echo "make rdense  -- compile and run all C and F90 examples for feast_dense interface "
	@echo "make rbanded -- compile and run all C and F90 examples for feast_banded interface "
	@echo "make rsparse -- compile and run all C and F90 examples for feast_sparse interface "
	@echo "************* Cleaning ************************** "
	@echo "make clean  -- clean all C and F90 examples"
	@echo "make cleanC -- clean all C examples"
	@echo "make cleanF -- clean all F90 examples"
	@echo "make clean_dense  -- clean all C and F90 examples for feast_dense interface "
	@echo "make clean_banded -- clean all C and F90 examples for feast_banded interface "
	@echo "make clean_sparse -- clean all C and F90 examples for feast_sparse interface "
	@echo "************* WARNING ************************** "
	@echo "!!!!Before use: Please correct accordingly compiler and libraries paths, change compiler options " 
	@echo " in file make.inc !!!!"
	@echo

#########################
dense: 
	cd ./C/1_dense ; make all ; cd $(FEASTROOT)/example/Fortran/1_dense; make all ; \

clean_dense:
	cd ./C/1_dense ; make clean; cd $(FEASTROOT)/example/Fortran/1_dense; make clean; \

banded:
	cd ./C/2_banded; make all ; cd $(FEASTROOT)/example/Fortran/2_banded; make all ; \

clean_banded:
	cd ./C/2_banded; make clean; cd $(FEASTROOT)/example/Fortran/2_banded; make clean; \

sparse:
	cd ./C/3_sparse; make all ; cd $(FEASTROOT)/example/Fortran/3_sparse; make all ; \
 
clean_sparse:
	cd ./C/3_sparse; make clean; cd $(FEASTROOT)/example/Fortran/3_sparse; make clean; \

all: dense banded sparse

allC:
	cd ./C/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all;\

allF:
	cd ./Fortran/1_dense; make all; cd ../2_banded; make all; cd ../3_sparse; make all;\

clean: clean_dense clean_banded clean_sparse

cleanC:
	cd ./C/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean;\

cleanF:
	cd ./Fortran/1_dense; make clean; cd ../2_banded; make clean; cd ../3_sparse; make clean;\


rdense: 
	cd ./C/1_dense ; make rall ; cd $(FEASTROOT)/example/Fortran/1_dense; make rall ; \

rbanded:
	cd ./C/2_banded; make rall ; cd $(FEASTROOT)/example/Fortran/2_banded; make rall ; \

rsparse:
	cd ./C/3_sparse; make rall ; cd $(FEASTROOT)/example/Fortran/3_sparse; make rall ; \


rall: rdense rbanded rsparse

rallC:
	cd ./C/1_dense; make rall; cd ../2_banded; make rall; cd ../3_sparse; make rall;\

rallF:
	cd ./Fortran/1_dense; make rall; cd ../2_banded; make rall; cd ../3_sparse; make rall;\


