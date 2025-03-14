SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .C .f90 .o

include ../make.common

API = $(shell pwd)

TARGET = libuptight_$(ARCH1).a

CSOURCES = UptInitException.C \
      	   UptSolverException.C 

MODULES = precision.mod \
	globals.mod \
	upt_param.mod \
	uptight.mod \
	savemofile.mod \
	sparse_matrix.mod \
	tb_ham.mod		

OBJECTS = $(CSOURCES:.C=.o)  

all: libuptight
	
tests: test_upt_api_c test_upt_api_f

.f90.o:
	$(FC) -c $(FFLAGS) $(OPENMP) -o $@ $< 
.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<
.C.o:
	$(CC) -c $(CFLAGS) -o $@ $<

libs:
	(cd ../lib_uptight; $(MAKE) uptight.a; cp exception_codes.h $(MODULES) $(API))

$(TARGET): upt_api.o    
	ar -r $@ $^ ../lib_uptight/*.o 



libuptight: $(TARGET) cpp_wrapped_binding 
		
upt_api.o: uptightAPICommon.o
uptightAPICommon.o: upt_param.mod
upt_param.mod: libs

c_binding: upt_api.f90
	./bind_fortran -f c \
            -m ./binding/map_upt_c.txt \
            -b ./binding/begin_c.txt -n -t \
            -e ./binding/end_c.txt \
            -s "$(FNAME_SUFFIX)" -p "$(FNAME_PREFIX)" \
            $(if $(FNAME_UPPERCASE),-u) $^ > uptight.h

cpp_wrapped_binding: upt_api.f90
	./bind_fortran -f c++ \
            -m ./binding/map_upt_cpp_wrapped.txt \
            -b ./binding/begin_cpp_wrapped.txt -t \
            -e ./binding/end_cpp_wrapped.txt \
            -s "$(FNAME_SUFFIX)" -p "$(FNAME_PREFIX)" \
            $(if $(FNAME_UPPERCASE),-u) $^ > uptight.h

test_upt_api_c: cpp_wrapped_binding test_upt_api_c.o $(TARGET) 	
	$(CC) $(LDFLAGS) -o $@ test_upt_api_c.o $(TARGET) $(LCFLAGS)

test_upt_api_f: test_upt_api_f.o $(TARGET)
	$(FC) $(LDFLAGS) -o $@ test_upt_api_f.o $(TARGET) $(LFLAGS)


.PHONY: clean
clean:
	rm -f *.o; rm -f *.mod *.a
	(cd ../lib_uptight; $(MAKE) clean)

tar_source:
	tar cvf uptight.tar ${SOURCES} Makefile




