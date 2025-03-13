include ./make.inc

ARCH:=none
#%if $(ARCH) == none
#	@echo "Usage: make ARCH=<arch> all; make ARCH=<arch> clean "
#	@echo "!!!! <arch> is your architecture type and the name of the directory ../lib"
#	@echo "where the FEAST library will be created !!!! "
#	@echo
#%endif

help:
	@echo "Usage: make ARCH=<arch>; make ARCH=<arch> clean "
	@echo "!!!! <arch> is your architecture type and the name of the directory in ../lib"
	@echo "where the FEAST library will be created !!!! "
	@echo

#==============================================
# Name of Libraries 
#==============================================
NAME	= libfeast.a    # current library 
NAME_dense = libfeast_dense.a
NAME_banded = libfeast_banded.a
NAME_sparse = libfeast_sparse.a
NAME_rci = libfeast_rci.a

  


#sources Fortran90
CODE90 =   feast_dense.f90 feast_banded.f90 feast_sparse.f90 feast_rci.f90 
CODE90_rci = feast_rci.f90 
CODE90_banded = feast_banded.f90 feast_rci.f90 
CODE90_dense =  feast_dense.f90  feast_rci.f90 
CODE90_sparse =  feast_sparse.f90  feast_rci.f90 



#sources C 
CODEC= f90_noruntime.c feast_rci_wc.c feast_dense_wc.c feast_banded_wc.o feast_sparse_wc.c
CODEC_rci=  f90_noruntime.c feast_rci_wc.c
CODEC_dense=  f90_noruntime.c feast_dense_wc.c feast_rci_wc.c
CODEC_banded=  f90_noruntime.c feast_banded_wc.c feast_rci_wc.c
CODEC_sparse=  f90_noruntime.c feast_sparse_wc.c feast_rci_wc.c

# source F77
#CODE77=

#################################

OBJECTS = $(CODEC:.c=.o) $(CODE90_cwrapper:.f90=.o) $(CODE90:.f90=.o) #$(CODE77:.f=.o)

OBJECTS_rci = $(CODEC_rci:.c=.o) $(CODE90_cwrapper:.f90=.o) $(CODE90_rci:.f90=.o) # $(CODE77:.f=.o)

OBJECTS_dense = $(CODEC_dense:.c=.o) $(CODE90_cwrapper:.f90=.o) $(CODE90_dense:.f90=.o) #  $(CODE77:.f=.o) 

OBJECTS_banded = $(CODEC_banded:.c=.o) $(CODE90_cwrapper:.f90=.o) $(CODE90_banded:.f90=.o)  # $(CODE77:.f=.o)

OBJECTS_sparse = $(CODEC_sparse:.c=.o) $(CODE90_cwrapper:.f90=.o) $(CODE90_sparse:.f90=.o)  # $(CODE77:.f=.o) 


# .SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: clean all

#.f90.o: $(CODE90)
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

#.f.o: $(CODE77)
#%.o: %.f
#	$(FC) $(FFLAGS)  -c $< -o $@

#.f.o: $(CODEC)
%.o: %.f
	$(CC) $(CFLAGS)  -c $< -o $@


all: $(NAME)

$(NAME): $(OBJECTS)
	ar r $(NAME) $(OBJECTS)
	ar r $(NAME_rci) $(OBJECTS_rci)
	ar r $(NAME_dense) $(OBJECTS_dense)
	ar r $(NAME_banded) $(OBJECTS_banded)
	ar r $(NAME_sparse) $(OBJECTS_sparse)
	if test -d ../lib/$(ARCH); then \
	echo "directory $(ARCH) already exists"; \
	else echo "create directory $(ARCH)" && mkdir ../lib/$(ARCH); fi;
	mv *.a ../lib/$(ARCH)
	#rm -f $(OBJECTS) 

clean: 
	rm -f $(OBJECTS) 
	rm ../lib/$(ARCH)/*.a
	rm -r ../lib/$(ARCH) 
touch:
	@touch $(CODE90)











