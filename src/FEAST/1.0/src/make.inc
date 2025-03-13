############################################################################
#  Compilers 
############################################################################ 

########## F90 compiler and flags
F90 = gfortran    
#ifort-11.1
F90FLAGS = -ffree-line-length-none 
#-O3

########## C compiler and flags
CC     = gcc
ICFLAGS = -O3 
CFLAGS += -fPIC
F90FLAGS += -fPIC
FFLAGS += -fPIC


