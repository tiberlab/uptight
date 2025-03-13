
FC = g95 
FLAGS = -W #-check all

OBJS = precision.o allocation.o prova.o 


target: prova

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FLAGS) -c $< 

prova: $(OBJS)
	$(FC) -o $@  $(OBJS)

clean:
	rm -rf $(OBJS) *.mod prova
