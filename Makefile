#FC = mpiifort -r8 -i8 -O3 -pg -g -mcmodel=medium
#FC = mpiifort -r8 -i8 -O3 -openmp -pg -g -mcmodel=medium
#FC = mpif90 -O3 -pg -g -mcmodel=medium
#FC = mpif90 -O3 -openmp -pg -g -mcmodel=medium
#FC = ifort -r8 -i8 -O3 -xT -mcmodel=medium -i-dynamic -vec-report0 -warn nounused
#FC = ifort -r8 -i8 -O3 -xT -mcmodel=medium -i-dynamic -vec-report0 -warn nounused -fpe0 -traceback
FC = ifort -r8 -i8 -O3 -openmp -mcmodel=medium -i-dynamic -vec-report0 -warn nounused -check all -fpe0 -traceback 

EXENAME =  3_exec

OBJS= modulex.o lica.o geom.o pcg.o two_phase.o 

$(EXENAME) : $(OBJS)
	$(FC) -o $(EXENAME) $(OBJS)

.f.o :
	$(FC)  -c $<

clean :
	rm -f *.o *.mod