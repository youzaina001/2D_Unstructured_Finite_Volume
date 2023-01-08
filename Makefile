PROGRAM = exec
OBJECTS = parameters.o cartesian_mesh.o flux.o solver.o high_order.o initial_condition.o boundary_conditions.o time_scheme.o vtk.o main.o
FC = gfortran
FFLAGS = -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all -finit-real=snan -fbacktrace

$(PROGRAM): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(PROGRAM) $(OBJECTS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)