FC = gfortran
OPT = -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all -finit-real=snan -fbacktrace
EXE = exec

$(EXE): parameters.o cartesian_mesh.o flux.o solver.o initial_condition.o boundary_conditions.o time_scheme.o vtk.o main.o
	gfortran -o $@ $^
%.o: %.f90
	gfortran -c $<

clean:
	rm -f *.o *.mod $(EXE)