FC = gfortran
OPT = -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all -finit-real=snan -fbacktrace
EXE = exec

$(EXE): mod_params.o data.o read_mesh.o mesh.o flux.o vtk.o main.o
	gfortran -o $@ $^
%.o: %.f90
	gfortran -c $<

clean:
	rm -f *.o *.mod *.vtk $(EXE)