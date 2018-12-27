CC=mpiifort

F90_FILES=precision.f90 henondim.f90 mpihenon.f90 
F90_OBJ=precision.o henondim.o mpihenon.o
APP=mpihenon
OPTIONS=-heap-arrays 100 -warn all -checkall
all:
	$(CC) $(OPTIONS)  -c $(F90_FILES)
	$(CC) $(OPTIONS)  $(F90_OBJ) -o $(APP)
#	$(CC) $(F90_FILES) -o $(F90_OBJ)
clean:
	rm -f *.o
	rm -f $(APP)
	rm -f *.mod
	rm -f *.out
	rm -f *.error


