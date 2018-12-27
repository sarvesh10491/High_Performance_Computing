CC=mpiifort

F90_MPI= jacobimpi.f90 
F90_OBJ_MPI=jacobimpi.o
APP_MPI=jacobimpi
OPTIONS=-heap-arrays 100 -traceback 
all:
	$(CC) $(OPTIONS)  -c $(F90_MPI)
	$(CC) $(OPTIONS)  $(F90_OBJ_MPI) -o $(APP_MPI)
clean:
	rm -f *.o
	rm -f $(APP_MPI)
	rm -f *.mod
	rm -f *.out
	rm -f *.error


