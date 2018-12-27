CC=ifort

F90_FILES=precision.f90 lsq.f90 henondim.f90 main.f90
F90_OBJ=precision.o lsq.o henondim.o main.o
APP=main
OPTIONS=-heap-arrays 100 -mkl -warn all -qopenmp # -check all  or -check bounds
	
EXPORT=$(shell export)
all:
	$(CC) $(OPTIONS)  -c $(F90_FILES) 
	$(CC) $(OPTIONS)  $(F90_OBJ) -o $(APP)

clean:
	rm -f *.o
	rm -f $(APP)
	rm -f *.mod


