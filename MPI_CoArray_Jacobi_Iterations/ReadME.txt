FORTRAN implementation of Jacobi iteration algorithm using MPI & COARRAY

Deployment :
+++++++++++++

On untar, it will produce two separate folders as below. Switch to respective directory for executing implementation.
1] mpi
2] coarray

MPI
=====
Contains following files :
jacobimpi.f90
job_jacobimpi
Makefile

To execute mpi implementation, run following commands :
----------------------------------------------------------
module load intel-mpi/2018x
sbatch job_jacobimpi


Output is stored in jacobimpi<JOB ID>.out files.

####################################################################

COARRAY
=========
Contains following files :
jacobicoarray.f90
Makefile

To execute coarrayimplementation, run following commands :
----------------------------------------------------------
module load intel-mpi/2018x
export FOR_COARRAY_NUM_IMAGES=16
make
./jacobicoarray


Output is displayed on console on completion of execution.

####################################################################

Authors
========
Sarvesh Patil (ASU ID : 1213353386)


License
========
This project is licensed under the ASU License developed for APM 525 course

