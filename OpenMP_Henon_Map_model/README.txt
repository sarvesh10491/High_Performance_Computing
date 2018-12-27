Getting Started

    These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 
    See deployment for notes on how to deploy the project on a live system.


Installing

Download below files in user directory on your machine/server


Deployment

   job_henon file contains following batch job related parameters that needs to be set to run respective batch job.
   #SBATCH -n <no. of cores>
   #SBATCH --ntasks-per-node <defines tasks per node>
   #SBATCH -J <current job name>
   #SBATCH -o <output_file_name>
   #SBATCH -e <error_file_name>

   export OMP_NUM_THREADS=<no. of threads>

   This job file also contains make commands which compiles the source files.
   module load intel/2018x 

   To run batch job, run following
   sbatch job_henon


Expected results

   On completion of batch job execution, output is attached in file specified in job description file.
  

Authors

Sarvesh Patil (ASU ID : 1213353386)
Shubham Nandanwankar(ASU ID :1213350370)

License

This project is licensed under the ASU License developed for APM 525 course

