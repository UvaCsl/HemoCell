#!/bin/bash
#SBATCH -t 00:10:00 #Ten minutes, adjust to need
#SBATCH -n 1       #number of needed processors, automatically calculate number of needed nodes
#SBATCH -p short    #partition allocation
module rm mpi fortran c
module load mpi/openmpi/2.0.1

srun ../performance_testing configs/config_1.xml 
