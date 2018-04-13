#!/bin/bash
#PBS -q normal
#PBS -P Personal
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -l walltime=5:00

#Load Modules
module load cmake
module load gcc/4.9.3
module load hdf5/1.8.16/gcc493/serial
module load openmpi/gcc493/1.10.2
module list

#Run this script from the directory containing the executable, preferably on the ~/scratch folder
cd ~/scratch/stretchMalaria/
mpirun -n 4 ./stretchMalaria config.xml

#compiling and processing should be done on the login node!
