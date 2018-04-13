#!/bin/bash
#PBS -q normal
#PBS -P Personal
#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l walltime=24:00:00

#Run this script from the directory containing the executable, preferably on the ~/scratch folder
mpirun ./stretchCell config.xml

#compiling and processing should be done on the login node!
