#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi
module load pre2019
module rm mpi fortran c
module load OpenMPI/3.0.0-GCC-6.4.0-2.28
module load cmake/3.7.2
module load hdf5/serial/intel/1.10.0-patch1
module load python/3.5.2
module list
