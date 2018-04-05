#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi
module rm mpi fortran c
module load mpi/openmpi/2.0.1
module load gcc/5.2.0
module load hdf5/ompi/gnu/1.8.17-ompi2.0.1
module load cmake/3.7.2
#module load python
module load python/3.5.0
module list
