#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

module load cmake
module load gcc/4.9
module unload mpi.ibm
module load mpi.intel
module load fftw
module load hdf5
module load java
