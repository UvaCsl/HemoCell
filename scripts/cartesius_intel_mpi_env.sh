#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi
module rm mpi fortran c
module load c/intel/18.0.1
module load mpi/impi/18.0.1
module load hdf5/impi/gnu
module load cmake/3.7.2
module load python
module list
export I_MPI_CXX=icpc
