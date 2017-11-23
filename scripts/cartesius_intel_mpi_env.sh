#!/bin/bash
[[ $_ != $0 ]] || (echo "Script is a subshell, this wont work, source it instead!"; exit)
module rm mpi fortran c
module load mpi/impi/17.0.5
module load gcc/5.2.0
module load hdf5/impi/gnu
module load cmake/3.7.2
module load python
module list
