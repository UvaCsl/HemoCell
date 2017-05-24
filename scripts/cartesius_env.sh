#!/bin/bash
[[ $_ != $0 ]] || (echo "Script is a subshell, this wont work, source it instead!"; exit)
module rm mpi fortran c
module load mpi/openmpi/1.10.2
module load gcc/5.2.0
module load hdf5/ompi/gnu
module load python
module list
