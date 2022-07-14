#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi


module load mpt
module load gcc
module load cmake
module load petsc/3.13.2-mpt
module load hdf5parallel/1.10.6-gcc8-mpt225
module list
