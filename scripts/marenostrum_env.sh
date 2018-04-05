#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi
module load gcc intel impi mkl hdf5 python/2.7.13 cmake
module list
