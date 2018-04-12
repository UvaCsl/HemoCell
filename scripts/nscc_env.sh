#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi
module load cmake
module load gcc/4.9.3
module load hdf5/1.8.16/gcc493/serial
module load openmpi/gcc493/1.10.2
module list
