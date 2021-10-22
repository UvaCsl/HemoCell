#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

# If module enviroment fails to init on Lisa uncomment next line
#. /sara/sw/modules/module/init/bash

module load 2021
module load CMake/3.20.1-GCCcore-10.3.0
module load GCC/10.3.0
module load HDF5/1.10.7-gompi-2021a
module load OpenMPI/4.1.1-GCC-10.3.0
module list
