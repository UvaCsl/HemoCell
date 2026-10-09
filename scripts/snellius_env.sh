#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

# If module enviroment fails to init on Lisa uncomment next line
#. /sara/sw/modules/module/init/bash

module load 2025
module load GCC/14.2.0
module load CMake/3.31.3-GCCcore-14.2.0
module load OpenMPI/5.0.7-GCC-14.2.0
module load HDF5/1.14.6-gompi-2025a
module list
