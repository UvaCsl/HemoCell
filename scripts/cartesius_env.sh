#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

# If module enviroment fails to init on Lisa uncomment next line
#. /sara/sw/modules/module/init/bash

module load 2020
module load CMake/3.16.4-GCCcore-9.3.0
module load GCC/9.3.0
module load HDF5/1.10.6-gompi-2020a
module load OpenMPI/4.0.3-GCC-9.3.0
module list
