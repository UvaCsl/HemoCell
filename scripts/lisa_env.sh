#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

# If module enviroment fails to init on Lisa uncomment next line 
#. /sara/sw/modules/module/init/bash

module load gcc/6.3.0
module load hdf5/intel/1.8.16-parallel
module load cmake/3.5.1
module load openmpi/gnu/
module list
