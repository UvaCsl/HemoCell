#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

# If module enviroment fails to init on Lisa uncomment next line 
#. /sara/sw/modules/module/init/bash
module load pre2019  #This is a temporary fix (24-2-2020)

module load pre2019

module load pre2019

module load GCC/7.3.0-2.30
module load hdf5
module load cmake/3.5.1
module load OpenMPI/3.1.1-GCC-7.3.0-2.30
module list
