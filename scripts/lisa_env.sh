#!/bin/bash

# If module enviroment fails to init on Lisa uncomment next line 
#. /sara/sw/modules/module/init/bash

module load gcc/4.8.1
module load hdf5/gnu
module load openmpi/gnu/1.8.2
module load cmake/3.0.1
module list
