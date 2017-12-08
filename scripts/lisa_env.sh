#!/bin/bash

# If module enviroment fails to init on Lisa uncomment next line 
#. /sara/sw/modules/module/init/bash

module load gcc/6.3.1
module load hdf5/gnu
module load openmpi/gnu/1.8.2
module load cmake/3.5.1
module list
