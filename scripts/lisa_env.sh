#!/bin/bash

# If module enviroment fails to init on Lisa uncomment next line 
#. /sara/sw/modules/module/init/bash

module load gcc/5.4.0
module load hdf5/gnu
module load openmpi/gnu/1.4.3
module load cmake/3.5.1
module list
