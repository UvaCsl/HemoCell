#!/bin/bash
#SBATCH -t 05:00 
#SBATCH -n 15      #number of needed processors, automatically calculate number of needed nodes
#SBATCH -p short    #partition allocation

#make sure there is some sane modules loaded
source ../../scripts/lisa_env.sh

#srun works better with the environment than mpirun
mpirun ./pipeflow config.xml
