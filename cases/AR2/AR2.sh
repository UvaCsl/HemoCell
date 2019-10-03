#!/bin/bash
#SBATCH -t 5-00:00:00
#SBATCH -n 45        # number of needed processors, automatically calculate number of needed nodes
#SBATCH -p normal    # partition allocation

#make sure there is some sane modules loaded
source ../../scripts/lisa_env.sh

#srun works better with the environment than mpirun
mpirun ./prog config.xml
