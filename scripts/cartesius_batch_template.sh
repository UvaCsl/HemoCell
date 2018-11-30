#!/bin/bash
#SBATCH -t 00:10:00 #Ten minutes, adjust to need
#SBATCH -n 40       #number of needed processors, automatically calculate number of needed nodes
#SBATCH -p short    #partition allocation
		#name     max   type description
		# normal    5-00:00:00    720       thin      the default partition
		# fat       5-00:00:00     16       fat       fat node partition
		# short       01:00:00    720       thin      for short jobs
		# gpu       5-00:00:00     48       gpu       gpu nodes for production runs
		# gpu_short   01:00:00     64       gpu       gpu nodes for test runs
		# staging   5-00:00:00   1 core     server    for accessing the archive and external systems
module rm mpi fortran c
module load OpenMPI
module load hdf5/serial/intel/1.10.0-patch1

#srun works better with the environment than mpirun
srun ./pipeflow tmp/checkpoint.xml 
