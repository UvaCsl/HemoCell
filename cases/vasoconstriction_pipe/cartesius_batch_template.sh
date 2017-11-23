#!/bin/bash
#SBATCH -t 100:00:00 #Ten minutes, adjust to need
#SBATCH -n 40       #number of needed processors, automatically calculate number of needed nodes
#SBATCH -p normal    #partition allocation
		#name     max   type description
		# normal    5-00:00:00    720       thin      the default partition
		# fat       5-00:00:00     16       fat       fat node partition
		# short       01:00:00    720       thin      for short jobs
		# gpu       5-00:00:00     48       gpu       gpu nodes for production runs
		# gpu_short   01:00:00     64       gpu       gpu nodes for test runs
		# staging   5-00:00:00   1 core     server    for accessing the archive and external systems
module rm mpi fortran c
module load mpi/openmpi/1.10.2
module load gcc/5.2.0

# This is necessary for runs over 128 nodes -> otherwise the communication queue overflows
export OMPI_MCA_btl_openib_receive_queues="X,128,256,192,128:X,2048,256,128,32:X,12288,256,128,32:X,65536,256,128,32"

#Place the actual command here
#Dont forget to place this file in the directory you want to run from
mpirun -n 40 ./vasoconstriction_pipe config.xml #dont forget to adjust -n !!
