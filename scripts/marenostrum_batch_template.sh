#!/bin/bash
#SBATCH --job-name=mpi
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=4
#SBATCH --qos=training
#SBATCH --time=00:04:00
srun ./pipeflow config.xml
