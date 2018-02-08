#!/bin/bash
[[ $_ != $0 ]] || (echo "Script is a subshell, this wont work, source it instead!"; exit)
module load gcc intel impi mkl hdf5 python/2.7.13 cmake
module list