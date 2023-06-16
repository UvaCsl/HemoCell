#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi


module load cmake
module load cray-hdf5
module load cray-python
module list

export CC=cc
export CXX=CC

# Comment out line 130 in CMakeLists.txt to avoid messy unused parameter warnings.

# compile by first running: cmake ..
# from hemocell/build directory
