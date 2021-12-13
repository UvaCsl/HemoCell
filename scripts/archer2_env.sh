#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi


module load cmake/3.21.3
module load cray-hdf5/1.12.0.7
module load cray-mpich/8.1.9
module load cray-python/3.9.4.1
module list

export CC=$(which clang)
export CXX=$(which clang++)
export LDFLAGS="-L${CRAYLIBS_X86_64} $(CC --cray-print-opts=libs) -lmpi"

# compile by first running: cmake ../ -DMPI_CXX_COMPILER=$(which CC) -DMPI_C_COMPILER=$(which cc) -DMPI_COMPILER_FLAGS="--cray-print-opts=all"
# from hemocell/build directory