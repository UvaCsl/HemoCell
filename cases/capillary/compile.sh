#!/bin/bash
trap "exit" INT

echo "=========== Building =========="
date


if [ ! -d "./build" ]; then
  echo "* Running CMake..."
  mkdir build
  cd build
  cmake ..
  cd ..
fi

echo "* Compiling..."
cd build; 
make -j 4; 

if [ $? -ne 0 ]
then
    echo make failed
    exit 1
fi

cd ..
rm -rf tmp

echo "* Running simulation..."
mpirun -n 4 capillary config.xml

if [ $? -ne 0 ]
then
    echo mpirun failed
    exit 1
fi
echo "* Post process..."
../../scripts/batchPostProcess.sh 
date
echo "=========== Done ==========="
