#!/bin/bash
trap "exit" INT

if [ ! -d "./build" ]; then
  echo "(Running CMake...)"
  mkdir build
  cd build
  cmake ..
  cd ..
fi

echo "=========== Executing a clean run =========="
date

echo "* Compiling..."
cd build; make -j 4 2>&1 >/dev/null | grep 'Error\|error\|\*\*\*'; cd ..

echo "* Removing previous results..."
rm output.log
rm results.tgz
rm -r ./tmp/

echo "* Setting up domain..."
timeout 6 ./pipeflow configHO.xml > init.log

echo "* Running simulation..."
mpirun -np 4 ./pipeflow tmp/checkpoint.xml > output.log

echo "* Creating output files..."
../../scripts/batchPostProcess.sh > /dev/null

echo "* Compressing output..."
tar -cvzf results.tgz ./tmp > /dev/null

date
echo "=========== Done ==========="
