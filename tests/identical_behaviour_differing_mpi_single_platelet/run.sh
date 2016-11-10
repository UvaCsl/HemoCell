#!/bin/bash
echo "======Building the binary======"
mkdir -p build
cd build
cmake ../ > /dev/null 2> /dev/null
make -j8 > /dev/null 2>/dev/null
make -j8 > /dev/null 2>/dev/null
cd ../
echo "======Setting up tests========="
rm -r case1 case2 2> /dev/null
mkdir -p case1 case2
cp cells.pos case1/cells.pos
cp cells.pos case2/cells.pos
cp configHO.xml case1/
cp configHO.xml case2/
cp pipeflow case1/
cp pipeflow case2/
cp tube.stl case1/
cp tube.stl case2/
echo "======Running Cases==========="
cd case1
./pipeflow configHO.xml 2> /dev/null  > /dev/null &
cd ../case2
timeout 5 ./pipeflow configHO.xml 2> /dev/null > /dev/null
mpirun -n 2 ./pipeflow tmp/checkpoint.xml 2> /dev/null > /dev/null &
cd ../
wait
echo "=====Checking Output=========="
difference=$(python compare_hdf5.py)
if [ -z "$difference" ]; then
	echo "No difference"
else
	echo "ERROR: similar case output different values, undefined behaviour!"
  echo "$difference"
	echo "Not cleaning up..."
  exit
fi
echo "=====Cleaning Up============="
rm -r case1 case2 build
rm pipeflow
