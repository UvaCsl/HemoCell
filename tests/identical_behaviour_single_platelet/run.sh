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
./pipeflow configHO.xml 2> /dev/null > /dev/null &
cd ../
wait
echo "=====Checking Output=========="
difference=$(find ./case1/ -exec basename {} \; |grep h5 | grep -v RBC | xargs -I filename -n1 h5diff -c case1/tmp/hdf5/filename case2/tmp/hdf5/filename |grep -v "Not comparable")
if [ -z $difference ]; then
	echo "No difference"
else
	echo "ERROR: similar case output different values, undefined behaviour!"
fi
echo "=====Cleaning Up============="
rm -r case1 case2 build
rm pipeflow
