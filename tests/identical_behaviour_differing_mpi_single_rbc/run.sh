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
timeout 10 ./pipeflow configHO.xml 2> /dev/null > /dev/null
mpirun -n 2 ./pipeflow tmp/checkpoint.xml 2> /dev/null > /dev/null &
cd ../
wait
echo "=====Merging h5files (crude)=="
find ./case1/ -exec basename {} \; |grep h5 | grep RBC |grep -v _Cell3D  | cut -d. -f1-3 |sort  |  xargs -I filename -n1 sh -c "h5merge -i case2/tmp/hdf5/filename.001.h5 -o case2/tmp/hdf5/filename.000.h5 pbcPosition 2>/dev/null > /dev/null || cp case2/tmp/hdf5/filename.001.h5 case2/tmp/hdf5/filename.000.h5"
echo "=====Checking Output=========="
difference=$(find ./case1/ -exec basename {} \; |grep h5 | grep RBC | grep -v _Cell | xargs -I filename -n1 h5diff -p 0.001 -c case1/tmp/hdf5/filename case2/tmp/hdf5/filename pbcPosition|grep -v "Not comparable")
difference1=$(find ./case1/ -exec basename {} \; |grep h5 | grep RBC | grep -v _Cell | xargs -I filename -n1 h5diff -p 0.001 -c case1/tmp/hdf5/filename case2/tmp/hdf5/filename force|grep -v "Not comparable")
difference2=$(find ./case1/ -exec basename {} \; |grep h5 | grep RBC | grep -v _Cell | xargs -I filename -n1 h5diff -p 0.001 -c case1/tmp/hdf5/filename case2/tmp/hdf5/filename velocity|grep -v "Not comparable")
difference="$difference$difference1$difference2"
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
