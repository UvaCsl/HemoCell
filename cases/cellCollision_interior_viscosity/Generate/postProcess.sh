#!/bin/bash

directory="cellCollision"
directoryCopy="cellCollision1"
CORES=16


for val in {1..$CORES}
do
echo $directory$val
cd $directory$val
python csvConcatentate.py
cd ..
done

mkdir outputCSV

for val in {1..$CORES}
do
cp $directory$val/RBC_HO_data.csv outputCSV/RBC_HO_data_$val.csv
cp $directory$val/PLT_data.csv outputCSV/PLT_data_$val.csv
done
