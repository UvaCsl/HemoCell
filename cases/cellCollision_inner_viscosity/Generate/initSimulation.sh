#!/bin/bash

directory="cellCollision"
CORES=32

for ((i = 1; i <= $CORES; i++))
do
echo "Creating $directory$i"
cp -r $directory $directory$i
done

echo "Change the orientations of RBC and PLT"

python generateRandom.py




