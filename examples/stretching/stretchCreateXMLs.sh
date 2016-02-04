#!/bin/bash
## This script creates several XML files for the validation of the elastic properties of a single RBC.
## It is suggested to keep the following folder structure and run the application inside these folders,
## redirecting the output to a .log file.
../../scripts/conficsion.py $1 --stretchForce 0 5 10 15 20 25 30 35 40 50 75 100 125 150 175 200
for i in *stretchForce*xml; do 
    mkdir ${i//.xml}; 
    mv ${i} ${i//.xml} 
done 