#!/bin/bash

# This script generates validation figures for the stretch cell example. A
# single RBC is subjected to an stretching force of different magnitudes for
# while the axial and transverse dimensions of the cell are tracked.
#
# The `stretchCell` example writes a brief log file containing the current
# iteration count and the axial, transverse dimensions as `stretch-*.log`. The
# script considers the following stretch values:

stretch_forces=(0 25 50 75 125 150 173 175)

# The resulting data is drawn with `gnuplot` if present on the system.
# Otherwise, a multitude of options could be used for visualisation.

executable="stretchCell"
[ -f "${executable}" ] || { echo "The executable: '${executable}' does not exist\
, please compile '${executable}' first."; exit 1; }

for force in ${stretch_forces[@]};
do
        sed -i "s#<stretchForce>.*#<stretchForce>${force}</stretchForce>#" config.xml
        ./stretchCell config.xml
done
mv stretch-*.log validation/

cd validation

echo "force iter axial transverse" > combined.txt
for force in ${stretch_forces[@]};
do
        echo -n "${force} " >> combined.txt
        tail -n1 "stretch-${force}.log" >> combined.txt
done

gnuplot validation.gpl

cd ..
