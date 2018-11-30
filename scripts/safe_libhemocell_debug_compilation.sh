#!/bin/bash
set -e
BASEDIR=$(dirname "$0")
HEMOCELL_DIR=${BASEDIR}/../
cd ${HEMOCELL_DIR}/build/hemocell_debug/

#only allow ONE process
exec 200>./lock
flock -n 200 || exit 0
reset -I ||
cmake ./

#Cmake and make are insane, use insanity to counter it!
parent=`ps -o ppid= -p $$`
parent1=`ps -o ppid= -p $parent`
parent2=`ps -o ppid= -p $parent1`
maybe_make=`ps -o command= $parent2`
echo $maybe_make
make `echo $maybe_make | cut -d" " -f2- | grep -v make`
