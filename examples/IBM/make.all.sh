#!/bin/bash

if [ "$1" == "" ]; then
export MAKEFILE=Makefile
else
export MAKEFILE=Makefile.$1
fi

for file in *.cpp; do
	echo make -f ${MAKEFILE} projectFiles=${file}
	make -f ${MAKEFILE} projectFiles=${file} || break
done


