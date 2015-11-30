#!/bin/bash
# using regular expressions to extract data
for i in */*.log; do 
	echo $(dirname $i | sed 's#.*stretchForce-\([^_]*\)_.*#\1#g') , $(cat $i|grep Stretch|tail -n1 | sed 's#^.*Iteration:\([^;]*\).*Stretch (\(.*\))#\1, \2#g') ; 
done

