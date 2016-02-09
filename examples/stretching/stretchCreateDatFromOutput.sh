#!/bin/bash
# Using regular expressions to extract data, provided that the standard output is forwarded to a log.
for i in */*.log; do 
	echo $(dirname $i | sed 's#.*stretchForce-\([^_]*\)_.*#\1#g') , $(cat $i|grep Stretch|tail -n1 | sed 's#^.*Iteration:\([^;]*\).*Stretch (\(.*\))#\1, \2#g') ; 
done

