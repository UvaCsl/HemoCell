#!/bin/bash

# This scripts performs the setup of external dependencies of hemocell,
# specifically the setup of `parmetis`
# (http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview).

tag="4.0.3"
target="http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-${tag}.tar.gz"

# removes the old installation if still present
if [ -d "./parmetis" ]; then
	rm -r ./parmetis
fi

# retrieve the specified version
if [ ! -e "parmetis.tar.gz" ]; then
	wget -qO parmetis.tar.gz $target || {
		echo "Error downloading parmetis, exiting..."
		exit 1
	}
fi

# unpack, rename, and compile metis and parmetis locally
tar -xzf parmetis.tar.gz
mv parmetis-* ./parmetis

for lib in "parmetis/metis" "parmetis"; do
	(
		cd "$lib" || {
			echo "${lib} not present"
			exit 1
		}
		make config prefix="${PWD}/build/local"
		make
		make install
	)
done
