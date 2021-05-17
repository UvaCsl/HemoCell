#!/bin/bash
trap "exit" INT

# This script invokes the compilation of the example present in the current
# directory.

example=${PWD}

if [ ! -d "../../build" ]; then
	(
		mkdir ../../build
		cd ../../build || {
			echo "The build directory should exist"
			exit 1
		}
		cmake ..
	)
fi

(
	cd ../../build || {
		echo "The build directory should exist"
		exit 1
	}
	cmake --build . --target "${example##*/}" --parallel "$(nproc)"
)
