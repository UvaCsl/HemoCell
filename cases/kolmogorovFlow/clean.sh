#!/bin/bash
echo "Warning! This will remove all files produced during the CMake build, except the binary!"
read -p "Are you sure? [y/n]" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
	echo '** Clearing palabos and hemocell compiled libraries...'
	find ../../build/. -type f ! -name 'CMakeLists.txt' -delete
	echo '** Clearing local build...'
	rm -r ./build
fi

