#!/bin/bash
trap "exit" INT

echo "=========== Building =========="
date

if [ ! -d "./build" ]; then
  echo "* Running CMake..."
  mkdir build
  cd build
  cmake ..
  cd ..
fi

echo "* Compiling..."
cd build; 
script -q -c "make -j 4 2>&1 >/dev/null | grep 'Error\|error\|\*\*\*'"; 
cd ..

date
echo "=========== Done ==========="
