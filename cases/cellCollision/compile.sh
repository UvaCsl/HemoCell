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
#cd build; make -j 4 2>&1 >/dev/null | grep 'Error\|error\|\*\*\*'; cd ..
cd build; make -j 12

date
echo "=========== Done ==========="
