#!/bin/bash
trap "exit" INT

# This script invokes the compilation of the example present in the current
# directory.

echo "=========== Building =========="
date

example=${PWD}

if [ ! -d "../../build" ]; then
  echo "* Running CMake..."
  mkdir ../../build
  cd ../../build || exit 1
  cmake ..
  cd "$example" || exit 1
fi

echo "* Compiling..."
cd ../../build || exit 1
cmake --build . --target "${example##*/}"
cd "$example" || exit 1

date
echo "=========== Done ==========="
