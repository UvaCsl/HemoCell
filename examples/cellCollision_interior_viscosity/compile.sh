#!/bin/bash
trap "exit" INT

# This script invokes the compilation of the example present in the current
# directory.
#
# Note: this example requires interior viscosity and the `INTERIOR_VISCOSITY`
# option should be enabled when compiling hemocell.

echo "=========== Building =========="
date

example=${PWD}

if [ ! -d "../../build" ]; then
  echo "* Running CMake..."
  mkdir ../../build
  cd ../../build || exit 1
  cmake .. -DINTERIOR_VISCOSITY=ON
  cd "$example" || exit 1
fi

echo "* Compiling..."
cd ../../build || exit 1
cmake --build . --target "${example##*/}"
cd "$example" || exit 1

date
echo "=========== Done ==========="
