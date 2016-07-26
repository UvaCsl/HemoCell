#!/bin/bash

getAbsFilename() {
  # $1 is the relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

patch $(getAbsFilename "./../palabos/src/offLattice/triangularSurfaceMesh.hh") ./triangularSurfaceMesh.hh.patch
patch $(getAbsFilename "./../palabos/src/particles/particleField3D.h") ./particleField3D.h.patch
patch $(getAbsFilename "./../palabos/src/particles/particleField3D.hh") ./particleField3D.hh.patch
