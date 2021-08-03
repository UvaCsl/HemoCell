#pragma once

#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"

class Wedge {
public:
  static std::tuple<unsigned, unsigned, unsigned> domain_size(unsigned resolution);
  static int geometry(plb::MultiBlockLattice3D<double, DESCRIPTOR> *&lattice, unsigned resolution);
  static plb::Array<double, 3> driving_force(double);
};
