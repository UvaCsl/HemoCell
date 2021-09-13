#pragma once

#include "helper/geometry.h"
#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"

class Bifurcation {
public:
  static std::tuple<unsigned, unsigned, unsigned> domain_size(unsigned resolution);

  static int geometry(plb::MultiBlockLattice3D<double, DESCRIPTOR> *&lattice,
                      unsigned resolution, double capillary_diameter);

  static plb::Array<double, 3> driving_force(double);
};
