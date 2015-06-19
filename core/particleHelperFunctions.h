#ifndef PARTICLE_HELPER_FUNCTIONALS_H
#define PARTICLE_HELPER_FUNCTIONALS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <map>
#include <algorithm>

using namespace plb;
using namespace std;

namespace plb {

template<typename T>
void serializeVector(HierarchicSerializer& serializer, std::vector<T> const& vec);


template<typename T>
std::vector<T> unserializeVector(HierarchicUnserializer& unserializer);


} //namespace plb


#include "particleHelperFunctions.hh"

#endif  // PARTICLE_HELPER_FUNCTIONALS_H

