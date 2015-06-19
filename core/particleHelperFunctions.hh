#ifndef PARTICLE_HELPER_FUNCTIONALS_HH
#define PARTICLE_HELPER_FUNCTIONALS_HH

#include "particleHelperFunctions.h"

namespace plb {

template<typename T>
void serializeVector(HierarchicSerializer& serializer, std::vector<T> const& vec)
{
    plint n = vec.size();
    serializer.addValue<plint>(n);
    for (int i = 0; i < n; ++i) {
        serializer.addValue<T>(vec[i]);
    }
}


template<typename T>
std::vector<T> unserializeVector(HierarchicUnserializer& unserializer)
{
    std::vector<T> vec;
    plint n;
    unserializer.readValue<plint>(n);
    for (int i = 0; i < n; ++i) {
        T value;
        unserializer.readValue<T>(value);
        vec.push_back(value);
    }
    return vec;
}

} //namespace plb


#endif  // PARTICLE_HELPER_FUNCTIONALS_HH

