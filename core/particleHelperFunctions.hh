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




template<typename T1, typename T2>
void serializeMap(HierarchicSerializer& serializer, std::map<T1,T2> const& dict)
{
    typename std::map<T1,T2>::const_iterator iter;
    serializer.addValue<plint>(dict.size());
    for (iter  = dict.begin(); iter != dict.end(); ++iter) {
        serializer.addValue<T1>(iter->first);
        serializer.addValue<T2>(iter->second);
    }


}


template<typename T1, typename T2>
std::map<T1,T2> unserializeMap(HierarchicUnserializer& unserializer)
{
    std::map<T1,T2> dict;
    plint n;
    unserializer.readValue<plint>(n);
    for (int i = 0; i < n; ++i) {
        T1 key;
        unserializer.readValue<T1>(key);
        T2 value;
        unserializer.readValue<T2>(value);
        dict[key] = value;
    }
    return dict;
}



void serializeString(HierarchicSerializer& serializer, std::string const& s)
{
	plint n = s.length();
    serializer.addValue<plint>(n);
	for (pluint i = 0; i < s.length(); ++i)
        serializer.addValue<int>( int(s[i]) );
}


std::string unserializeString(HierarchicUnserializer& unserializer)
{
    plint n;
    int c;
    unserializer.readValue<plint>(n);
    std::string s(n, ' ');
	for (int i = 0; i < n; ++i) {
        unserializer.readValue<int>(c);
        s[i] = char(c);
	}
    return s;
}


} //namespace plb


#endif  // PARTICLE_HELPER_FUNCTIONALS_HH

