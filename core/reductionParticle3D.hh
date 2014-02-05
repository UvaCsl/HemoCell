/* This file is part of the Palabos library.
 * Copyright (C) 2009, 2010 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef REDUCTION_PARTICLE_3D_HH
#define REDUCTION_PARTICLE_3D_HH

#include "reductionParticle3D.h"

namespace plb {

/* *************** class ReductionParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
ReductionParticle3D<T,Descriptor>::ReductionParticle3D()
    : Particle3D<T,Descriptor>(), cellId(-1), processor(-1)
{ }

template<typename T, template<typename U> class Descriptor>
ReductionParticle3D<T,Descriptor>::ReductionParticle3D(plint tag_, Array<T,3> const& position)
    :  Particle3D<T,Descriptor>(tag_, position), cellId(tag_), processor(this->getMpiProcessor())
{ }


template<typename T, template<typename U> class Descriptor>
ReductionParticle3D<T,Descriptor>::ReductionParticle3D(plint tag_, Array<T,3> const& position,
        plint cellId_, plint processor_, plint nParticles_,
        std::map<plint, T > const& quantities1D_,
        std::map<plint, Array<T,3> > const& quantities3D_,
        std::map<plint, std::vector<T> > const& quantitiesND_)
    :  Particle3D<T,Descriptor>(tag_, position), cellId(cellId_), processor(processor_), nParticles(nParticles_),
       quantities1D(quantities1D_), quantities3D(quantities3D_),quantitiesND(quantitiesND_)
{ }


template<typename T, template<typename U> class Descriptor>
void ReductionParticle3D<T,Descriptor>::advance() {
//    this->getPosition() = quantities3D[13]; // 13 = CCR_NO_PBC_POSITION_MEAN, but it is unnecessary, since these particles are temporary
}

template<typename T, template<typename U> class Descriptor>
int ReductionParticle3D<T,Descriptor>::getId() const {
    return id;
}

//template<typename T1, typename T2>
//class MapIterator {
//public:
//    MapIterator(map<T1, T2> const& mymap_)
//        : mymap(mymap_) { myiterator = mymap.begin(); };
//
//private:
//    map<T1, T2> const& mymap;
//    map<T1, T2>::iterator myiterator;
//};
//std::pair<T1, T2> iterateMap(map<T1, T2> mymap) {
//
//}

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



template<typename T, template<typename U> class Descriptor>
void ReductionParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
//    serializer.addValues<T,3>(force);
//    serializer.addValue<T>(E_repulsive);
    Particle3D<T,Descriptor>::serialize(serializer);
    serializer.addValue<plint>(cellId);
    serializer.addValue<plint>(processor);
    serializer.addValue<plint>(nParticles);

    typename std::map<plint, T >::const_iterator iter1D;
    typename std::map<plint, Array<T,3> >::const_iterator iter3D;
    typename std::map<plint, std::vector<T> >::const_iterator iterND;


    serializer.addValue<plint>(quantities1D.size());
    for (iter1D  = quantities1D.begin(); iter1D != quantities1D.end(); ++iter1D) {
        serializer.addValue<plint>(iter1D->first);
        serializer.addValue<T>(iter1D->second);
    }
    serializer.addValue<plint>(quantities3D.size());
    for (iter3D  = quantities3D.begin(); iter3D != quantities3D.end(); ++iter3D) {
        serializer.addValue<plint>(iter3D->first);
        serializer.addValues<T,3>(iter3D->second);
    }
    serializer.addValue<plint>(quantitiesND.size());
    for (iterND  = quantitiesND.begin(); iterND != quantitiesND.end(); ++iterND) {
        serializer.addValue<plint>(iterND->first);
        serializeVector<T>(serializer, iterND->second);
    }

}

template<typename T, template<typename U> class Descriptor>
void ReductionParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
//    unserializer.readValues<T,3>(force);
//    unserializer.readValue<T>(E_repulsive);
    quantities1D.clear();      quantities3D.clear();     quantitiesND.clear();

    Particle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue<plint>(cellId);
    unserializer.readValue<plint>(processor);
    unserializer.readValue<plint>(nParticles);


    plint size1D, size3D, sizeND;
    unserializer.readValue<plint>(size1D);
    for (int id = 0; id < size1D; ++id) {
        plint qId;
        unserializer.readValue<plint>(qId);
        unserializer.readValue<T>(quantities1D[qId]);
    }

    unserializer.readValue<plint>(size3D);
    for (int id = 0; id < size3D; ++id) {
        plint qId;
        unserializer.readValue<plint>(qId);
        unserializer.readValues<T,3>( (quantities3D[qId]) );
    }

    unserializer.readValue<plint>(sizeND);
    for (int id = 0; id < sizeND; ++id) {
        plint qId;
        unserializer.readValue<plint>(qId);
        quantitiesND[qId] = unserializeVector<T>(unserializer);
    }
}


template<typename T, template<typename U> class Descriptor>
ReductionParticle3D<T,Descriptor>* ReductionParticle3D<T,Descriptor>::clone() const {
    return new ReductionParticle3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
bool ReductionParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}

template<typename T, template<typename U> class Descriptor>
std::string ReductionParticle3D<T,Descriptor>::getVectorName(plint whichVector) const {
    return "empty";
}

template<typename T, template<typename U> class Descriptor>
plint ReductionParticle3D<T,Descriptor>::getVectorsNumber() const {
        return 0;
}

/* Same for scalars */
template<typename T, template<typename U> class Descriptor>
bool ReductionParticle3D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    return Particle3D<T,Descriptor>::getScalar(whichScalar, scalar);
}


template<typename T, template<typename U> class Descriptor>
std::string ReductionParticle3D<T,Descriptor>::getScalarName(plint whichScalar) const {
    return "empty";
}


template<typename T, template<typename U> class Descriptor>
plint ReductionParticle3D<T,Descriptor>::getScalarsNumber() const {
        return 0;
}

}  // namespace plb

#endif  // REDUCTION_PARTICLE_3D_HH
