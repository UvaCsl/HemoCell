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

#ifndef IMMERSED_WALL_PARTICLE_3D_HH
#define IMMERSED_WALL_PARTICLE_3D_HH

#include "core/globalDefs.h"
#include "immersedCellParticle3D.h"

namespace plb {

/* *************** class ImmersedCellParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int ImmersedCellParticle3D<T,Descriptor>::id =
        meta::registerImmersedCellParticle3D<T,Descriptor,ImmersedCellParticle3D<T,Descriptor> >("ImmersedCell");

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D()
    : v(T(),T(),T()),
      vHalfTime(T(),T(),T()),
      a(T(),T(),T()), force(T(),T(),T()), vPrevious(T(),T(),T()), cellId(-1)
{ }

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D (
        plint tag_, Array<T,3> const& position, plint cellId_ )
    : Particle3D<T,Descriptor>(tag_, position), 
      v(T(),T(),T()),
      vHalfTime(T(),T(),T()),
      a(T(),T(),T()),
      force(T(),T(),T()),
      vPrevious(T(),T(),T()),
      cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D (
        plint tag_, Array<T,3> const& position,
        Array<T,3> const& v_, Array<T,3> const& vHalfTime_,
        Array<T,3> const& a_, Array<T,3> const& force_,  Array<T,3> const& vPrevious_, plint cellId_ )
    : Particle3D<T,Descriptor>(tag_, position),
      v(v_),
      vHalfTime(vHalfTime_),
      a(a_),
      force(force_),
      vPrevious(vPrevious_),
      cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::advance() {
// Velocity Verlet
//    vHalfTime = v + (T)0.5*a;
//    this->getPosition() += vHalfTime;
// Adams-Bashforth update scheme
    this->getPosition() += 1.5*v - 0.5*vPrevious;
    vPrevious = v;
// Euler update scheme
//    this->getPosition() += v;
}

template<typename T, template<typename U> class Descriptor>
int ImmersedCellParticle3D<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::reset(Array<T,3> const& position_)
{
    Particle3D<T,Descriptor>::reset(position_);
    v.resetToZero();
    vHalfTime.resetToZero();
    a.resetToZero();
    force.resetToZero();
    vPrevious.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<T,Descriptor>::serialize(serializer);
    serializer.addValues<T,3>(v);
    serializer.addValues<T,3>(vHalfTime);
    serializer.addValues<T,3>(a);
    serializer.addValues<T,3>(force);
    serializer.addValues<T,3>(vPrevious);
    serializer.addValue<plint>(cellId);
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValues<T,3>(v);
    unserializer.readValues<T,3>(vHalfTime);
    unserializer.readValues<T,3>(a);
    unserializer.readValues<T,3>(force);
    unserializer.readValues<T,3>(vPrevious);
    unserializer.readValue<plint>(cellId);
}

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>* ImmersedCellParticle3D<T,Descriptor>::clone() const {
    return new ImmersedCellParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
bool ImmersedCellParticle3D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    if (whichScalar==0) {
        scalar = get_cellId();
        return true;
    }
    return Particle3D<T,Descriptor>::getScalar(whichScalar, scalar);
}

template<typename T, template<typename U> class Descriptor>
bool ImmersedCellParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    if (whichVector==0) {
        vector = get_v();
        return true;
    } else if (whichVector==1) {
        vector = get_a();
        return true;
    } else if (whichVector==2) {
        vector = get_vHalfTime();
        return true;
    } else if (whichVector==3) {
        vector = get_force();
        return true;
    } else if (whichVector==4) {
        vector = get_vPrevious();
        return true;
    }
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}

}  // namespace plb

#endif  // IMMERSED_WALL_PARTICLE_3D_HH
