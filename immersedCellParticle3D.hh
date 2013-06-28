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

#include "palabos3D.h"
#include "palabos3D.hh"
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
      a(T(),T(),T()), force(T(),T(),T()), vPrevious(T(),T(),T()),
      f_wlc(T(),T(),T()), f_bending(T(),T(),T()), f_volume(T(),T(),T()), f_surface(T(),T(),T()), f_shear(T(),T(),T()), f_viscosity(T(),T(),T()),
      stress(T(),T(),T()), E_bending(T(),T(),T()),
      cellId(-1)
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
      f_wlc(T(),T(),T()), f_bending(T(),T(),T()), f_volume(T(),T(),T()), f_surface(T(),T(),T()), f_shear(T(),T(),T()), f_viscosity(T(),T(),T()),
      stress(T(),T(),T()), E_bending(T(),T(),T()),
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
      f_wlc(T(),T(),T()), f_bending(T(),T(),T()), f_volume(T(),T(),T()), f_surface(T(),T(),T()), f_shear(T(),T(),T()), f_viscosity(T(),T(),T()),
      stress(T(),T(),T()), E_bending(T(),T(),T()),
      cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D (
        plint tag_, Array<T,3> const& position,
        Array<T,3> const& v_, Array<T,3> const& vHalfTime_,
        Array<T,3> const& a_, Array<T,3> const& force_,  Array<T,3> const& vPrevious_,
        Array<T,3> const& f_wlc_, Array<T,3> const& f_bending_, Array<T,3> const& f_volume_, Array<T,3> const& f_surface_, Array<T,3> const& f_shear_, Array<T,3> const& f_viscosity_,
        Array<T,3> const& stress_, Array<T,3> const& E_bending_,
        plint cellId_ )
    : Particle3D<T,Descriptor>(tag_, position),
      v(v_),
      vHalfTime(vHalfTime_),
      a(a_),
      force(force_),
      vPrevious(vPrevious_),
      f_wlc(f_wlc_), f_bending(f_bending_), f_volume(f_volume_), f_surface(f_surface_), f_shear(f_shear_), f_viscosity(f_viscosity_),
      stress(stress_), E_bending(E_bending_),
      cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::advance() {
// No fluid interaction
//    v += force;
//    this->getPosition() += v + 0.5*force;
// Velocity Verlet
//    vHalfTime = v + (T)0.5*force;
//    this->getPosition() += vHalfTime;
// Adams-Bashforth update scheme
//    this->getPosition() += 1.5*v - 0.5*vPrevious;
//    vPrevious = v;
// Euler update scheme
    this->getPosition() += v;
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

    f_wlc.resetToZero();
    f_bending.resetToZero();
    f_volume.resetToZero();
    f_surface.resetToZero();
    f_shear.resetToZero();
    f_viscosity.resetToZero();
    stress.resetToZero();
    E_bending.resetToZero();

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

    serializer.addValues<T,3>(f_wlc);
    serializer.addValues<T,3>(f_bending);
    serializer.addValues<T,3>(f_volume);
    serializer.addValues<T,3>(f_surface);
    serializer.addValues<T,3>(f_shear);
    serializer.addValues<T,3>(f_viscosity);
    serializer.addValues<T,3>(stress);
    serializer.addValues<T,3>(E_bending);

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

    unserializer.readValues<T,3>(f_wlc);
    unserializer.readValues<T,3>(f_bending);
    unserializer.readValues<T,3>(f_volume);
    unserializer.readValues<T,3>(f_surface);
    unserializer.readValues<T,3>(f_shear);
    unserializer.readValues<T,3>(f_viscosity);
    unserializer.readValues<T,3>(stress);
    unserializer.readValues<T,3>(E_bending);

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
        vector = get_vHalfTime();
        return true;
    } else if (whichVector==2) {
        vector = get_a();
        return true;
    } else if (whichVector==3) {
        vector = get_force();
        return true;
    } else if (whichVector==4) {
        vector = get_vPrevious();
        return true;
    } else if (whichVector==5) {
        vector = get_f_wlc();
        return true;
    } else if (whichVector==6) {
        vector = get_f_bending();
        return true;
    } else if (whichVector==7) {
        vector = get_f_volume();
        return true;
    } else if (whichVector==8) {
        vector = get_f_surface();
        return true;
    } else if (whichVector==9) {
        vector = get_f_shear();
        return true;
    } else if (whichVector==10) {
        vector = get_f_viscosity();
        return true;
    } else if (whichVector==11) {
        vector = get_stress();
        return true;
    } else if (whichVector==12) {
        vector = get_E_bending();
        return true;
    }
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}

template<typename T, template<typename U> class Descriptor>
std::string ImmersedCellParticle3D<T,Descriptor>::getVectorName(plint whichVector) const {
    if (whichVector==0) {
        return "velocity";
    } else if (whichVector==1) {
        return "vHalfTime";
    } else if (whichVector==2) {
        return "acceleration";
    } else if (whichVector==3) {
        return "force";
    } else if (whichVector==4) {
        return "vPrevious";
    } else if (whichVector==5) {
        return "f_wlc";
    } else if (whichVector==6) {
        return "f_bending";
    } else if (whichVector==7) {
        return "f_volume";
    } else if (whichVector==8) {
        return "f_surface";
    } else if (whichVector==9) {
        return "f_shear";
    } else if (whichVector==10) {
        return "f_viscosity";
    } else if (whichVector==11) {
        return "stress";
    } else if (whichVector==12) {
        return "E_bending";
    }
    return "empty";
}

template<typename T, template<typename U> class Descriptor>
plint ImmersedCellParticle3D<T,Descriptor>::getVectorsNumber() const {
        return 13;
}

}  // namespace plb

#endif  // IMMERSED_WALL_PARTICLE_3D_HH
