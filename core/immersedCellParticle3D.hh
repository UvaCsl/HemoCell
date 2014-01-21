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

#include "immersedCellParticle3D.h"

namespace plb {

/* *************** class ImmersedCellParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int ImmersedCellParticle3D<T,Descriptor>::id =
        meta::registerImmersedCellParticle3D<T,Descriptor,ImmersedCellParticle3D<T,Descriptor> >("ImmersedCell");

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D()
    : v(T(),T(),T()),
      pbcPosition(this->getPosition()),
      a(T(),T(),T()), force(T(),T(),T()), vPrevious(T(),T(),T()),
      f_wlc(T(),T(),T()), f_bending(T(),T(),T()), f_volume(T(),T(),T()),
      f_surface(T(),T(),T()), f_shear(T(),T(),T()), f_viscosity(T(),T(),T()),
      f_repulsive(T(),T(),T()),
      stress(T(),T(),T()),
      E_other(T()),
      E_inPlane(T()), E_bending(T()), E_area(T()),  E_volume(T()),
      E_repulsive(T()),
      processor(0), cellId(-1)
{ }

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D (
        plint tag_, Array<T,3> const& position, plint cellId_ )
    : Particle3D<T,Descriptor>(tag_, position), 
      v(T(),T(),T()),
      pbcPosition(position),
      a(T(),T(),T()),
      force(T(),T(),T()),
      vPrevious(T(),T(),T()),
      f_wlc(T(),T(),T()), f_bending(T(),T(),T()), f_volume(T(),T(),T()), f_surface(T(),T(),T()), f_shear(T(),T(),T()),
      f_viscosity(T(),T(),T()), f_repulsive(T(),T(),T()),
      stress(T(),T(),T()),
      E_other(T()),
      E_inPlane(T()), E_bending(T()), E_area(T()),  E_volume(T()), E_repulsive(T()),
      processor(0), cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D (
        plint tag_, Array<T,3> const& position,
        Array<T,3> const& v_, Array<T,3> const& pbcPosition_,
        Array<T,3> const& a_, Array<T,3> const& force_,  Array<T,3> const& vPrevious_, plint cellId_ )
    : Particle3D<T,Descriptor>(tag_, position),
      v(v_),
      pbcPosition(pbcPosition_),
      a(a_),
      force(force_),
      vPrevious(vPrevious_),
      f_wlc(T(),T(),T()), f_bending(T(),T(),T()), f_volume(T(),T(),T()), f_surface(T(),T(),T()), f_shear(T(),T(),T()), f_viscosity(T(),T(),T()),
      f_repulsive(T(),T(),T()),
      stress(T(),T(),T()),
      E_other(T()),
      E_inPlane(T()), E_bending(T()), E_area(T()),  E_volume(T()), E_repulsive(T()),
      processor(0), cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>::ImmersedCellParticle3D (
        plint tag_, Array<T,3> const& position,
        Array<T,3> const& v_, Array<T,3> const& pbcPosition_,
        Array<T,3> const& a_, Array<T,3> const& force_,  Array<T,3> const& vPrevious_,
        Array<T,3> const& f_wlc_, Array<T,3> const& f_bending_, Array<T,3> const& f_volume_,
        Array<T,3> const& f_surface_, Array<T,3> const& f_shear_, Array<T,3> const& f_viscosity_,
        Array<T,3> const& f_repulsive_,
        Array<T,3> const& stress_,
        T const& E_other_,
        T const& E_inPlane_, T const& E_bending_,
        T const& E_area_, T const& E_volume_, T const& E_repulsive_,
        plint processor_, plint cellId_ )
    : Particle3D<T,Descriptor>(tag_, position),
      v(v_),
      pbcPosition(pbcPosition_),
      a(a_),
      force(force_),
      vPrevious(vPrevious_),
      f_wlc(f_wlc_), f_bending(f_bending_), f_volume(f_volume_),
      f_surface(f_surface_), f_shear(f_shear_), f_viscosity(f_viscosity_), f_repulsive(f_repulsive_),
      stress(stress_),
      E_other(E_other_),
      E_inPlane(E_inPlane_), E_bending(E_bending_),
      E_area(E_area_), E_volume(E_volume_), E_repulsive(E_repulsive_),
      processor(processor_), cellId(cellId_)
{ }

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::advance() {
// No fluid interaction
//    v += force;
//    this->getPosition() += v + 0.5*force;
// Velocity Verlet
//    pbcPosition = v + (T)0.5*force;
//    this->getPosition() += pbcPosition;
// Adams-Bashforth update scheme
//    this->getPosition() += 1.5*v - 0.5*vPrevious;
//    vPrevious = v;
// Euler update scheme
    this->getPosition() += vPrevious;
    pbcPosition += vPrevious;
    vPrevious.resetToZero();
    processor = this->getMpiProcessor();
}

template<typename T, template<typename U> class Descriptor>
int ImmersedCellParticle3D<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::reset(Array<T,3> const& position_, Array<T,3> const& velocity_, bool allVariables) {
        Particle3D<T,Descriptor>::reset(position_);
        if (allVariables) {
            pbcPosition = position_;
        }

        v = velocity_;
        vPrevious = velocity_;

        a.resetToZero();
        force.resetToZero();

        f_wlc.resetToZero();
        f_bending.resetToZero();
        f_volume.resetToZero();
        f_surface.resetToZero();
        f_shear.resetToZero();
        f_viscosity.resetToZero();
        f_repulsive.resetToZero();
        stress.resetToZero();

        E_other = T();
        E_inPlane = T();
        E_bending = T();
        E_area = T();
        E_volume = T();
        E_repulsive = T();

        processor = this->getMpiProcessor();
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::reset(Array<T,3> const& position_)
{
        reset(position_, Array<T,3>(0.,0.,0.));
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<T,Descriptor>::serialize(serializer);
    serializer.addValues<T,3>(v);
    serializer.addValues<T,3>(pbcPosition);
    serializer.addValues<T,3>(a);
    serializer.addValues<T,3>(force);
    serializer.addValues<T,3>(vPrevious);

    serializer.addValues<T,3>(f_wlc);
    serializer.addValues<T,3>(f_bending);
    serializer.addValues<T,3>(f_volume);
    serializer.addValues<T,3>(f_surface);
    serializer.addValues<T,3>(f_shear);
    serializer.addValues<T,3>(f_viscosity);
    serializer.addValues<T,3>(f_repulsive);
    serializer.addValues<T,3>(stress);

    serializer.addValue<T>(E_other);
    serializer.addValue<T>(E_inPlane);
    serializer.addValue<T>(E_bending);
    serializer.addValue<T>(E_area);
    serializer.addValue<T>(E_volume);
    serializer.addValue<T>(E_repulsive);

    serializer.addValue<plint>(processor);
    serializer.addValue<plint>(cellId);
}

template<typename T, template<typename U> class Descriptor>
void ImmersedCellParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValues<T,3>(v);
    unserializer.readValues<T,3>(pbcPosition);
    unserializer.readValues<T,3>(a);
    unserializer.readValues<T,3>(force);
    unserializer.readValues<T,3>(vPrevious);

    unserializer.readValues<T,3>(f_wlc);
    unserializer.readValues<T,3>(f_bending);
    unserializer.readValues<T,3>(f_volume);
    unserializer.readValues<T,3>(f_surface);
    unserializer.readValues<T,3>(f_shear);
    unserializer.readValues<T,3>(f_viscosity);
    unserializer.readValues<T,3>(f_repulsive);
    unserializer.readValues<T,3>(stress);

    unserializer.readValue<T>(E_other);
    unserializer.readValue<T>(E_inPlane);
    unserializer.readValue<T>(E_bending);
    unserializer.readValue<T>(E_area);
    unserializer.readValue<T>(E_volume);
    unserializer.readValue<T>(E_repulsive);

    unserializer.readValue<plint>(processor);
    unserializer.readValue<plint>(cellId);
}


template<typename T, template<typename U> class Descriptor>
ImmersedCellParticle3D<T,Descriptor>* ImmersedCellParticle3D<T,Descriptor>::clone() const {
    return new ImmersedCellParticle3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
bool ImmersedCellParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    if (whichVector==0) {
        vector = get_v();
        return true;
    } else if (whichVector==1) {
        vector = get_pbcPosition();
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
        vector = get_f_repulsive();
        return true;
    } else if (whichVector==12) {
        vector = get_stress();
        return true;
    }
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}

template<typename T, template<typename U> class Descriptor>
std::string ImmersedCellParticle3D<T,Descriptor>::getVectorName(plint whichVector) const {
    if (whichVector==0) {
        return "velocity";
    } else if (whichVector==1) {
        return "pbcPosition";
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
        return "f_repulsive";
    } else if (whichVector==12) {
        return "stress";
    }
    return "empty";
}

template<typename T, template<typename U> class Descriptor>
plint ImmersedCellParticle3D<T,Descriptor>::getVectorsNumber() const {
        return 13;
}

/* Same for scalars */
template<typename T, template<typename U> class Descriptor>
bool ImmersedCellParticle3D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    if (whichScalar==0) {
    	scalar = T(this->getTag());
        return true;
    } else if (whichScalar==1) {
        scalar = T(get_cellId());
        return true;
    } else if (whichScalar==2) {
        scalar = T(get_processor());
        return true;
    } else if (whichScalar==3) {
        scalar = get_E_total();
        return true;
    } else if (whichScalar==4) {
        scalar = T(get_E_inPlane());
        return true;
    } else if (whichScalar==5) {
        scalar = T(get_E_bending());
        return true;
    } else if (whichScalar==6) {
        scalar = T(get_E_area());
        return true;
    } else if (whichScalar==7) {
        scalar = T(get_E_volume());
        return true;
    } else if (whichScalar==8) {
        scalar = T(get_E_repulsive());
        return true;
    } else if (whichScalar==9) {
        scalar = T(get_E_other());
        return true;
    }
    return Particle3D<T,Descriptor>::getScalar(whichScalar, scalar);
}


template<typename T, template<typename U> class Descriptor>
std::string ImmersedCellParticle3D<T,Descriptor>::getScalarName(plint whichScalar) const {
    if (whichScalar==0) {
        return "tag";
    } else if (whichScalar==1) {
        return "processor";
    } else if (whichScalar==2) {
        return "processor";
    } else if (whichScalar==3) {
        return "E_total";
    } else if (whichScalar==4) {
        return "E_inPlane";
    } else if (whichScalar==5) {
        return "E_bending";
    } else if (whichScalar==6) {
        return "E_area";
    } else if (whichScalar==7) {
        return "E_volume";
    } else if (whichScalar==8) {
        return "E_repulsive";
    } else if (whichScalar==9) {
        return "E_other";
    }
    return "empty";
}


template<typename T, template<typename U> class Descriptor>
plint ImmersedCellParticle3D<T,Descriptor>::getScalarsNumber() const {
        return 10;
}

}  // namespace plb

#endif  // IMMERSED_WALL_PARTICLE_3D_HH
