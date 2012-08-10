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

#ifndef IMMERSED_WALL_PARTICLE_3D_H
#define IMMERSED_WALL_PARTICLE_3D_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "particles/particle3D.h"
#include "particles/particleIdentifiers3D.h"
#include "atomicBlock/blockLattice3D.h"
#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ImmersedWallParticle3D : public Particle3D<T,Descriptor> {
public:
    ImmersedWallParticle3D();
    ImmersedWallParticle3D( plint tag_, Array<T,3> const& position, plint cellId_ = -1 );
    ImmersedWallParticle3D( plint tag_, Array<T,3> const& position,
                          Array<T,3> const& v_, Array<T,3> const& vHalfTime_,
                            Array<T,3> const& a_, Array<T,3> const& force_, Array<T,3> const& vPrevious_,
                            plint cellId_ = -1);
    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.) { }
    virtual void fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling=1.) { }
    /// Implements "steps 1 and 2" of the Verlet algorithm: given
    ///   x(t), v(t), and a(t), it computes v(t+1/2) and x(t+1).
    virtual void advance();
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    virtual int getId() const;
    virtual void reset(Array<T,3> const& position);
    virtual ImmersedWallParticle3D<T,Descriptor>* clone() const;
    /// Return the cellId through a generic interface (vector id=0).
    virtual bool getScalar(plint whichScalar, T& scalar) const;
    /// Return the velocity, acceleration or vHalfTime through a generic interface (vector id=0,1,2).
    virtual bool getVector(plint whichVector, Array<T,3>& vector) const;
    Array<T,3> const& get_v() const { return v; }
    Array<T,3> const& get_vHalfTime() const { return vHalfTime; }
    Array<T,3> const& get_vPrevious() const { return vPrevious; }
    Array<T,3> const& get_a() const { return a; }
    Array<T,3> const& get_force() const { return force; }
    plint const& get_cellId() const { return cellId; }
    Array<T,3>& get_v() { return v; }
    Array<T,3>& get_vHalfTime() { return vHalfTime; }
    Array<T,3>& get_vPrevious() { return vPrevious; }
    Array<T,3>& get_a() { return a; }
    Array<T,3>& get_force() { return force; }
    plint& get_cellId() { return cellId; }
private:
    Array<T,3> v, vHalfTime, a, force, vPrevious;
    plint cellId;
    static int id;
};

namespace meta {

template<typename T, template<typename U> class Descriptor>
ParticleRegistration3D<T,Descriptor>& particleRegistration3D();


template< typename T,
          template<typename U> class Descriptor,
          class ImmersedWallParticle >
class ImmersedWallParticleGenerator3D : public ParticleGenerator3D<T,Descriptor>
{
    virtual Particle3D<T,Descriptor>* generate (
            HierarchicUnserializer& unserializer ) const
    {
        // tag, position, scalars, vectors.
        plint tag;
        unserializer.readValue(tag);
        Array<T,3> position;
        unserializer.readValues<T,3>(position);
        Array<T,3> v, vHalfTime, vPrevious, a, force;
        unserializer.readValues<T,3>(v);
        unserializer.readValues<T,3>(vHalfTime);
        unserializer.readValues<T,3>(a);
        unserializer.readValues<T,3>(force);
        unserializer.readValues<T,3>(vPrevious);
        plint cellId;
        unserializer.readValue(cellId);
        return new ImmersedWallParticle(tag, position, v, vHalfTime, vHalfTime, a, force, cellId);
    }
};


template< typename T,
          template<typename U> class Descriptor,
          class ImmersedWallParticle >
int registerImmersedWallParticle3D(std::string name) {
    return particleRegistration3D<T,Descriptor>().announce (
               name, new ImmersedWallParticleGenerator3D<T,Descriptor,ImmersedWallParticle> );
}

}  // namespace meta

}  // namespace plb

#endif  // IMMERSED_WALL_PARTICLE_3D_H

