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

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ImmersedCellParticle3D : public Particle3D<T,Descriptor> {
public:
    ImmersedCellParticle3D();
    ImmersedCellParticle3D( plint tag_, Array<T,3> const& position, plint cellId_ = -1 );
    ImmersedCellParticle3D( plint tag_, Array<T,3> const& position,
                          Array<T,3> const& v_, Array<T,3> const& pbcPosition_,
                            Array<T,3> const& a_, Array<T,3> const& force_, Array<T,3> const& vPrevious_,
                            plint cellId_ = -1);
    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.) { }
    virtual void rhoBarJtoParticle(NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling=1.) { }
    virtual void fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling=1.) { }
    /// Implements Euler integration with velocity alone.
    virtual void advance();
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    virtual int getId() const;
    virtual void reset(Array<T,3> const& position, Array<T,3> const& velocity_);
    virtual void reset(Array<T,3> const& position);
    virtual ImmersedCellParticle3D<T,Descriptor>* clone() const;
    /// Return the cellId through a generic interface (vector id=0).
    virtual bool getScalar(plint whichScalar, T& scalar) const;
    std::string getScalarName(plint whichScalar) const;
    plint getScalarsNumber() const ;
    /// Return the velocity, acceleration or pbcPosition through a generic interface (vector id=0,1,2).
    virtual bool getVector(plint whichVector, Array<T,3>& vector) const;
    std::string getVectorName(plint whichVector) const;
    plint getVectorsNumber() const ;
    Array<T,3> const& get_v() const { return v; }
    Array<T,3> const& get_pbcPosition() const { return pbcPosition; }
    Array<T,3> const& get_vPrevious() const { return vPrevious; }
    Array<T,3> const& get_a() const { return a; }
    Array<T,3> const& get_force() const { return force; }
    Array<T,3>& get_v() { return v; }
    Array<T,3>& get_pbcPosition() { return pbcPosition; }
    Array<T,3>& get_vPrevious() { return vPrevious; }
    Array<T,3>& get_a() { return a; }
    Array<T,3>& get_force() { return force; }
    plint const& get_cellId() const { return cellId; }
    plint& get_cellId() { return cellId; }
    plint const& get_processor() const { return processor; }
    plint& get_processor() { return processor; }
    plint getMpiProcessor() {
    	int myrank = 0;
#ifdef PLB_MPI_PARALLEL
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    	return plint(myrank); }
private:
    Array<T,3> v, pbcPosition, a, force, vPrevious;
    static int id;
private:
    Array<T,3> f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity;
    Array<T,3> stress;
    T E_other, E_inPlane, E_bending, E_area,  E_volume;

private:
    plint processor;
    plint cellId;
public:
    ImmersedCellParticle3D (
            plint tag_, Array<T,3> const& position,
            Array<T,3> const& v_, Array<T,3> const& pbcPosition_,
            Array<T,3> const& a_, Array<T,3> const& force_,  Array<T,3> const& vPrevious_,
            Array<T,3> const& f_wlc_, Array<T,3> const& f_bending_, Array<T,3> const& f_volume_, Array<T,3> const& f_surface_, Array<T,3> const& f_shear_, Array<T,3> const& f_viscosity_,
            Array<T,3> const& stress_,
            T const& E_other_,
            T const& E_inPlane_, T const& E_bending_,
            T const& E_area_, T const& E_volume_,
            plint processor_, plint cellId_ );

    Array<T,3> const& get_f_wlc() const { return f_wlc; }
    Array<T,3> const& get_f_bending() const { return f_bending; }
    Array<T,3> const& get_f_volume() const { return f_volume; }
    Array<T,3> const& get_f_surface() const { return f_surface; }
    Array<T,3> const& get_f_shear() const { return f_shear; }
    Array<T,3> const& get_f_viscosity() const { return f_viscosity; }
    Array<T,3> const& get_stress() const { return stress; }

    T const& get_E_other() const { return E_other; }
    T const& get_E_inPlane() const { return E_inPlane; }
    T const& get_E_bending() const { return E_bending; }
    T const& get_E_area() const { return E_area; }
    T const& get_E_volume() const { return E_volume; }

    Array<T,3>& get_f_wlc() { return f_wlc; }
    Array<T,3>& get_f_bending() { return f_bending; }
    Array<T,3>& get_f_volume() { return f_volume; }
    Array<T,3>& get_f_surface() { return f_surface; }
    Array<T,3>& get_f_shear() { return f_shear; }
    Array<T,3>& get_f_viscosity() { return f_viscosity; }
    Array<T,3>& get_stress() { return stress; }

    T& get_E_other() { return E_other; }
    T& get_E_inPlane() { return E_inPlane; }
    T& get_E_bending() { return E_bending; }
    T& get_E_area() { return E_area; }
    T& get_E_volume() { return E_volume; }

    T const get_E_total() const { return (E_other + E_inPlane + E_bending + E_area + E_volume);}

};

namespace meta {

template<typename T, template<typename U> class Descriptor>
ParticleRegistration3D<T,Descriptor>& particleRegistration3D();


template< typename T,
          template<typename U> class Descriptor,
          class ImmersedCellParticle >
class ImmersedCellParticleGenerator3D : public ParticleGenerator3D<T,Descriptor>
{
    virtual Particle3D<T,Descriptor>* generate (
            HierarchicUnserializer& unserializer ) const
    {
        // tag, position, scalars, vectors.
        plint tag;
        Array<T,3> position;
        Array<T,3> v, pbcPosition, a, force, vPrevious;
        Array<T,3> f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity, stress;
        T E_other;
        T E_inPlane, E_bending, E_area,  E_volume;
        plint processor, cellId;

        unserializer.readValue(tag);
        unserializer.readValues<T,3>(position);
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
        unserializer.readValues<T,3>(stress);

        unserializer.readValue(E_other);
        unserializer.readValue(E_inPlane);
        unserializer.readValue(E_bending);
        unserializer.readValue(E_area);
        unserializer.readValue(E_volume);

        unserializer.readValue(processor);
        unserializer.readValue(cellId);

        return new ImmersedCellParticle(tag, position, v, pbcPosition, a, force, vPrevious,
                f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity, stress,
                E_other, E_inPlane, E_bending, E_area,E_volume,
                processor, cellId);
    }
};


template< typename T,
          template<typename U> class Descriptor,
          class ImmersedCellParticle >
int registerImmersedCellParticle3D(std::string name) {
    return particleRegistration3D<T,Descriptor>().announce (
               name, new ImmersedCellParticleGenerator3D<T,Descriptor,ImmersedCellParticle> );
}

}  // namespace meta

}  // namespace plb

#include "immersedCellParticle3D.hh"

#endif  // IMMERSED_WALL_PARTICLE_3D_H

