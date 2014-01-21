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

#ifndef REDUCTION_PARTICLE_3D_H
#define REDUCTION_PARTICLE_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <map>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ReductionParticle3D : public Particle3D<T,Descriptor> {
public:
    ReductionParticle3D();
//    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.) { }
//    virtual void rhoBarJtoParticle(NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling=1.) { }
//    virtual void fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling=1.) { }
    /// Implements Euler integration with velocity alone.
    virtual void advance();
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    virtual ReductionParticle3D<T,Descriptor>* clone() const;
    /// Return the cellId through a generic interface (vector id=0).

    plint getMpiProcessor() {
    	int myrank = 0;
#ifdef PLB_MPI_PARALLEL
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    	return plint(myrank); }
private:
    static int id;
    plint processor;
    plint cellId;
    plint nParticles;
    std::map<plint, T > quantities1D; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, Array<T,3> > quantities3D;
    std::map<plint, std::vector<T> > quantitiesND;
public:
    virtual int getId() const;

    plint const& get_cellId() const { return cellId; }
    plint& get_cellId() { return cellId; }
    plint const& get_processor() const { return processor; }
    plint& get_processor() { return processor; }
    plint const& get_nParticles() const { return nParticles; }
    plint& get_nParticles() { return nParticles; }

    virtual bool getScalar(plint whichScalar, T& scalar) const;
    std::string getScalarName(plint whichScalar) const;
    plint getScalarsNumber() const ;

    virtual bool getVector(plint whichVector, Array<T,3>& vector) const;
    std::string getVectorName(plint whichVector) const;
    plint getVectorsNumber() const ;

};

namespace meta {

template<typename T, template<typename U> class Descriptor>
ParticleRegistration3D<T,Descriptor>& particleRegistration3D();


template< typename T,
          template<typename U> class Descriptor,
          class ReductionParticle >
class ReductionParticleGenerator3D : public ParticleGenerator3D<T,Descriptor>
{
    virtual Particle3D<T,Descriptor>* generate (
            HierarchicUnserializer& unserializer ) const
    {
        // tag, position, scalars, vectors.

        plint processor, cellId;
        std::map<plint, T > quantities1D; // quantities1D[CCR_VOLUME] = CELL_VOLUME
        std::map<plint, Array<T,3> > quantities3D;
        std::map<plint, std::vector<T> > quantitiesND;

        unserializer.readValue<plint>(processor);
        unserializer.readValue<plint>(cellId);

        plint size1D, size3D, sizeND, qId;
        unserializer.readValue<plint>(size1D);
        for (int id = 0; id < size1D; ++id) {
            unserializer.readValue<plint>(qId);
            unserializer.readValue<T>(quantities1D[qId]);
        }

        unserializer.readValue<plint>(size3D);
        for (int id = 0; id < size3D; ++id) {
            unserializer.readValue<plint>(qId);
            unserializer.readValue<T,3>(quantities3D[qId]);
        }

        unserializer.readValue<plint>(sizeND);
        for (int id = 0; id < sizeND; ++id) {
            unserializer.readValue<plint>(qId);
            unserializer.readValue<T,3>(quantitiesND[qId]);
        }

        return new ReductionParticle();
    }
};


template< typename T,
          template<typename U> class Descriptor,
          class ReductionParticle >
int registerReductionParticle3D(std::string name) {
    return particleRegistration3D<T,Descriptor>().announce (
               name, new ReductionParticleGenerator3D<T,Descriptor,ReductionParticle> );
}

}  // namespace meta

}  // namespace plb

#include "reductionParticle3D.hh"

#endif  // REDUCTION_PARTICLE_3D_H

