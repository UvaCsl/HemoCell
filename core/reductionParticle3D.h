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
#include "cellReductionTypes.h"
using namespace plb;
using namespace std;

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ReductionParticle3D : public Particle3D<T,Descriptor> {
public:
    ReductionParticle3D();
    ReductionParticle3D(plint tag_, Array<T,3> const& position);
    ReductionParticle3D(plint tag_, Array<T,3> const& position,
                plint cellId_, plint processor_, plint nParticles_,
                std::map<plint, T > const& quantities1D_,
                std::map<plint, Array<T,3> > const& quantities3D_,
                std::map<plint, std::vector<T> > const& quantitiesND_);

    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.) { }
    virtual void rhoBarJtoParticle(NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling=1.) { }
    virtual void fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling=1.) { }
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
    static int id ;
    plint cellId;
    plint processor;
    plint nParticles;
    std::map<plint, T > quantities1D; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, Array<T,3> > quantities3D;
    std::map<plint, std::vector<T> > quantitiesND;
public:
    void insert(plint ccrId, T value) { quantities1D[ccrId] = value; }
    void insert(plint ccrId, Array<T,3> value) { quantities3D[ccrId] = value; }
    void insert(plint ccrId, std::vector<T> value) { quantitiesND[ccrId] = value; }
    T const& get1D(plint ccrId) { return quantities1D[ccrId]; }
    Array<T,3> const& get3D(plint ccrId) { return quantities3D[ccrId]; }
    std::vector<T> const& getND(plint ccrId) { return quantitiesND[ccrId]; }

    std::map<plint, T >& getQuantities1D() { return quantities1D; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, Array<T,3> >& getQuantities3D() { return quantities3D; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME
    std::map<plint, std::vector<T> >& getQuantitiesND() { return quantitiesND; }; // quantities1D[CCR_VOLUME] = CELL_VOLUME

public:
    virtual int getId() const;

    plint const& get_cellId() const { return cellId; }
    plint& get_cellId() { return cellId; }
    plint const& get_processor() const { return processor; }
    plint& get_processor() { return processor; }
    plint const& get_nParticles() const { return nParticles; }
    plint& get_nParticles() { return nParticles; }

    T const& getVolume() const { return quantities1D[CCR_VOLUME]; }
    T& getVolume() { return quantities1D[CCR_VOLUME]; }

    T const& getSurface() const { return quantities1D[CCR_SURFACE]; }
    T& getSurface() { return quantities1D[CCR_SURFACE]; }

    virtual bool getScalar(plint whichScalar, T& scalar) const;
    std::string getScalarName(plint whichScalar) const;
    plint getScalarsNumber() const ;

    virtual bool getVector(plint whichVector, Array<T,3>& vector) const;
    std::string getVectorName(plint whichVector) const;
    plint getVectorsNumber() const ;

};


template<typename T, template<typename U> class Descriptor>
int ReductionParticle3D<T,Descriptor>::id = meta::registerGenericParticle3D<T,Descriptor,ReductionParticle3D<T,Descriptor> >("ReductionParticle3D");

template<typename T>
void serializeVector(HierarchicSerializer& serializer, std::vector<T> const& vec);

template<typename T>
std::vector<T> unserializeVector(HierarchicUnserializer& unserializer);


}  // namespace plb

#include "reductionParticle3D.hh"

#endif  // REDUCTION_PARTICLE_3D_H

