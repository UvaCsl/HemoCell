#ifndef BOND_PARTICLE_3D_H
#define BOND_PARTICLE_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "surfaceParticle3D.h"
#include "particleHelperFunctions.h"
#include <string>
#include <map>

using namespace plb;

namespace trombocit {



template<typename T, template<typename U> class Descriptor>
class BondParticle3D : public Particle3D<T,Descriptor> {
public:
    BondParticle3D() :
        Particle3D<T,Descriptor>(),
        processor(getMpiProcessor()),
              r(0), eij(Array<T,3>(0,0,0)), bondTime(0)
        {
            for (int var = 0; var < 2; ++var) {
                particles[var] = NULL;
                positions[var] = Array<T,3>(0,0,0);
                processors[var] = getMpiProcessor();
                particles[var] = NULL;
                cellId[var] = 0;
                vertexId[var] = 0;
            }
            this->setTag(-1);
        } ;


    BondParticle3D(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r_, Array<T,3> eij_, plint tag_, std::string uid_) :
        processor(getMpiProcessor()), bondTime(0), uid(uid_)
    {
        positions[0] = positions[1] = Array<T,3>(0,0,0);
        this->setTag(tag_);
        update(p0, p1, r_, eij_);
    };


    virtual ~BondParticle3D() { };

    BondParticle3D(BondParticle3D<T,Descriptor> const& rhs)
    : Particle3D<T,Descriptor>(rhs.getTag(), rhs.getPosition()),
      processor(rhs.processor),
      r(rhs.r), eij(rhs.eij), bondTime(rhs.bondTime), uid(rhs.uid)
      {
        for (int var = 0; var < 2; ++var) {
            particles[var] = rhs.particles[var];
            positions[var] = rhs.positions[var];
            processors[var] = rhs.processors[var];
            cellId[var] = rhs.cellId[var];
            vertexId[var] = rhs.vertexId[var];
        }
    };

    virtual BondParticle3D<T,Descriptor>* clone() const;

    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.) { }
    virtual void rhoBarJtoParticle(NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling=1.) { }
    virtual void fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling=1.) { }
    virtual void advance();
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);


    virtual int getId() const { return id; } ;
    plint const& get_processor() const { return processor; }
    plint& get_processor() { return processor; }
    plint getMpiProcessor() {
        int myrank = 0;
#ifdef PLB_MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
        return plint(myrank);
    }
    bool isInTheSameDomain(plint mpiProcessor) {
        return (mpiProcessor == getMpiProcessor());
    }
    bool isInTheSameDomain(Particle3D<T,Descriptor> * p0) {
        if (p0==NULL) { return false; }
        return (castParticleToICP3D(p0)->get_processor() == getMpiProcessor());
    }


public:
    void updateFromBondParticle(Particle3D<T,Descriptor> * p0) {
    	BondParticle3D<T,Descriptor> * bp = dynamic_cast<BondParticle3D<T,Descriptor>*>(p0);
    	this->update(bp->particles[0], bp->particles[1], bp->r, bp->eij);
    }

    void update(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r_, Array<T,3> eij_) {
    	dx = Array<T,3>(0,0,0);
    	if (p0 != NULL) {
    		particles[0] = p0;
            positions[0] = p0->getPosition();
            processors[0] = castParticleToICP3D(p0)->getMpiProcessor();
            cellId[0] = castParticleToICP3D(p0)->get_cellId();
            vertexId[0] = castParticleToICP3D(p0)->getVertexId();
            dx += castParticleToICP3D(p0)->get_v() * 0.5;

    	}
    	if (p1 != NULL) {
    		particles[1] = p1;
            positions[1] = p1->getPosition();
            processors[1] = castParticleToICP3D(p1)->getMpiProcessor();
            cellId[1] = castParticleToICP3D(p1)->get_cellId();
            vertexId[1] = castParticleToICP3D(p1)->getVertexId();
            dx += castParticleToICP3D(p1)->get_v()  * 0.5 ;
    	}
        eij = eij_; //positions[0] - positions[1];
        r = r_; // norm(eij); eij = eij * (1.0/r);
        this->getPosition() = (positions[0] + positions[1]) * 0.5;
    }

    bool applyForce(Array<T,3> force) {
        if (particles[0] != NULL) { castParticleToICP3D(particles[0])->get_force() -= force; };
        if (particles[1] != NULL) { castParticleToICP3D(particles[1])->get_force() += force; };
        return true;
    };

    std::pair<plint,plint> getCellIdVertexIdPair(plint pid) {
    	PLB_ASSERT(pid<2);
        return std::make_pair<plint,plint>(cellId[pid], vertexId[pid]);
    }

    T & get_r() { return r; }
    Array<T,3> & get_eij() { return eij; }
    Array<T,3> get_rVector() { return r*eij; }

    Array<T,3> getPositions(plint pid) { PLB_ASSERT(pid<2); return positions[pid]; };
//    Array<T,3> getVelocities(plint pid) { PLB_ASSERT(pid<2); return velocities[pid]; };
    plint getProcessors(plint pid) { PLB_ASSERT(pid<2); return processors[pid]; };
    plint getCellIds(plint pid) { PLB_ASSERT(pid<2); return cellId[pid]; };
    int getBondTime() { return bondTime; }

    Particle3D<T,Descriptor>* getParticle(plint pid) {
    	PLB_ASSERT(pid<2);
    	return particles[pid];
    } ;

    std::string getUID() { return uid; }
private:
    static int id;
    plint processor;
private:
    T r, bondTime;
    plint cellId[2];
    plint vertexId[2];
    Array<T,3> dx;
    Array<T,3> eij;
    Array<T,3> positions[2];
//    Array<T,3> velocities[2];
    std::string uid;
    plint processors[2];
    Particle3D<T,Descriptor>* particles[2];
};

template<typename T, template<typename U> class Descriptor>
int BondParticle3D<T,Descriptor>::id = meta::registerGenericParticle3D<T,Descriptor,BondParticle3D<T,Descriptor> >("BondParticle3D");

template<typename T, template<typename U> class Descriptor>
BondParticle3D<T,Descriptor>* castParticle3DToBondParticle3D(Particle3D<T,Descriptor>* particle) {
    return  dynamic_cast<BondParticle3D<T,Descriptor>*> (particle);
}


#include "BondParticle3D.hh"


} // namespace trombocit


#endif  // BOND_PARTICLE_3D_H
