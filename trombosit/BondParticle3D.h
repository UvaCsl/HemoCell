#ifndef BOND_PARTICLE_3D_H
#define BOND_PARTICLE_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include <string>

using namespace plb;

namespace trombocit {

std::string createBondId(plint cellId0, plint cellId1, plint vertexId0, plint vertexId1) {
    std::string ret;
    std::ostringstream stm0, stm1;
    stm0 << cellId0 << "-" << vertexId0;
    stm1 << cellId1 << "-" << vertexId1;
    if (cellId0 < cellId1) { ret = stm0.str() + "_" + stm1.str(); }
    else { ret = stm1.str() + "_" + stm0.str();  }
    return ret;
}


template<typename T, template<typename U> class Descriptor>
std::string createBondId(Particle3D<T,Descriptor>* particle0, Particle3D<T,Descriptor>* particle1) {
    ImmersedCellParticle3D<T,Descriptor>* p0 = castParticleToICP3D(particle0);
    ImmersedCellParticle3D<T,Descriptor>* p1 = castParticleToICP3D(particle1);

    return createBondId(p0->get_cellId(), p1->get_cellId(),  p0->getVertexId(), p1->getVertexId());
}




template<typename T, template<typename U> class Descriptor>
class BondParticle3D : public Particle3D<T,Descriptor> {
public:
    BondParticle3D() :
        processor(getMpiProcessor()), bondType(0),
              r(0), eij(Array<T,3>(0,0,0))
        {
            for (int var = 0; var < 2; ++var) {
                positions[var] = Array<T,3>(0,0,0);
                velocities[var] = Array<T,3>(0,0,0);
                processors[var] = getMpiProcessor();
                particles[var] = NULL;
            }
        } ;
    virtual ~BondParticle3D() {     };

    BondParticle3D(BondParticle3D<T,Descriptor> const& rhs)
    : processor(rhs.processor), bondType(rhs.bondType),
      r(rhs.r), eij(rhs.eij)
      {
        for (int var = 0; var < 2; ++var) {
            positions[var] = rhs.positions[var];
            velocities[var] = rhs.velocities[var];
            processors[var] = rhs.processors[var];
            particles[var] = rhs.particles[var];
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
    void update(T r, Array<T,3> eij) {
        this->r=r;
        this->eij=eij;
    }
    void update(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, bool updateDistance=true) {
        processors[0] = castParticleToICP3D(p0)->get_processor();
        processors[1] = castParticleToICP3D(p1)->get_processor();
        if ( isInTheSameDomain(p0) ) { particles[0] = p0; }
        else { particles[0] = NULL; }
        if ( isInTheSameDomain(p1) ) { particles[1] = p1; }
        else { particles[1] = NULL; }
    }
    void update(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
        this->update(p0,p1);
        this->update(r,eij);
    }
private:
    static int id;
    plint processor;
    plint bondType;

private:
    T r;
    Array<T,3> eij;
    Array<T,3> positions[2];
    Array<T,3> velocities[2];
    plint processors[2];
    std::map<plint, T> bondTypeMaps[2];
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
