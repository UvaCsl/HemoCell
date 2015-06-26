#ifndef BOND_FIELD_3D_H
#define BOND_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellCellForces3D.h"
#include <map>
#include <string>
#include <sstream>

using namespace plb;


namespace trombocit {


template<typename T, template<typename U> class Descriptor>
class BondField3D {
public:
    BondField3D(CellField3D<T, Descriptor> & cellField1, CellField3D<T, Descriptor> & cellField2) {
        MultiBlockManagement3D const& particleManagement(cellField1.getParticleField3D().getMultiBlockManagement());
        BondParticles3D = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
                particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
        BondParticles3D->periodicity().toggleAll(true);
        BondParticles3D->toggleInternalStatistics(false);

        particleParticleBondArg.push_back( &(cellField1.getParticleField3D()) );
        particleParticleBondArg.push_back( &(cellField2.getParticleField3D()) );
        particleParticleBondArg.push_back( BondParticles3D );
    } ;

    BondField3D(CellField3D<T, Descriptor> & cellField1, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField2) {
        MultiBlockManagement3D const& particleManagement(cellField1.getParticleField3D().getMultiBlockManagement());
        BondParticles3D = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
                particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
        BondParticles3D->periodicity().toggleAll(true);
        BondParticles3D->toggleInternalStatistics(false);

        particleParticleBondArg.push_back( &(cellField1.getParticleField3D()) );
        particleParticleBondArg.push_back( &particleField2 );
        particleParticleBondArg.push_back( BondParticles3D );
    } ;

    BondField3D(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField1, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField2) {
        MultiBlockManagement3D const& particleManagement(particleField1.getMultiBlockManagement());
        BondParticles3D = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
                particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
        BondParticles3D->periodicity().toggleAll(true);
        BondParticles3D->toggleInternalStatistics(false);

        particleParticleBondArg.push_back( &particleField1 );
        particleParticleBondArg.push_back( &particleField2 );
        particleParticleBondArg.push_back( BondParticles3D );
    } ;



    BondField3D(BondField3D<T, Descriptor> & rhs)
        : BondParticles3D(rhs.particleParticleBondArg),
          particleParticleBondArg(rhs.particleParticleBondArg) {} ;
    virtual ~BondField3D() { delete BondParticles3D; } ;


public:
    // Bond doesn't exist, particle is not saturated and for SameCellFields, cellId is not the same.
    virtual bool isNewBondPossible(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij)=0;
    // Insert BondParticle, update data structures
    virtual bool createBond(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij)=0;
    // Change BondParticle getTag() to -1, remove from data structures.
    virtual void breakBonds()=0;
    // Scan all bonds and apply force (or other dynamics).
    virtual bool applyBondDynamics()=0;

public:
    std::vector<MultiBlock3D*> & getParticleParticleBondArg()  { return particleParticleBondArg; }
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & getBondParticles3D()  { return *BondParticles3D; }
    std::map<plint, BondType<T,Descriptor> &> & getBondTypes() { return bondTypes; } ;

private:
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> >* BondParticles3D;
    std::vector<MultiBlock3D*> particleParticleBondArg;
    std::map<plint, BondType<T,Descriptor> &> bondTypes;
    std::map<plint, std::string > bondTypeUIDs;
};



// This function object defines the force between two LSPs of different CellField3D, once their LSPs are in proximity.
// It it to be used as an argument to ApplyProximityDynamics3D
template<typename T, template<typename U> class Descriptor>
class BondProximityDynamics3D {
public:
    BondProximityDynamics3D (BondField3D<T, Descriptor> & bondField_) : bondField(bondField_) { };
    BondProximityDynamics3D (BondProximityDynamics3D<T,Descriptor> const& rhs)  : bondField(rhs.bondField) { };
    virtual ~BondProximityDynamics3D () {};
    virtual bool operator()(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
        bool conditionsMet = bondField.isNewBondPossible(p0, p1, r, eij);
        if (conditionsMet) { bondField.createBond(p0, p1, r, eij); }
        return conditionsMet;
    }

    virtual void open(Box3D domain, std::vector<AtomicBlock3D*> fields) { };

    virtual void close(Box3D domain, std::vector<AtomicBlock3D*> fields) {
        ParticleField3D<T,Descriptor>& bondParticleField = // 3rd field is the bondparticleField.
            *dynamic_cast<ParticleField3D<T,Descriptor>*>(fields[2]);
        bondField.applyBondDynamics();
        bondField.breakBonds();
        bondParticleField.removeParticles(bondParticleField.getBoundingBox(), -1);
        bondParticleField.advanceParticles(bondParticleField.getBoundingBox(), -1);
    };
    // conditionsAreMet is not necessary here.
    virtual bool conditionsAreMet(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) { return true; }
private:
    BondField3D<T, Descriptor> & bondField;
};




template<typename T, template<typename U> class Descriptor>
class SameCellFieldBondField3D : BondField3D<T,Descriptor> {
public:
    SameCellFieldBondField3D(CellField3D<T, Descriptor> & cellField1) ;
    SameCellFieldBondField3D(SameCellFieldBondField3D<T, Descriptor> & rhs) ;
    virtual ~SameCellFieldBondField3D() { } ;

public:
    // Bond doesn't exist, particle is not saturated and for SameCellFields, cellId is not the same.
    virtual bool isNewBondPossible(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij)=0;
    // Insert BondParticle, update data structures
    virtual bool createBond(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij)=0;
    // Change BondParticle getTag() to -1, remove from data structures.
    virtual void breakBonds()=0;
    // Scan all bonds and apply force (or other dynamics).
    virtual bool applyBondDynamics()=0;
    // Get unique id for bond
    virtual std::string getUID(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1) {
        ImmersedCellParticle3D<T,Descriptor>* sp0 = castParticleToICP3D(p0);
        ImmersedCellParticle3D<T,Descriptor>* sp1 = castParticleToICP3D(p1);
        plint cellId0=sp0->get_cellId(), cellId1=sp1->get_cellId();
        plint vertexId0=sp0->getVertexId(), vertexId1=sp1->getVertexId();

        std::string ret;
        std::ostringstream stm0, stm1;
        stm0 << cellId0 << "-" << vertexId0;
        stm1 << cellId1 << "-" << vertexId1;
        if (cellId0 < cellId1) { ret = stm0.str() + "_" + stm1.str(); }
        else { ret = stm1.str() + "_" + stm0.str();  }
        return ret;
    }

private:
    std::map<std::string, Particle3D<T,Descriptor> *> BondDictionary;

};




} // namespace trombocit


#include "BondField3D.hh"


#endif  // BOND_FIELD_3D_HH
