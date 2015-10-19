#ifndef BOND_FIELD_3D_H
#define BOND_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellCellForces3D.h"
#include "BondFunctionals3D.h"
#include <map>
#include <set>
#include <string>
#include <sstream>

using namespace plb;


namespace trombocit {


template<typename T, template<typename U> class Descriptor>
class BondField3D {
public:
	// Single CellField
    BondField3D(CellField3D<T, Descriptor> & cellField, BondType<T,Descriptor> & bondType_) : bondType(bondType_) {
        MultiBlockManagement3D const& particleManagement(cellField.getParticleField3D().getMultiBlockManagement());
        BondParticles3D = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
                particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
        BondParticles3D->periodicity().toggleAll(true);
        BondParticles3D->toggleInternalStatistics(false);

        particleParticleBondArg.push_back( &(cellField.getParticleField3D()) );
        particleParticleBondArg.push_back( &(cellField.getParticleField3D()) );
        particleParticleBondArg.push_back( BondParticles3D );
    } ;

    // Two CellFields
    BondField3D(CellField3D<T, Descriptor> & cellField1, CellField3D<T, Descriptor> & cellField2, BondType<T,Descriptor> & bondType_) : bondType(bondType_) {
        MultiBlockManagement3D const& particleManagement(cellField1.getParticleField3D().getMultiBlockManagement());
        BondParticles3D = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
                particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
        BondParticles3D->periodicity().toggleAll(true);
        BondParticles3D->toggleInternalStatistics(false);

        particleParticleBondArg.push_back( &(cellField1.getParticleField3D()) );
        particleParticleBondArg.push_back( &(cellField2.getParticleField3D()) );
        particleParticleBondArg.push_back( BondParticles3D );
    } ;

    // CellField and ParticleField
    BondField3D(CellField3D<T, Descriptor> & cellField1, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField2,
    		BondType<T,Descriptor> & bondType_) : bondType(bondType_) {
        MultiBlockManagement3D const& particleManagement(cellField1.getParticleField3D().getMultiBlockManagement());
        BondParticles3D = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
                particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
        BondParticles3D->periodicity().toggleAll(true);
        BondParticles3D->toggleInternalStatistics(false);

        particleParticleBondArg.push_back( &(cellField1.getParticleField3D()) );
        particleParticleBondArg.push_back( &particleField2 );
        particleParticleBondArg.push_back( BondParticles3D );
    } ;

    // Two ParticleFields
    BondField3D(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField1, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField2,
    		BondType<T,Descriptor> & bondType_) : bondType(bondType_) {
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
          particleParticleBondArg(rhs.particleParticleBondArg),
          bondType(rhs.bondType){} ;
    virtual ~BondField3D() { delete BondParticles3D; } ;


public:
    // Bond doesn't exist, particle is not saturated and for SameCellFields, cellId is not the same.
    virtual bool isBondPossible(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
    	std::string uid = bondType.getUID(p0, p1) ;
    	return (bondUIDs.count(uid) == 0) and bondType.isBondPossible(p0, p1, r, eij);
    }
    // Insert BondParticle, update data structures
    virtual BondParticle3D<T,Descriptor> * createBondParticle(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
    	bool created = bondType.createBond(p0, p1, r, eij);
    	PLB_ASSERT(created);
		std::string uid = bondType.getUID(p0, p1) ;
		BondParticle3D<T,Descriptor> * bp = new BondParticle3D<T,Descriptor>(p0, p1, r, eij, uid);
		return bp;
    }
public:
    std::vector<MultiBlock3D*> & getParticleParticleBondArg()  { return particleParticleBondArg; }
    Box3D getBoundingBox()  { return BondParticles3D->getBoundingBox(); }
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & getBondParticles3D()  { return *BondParticles3D; }
    std::set<std::string> & getBondUIDs() { return bondUIDs; } ;
    BondType<T,Descriptor> & getBondType() { return bondType; } ;

private:
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> >* BondParticles3D;
    std::vector<MultiBlock3D*> particleParticleBondArg;
    BondType<T,Descriptor> & bondType;
    std::set<std::string> bondUIDs;
};



// This function object defines the force between two LSPs of different CellField3D, once their LSPs are in proximity.
// It it to be used as an argument to ApplyProximityDynamics3D
template<typename T, template<typename U> class Descriptor>
class BondProximityDynamics3D : public ProximityDynamics3D<T,Descriptor>  {
public:
    BondProximityDynamics3D (BondField3D<T, Descriptor> & bondField_) : bondField(bondField_), bondParticleField(0) { };
    BondProximityDynamics3D (BondProximityDynamics3D<T,Descriptor> const& rhs)  : bondField(rhs.bondField), bondParticleField(rhs.bondParticleField)  { };
    virtual ~BondProximityDynamics3D () {};
    virtual bool operator()(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
        bool conditionsMet = bondField.isBondPossible(p0, p1, r, eij);
        if (conditionsMet) {
        	BondParticle3D<T,Descriptor> * bp = bondField.createBondParticle(p0, p1, r, eij);
        	bondParticleField->addParticle(bondParticleField->getBoundingBox(), bp);
        }
        return conditionsMet;
    }

    virtual void open(Box3D domain, std::vector<AtomicBlock3D*> fields) {
        bondParticleField = dynamic_cast<ParticleField3D<T,Descriptor>*>(fields[2]);
        std::vector<Particle3D<T,Descriptor>*> bondParticles;
        bondParticleField->findParticles(bondParticleField->getBoundingBox(), bondParticles);
        bondField.getBondUIDs().clear();
        for (pluint iParticle=0; iParticle<bondParticles.size(); ++iParticle) {
            std::string uid = castParticle3DToBondParticle3D(bondParticles[iParticle])->getUID() ;
            bondField.getBondUIDs().insert(uid);
        }
    };

    virtual void close(Box3D domain, std::vector<AtomicBlock3D*> fields) {
    };
    // conditionsAreMet is not necessary here.
    virtual int conditionsAreMet(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) { return bondField.isBondPossible(p0, p1, r, eij); }
private:
    BondField3D<T, Descriptor> & bondField;
    ParticleField3D<T,Descriptor> * bondParticleField;
};



template<typename T, template<typename U> class Descriptor>
class BondFieldWrapper3D {
public:
	BondFieldWrapper3D(CellField3D<T, Descriptor> & cellField, BondType<T,Descriptor> & bondType_) :
		bondField(cellField, bondType_) {    } ;

    // Two CellFields
	BondFieldWrapper3D(CellField3D<T, Descriptor> & cellField1, CellField3D<T, Descriptor> & cellField2, BondType<T,Descriptor> & bondType_) :
		bondField(cellField1, cellField2, bondType_) {    } ;

    // CellField and ParticleField
	BondFieldWrapper3D(CellField3D<T, Descriptor> & cellField1, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField2,
    		BondType<T,Descriptor> & bondType_) :
		bondField(cellField1, particleField2, bondType_) {    } ;

    // Two ParticleFields
	BondFieldWrapper3D(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField1, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField2,
    		BondType<T,Descriptor> & bondType_) :
		bondField(particleField1, particleField2, bondType_) {    } ;


	BondFieldWrapper3D(BondField3D<T, Descriptor> & bondField_)
        : bondField(bondField_) {} ;

	BondFieldWrapper3D(BondFieldWrapper3D<T, Descriptor> & rhs)
        : bondField(rhs.bondField) {} ;

    virtual ~BondFieldWrapper3D() { } ;

    void update() {
        BondProximityDynamics3D<T, Descriptor> bpd(bondField);
        T breakDistance = bondField.getBondType().getBreakDistance();
        applyProcessingFunctional (
            new ApplyProximityDynamics3D<T,Descriptor>(bpd, breakDistance),
            bondField.getBoundingBox(), bondField.getParticleParticleBondArg() );
        applyProcessingFunctional (
            new UpdateBondParticles3D<T,Descriptor>(bondField),
            bondField.getBoundingBox(), bondField.getParticleParticleBondArg() );
    }
private:
    BondField3D<T, Descriptor> & bondField;
};

} // namespace trombocit


#include "BondField3D.hh"


#endif  // BOND_FIELD_3D_HH
