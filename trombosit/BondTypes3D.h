#ifndef BOND_TYPES_3D_H
#define BOND_TYPES_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellCellForces3D.h"
#include <map>
#include <string>
#include <iostream>     // std::cout, std::right, std::endl
#include <vector>

using namespace plb;


namespace trombocit {


template<typename T, template<typename U> class Descriptor>
class BondType {
public:
	// areSameCellType Means if the two edges of the bond belong to the same type --> PLT<->PLT is true, PLT<->Wall is false
	// This is used to distinguish if interactions are possible or not;
    BondType(T r_create_, T r_break_, bool areSameCellType_)
        : r_create(r_create_), r_break(r_break_), bondTypeId(globalBondTypeId++), areSameCellType(areSameCellType_)
          {} ;
    virtual ~BondType() {} ;
public:
    T & getBreakDistance() { return r_break; } ;
    T & getCreateDistance() { return r_create; } ;
    plint getBondTypeId() { return bondTypeId; }
    void setBondTypeId(plint bondTypeId_) { bondTypeId = bondTypeId_; }
    bool getAreSameCellType() { return areSameCellType; }
    void getAreSameCellType(bool aSCT) { areSameCellType = aSCT; }
    // Particle p0 and p1 should exist in the same domain (even in envelopes) for a bond to be created.
    //    Hence p0 and p1 always exist and the envelopes should be large enough (>r_create) for this assumption to hold.
    bool isBondPossible(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
    	if (r > this->getCreateDistance() ) { return false;}
    	// If they are the same type they should: (1) have different cellId (2) be considered only once, hence only when cId1>cId2;
    	if (areSameCellType and (castParticleToICP3D(p0)->get_cellId() <= castParticleToICP3D(p1)->get_cellId())) { return false; }
    	return true;
    }
    // Create bond can also perform other function in bondParticle, like change saturation etc.
    virtual bool createBond(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
        return isBondPossible(p0, p1, r, eij) ;
    }
    // Break bond can also perform other functions in bondParticle, like change saturation etc.
    virtual bool breakBond(Particle3D<T,Descriptor> * bondParticle) {
        bool breakB = castParticle3DToBondParticle3D(bondParticle)->get_r() > this->getBreakDistance();
        if (breakB) {  bondParticle->setTag(-1); }  // Prepare bondParticle for deletion
        return breakB;
    }
    // UID should be unique and it is suggested to include the bondTypeId.
    virtual std::string getUID(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1) {
        ImmersedCellParticle3D<T,Descriptor>* sp0 = castParticleToICP3D(p0);
        ImmersedCellParticle3D<T,Descriptor>* sp1 = castParticleToICP3D(p1);
        plint cellId0=sp0->get_cellId(), cellId1=sp1->get_cellId();
        plint vertexId0=sp0->getVertexId(), vertexId1=sp1->getVertexId();
        std::string ret;
        std::ostringstream ss, stm0, stm1;
        ss << bondTypeId ;
        stm0 << cellId0 << "-" << vertexId0;
        stm1 << cellId1 << "-" << vertexId1;
        if (areSameCellType and cellId0 > cellId1) { ret = ss.str() + "_" + stm1.str() + "_" + stm0.str(); }
        else { ret = ss.str() + "_" + stm0.str() + "_" + stm1.str();  }
    	return ret;
    } ;
    // Force will be applied via the bondParticle.
    virtual void applyForce(Particle3D<T,Descriptor> * bondParticle)=0;
private:
    T r_create, r_break;
    static plint globalBondTypeId;
    plint bondTypeId;
    bool areSameCellType; // Means if the two edges of the bond belong to the same type --> PLT<->PLT is true, PLT<->Wall is false
};

template<typename T, template<typename U> class Descriptor>
plint BondType<T,Descriptor>::globalBondTypeId=0;



template<typename T, template<typename U> class Descriptor>
class SimpleUnsaturatedBond: public BondType<T,Descriptor> {
public:
    /*
     * Example class, applying force to the attached particles.
     */
    SimpleUnsaturatedBond(CellCellForce3D<T> & forceType_, T r_create_, T r_break_, bool areSameCellType_=false)
       : forceType(forceType_), BondType<T,Descriptor>(r_create_, r_break_, areSameCellType_)  {
        pcout << "SimpleUnsaturatedBond, typeId " << this->getBondTypeId() << std::endl;
    };
    virtual ~SimpleUnsaturatedBond() {};
public:
    virtual void applyForce(Particle3D<T,Descriptor> * bondParticle) {
        BondParticle3D<T,Descriptor> * bondP=castParticle3DToBondParticle3D(bondParticle);
        T r = bondP->get_r();
        Array<T,3> eij = bondP->get_eij();
        bondP->applyForce( forceType(r, eij) );
    }
private:
    CellCellForce3D<T> & forceType;
};


template<typename T, template<typename U> class Descriptor>
class SimpleAsymmetricSaturatedBond: public BondType<T,Descriptor> {
public:
    /*
     * This is to be used when there are two species (eg. platelet-ECM) with different receptor strengths.
     * This class takes care of the Saturation, as well.
     * TODO: Scaling force w/ Saturation
     */
    SimpleAsymmetricSaturatedBond(CellCellForce3D<T> & forceType_, T r_create_, T r_break_,
    		T deltaSaturation0_, T deltaSaturation1_,
    		T maxSaturation0_, T maxSaturation1_,
    		bool areSameCellType_=false)
       : forceType(forceType_), BondType<T,Descriptor>(r_create_, r_break_, areSameCellType_) {
        deltaSaturation.push_back(deltaSaturation0_);
        deltaSaturation.push_back(deltaSaturation1_);
        maxSaturation.push_back(maxSaturation0_);
        maxSaturation.push_back(maxSaturation1_);
        pcout << "SimpleAsymmetricSaturatedBond, typeId " << this->getBondTypeId() << std::endl;
    };
    virtual ~SimpleAsymmetricSaturatedBond() {};
public:
    virtual void applyForce(Particle3D<T,Descriptor> * bondParticle) {
        BondParticle3D<T,Descriptor> * bondP=castParticle3DToBondParticle3D(bondParticle);
        T r = bondP->get_r();
        Array<T,3> eij = bondP->get_eij();
        bondP->applyForce( forceType(r, eij) );
    }
    virtual bool createBond(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
    	if (not isBondPossible(p0, p1, r, eij)) { return false; }
        bool p0Saturated = castParticleToICP3D(p0)->getBondTypeSaturation(this->getBondTypeId())+deltaSaturation[0] > maxSaturation[0];
        if ( p0Saturated ) { return false; }
        bool p1Saturated = castParticleToICP3D(p1)->getBondTypeSaturation(this->getBondTypeId())+deltaSaturation[1] > maxSaturation[1];
        if ( p1Saturated ) { return false; }
		increaseSaturation(p0, deltaSaturation[0]);
		increaseSaturation(p1, deltaSaturation[1]);
        return true;

    }
    virtual bool breakBond(Particle3D<T,Descriptor> * bondParticle) {
    	T r = castParticle3DToBondParticle3D(bondParticle)->get_r();
        bool distanceLarger = r > this->getBreakDistance();
        if (distanceLarger) {
        	Particle3D<T,Descriptor> * p0 = castParticle3DToBondParticle3D(bondParticle)->getParticle(0);
        	Particle3D<T,Descriptor> * p1 = castParticle3DToBondParticle3D(bondParticle)->getParticle(1);
            increaseSaturation(p0, -1.0 * deltaSaturation[0]);
            increaseSaturation(p1, -1.0 * deltaSaturation[1]);
            bondParticle->setTag(-1); // Prepare bondParticle for deletion
        }
        return distanceLarger;
    }
private:
    bool increaseSaturation(Particle3D<T,Descriptor> * pp0, T dSaturation) {
        if (pp0 == NULL) {
        	pcout << "pp0 failure;" << std::endl;
        	return false;
        }
        castParticleToICP3D(pp0)->getBondTypeSaturation( this->getBondTypeId() ) = castParticleToICP3D(pp0)->getBondTypeSaturation( this->getBondTypeId() ) + dSaturation;
        return true;
    };
    bool resetSaturation(Particle3D<T,Descriptor> * pp0) {
        if (pp0 == NULL) { return false; }
        castParticleToICP3D(pp0)->getBondTypeSaturation( this->getBondTypeId() ) = 0;
        return true;
    };
    CellCellForce3D<T> & forceType;
    std::vector<T> deltaSaturation, maxSaturation;
};



template<typename T, template<typename U> class Descriptor>
class SimpleSaturatedBond: public SimpleAsymmetricSaturatedBond<T,Descriptor> {
public:
    /*
     * SimpleSaturatedBond is basically the symmetric case of SimpleAsymmetricSaturatedBond
     */
    SimpleSaturatedBond(CellCellForce3D<T> & forceType_, T r_create_, T r_break_, T deltaSaturation_, T maxSaturation_, bool areSameCellType_=false)
       : SimpleAsymmetricSaturatedBond<T,Descriptor>(forceType_, r_create_, r_break_, deltaSaturation_, deltaSaturation_, maxSaturation_, maxSaturation_, areSameCellType_)
     {
        pcout << "SimpleSaturatedBond, typeId " << this->getBondTypeId() << std::endl;
    };
    virtual ~SimpleSaturatedBond() {};
private:
};




} // namespace trombocit


#include "BondTypes3D.hh"


#endif  // BOND_TYPES_3D_H
