#ifndef BOND_TYPES_3D_H
#define BOND_TYPES_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellCellForces3D.h"
#include <map>
#include <string>
#include <vector>

using namespace plb;


namespace trombocit {


template<typename T, template<typename U> class Descriptor>
class BondType {
public:
    BondType(T r_create_, T r_break_)
        : r_create(r_create_), r_break(r_break_), bondTypeId(globalBondTypeId++) {} ;
    virtual ~BondType() {} ;
public:
    T & getBreakDistance() { return r_create; } ;
    T & getCreateDistance() { return r_break; } ;
    plint getBondTypeId() { return bondTypeId; }
    // Particle p0 and p1 should exist in the same domain (even in envelopes) for a bond to be created.
    //    Hence p0 and p1 always exist and the envelopes should be large enough (>r_create) for this assumption to hold.
    virtual bool createBond(Particle3D<T,Descriptor> * p0, Particle3D<T,Descriptor> * p1, T r, Array<T,3> eij) {
        return r < this->getCreateDistance();
    }
    // Force will be applied via the bondParticle.
    virtual void applyForce(Particle3D<T,Descriptor> * bondParticle)=0;
    // Break bond can also perform other function in bondParticle, like change the type of bond (break, change etc.)
    virtual bool breakBond(Particle3D<T,Descriptor> * bondParticle) {
        return castParticle3DToBondParticle3D(bondParticle)->get_r() > this->getBreakDistance();
    }
    // UID should be unique and it is suggested to include the bondTypeId.
    virtual std::string getUID() {
        return "";
    };
private:
    T r_create, r_break;
    static plint globalBondTypeId;
    plint bondTypeId;
};

template<typename T, template<typename U> class Descriptor>
plint BondType<T,Descriptor>::globalBondTypeId=0;



template<typename T, template<typename U> class Descriptor>
class SimpleUnsaturatedBond: public BondType<T,Descriptor> {
public:
    /*
     * Example class, applying force to the attached particles.
     */
    SimpleUnsaturatedBond(CellCellForce3D<T> & forceType_, T r_create_, T r_break_)
       : forceType(forceType_), BondType<T,Descriptor>(r_create_, r_break_)  {
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
     * This is to be used when there are two species (eg. platelet-ECM) with different receptor strengths
     */
    SimpleAsymmetricSaturatedBond(CellCellForce3D<T> & forceType_, T r_create_, T r_break_, T deltaSaturation0_, T deltaSaturation1_, T maxSaturation0_, T maxSaturation1_)
       : forceType(forceType_), BondType<T,Descriptor>(r_create_, r_break_) {
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
    virtual bool createBond(Particle3D<T,Descriptor> * pp0, Particle3D<T,Descriptor> * pp1, T r, Array<T,3> eij) {
        bool distanceSmaller = r < this->getCreateDistance();
        ImmersedCellParticle3D<T,Descriptor> * p0 = castParticleToICP3D(pp0);
        ImmersedCellParticle3D<T,Descriptor> * p1 = castParticleToICP3D(pp1);
        bool p0Unsaturated = p0->getBondTypeSaturation(this->getBondTypeId())+deltaSaturation[0] <= maxSaturation[0];
        bool p1Unsaturated = p1->getBondTypeSaturation(this->getBondTypeId())+deltaSaturation[1] <= maxSaturation[1];
        bool canCreate = distanceSmaller and p0Unsaturated and p1Unsaturated;
        if (canCreate) {
            p0->getBondTypeSaturation(this->getBondTypeId()) += deltaSaturation[0];
            p1->getBondTypeSaturation(this->getBondTypeId()) += deltaSaturation[1];
        }
        return canCreate;
    }
    virtual bool breakBond(Particle3D<T,Descriptor> * bondParticle) {
        bool distanceLarger = castParticle3DToBondParticle3D(bondParticle)->get_r() > this->getBreakDistance();
        if (distanceLarger) {
            castParticle3DToBondParticle3D(bondParticle)->increaseSaturation(-deltaSaturation[0], 0);
            castParticle3DToBondParticle3D(bondParticle)->increaseSaturation(-deltaSaturation[1], 1);
        }
        return distanceLarger;
    }
private:
    CellCellForce3D<T> & forceType;
    std::vector<T> deltaSaturation, maxSaturation;
};



template<typename T, template<typename U> class Descriptor>
class SimpleSaturatedBond: public SimpleAsymmetricSaturatedBond<T,Descriptor> {
public:
    /*
     * SimpleSaturatedBond is basically the symmetric case of SimpleAsymmetricSaturatedBond
     */
    SimpleSaturatedBond(CellCellForce3D<T> & forceType_, T r_create_, T r_break_, T deltaSaturation_, T maxSaturation_)
       : SimpleAsymmetricSaturatedBond<T,Descriptor>(forceType_, r_create_, r_break_, deltaSaturation_, deltaSaturation_, maxSaturation_, maxSaturation_)
     {
        pcout << "SimpleSaturatedBond, typeId " << this->getBondTypeId() << std::endl;
    };
    virtual ~SimpleSaturatedBond() {};
private:
};




} // namespace trombocit


#include "BondTypes3D.hh"


#endif  // BOND_TYPES_3D_H
