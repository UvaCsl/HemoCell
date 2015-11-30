#ifndef PROXIMITY_DYNAMICS_H
#define PROXIMITY_DYNAMICS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellCellForces3D.h"

#include <math.h>
using namespace std;
using namespace plb;

template<typename T>
class CellCellForce3D;

template<typename T, template<typename U> class Descriptor>
class CellField3D;

// This function object defines the force between two LSPs of different CellField3D, once their LSPs are in proximity.
// It it to be used as an argument to ApplyProximityDynamics3D
template<typename T, template<typename U> class Descriptor>
class ProximityDynamics3D {
public:
    ProximityDynamics3D () { };
    ProximityDynamics3D (ProximityDynamics3D<T,Descriptor> const& rhs) { }
    virtual ~ProximityDynamics3D () {} ;
    virtual bool operator()(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij) =0;
    virtual void open(Box3D domain, std::vector<AtomicBlock3D*> fields) = 0;
    virtual void close(Box3D domain, std::vector<AtomicBlock3D*> fields) = 0;
    virtual bool conditionsAreMet(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij) =0;
private:
};


/* ******** ApplyProximityDynamics3D *********************************** */

template<typename T, template<typename U> class Descriptor>
/* This class changes and syncs every field it receives */
class ApplyProximityDynamics3D : public BoxProcessingFunctional3D
{
public:
    /* ProximityDynamics3D<T,Descriptor> proximityDynamics
        * should have the following methods:
            * actually performing the proximityDynamics. The checking of additional criteria should be made there.
                bool operator()(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij);
            * opening the dynamics, in case it is necessary.
                bool open(Box3D domain, std::vector<AtomicBlock3D*> fields);
            * closing the dynamics, in case it is necessary.
                bool close(Box3D domain, std::vector<AtomicBlock3D*> fields);
		ProximityDynamics3D<T,Descriptor> * proximityDynamics is deleted afterwards;
    */
    ApplyProximityDynamics3D (ProximityDynamics3D<T,Descriptor> & proximityDynamics_, T cutoffRadius_);
    virtual ~ApplyProximityDynamics3D();
    ApplyProximityDynamics3D(ApplyProximityDynamics3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ApplyProximityDynamics3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
public:
    /*
     All of the following conditions have to be met, in order for the force to be applied:
           * Their distance has to be less that the cutoffRadius.
     */
    bool checkDistance(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T& r, Array<T,3>& eij);
private:
    ProximityDynamics3D<T,Descriptor> & proximityDynamics;
    T cutoffRadius;
};




// This function object defines the force between two LSPs of different CellField3D, once their LSPs are in proximity.
// It it to be used as an argument to ApplyProximityDynamics3D
template<typename T, template<typename U> class Descriptor>
class ProximityCellCellForce3D : public ProximityDynamics3D<T,Descriptor> {
public:
    ProximityCellCellForce3D (CellCellForce3D<T> & forceType_, T cutoffRadius_)
        : forceType(forceType_), cutoffRadius(cutoffRadius_) { }
    ProximityCellCellForce3D (ProximityCellCellForce3D<T,Descriptor> const& rhs)
        : forceType(rhs.forceType), cutoffRadius(rhs.cutoffRadius) { }
    virtual ~ProximityCellCellForce3D () { } ;
    virtual bool operator()(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij) {
        bool conditionsMet = conditionsAreMet(p1,p2,r,eij);
        if (conditionsMet) {
            Array<T,3> force = forceType(r, eij);
            castParticleToICP3D(p1)->get_force() += force;
            castParticleToICP3D(p2)->get_force() -= force;
        }
        return conditionsMet>0;
    }
    virtual void open(Box3D domain, std::vector<AtomicBlock3D*> fields) { };
    virtual void close(Box3D domain, std::vector<AtomicBlock3D*> fields) { };
    virtual bool conditionsAreMet(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij) {
        return r < cutoffRadius;
    };
private:
    CellCellForce3D<T> & forceType;
    T cutoffRadius;
};


// This function object defines the force between two LSPs of the same CellField3D, once the LSPs of two different cells are in proximity.
// It it to be used as an argument to ApplyProximityDynamics3D
template<typename T, template<typename U> class Descriptor>
class ProximitySameCellFieldForce3D : public ProximityDynamics3D<T,Descriptor>  {
public:
    ProximitySameCellFieldForce3D (CellCellForce3D<T> & forceType_, T cutoffRadius_)
        : forceType(forceType_), cutoffRadius(cutoffRadius_) {};
    ProximitySameCellFieldForce3D (ProximitySameCellFieldForce3D<T,Descriptor> const& rhs)
        : forceType(rhs.forceType), cutoffRadius(rhs.cutoffRadius) { }
    virtual ~ProximitySameCellFieldForce3D () { };
    virtual bool operator()(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij) {
        bool conditionsMet = conditionsAreMet(p1,p2,r,eij);
        if (conditionsMet) {
            Array<T,3> force = forceType(r, eij);
            castParticleToICP3D(p1)->get_force() += force;
            castParticleToICP3D(p2)->get_force() -= force;
        }
        return conditionsMet>0;
    }
    virtual void open(Box3D domain, std::vector<AtomicBlock3D*> fields) { };
    virtual void close(Box3D domain, std::vector<AtomicBlock3D*> fields) {};
    virtual bool conditionsAreMet(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> eij) {
        if (r > cutoffRadius) { return false; }
        SurfaceParticle3D<T,Descriptor>* cParticle = castParticleToICP3D(p1);
        SurfaceParticle3D<T,Descriptor>* nParticle = castParticleToICP3D(p2);
        // If they belong to the same cell, don't do anything;
        if (cParticle->get_cellId() <= nParticle->get_cellId()) { return false; }
        return true;
    };
private:
    CellCellForce3D<T> & forceType;
    T cutoffRadius;
};


template<typename T, template<typename U> class Descriptor>
void applySameCellFieldForces(CellField3D<T, Descriptor> & cellField, CellCellForce3D<T> & forceType, T cutoffRadius);



#include "proximityDynamics3D.hh"

#endif  // PROXIMITY_DYNAMICS_H

