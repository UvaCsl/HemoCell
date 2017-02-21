#ifndef PROXIMITY_DYNAMICS_HH
#define PROXIMITY_DYNAMICS_HH

#include "proximityDynamics3D.h"




template<typename T, template<typename U> class Descriptor>
void applySameCellFieldForces(HemoCellField & cellField, CellCellForce3D<T> & forceType, T cutoffRadius) {
    std::vector<MultiBlock3D*> particleFieldArg;
    particleFieldArg.push_back(cellField.getParticleArg()[0]);
    particleFieldArg.push_back(cellField.getParticleArg()[0]);
    ProximitySameCellFieldForce3D<T, Descriptor> pscd(forceType, cutoffRadius);
    applyProcessingFunctional (
        new ApplyProximityDynamics3D<T,Descriptor>(pscd, cutoffRadius),
        particleFieldArg[0]->getBoundingBox(), particleFieldArg );
}




/* ******** ApplyProximityDynamics3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ApplyProximityDynamics3D<T,Descriptor>::ApplyProximityDynamics3D (ProximityDynamics3D<T,Descriptor> & proximityDynamics_, T cutoffRadius_)
: proximityDynamics(proximityDynamics_), cutoffRadius(cutoffRadius_) { }

template<typename T, template<typename U> class Descriptor>
ApplyProximityDynamics3D<T,Descriptor>::~ApplyProximityDynamics3D() { }

template<typename T, template<typename U> class Descriptor>
ApplyProximityDynamics3D<T,Descriptor>::ApplyProximityDynamics3D (ApplyProximityDynamics3D<T,Descriptor> const& rhs)
: proximityDynamics(rhs.proximityDynamics), cutoffRadius(rhs.cutoffRadius) { }

template<typename T, template<typename U> class Descriptor>
bool ApplyProximityDynamics3D<T,Descriptor>::checkDistance (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij)
{
    eij = p1->getPosition() - p2->getPosition();
    r = norm(eij);
    if (r > cutoffRadius) { return false; }
    eij = eij * (1.0/r);
    return true;
}

template<typename T, template<typename U> class Descriptor>
void ApplyProximityDynamics3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()>=2 );
    ParticleField3D<T,Descriptor>& particleField1 =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& particleField2 =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);
    Box3D extendedDomain1 = particleField1.getBoundingBox();
    Box3D extendedDomain2 = particleField2.getBoundingBox();

    Dot3D offset = computeRelativeDisplacement(particleField1, particleField2);

    T r;
    Array<T,3> eij;
    plint dR = ceil(cutoffRadius);

    std::vector<Particle3D<T,Descriptor>*> currentParticles, neighboringParticles;
    proximityDynamics.open(domain, blocks);
    for (plint iX=extendedDomain1.x0; iX<=extendedDomain1.x1; ++iX) {
        for (plint iY=extendedDomain1.y0; iY<=extendedDomain1.y1; ++iY) {
            for (plint iZ=extendedDomain1.z0; iZ<=extendedDomain1.z1; ++iZ) {
                currentParticles.clear(); neighboringParticles.clear();

                plint x0,x1,y0,y1,z0,z1;
                plint iXo=iX-offset.x, iYo=iY-offset.y, iZo=iZ-offset.z;
                x0 = max(extendedDomain2.x0, iXo-dR); x1 = min(extendedDomain2.x1, iXo+dR);
                y0 = max(extendedDomain2.y0, iYo-dR); y1 = min(extendedDomain2.y1, iYo+dR);
                z0 = max(extendedDomain2.z0, iZo-dR); z1 = min(extendedDomain2.z1, iZo+dR);

                Box3D currentBox(iX,iX,iY,iY,iZ,iZ);
                Box3D neighboringBoxes(x0, x1, y0, y1, z0, z1);

                particleField1.findParticles(currentBox, currentParticles);
                if (currentParticles.size() == 0) { continue; }

                particleField2.findParticles(neighboringBoxes, neighboringParticles);
                if (neighboringParticles.size() == 0) { continue; }

                for (pluint cP=0; cP<currentParticles.size(); ++cP) {
                    for (pluint nP=0; nP<neighboringParticles.size(); ++nP) {
                        if (checkDistance(currentParticles[cP], neighboringParticles[nP], r, eij)) {
                            proximityDynamics(currentParticles[cP], neighboringParticles[nP], r, eij);
                        }
                    }
                }
            }
        }
    }
    proximityDynamics.close(domain, blocks);
}

template<typename T, template<typename U> class Descriptor>
ApplyProximityDynamics3D<T,Descriptor>*
    ApplyProximityDynamics3D<T,Descriptor>::clone() const
{
    return new ApplyProximityDynamics3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ApplyProximityDynamics3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (int i = 0; i < isWritten.size(); ++i)
        { isWritten[i] = true; }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ApplyProximityDynamics3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk; // Takes the whole domain in the processGenericBlocks
}

template<typename T, template<typename U> class Descriptor>
void ApplyProximityDynamics3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (int i = 0; i < modified.size(); ++i)
        { modified[i] = modif::allVariables; }
}




#endif  // PROXIMITY_DYNAMICS_HH

