#ifndef CELL_CELL_FORCES_3D_HH
#define CELL_CELL_FORCES_3D_HH

#include "cellCellForces3D.h"



/* ******** ComputeCellCellForces3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeCellCellForces3D<T,Descriptor>::ComputeCellCellForces3D (CellCellForce3D<T> & calcForce_, T cutoffRadius_)
: calcForce(calcForce_), cutoffRadius(cutoffRadius_) { }

template<typename T, template<typename U> class Descriptor>
ComputeCellCellForces3D<T,Descriptor>::~ComputeCellCellForces3D() {
//    std::cout <<" ~ComputeCellCellForces3D() " << global::mpi().getRank() << std::endl;
}

template<typename T, template<typename U> class Descriptor>
ComputeCellCellForces3D<T,Descriptor>::ComputeCellCellForces3D ( ComputeCellCellForces3D<T,Descriptor> const& rhs)
: calcForce(rhs.calcForce), cutoffRadius(rhs.cutoffRadius) { }

template<typename T, template<typename U> class Descriptor>
bool ComputeCellCellForces3D<T,Descriptor>::conditionsAreMet (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij)
{
    SurfaceParticle3D<T,Descriptor>* cParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p1);
    SurfaceParticle3D<T,Descriptor>* nParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p2);
    if (cParticle->get_cellId() == nParticle->get_cellId()) { return false; }
    if (cParticle->getTag() < nParticle->getTag()) { return false; }
    eij = cParticle->getPosition() - nParticle->getPosition();  // TODO: when force is positive, it will lead to adhesion
    r = norm(eij);
    if (r > cutoffRadius) { return false; }
    eij = eij * (1.0/r);
    return true;
}


template<typename T, template<typename U> class Descriptor>
void ComputeCellCellForces3D<T,Descriptor>::applyForce (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> & eij)
{
        SurfaceParticle3D<T,Descriptor>* cParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p1);
        SurfaceParticle3D<T,Descriptor>* nParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p2);
        Array<T,3> force = calcForce(r, eij);
        cParticle->get_force() += force;
        nParticle->get_force() -= force;

#ifdef PLB_DEBUG // Less Calculations
        cParticle->get_f_repulsive() += force;
        nParticle->get_f_repulsive() -= force;
        T potential = calcForce.calculatePotential(r,eij);
        nParticle->get_E_repulsive() += potential;
        nParticle->get_E_repulsive() += potential;
#endif
}


template<typename T, template<typename U> class Descriptor>
void ComputeCellCellForces3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    Box3D extendedDomain = particleField.getBoundingBox();
    T r;
    Array<T,3> eij;
    plint dR = ceil(cutoffRadius);

    std::vector<Particle3D<T,Descriptor>*> currentParticles, neighboringParticles;
    for (plint iX=extendedDomain.x0; iX<=extendedDomain.x1; ++iX) {
        for (plint iY=extendedDomain.y0; iY<=extendedDomain.y1; ++iY) {
            for (plint iZ=extendedDomain.z0; iZ<=extendedDomain.z1; ++iZ) {
                currentParticles.clear(); neighboringParticles.clear();

                Box3D currentBox(iX,iX,iY,iY,iZ,iZ);
                plint x0, x1, y0, y1, z0, z1;
                x0 = max(extendedDomain.x0, iX-dR); x1 = min(extendedDomain.x1, iX+dR);
                y0 = max(extendedDomain.y0, iY-dR); y1 = min(extendedDomain.y1, iY+dR);
                z0 = max(extendedDomain.z0, iZ-dR); z1 = min(extendedDomain.z1, iZ+dR);

                Box3D neighboringBoxes(x0, x1, y0, y1, z0, z1);
                particleField.findParticles(currentBox, currentParticles);
                if (currentParticles.size() == 0) { continue; }
                particleField.findParticles(neighboringBoxes, neighboringParticles);
                for (pluint cP=0; cP<currentParticles.size(); ++cP) {
                    for (pluint nP=0; nP<neighboringParticles.size(); ++nP) {
                        if (conditionsAreMet(currentParticles[cP], neighboringParticles[nP], r, eij)) {
                            applyForce(currentParticles[cP], neighboringParticles[nP], r, eij);
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ComputeCellCellForces3D<T,Descriptor>*
    ComputeCellCellForces3D<T,Descriptor>::clone() const
{
    return new ComputeCellCellForces3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellCellForces3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeCellCellForces3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellCellForces3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}




/* ******** ComputeDifferentCellForces3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeDifferentCellForces3D<T,Descriptor>::ComputeDifferentCellForces3D (CellCellForce3D<T> & calcForce_, T cutoffRadius_)
: calcForce(calcForce_), cutoffRadius(cutoffRadius_) { }

template<typename T, template<typename U> class Descriptor>
ComputeDifferentCellForces3D<T,Descriptor>::~ComputeDifferentCellForces3D() { }

template<typename T, template<typename U> class Descriptor>
ComputeDifferentCellForces3D<T,Descriptor>::ComputeDifferentCellForces3D ( ComputeDifferentCellForces3D<T,Descriptor> const& rhs)
: calcForce(rhs.calcForce), cutoffRadius(rhs.cutoffRadius) { }

template<typename T, template<typename U> class Descriptor>
bool ComputeDifferentCellForces3D<T,Descriptor>::conditionsAreMet (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij)
{
    eij = p1->getPosition() - p2->getPosition();
    r = norm(eij);
    if (r > cutoffRadius) { return false; }
    eij = eij * (1.0/r);
    return true;
}


template<typename T, template<typename U> class Descriptor>
void ComputeDifferentCellForces3D<T,Descriptor>::applyForce (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij)
{
        SurfaceParticle3D<T,Descriptor>* cParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p1);
        SurfaceParticle3D<T,Descriptor>* nParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p2);
        Array<T,3> force = calcForce(r, eij);
        cParticle->get_force() += force;
        nParticle->get_force() -= force;
#ifdef PLB_DEBUG // Less Calculations
        cParticle->get_f_repulsive() -= force;
        nParticle->get_f_repulsive() += force;
        T potential = calcForce.calculatePotential(r,eij);
        cParticle->get_E_repulsive() -= potential;
        nParticle->get_E_repulsive() -= potential;
#endif
}


template<typename T, template<typename U> class Descriptor>
void ComputeDifferentCellForces3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField1 =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& particleField2 =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);

    T r;
    Array<T,3> eij;
    plint dR = ceil(cutoffRadius);

    std::vector<Particle3D<T,Descriptor>*> currentParticles, neighboringParticles;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                currentParticles.clear(); neighboringParticles.clear();

                Box3D currentBox(iX,iX,iY,iY,iZ,iZ);
                Box3D neighboringBoxes(iX-dR,iX+dR,iY-dR,iY+dR,iZ-dR,iZ+dR);
                particleField1.findParticles(currentBox, currentParticles);
                if (currentParticles.size() == 0) { continue; }
                particleField2.findParticles(neighboringBoxes, neighboringParticles);
                if (neighboringParticles.size() == 0) { continue; }

                for (pluint cP=0; cP<currentParticles.size(); ++cP) {
                    for (pluint nP=0; nP<neighboringParticles.size(); ++nP) {
                        if (conditionsAreMet(currentParticles[cP], neighboringParticles[nP], r, eij)) {
                            applyForce(currentParticles[cP], neighboringParticles[nP], r, eij);
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ComputeDifferentCellForces3D<T,Descriptor>*
    ComputeDifferentCellForces3D<T,Descriptor>::clone() const
{
    return new ComputeDifferentCellForces3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeDifferentCellForces3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeDifferentCellForces3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeDifferentCellForces3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}


/* ******** ComputeWallCellForces3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeWallCellForces3D<T,Descriptor>::ComputeWallCellForces3D (CellCellForce3D<T> & calcForce_, T cutoffRadius_)
: calcForce(calcForce_), cutoffRadius(cutoffRadius_) { }

template<typename T, template<typename U> class Descriptor>
ComputeWallCellForces3D<T,Descriptor>::~ComputeWallCellForces3D() { }

template<typename T, template<typename U> class Descriptor>
ComputeWallCellForces3D<T,Descriptor>::ComputeWallCellForces3D ( ComputeWallCellForces3D<T,Descriptor> const& rhs)
: calcForce(rhs.calcForce), cutoffRadius(rhs.cutoffRadius) { }

template<typename T, template<typename U> class Descriptor>
bool ComputeWallCellForces3D<T,Descriptor>::conditionsAreMet (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij)
{
    eij = p2->getPosition() - p1->getPosition(); //NOTE: force < 0: adhesion;  force > 0: repulsion
    r = norm(eij);
    if (r > cutoffRadius) { return false; }
    eij = eij * (1.0/r);
    return true;
}


template<typename T, template<typename U> class Descriptor>
void ComputeWallCellForces3D<T,Descriptor>::applyForce (
        Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij)
{
        SurfaceParticle3D<T,Descriptor>* nParticle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (p2);
        Array<T,3> force = calcForce(r, eij);
        nParticle->get_force() += force;
        //cout << "(ComputeWallCellForces3D) r=" << r << ", norm-force=" << norm(force) << std::endl;
#ifdef PLB_DEBUG // Less Calculations
        nParticle->get_f_repulsive() += force;
        T potential = calcForce.calculatePotential(r,eij);
        nParticle->get_E_repulsive() -= potential;
#endif
}


template<typename T, template<typename U> class Descriptor>
void ComputeWallCellForces3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& wallParticleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);

    T r;
    Array<T,3> eij;
    plint dR = ceil(cutoffRadius);

    std::vector<Particle3D<T,Descriptor>*> currentWallParticles, neighboringCellParticles;

    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                currentWallParticles.clear(); neighboringCellParticles.clear();

                Box3D currentBox(iX,iX,iY,iY,iZ,iZ);
                Box3D neighboringBoxes(iX-dR,iX+dR,iY-dR,iY+dR,iZ-dR,iZ+dR);
                wallParticleField.findParticles(currentBox, currentWallParticles);
                if (currentWallParticles.size() == 0) { continue; }
                particleField.findParticles(neighboringBoxes, neighboringCellParticles);
                if (neighboringCellParticles.size() == 0) { continue; }

                for (pluint cP=0; cP<currentWallParticles.size(); ++cP) {
                    for (pluint nP=0; nP<neighboringCellParticles.size(); ++nP) {
                        if (conditionsAreMet(currentWallParticles[cP], neighboringCellParticles[nP], r, eij)) {
                            applyForce(currentWallParticles[cP], neighboringCellParticles[nP], r, eij);
                        }
                    }
                }
            }
        }
    }
    
    /*
    Dot3D offset = wallParticleField.getLocation();
    Box3D rdomain(domain.x0 + offset.x, domain.x1 + offset.x, domain.y0 + offset.y, domain.y1 + offset.y, domain.z0 + offset.z, domain.z1 + offset.z);

    // TODO: this does not change -> should be turned into a static lookup list
    wallParticleField.findParticles(rdomain, currentWallParticles);
    if (currentWallParticles.size() == 0) { return; }

    for (pluint cP=0; cP<currentWallParticles.size(); ++cP) {
        
        neighboringCellParticles.clear();

        Array<T,3> position = currentWallParticles[cP]->getPosition();
        plint iX = util::roundToInt(position[0]);//-offset.x;
        plint iY = util::roundToInt(position[1]);//-offset.y;
        plint iZ = util::roundToInt(position[2]);//-offset.z;

        Box3D neighboringBoxes(iX-dR,iX+dR,iY-dR,iY+dR,iZ-dR,iZ+dR);
        
        //cout << "WP position: " << position[0] << " " << position[1] << " " << position[2] << " Domain: " << domain.x0 << " " << domain.x1 << " " << domain.y0 << " " << domain.y1 << " " << domain.z0 << " " << domain.z1 << " wall part offset: "<< offset.x << " " << offset.y << " " <<  offset.z <<endl;
        
        particleField.findParticles(neighboringBoxes, neighboringCellParticles);
        
        if (neighboringCellParticles.size() == 0) { continue; }

        for (pluint nP=0; nP<neighboringCellParticles.size(); ++nP) {
            if (conditionsAreMet(currentWallParticles[cP], neighboringCellParticles[nP], r, eij)) {
                applyForce(currentWallParticles[cP], neighboringCellParticles[nP], r, eij);
            }
        }
    }
    */

}

template<typename T, template<typename U> class Descriptor>
ComputeWallCellForces3D<T,Descriptor>*
    ComputeWallCellForces3D<T,Descriptor>::clone() const
{
    return new ComputeWallCellForces3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeWallCellForces3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;  // Wall Particle field.
    isWritten[1] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeWallCellForces3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeWallCellForces3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Wall Particle field.
    modified[1] = modif::allVariables; // Particle field.
}


/* ============== CellCellForces =================== */

template<typename T>
Array<T,3> CellCellForce3D<T>::operator() (T r, Array<T,3> & eij) {
    return calculateForce (r, eij);
}


template<typename T>
Array<T,3> CellCellForce3D<T>::operator() (Array<T,3> & x1, Array<T,3> & x2) {
    Array<T,3> eij = (x1 - x2)*1.0;
    T r = norm(eij);
    eij = eij * (1.0/r);
    return calculateForce (r, eij);
}

template<typename T>
T CellCellForce3D<T>::calculatePotential (T r, Array<T,3> & eij) {
    return calculatePotential (r);
}




#endif  // CELL_CELL_FORCES_3D_HH

