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

#ifndef CELL_STRETCHING_FORCES_3D_HH
#define CELL_STRETCHING_FORCES_3D_HH

#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "immersedCellParticleFunctional3D.h"
#include "immersedCellParticle3D.h"
#include "immersedBoundaryMethod3D.h"
#include <map>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include "cellStretchingForces3D.h"
#include "immersedCellsReductions.hh"

namespace plb {


/* ******** ApplyStretchingForce3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ApplyStretchingForce3D<T,Descriptor>::ApplyStretchingForce3D (
        std::vector<plint> const& outerLeftTags_, std::vector<plint> const& outerRightTags_,
        Array<T,3> const& stretchingForce_, T cellDensity_)
    : outerLeftTags(outerLeftTags_),
      outerRightTags(outerRightTags_),
      stretchingForce(stretchingForce_),
      cellDensity(cellDensity_)
{ }


template<typename T, template<typename U> class Descriptor>
ApplyStretchingForce3D<T,Descriptor>::ApplyStretchingForce3D (
        ApplyStretchingForce3D<T,Descriptor> const& rhs)
    : outerLeftTags(rhs.outerLeftTags),
      outerRightTags(rhs.outerRightTags),
      stretchingForce(rhs.stretchingForce),
      cellDensity(rhs.cellDensity)
{ }


template<typename T, template<typename U> class Descriptor>
void ApplyStretchingForce3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);
    map<plint, ImmersedCellParticle3D<T,Descriptor>*> tagToParticle;
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (found[iParticle]);
        tagToParticle[particle->getTag()] = particle;
    }
    plint numOuterLeftTags = outerLeftTags.size();
    plint numOuterRightTags = outerRightTags.size();
    for(plint it = 0; it < numOuterLeftTags; it++) {
        plint tag = outerLeftTags[it];
        if (tagToParticle.find(tag) != tagToParticle.end()) {
            ImmersedCellParticle3D<T,Descriptor>* particle = (tagToParticle[tag]);
            //particle->get_a() += -(stretchingForce * 1.0/numOuterLeftTags) * 1.0/cellDensity;
            particle->get_force() = particle->get_force() -stretchingForce * (1.0/numOuterLeftTags);
        } else pcout << "ImmerseCellParticle3D not found! Something is wrong here!" << std::endl;
    }
    for(plint it = 0; it < numOuterRightTags; it++) {
        plint tag = outerRightTags[it];
        if (tagToParticle.find(tag) != tagToParticle.end()) {
            ImmersedCellParticle3D<T,Descriptor>* particle = (tagToParticle[tag]);
            particle->get_a() = particle->get_a() + stretchingForce * (1.0/numOuterRightTags)/cellDensity;
            particle->get_force() = particle->get_force() + stretchingForce * (1.0/numOuterRightTags);
        } else pcout << "ImmerseCellParticle3D not found! Something is wrong here!" << std::endl;
    }
//    pcout << (stretchingForce * (1.0/numOuterRightTags))[0] <<", " << (stretchingForce * (1.0/numOuterRightTags))[1]
   //<<", " << (stretchingForce * (1.0/numOuterRightTags))[2] << std::endl;

}

template<typename T, template<typename U> class Descriptor>
ApplyStretchingForce3D<T,Descriptor>*
    ApplyStretchingForce3D<T,Descriptor>::clone() const
{
    return new ApplyStretchingForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ApplyStretchingForce3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ApplyStretchingForce3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ApplyStretchingForce3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::allVariables; // Particle field.
}


/* ******** FindTagsOfLateralCellParticles3D *********************************** */
template<typename T, template<typename U> class Descriptor>
bool compareParticlesInX (Particle3D<T,Descriptor>* iParticle, Particle3D<T,Descriptor>* jParticle) {
    T iX = iParticle->getPosition()[0];
    T jX = jParticle->getPosition()[0];
    return (iX<jX);
}


template<typename T, template<typename U> class Descriptor>
bool compareParticlesInY (Particle3D<T,Descriptor>* iParticle, Particle3D<T,Descriptor>* jParticle) {
    T iY = iParticle->getPosition()[1];
    T jY = jParticle->getPosition()[1];
    return (iY<jY);
}


template<typename T, template<typename U> class Descriptor>
bool compareParticlesInZ (Particle3D<T,Descriptor>* iParticle, Particle3D<T,Descriptor>* jParticle) {
    T iZ = iParticle->getPosition()[2];
    T jZ = jParticle->getPosition()[2];
    return (iZ<jZ);
}


template<typename T, template<typename U> class Descriptor>
FindTagsOfLateralCellParticles3D<T,Descriptor>::FindTagsOfLateralCellParticles3D
       (plint numParticlesPerSide_, std::vector<plint> * outerLeftTags_,
               std::vector<plint> * outerRightTags_, pluint direction_=0)
    : numParticlesPerSide(numParticlesPerSide_),
      outerLeftTags(outerLeftTags_),
      outerRightTags(outerRightTags_),
      direction(direction_)
{ }


template<typename T, template<typename U> class Descriptor>
FindTagsOfLateralCellParticles3D<T,Descriptor>::FindTagsOfLateralCellParticles3D (
        FindTagsOfLateralCellParticles3D<T,Descriptor> const& rhs)
    : numParticlesPerSide(rhs.numParticlesPerSide),
      outerLeftTags(rhs.outerLeftTags),
      outerRightTags(rhs.outerRightTags),
      direction(rhs.direction)
{ }


template<typename T, template<typename U> class Descriptor>
void FindTagsOfLateralCellParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
            PLB_PRECONDITION( blocks.size()==1 );
            ParticleField3D<T,Descriptor>& particleField =
                *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

            std::vector<Particle3D<T,Descriptor>*> found;
            particleField.findParticles(domain, found);
            if (direction == 0) {
                std::sort(found.begin(), found.end(), compareParticlesInX<T,Descriptor>);
            } else if (direction == 1) {
                std::sort(found.begin(), found.end(), compareParticlesInY<T,Descriptor>);
            } else if (direction == 2) {
                std::sort(found.begin(), found.end(), compareParticlesInZ<T,Descriptor>);
            }
            plint numParticles = found.size();
            for (plint iP = 0; iP < numParticlesPerSide; ++iP) {
                outerLeftTags->push_back(found[iP]->getTag());
                outerRightTags->push_back(found[numParticles - iP - 1]->getTag());
            }
}

template<typename T, template<typename U> class Descriptor>
FindTagsOfLateralCellParticles3D<T,Descriptor>*
FindTagsOfLateralCellParticles3D<T,Descriptor>::clone() const
{
    return new FindTagsOfLateralCellParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FindTagsOfLateralCellParticles3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT FindTagsOfLateralCellParticles3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void FindTagsOfLateralCellParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}


}  // namespace plb

#endif  // CELL_STRETCHING_FORCES_3D_HH

