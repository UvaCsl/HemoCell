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

#ifndef CELL_STRETCHING_FORCES_3D_H
#define CELL_STRETCHING_FORCES_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "immersedCellParticle3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "immersedBoundaryMethod3D.h"
#include <map>

#include "cellStretchingForces3D.hh"

namespace plb {


template<typename T, template<typename U> class Descriptor>
class ApplyStretchingForce3D : public BoxProcessingFunctional3D
{
public:
    ApplyStretchingForce3D (std::vector<plint> const& outerLeftTags_, std::vector<plint> const& outerRightTags_,
                            Array<T,3> const& stretchingForce_, T cellDensity_);
    ~ApplyStretchingForce3D() {} ;
    ApplyStretchingForce3D(ApplyStretchingForce3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ApplyStretchingForce3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> const& outerLeftTags;
    std::vector<plint> const& outerRightTags;
    Array<T,3> const& stretchingForce;
    T const& cellDensity;
};


template<typename T, template<typename U> class Descriptor>
class MeasureCellStretchDeformation3D : public BoxProcessingFunctional3D
{
public:
    MeasureCellStretchDeformation3D (std::vector<std::vector<plint>*> const& tags_, std::vector<T> * meanDeformation_);
    ~MeasureCellStretchDeformation3D() {} ;
    MeasureCellStretchDeformation3D(MeasureCellStretchDeformation3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MeasureCellStretchDeformation3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<std::vector<plint>*> tags;
    std::vector<T> * meanDeformation;
};


template<typename T, template<typename U> class Descriptor>
bool compareParticlesInX (Particle3D<T,Descriptor>* iParticle, Particle3D<T,Descriptor>* jParticle) ;


template<typename T, template<typename U> class Descriptor>
bool compareParticlesInY (Particle3D<T,Descriptor>* iParticle, Particle3D<T,Descriptor>* jParticle) ;


template<typename T, template<typename U> class Descriptor>
bool compareParticlesInZ (Particle3D<T,Descriptor>* iParticle, Particle3D<T,Descriptor>* jParticle) ;

template<typename T, template<typename U> class Descriptor>
class FindTagsOfLateralCellParticles3D : public BoxProcessingFunctional3D
{
public:
    FindTagsOfLateralCellParticles3D (plint numParticlesPerSide_, std::vector<plint> * outerLeftTags_, std::vector<plint> * outerRightTags_, pluint direction_);
    ~FindTagsOfLateralCellParticles3D() {} ;
    FindTagsOfLateralCellParticles3D(FindTagsOfLateralCellParticles3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FindTagsOfLateralCellParticles3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint const& numParticlesPerSide;
    std::vector<plint> * outerLeftTags;
    std::vector<plint> * outerRightTags;
    pluint direction;
};


}  // namespace plb

#endif  // CELL_STRETCHING_FORCES_3D_H

