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

#ifndef IMMERSED_WALL_PARTICLE_FUNCTIONAL_3D_H
#define IMMERSED_WALL_PARTICLE_FUNCTIONAL_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <map>

#include "immersedCellParticle3D.h"
#include "cellModel3D.h"
#include "immersedBoundaryMethod3D.h"
#include "meshGeneratingFunctions.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class CreateImmersedCellParticle3D : public BoxProcessingFunctional3D
{
public:
    CreateImmersedCellParticle3D (
            TriangleBoundary3D<T> const& triangleBoundary_ );
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CreateImmersedCellParticle3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangleBoundary3D<T> const& triangleBoundary;
};

template<typename T, template<typename U> class Descriptor>
class CreateTaggedImmersedCellParticle3D : public BoxProcessingFunctional3D
{
public:
    CreateTaggedImmersedCellParticle3D (
            TriangleBoundary3D<T> const& triangleBoundary_, plint tag_, plint numPartsPerTag_ );
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CreateTaggedImmersedCellParticle3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    TriangleBoundary3D<T> const& triangleBoundary;
    plint tag, numPartsPerTag;
};

template<typename T, template<typename U> class Descriptor>
class FluidVelocityToImmersedCell3D : public BoxProcessingFunctional3D
{
public:
    FluidVelocityToImmersedCell3D (plint ibmKernel_=4 );
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FluidVelocityToImmersedCell3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint ibmKernel;
};

template<typename T, template<typename U> class Descriptor>
class ForceToFluid3D : public BoxProcessingFunctional3D
{
public:
    ForceToFluid3D (plint ibmKernel_=4 );
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ForceToFluid3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint ibmKernel;
};

/// Count the number of particles, no matter which cellId, found inside the domain.
template<typename T, template<typename U> class Descriptor>
class CountTaggedParticlesFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CountTaggedParticlesFunctional3D(plint tag_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountTaggedParticlesFunctional3D<T,Descriptor>* clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint numParticlesId, tag;
};

/// Remove all particles of a certain tag from a given domain.
template<typename T, template<typename U> class Descriptor>
class AbsorbTaggedParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    AbsorbTaggedParticlesFunctional3D(plint tag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    /// Argument: Particle-field.
    virtual AbsorbTaggedParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint tag;
};

template<typename T, template<typename U> class Descriptor>
class GetTaggedParticleVelocity3D : public BoxProcessingFunctional3D
{
public:
    GetTaggedParticleVelocity3D(plint tag_);
    /// Arguments: [0] Orignal field (point-particles), [1] Cloned field (visual particles).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual GetTaggedParticleVelocity3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint tag;
};

template<typename T, template<typename U> class Descriptor>
class ComputeImmersedElasticForce3D : public BoxProcessingFunctional3D
{
public:
    ComputeImmersedElasticForce3D (
            TriangleBoundary3D<T> const& triangleBoundary_,
            RBCModel3D<T>* cellModel_,
            std::vector<T> const& cellsVolume_, std::vector<T> const& cellsSurface_);
    ~ComputeImmersedElasticForce3D();
    ComputeImmersedElasticForce3D(ComputeImmersedElasticForce3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeImmersedElasticForce3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangleBoundary3D<T> const& triangleBoundary;
    RBCModel3D<T>* cellModel;
    std::vector<T> const& cellsVolume;
    std::vector<T> const& cellsSurface;
};

/// Requirement: particles must be of type point-particle.
template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    getParticlePosAndVelocity (
            MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& originalParticles, plint tag );



}  // namespace plb

#include "immersedCellParticleFunctional3D.hh"
#endif  // IMMERSED_WALL_PARTICLE_FUNCTIONAL_3D_H

