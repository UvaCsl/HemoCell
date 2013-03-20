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

#ifndef IMMERSED_WALL_PARTICLE_FUNCTIONAL_3D_HH
#define IMMERSED_WALL_PARTICLE_FUNCTIONAL_3D_HH

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

namespace plb {

/* ******** CreateImmersedCellParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CreateImmersedCellParticle3D<T,Descriptor>::CreateImmersedCellParticle3D (
        TriangleBoundary3D<T> const& triangleBoundary_ )
    : triangleBoundary(triangleBoundary_)
{ }

template<typename T, template<typename U> class Descriptor>
void CreateImmersedCellParticle3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    TriangularSurfaceMesh<T> const& mesh = triangleBoundary.getMesh();
    for (plint iVertex=0; iVertex<mesh.getNumVertices(); ++iVertex) {
        Array<T,3> vertex(mesh.getVertex(iVertex));
        ImmersedCellParticle3D<T,Descriptor>* particle
            = new ImmersedCellParticle3D<T,Descriptor>(iVertex, vertex);
        particleField.addParticle(domain, particle);
    }
}

template<typename T, template<typename U> class Descriptor>
CreateImmersedCellParticle3D<T,Descriptor>* CreateImmersedCellParticle3D<T,Descriptor>::clone() const {
    return new CreateImmersedCellParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CreateImmersedCellParticle3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}

/* ******** CreateTaggedImmersedCellParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CreateTaggedImmersedCellParticle3D<T,Descriptor>::CreateTaggedImmersedCellParticle3D (
        TriangleBoundary3D<T> const& triangleBoundary_, plint tag_, plint numPartsPerTag_ )
    : triangleBoundary(triangleBoundary_), tag(tag_), numPartsPerTag(numPartsPerTag_)
{ }

template<typename T, template<typename U> class Descriptor>
void CreateTaggedImmersedCellParticle3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    TriangularSurfaceMesh<T> const& mesh = triangleBoundary.getMesh();
    for (plint iVertex= tag * numPartsPerTag; iVertex < (tag+1) * numPartsPerTag; ++iVertex) {
        Array<T,3> vertex(mesh.getVertex(iVertex));
        ImmersedCellParticle3D<T,Descriptor>* particle
            = new ImmersedCellParticle3D<T,Descriptor>(iVertex, vertex, tag);
            
        particleField.addParticle(domain, particle);

//        pcout << iVertex << " (" << vertex[0] << " " << vertex[1] << " " << vertex[2] << ") " << tag << ", " << particle->get_cellId() << " " << particle->getTag() << std::endl;
    }
//     pcout << std::endl;
    
//     std::vector<Particle3D<T,Descriptor>*> particles;
//     particleField.findParticles(domain, particles);
//     for (plint iVertex= tag * numPartsPerTag; iVertex < (tag+1) * numPartsPerTag; ++iVertex) {
//         pcout << "particles[" << iVertex << "].cellId ="  <<
//         (dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iVertex]))->get_cellId() << ", real iD = " <<
//         (dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iVertex]))->getTag() << ", tag = " << tag << "; " ;
//         pcout << std::endl;
//     }
//     pcout << std::endl;
}

template<typename T, template<typename U> class Descriptor>
CreateTaggedImmersedCellParticle3D<T,Descriptor>* CreateTaggedImmersedCellParticle3D<T,Descriptor>::clone() const {
    return new CreateTaggedImmersedCellParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CreateTaggedImmersedCellParticle3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}


/* ******** FluidVelocityToImmersedCell3D *********************************** */


template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedCell3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    std::vector<Cell<T,Descriptor>*> cells;
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        Array<T,3> velocity; velocity.resetToZero();
        interpolationCoefficientsPhi2(fluid, position, cellPos, weights);
        particle->get_v().resetToZero();
        for (pluint iCell=0; iCell < weights.size(); ++iCell) {
            velocity.resetToZero();
            fluid.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z).computeVelocity(velocity);
            particle->get_v() += weights[iCell] * velocity;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FluidVelocityToImmersedCell3D<T,Descriptor>* FluidVelocityToImmersedCell3D<T,Descriptor>::clone() const {
    return new FluidVelocityToImmersedCell3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedCell3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing; // Fluid field.
}


/* ******** ForceToFluid3D *********************************** */


template<typename T, template<typename U> class Descriptor>
void ForceToFluid3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    Cell<T,Descriptor>* cell;
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        interpolationCoefficientsPhi2(fluid, position, cellPos, weights);

        Array<T,3> elasticForce = particle->get_force();
        // pcout << "elastic force: (" << elasticForce[0] << ", "<< elasticForce[1] << ", "<< elasticForce[2] << ")\n";
        for (pluint iCell = 0; iCell < weights.size(); ++iCell) {
            cell = &fluid.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z);
            T *locForce = cell->getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
            for (pluint iA = 0; iA < 3; ++iA) {
                locForce[iA] += weights[iCell]*elasticForce[iA];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ForceToFluid3D<T,Descriptor>* ForceToFluid3D<T,Descriptor>::clone() const {
    return new ForceToFluid3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ForceToFluid3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
    modified[1] = modif::staticVariables; // Fluid field.
}


/* ******** CountTaggedParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CountTaggedParticlesFunctional3D<T,Descriptor>::CountTaggedParticlesFunctional3D(plint tag_)
    : numParticlesId(this->getStatistics().subscribeIntSum()), tag(tag_)
{ }

template<typename T, template<typename U> class Descriptor>
void CountTaggedParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    plint numParts = 0;
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
            
        if (particle->get_cellId() == tag) {
//         if (particle->getTag() == tag) {
            ++numParts;
        }
    }
    this->getStatistics().gatherIntSum(numParticlesId, numParts);
}

template<typename T, template<typename U> class Descriptor>
CountTaggedParticlesFunctional3D<T,Descriptor>* CountTaggedParticlesFunctional3D<T,Descriptor>::clone() const {
    return new CountTaggedParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountTaggedParticlesFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
plint CountTaggedParticlesFunctional3D<T,Descriptor>::getNumParticles() const {
    return this->getStatistics().getIntSum(numParticlesId);
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
plint countParticles (
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, plint tag )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    CountTaggedParticlesFunctional3D<T,Descriptor> functional(tag);
    applyProcessingFunctional(functional, domain, particleArg);
    return functional.getNumParticles();
}

/* ******** AbsorbTaggedParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
AbsorbTaggedParticlesFunctional3D<T,Descriptor>::
    AbsorbTaggedParticlesFunctional3D(plint tag_) : tag(tag_)
{ }

template<typename T, template<typename U> class Descriptor>
void AbsorbTaggedParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<plint> cellIds;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iP];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);

        if (particle->get_cellId() == tag) {
            cellIds.push_back(particle->getTag());
        }
    }

//     pcout << cellIds.size() << std::endl;
    for (pluint iP = 0; iP < cellIds.size(); ++iP) {
        particleField.removeParticles(domain,cellIds[iP]);
    }

//     particleField.removeParticles(domain,tag);
}

template<typename T, template<typename U> class Descriptor>
AbsorbTaggedParticlesFunctional3D<T,Descriptor>* AbsorbTaggedParticlesFunctional3D<T,Descriptor>::clone() const {
    return new AbsorbTaggedParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AbsorbTaggedParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** GetTaggedParticleVelocity3D *********************************** */

template<typename T, template<typename U> class Descriptor>
GetTaggedParticleVelocity3D<T,Descriptor>::GetTaggedParticleVelocity3D(plint tag_) : tag(tag_)
{ }

template<typename T, template<typename U> class Descriptor>
void GetTaggedParticleVelocity3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& originalField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& clonedField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> found;
    originalField.findParticles(domain, found);

    Dot3D offset = computeRelativeDisplacement(originalField, clonedField);
    Box3D clonedDomain(domain.shift(offset.x, offset.y, offset.z));

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (found[iParticle]);

//         if (particle->getTag() == tag) {
        if (particle->get_cellId() == tag) {
            Particle3D<T,Descriptor>& originalParticle = *found[iParticle];
            
            Particle3D<T,Descriptor>* clonedParticle = new VisualParticle3D<T,Descriptor> (
                    originalParticle.getTag(), originalParticle.getPosition() );
            std::vector<Array<T,3> > vectors;
            Array<T,3> velocity;
    #ifdef PLB_DEBUG
            bool ok =
    #endif
                originalParticle.getVector(0, velocity);
            PLB_ASSERT( ok );
            vectors.push_back(velocity);
            clonedParticle->setVectors(vectors),
            clonedField.addParticle(clonedDomain, clonedParticle);
        }
        
    }
}

template<typename T, template<typename U> class Descriptor>
GetTaggedParticleVelocity3D<T,Descriptor>*
    GetTaggedParticleVelocity3D<T,Descriptor>::clone() const
{
    return new GetTaggedParticleVelocity3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void GetTaggedParticleVelocity3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Original field.
    modified[1] = modif::dynamicVariables;  // Cloned field.
}


/* ******** ComputeImmersedElasticForce3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeImmersedElasticForce3D<T,Descriptor>::ComputeImmersedElasticForce3D (
        TriangleBoundary3D<T> const& triangleBoundary_,
        CellModel3D<T>* cellModel_,
        std::vector<T> const& cellsVolume_, std::vector<T> const& cellsSurface_)
    : triangleBoundary(triangleBoundary_),
      cellModel(cellModel_),
      cellsVolume(cellsVolume_),
      cellsSurface(cellsSurface_)
{ }

template<typename T, template<typename U> class Descriptor>
ComputeImmersedElasticForce3D<T,Descriptor>::~ComputeImmersedElasticForce3D()
{
    delete cellModel;
}

template<typename T, template<typename U> class Descriptor>
ComputeImmersedElasticForce3D<T,Descriptor>::ComputeImmersedElasticForce3D (
            ComputeImmersedElasticForce3D<T,Descriptor> const& rhs)
    : triangleBoundary(rhs.triangleBoundary),
      cellModel(rhs.cellModel->clone()),
      cellsVolume(rhs.cellsVolume),
      cellsSurface(rhs.cellsSurface)

{ }


template<typename T, template<typename U> class Descriptor>
void ComputeImmersedElasticForce3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);

    std::map< plint, Array<T,3> > particleVelocity;
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getTag();
        particleVelocity[vertexId] = particle->get_v();
    }
    // T eqArea = cellModel->getEquilibriumTriangleArea();
    Array<T,3> sforce; sforce.resetToZero();
    Array<T,3> sf_wlc, sf_bending, sf_volume, sf_surface, sf_shear, sf_viscosity;
    sf_wlc.resetToZero(); sf_bending.resetToZero(); sf_volume.resetToZero(); sf_surface.resetToZero();
    sf_shear.resetToZero(); sf_viscosity.resetToZero();
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getTag();
        T iSurface = T(); // one third of the sum of the areas of all triangles that share the given vertex.
        plint cellId = particle->get_cellId();
        if (!isRigid(triangleBoundary.getVertexProperty(vertexId))) {
            Array<T,3> elasticForce = cellModel->computeElasticForce (
                    triangleBoundary, vertexId );
            Array<T,3> f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity;
            f_wlc.resetToZero(); f_bending.resetToZero(); f_volume.resetToZero();
            f_surface.resetToZero(); f_shear.resetToZero(); f_viscosity.resetToZero();
            Array<T,3> cellForce = cellModel->computeCellForce (
                    triangleBoundary, cellsVolume[cellId], cellsSurface[cellId], iSurface, particleVelocity, vertexId,
                    f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity);
            particle->get_f_wlc() = f_wlc;
            particle->get_f_bending() = f_bending;
            particle->get_f_volume() = f_volume;
            particle->get_f_surface() = f_surface;
            particle->get_f_shear() = f_shear;
            particle->get_f_viscosity() = f_viscosity;

//            pcout << "f_wlc (" <<
//                    particle->get_f_wlc()[0] << ", " <<
//                    particle->get_f_wlc()[1] << ", " <<
//                    particle->get_f_wlc()[2] << ") " << std::endl;
//            T neF=norm(elasticForce), ncF=norm(cellForce);
//            pcout << "forces: "<< neF << " " << ncF << " " << neF/ncF  << ", " << norm(elasticForce-cellForce) <<
//                    " (" << (elasticForce/neF)[0] << ", "  << (cellForce/ncF)[0] << "), " <<
//                    " (" << (elasticForce/neF)[1] << ", "  << (cellForce/ncF)[1] << "), " <<
//                    " (" << (elasticForce/neF)[2] << ", "  << (cellForce/ncF)[2] << "), " <<
//                    std::endl;
            Array<T,3> force, acc; force.resetToZero(); acc.resetToZero();
            force = elasticForce + cellForce;
            // force *= iSurface*1.0/eqArea; // same as triangleBoundary.getMesh().computeVertexArea(vertexId);
            acc = force*1.0 / cellModel->getDensity();
            particle->get_a() = acc;
            particle->get_force() = force;



            sforce += force;
            sf_wlc += f_wlc;
            sf_bending += f_bending;
            sf_volume += f_volume;
            sf_surface += f_surface;
            sf_shear += f_shear;
            sf_viscosity += f_viscosity;
        }
    }
    if (fabs(sforce[0]) + fabs(sforce[1]) + fabs(sforce[2]) > 1e-5) {
        pcout << "sforce (" << sforce[0] << ", " << sforce[1] << ", " << sforce[2] << ") " << std::endl;
        pcout << "sf_wlc (" << sf_wlc[0] << ", " << sf_wlc[1] << ", " << sf_wlc[2] << ") " << std::endl;
        pcout << "sf_bending (" << sf_bending[0] << ", " << sf_bending[1] << ", " << sf_bending[2] << ") " << std::endl;
        pcout << "sf_volume (" << sf_volume[0] << ", " << sf_volume[1] << ", " << sf_volume[2] << ") " << std::endl;
        pcout << "sf_surface (" << sf_surface[0] << ", " << sf_surface[1] << ", " << sf_surface[2] << ") " << std::endl;
        pcout << "sf_shear (" << sf_shear[0] << ", " << sf_shear[1] << ", " << sf_shear[2] << ") " << std::endl;
        pcout << "sf_viscosity (" << sf_viscosity[0] << ", " << sf_viscosity[1] << ", " << sf_viscosity[2] << ") " << std::endl;
    }
}

template<typename T, template<typename U> class Descriptor>
ComputeImmersedElasticForce3D<T,Descriptor>*
    ComputeImmersedElasticForce3D<T,Descriptor>::clone() const
{
    return new ComputeImmersedElasticForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeImmersedElasticForce3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeImmersedElasticForce3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeImmersedElasticForce3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}



template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    getParticlePosAndVelocity (
            MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& originalParticles, plint tag )
{
    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particles( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       originalParticles.getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );

    std::vector<MultiBlock3D*> particleParticleArg;
    particleParticleArg.push_back(&originalParticles);
    particleParticleArg.push_back(particles.get());
    applyProcessingFunctional (
            new GetTaggedParticleVelocity3D<T,Descriptor>(tag), particles->getBoundingBox(), particleParticleArg );

    return particles;
}

}  // namespace plb

#endif  // IMMERSED_WALL_PARTICLE_FUNCTIONAL_3D_HH

