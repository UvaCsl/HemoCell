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

#include "immersedWallParticleFunctional3D.h"
#include "immersedWallParticle3D.h"
#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

/* ******** CreateImmersedWallParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CreateImmersedWallParticle3D<T,Descriptor>::CreateImmersedWallParticle3D (
        TriangleBoundary3D<T> const& triangleBoundary_ )
    : triangleBoundary(triangleBoundary_)
{ }

template<typename T, template<typename U> class Descriptor>
void CreateImmersedWallParticle3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    TriangularSurfaceMesh<T> const& mesh = triangleBoundary.getMesh();
    for (plint iVertex=0; iVertex<mesh.getNumVertices(); ++iVertex) {
        Array<T,3> vertex(mesh.getVertex(iVertex));
        ImmersedWallParticle3D<T,Descriptor>* particle
            = new ImmersedWallParticle3D<T,Descriptor>(iVertex, vertex);
        particleField.addParticle(domain, particle);
    }
}

template<typename T, template<typename U> class Descriptor>
CreateImmersedWallParticle3D<T,Descriptor>* CreateImmersedWallParticle3D<T,Descriptor>::clone() const {
    return new CreateImmersedWallParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CreateImmersedWallParticle3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}

/* ******** CreateTaggedImmersedWallParticle3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CreateTaggedImmersedWallParticle3D<T,Descriptor>::CreateTaggedImmersedWallParticle3D (
        TriangleBoundary3D<T> const& triangleBoundary_, plint tag_, plint numPartsPerTag_ )
    : triangleBoundary(triangleBoundary_), tag(tag_), numPartsPerTag(numPartsPerTag_)
{ }

template<typename T, template<typename U> class Descriptor>
void CreateTaggedImmersedWallParticle3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    TriangularSurfaceMesh<T> const& mesh = triangleBoundary.getMesh();
    for (plint iVertex= tag * numPartsPerTag; iVertex < (tag+1) * numPartsPerTag; ++iVertex) {
        Array<T,3> vertex(mesh.getVertex(iVertex));
        ImmersedWallParticle3D<T,Descriptor>* particle
            = new ImmersedWallParticle3D<T,Descriptor>(iVertex, vertex, tag);
            
        particleField.addParticle(domain, particle);
    }
}

template<typename T, template<typename U> class Descriptor>
CreateTaggedImmersedWallParticle3D<T,Descriptor>* CreateTaggedImmersedWallParticle3D<T,Descriptor>::clone() const {
    return new CreateTaggedImmersedWallParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CreateTaggedImmersedWallParticle3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}


/* ******** FluidVelocityToImmersedWall3D *********************************** */


template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedWall3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    std::vector<Cell<T,Descriptor>*> cells(8);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iParticle];
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        linearInterpolationCoefficients(fluid, position, cellPos, weights);

        // Use copy constructor in order to initialize dynamics object.
        Cell<T,Descriptor>* cellOnVertex;
        for (plint iCell=0; iCell<8; ++iCell) {
            cells[iCell] = &fluid.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z);
        }
        cellOnVertex = new Cell<T,Descriptor>(*cells[0]);
        for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
            (*cellOnVertex)[iPop] =
                weights[0]*(*cells[0])[iPop] + weights[1]*(*cells[1])[iPop] + weights[2]*(*cells[2])[iPop] +
                weights[3]*(*cells[3])[iPop] + weights[4]*(*cells[4])[iPop] + weights[5]*(*cells[5])[iPop] +
                weights[6]*(*cells[6])[iPop] + weights[7]*(*cells[7])[iPop];
        }
        cellOnVertex->computeVelocity(particle->get_v());
        delete cellOnVertex;
    }
}

template<typename T, template<typename U> class Descriptor>
FluidVelocityToImmersedWall3D<T,Descriptor>* FluidVelocityToImmersedWall3D<T,Descriptor>::clone() const {
    return new FluidVelocityToImmersedWall3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedWall3D<T,Descriptor>::getTypeOfModification (
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
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    std::vector<Cell<T,Descriptor>*> cells(8);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iParticle];
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        linearInterpolationCoefficients(fluid, position, cellPos, weights);
        
        Array<T,3> elasticForce = particle->get_force();

        // Use copy constructor in order to initialize dynamics object.
        for (plint iCell=0; iCell<8; ++iCell) {
            cells[iCell] = &fluid.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z);
        }
        for (plint iCell = 0; iCell < 8; ++iCell) {
            T *locForce = cells[iCell]->getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
            for (plint iA = 0; iA < 3; ++iA) {
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
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
            
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


/* ******** ComputeCellVolumeParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::ComputeCellVolumeParticlesFunctional3D(
		TriangleBoundary3D<T> const& triangleBoundary_, plint numberOfCells_)
    : numberOfCells(numberOfCells_), triangleBoundary(triangleBoundary_)
{ }

template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
	PLB_PRECONDITION( blocks.size()==1 );
	std::map<plint, T> tagsVolume;
	typename std::map<plint, T>::iterator tagsVolumeIterator;

	plint meshID = 1, cellId;
	triangleBoundary.pushSelect(0,meshID);
	std::vector<Particle3D<T,Descriptor>*> found;
	TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
	T triangleVolume;
	Array<T,3> areaTimesNormal;
    Array<T,3> v0, v1, v2, tmp;

    for (plint iTriangle=0; iTriangle<triangleMesh.getNumTriangles(); ++iTriangle) {
    	v0 = triangleMesh.getVertex(iTriangle, 0);
    	v1 = triangleMesh.getVertex(iTriangle, 1);
    	v2 = triangleMesh.getVertex(iTriangle, 2);
    	crossProduct(v2, v1, tmp);
    	cellId = triangleMesh.getVertexId(iTriangle, 0);
		tagsVolume[cellId] += VectorTemplate<T,Descriptor>::scalarProduct(tmp,v0)/(T) 6.0;
	}
	triangleBoundary.popSelect();

	for (pluint i=0; i< (pluint) numberOfCells+1; ++i) {
		pcout << this->getStatistics().subscribeSum() << std::endl;
	}
	pcout << "done" << std::endl;
	for ( tagsVolumeIterator=tagsVolume.begin() ; tagsVolumeIterator != tagsVolume.end(); tagsVolumeIterator++ ) {
		cellId = (*tagsVolumeIterator).first;
		pcout << cellId << std::endl;
		triangleVolume = (*tagsVolumeIterator).second;
		this->getStatistics().gatherSum(cellId, (T) triangleVolume);
	}
}

//template<typename T, template<typename U> class Descriptor>
//void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
//        Box3D domain, std::vector<AtomicBlock3D*> blocks )
//{
//	PLB_PRECONDITION( blocks.size()==1 );
//	ParticleField3D<T,Descriptor>& particleField =
//		*dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
//
//	std::vector<Particle3D<T,Descriptor>*> found;
//	particleField.findParticles(domain, found);
//
//	T totalVolume = (T) 0.0;
//	plint meshID = 1;
//	for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
//		Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
//		ImmersedWallParticle3D<T,Descriptor>* particle =
//			dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
//		if (particle->get_cellId() == cellId) {
//			plint vertexId = particle->getTag();
//
//			triangleBoundary.pushSelect(0,meshID);
//			Array<T,3> centralVertex = triangleBoundary.getMesh().getVertex(vertexId);
//			std::vector<plint> neighbors = triangleBoundary.getMesh().getNeighborVertexIds(vertexId);
//			for (pluint iA = 0; iA < neighbors.size()-1; ++iA) {
//				Array<T,3> vertexOne = triangleBoundary.getMesh().getVertex(neighbors[iA]);
//				Array<T,3> vertexTwo = triangleBoundary.getMesh().getVertex(neighbors[iA+1]);
//
//				Array<T,3> baryCenter = (centralVertex + vertexOne + vertexTwo)/(T)3;
//				Array<T,3> normal = triangleBoundary.getMesh().computeTriangleNormal(vertexId,neighbors[iA],neighbors[iA+1]);
//				T surface = triangleBoundary.getMesh().computeTriangleArea(vertexId,neighbors[iA],neighbors[iA+1]);
//
//				T volumeFraction =
//						surface * VectorTemplate<T,Descriptor>::scalarProduct(normal,baryCenter) / (T)3;
//				volumeFraction /= (T)3; // all the volumes are computed three times
//
//				totalVolume += volumeFraction;
//			}
//			triangleBoundary.popSelect();
//		}
//	}
//
//    this->getStatistics().gatherSum(cellVolumeId, (T) totalVolume);
//}

template<typename T, template<typename U> class Descriptor>
ComputeCellVolumeParticlesFunctional3D<T,Descriptor>* ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::clone() const {
    return new ComputeCellVolumeParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}


template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::getCellVolumeArray(std::vector<T>& cellVolumes) const {
	for (pluint i=0; i< (pluint) numberOfCells; ++i)
		cellVolumes.push_back(this->getStatistics().getSum(i));
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
    std::vector<plint> tags;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iP];
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);

        if (particle->get_cellId() == tag) {
            tags.push_back(particle->getTag());
        }
    }

//     pcout << tags.size() << std::endl;
    for (pluint iP = 0; iP < tags.size(); ++iP) {
        particleField.removeParticles(domain,tags[iP]);
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
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (found[iParticle]);

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
        ShellModel3D<T>* shellModel_, int meshID_ )
    : triangleBoundary(triangleBoundary_),
      shellModel(shellModel_), meshID(meshID_)
{ }

template<typename T, template<typename U> class Descriptor>
ComputeImmersedElasticForce3D<T,Descriptor>::~ComputeImmersedElasticForce3D()
{
    delete shellModel;
}

template<typename T, template<typename U> class Descriptor>
ComputeImmersedElasticForce3D<T,Descriptor>::ComputeImmersedElasticForce3D (
            ComputeImmersedElasticForce3D<T,Descriptor> const& rhs)
    : triangleBoundary(rhs.triangleBoundary),
      shellModel(rhs.shellModel->clone()), meshID(rhs.meshID)
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

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getTag();

        if (!isRigid(triangleBoundary.getVertexProperty(vertexId))) {
            Array<T,3> elasticForce = shellModel->computeElasticForce (
                    triangleBoundary, vertexId );
            triangleBoundary.pushSelect(0,meshID);  // Open, dynamic mesh
            T mass = shellModel->getDensity() * triangleBoundary.getMesh().computeVertexArea(vertexId);
            triangleBoundary.popSelect();
            particle->get_a() += elasticForce/mass;
            particle->get_force() += elasticForce;
        }
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

