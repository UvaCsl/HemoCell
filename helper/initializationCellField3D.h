#ifndef INITIALIZATION_CELL_FIELD_3D_H
#define INITIALIZATION_CELL_FIELD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellField3D.h"
#include "meshMetrics.h"
#include <stdlib.h>     /* srand, rand */


using namespace std;
using namespace plb;


double guessRandomNumber() {
    return rand()*1.0/RAND_MAX;
}

// Rotates the mesh in a random angle.
// Mesh is supposed to be centered in (0,0,0)
template<typename T>
void meshRandomRotation (TriangularSurfaceMesh<T> * mesh) {
    const T pi = 4.*atan(1.);
    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    mesh->translate(-1.0 * meshCenter);
    mesh->rotate(guessRandomNumber() * 2 * pi, guessRandomNumber() * pi, guessRandomNumber() * 2 * pi);

    mesh->computeBoundingBox (xRange, yRange, zRange);
    mesh->translate(Array<T,3>(-xRange[0], -yRange[0], -zRange[0]));
}


template<typename T, template<typename U> class Descriptor>
void positionCellInParticleField(ParticleField3D<T,Descriptor>& particleField, BlockLattice3D<T,Descriptor>& fluid,
                                            TriangularSurfaceMesh<T> * mesh, Array<T,3> startingPoint, plint cellId) {
    plint nVertices=mesh->getNumVertices();
    Dot3D relativeDisplacement = computeRelativeDisplacement(fluid, particleField);
    Dot3D fluidLocationDot3D = fluid.getLocation();
    plint iX, iY, iZ;
    plint maxNx = fluid.getNx(), maxNy = fluid.getNy(), maxNz = fluid.getNz();
    for (plint iVertex=0; iVertex < nVertices; ++iVertex) {
        Array<T,3> vertex = startingPoint + mesh->getVertex(iVertex);
        particleField.computeGridPosition (vertex, iX, iY, iZ);
        Dot3D fluidDomainCell = Dot3D(iX, iY, iZ) - relativeDisplacement;
        particleField.addParticle(particleField.getBoundingBox(), new ImmersedCellParticle3D<T,Descriptor>(vertex, cellId, iVertex) );
        bool isInsideFluidDomain = (fluidDomainCell.x >= 0 and fluidDomainCell.y >= 0 and fluidDomainCell.z >= 0) and
        		(fluidDomainCell.x < maxNx and fluidDomainCell.y < maxNy and fluidDomainCell.z < maxNz);
        if (isInsideFluidDomain and fluid.get(fluidDomainCell.x, fluidDomainCell.y, fluidDomainCell.z).getDynamics().isBoundary()) {
        	break;
        }
        else {
//            particleField.addParticle(particleField.getBoundingBox(), new ImmersedCellParticle3D<T,Descriptor>(vertex, cellId, iVertex) );
        }
    }
}


template<typename T, template<typename U> class Descriptor>
plint deleteIncompleteCells(ParticleField3D<T,Descriptor>& particleField, BlockLattice3D<T,Descriptor>& fluid, Box3D domain, plint nVertices) {
    // DELETE CELLS THAT ARE NOT WHOLE
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::map<plint, plint> cellIdToNumberOfVertices;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iP]);
        cellIdToNumberOfVertices[particle->get_cellId()]++;
    }
    plint cellsDeleted=0;
    typename std::map<plint, plint >::iterator itrt;
    for (itrt  = cellIdToNumberOfVertices.begin(); itrt != cellIdToNumberOfVertices.end(); ++itrt) {
        if (itrt->second < nVertices) {
            cellsDeleted++;
            particleField.removeParticles(particleField.getBoundingBox(), itrt->first);
        }
    }
    return cellsDeleted;
}


template<typename T, template<typename U> class Descriptor>
plint deleteIncompleteCells(ParticleField3D<T,Descriptor>& particleField, BlockLattice3D<T,Descriptor>& fluid, plint nVertices) {
    return deleteIncompleteCells(particleField, fluid, particleField.getBoundingBox(), nVertices);
}


template<typename T, template<typename U> class Descriptor>
class PositionCellParticles3D : public BoxProcessingFunctional3D
{
public:
    PositionCellParticles3D (
            TriangularSurfaceMesh<T> const& elementaryMesh_, std::vector<Array<T,3> > const& cellOrigins_) :
    elementaryMesh(elementaryMesh_), cellOrigins(cellOrigins_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual PositionCellParticles3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    TriangularSurfaceMesh<T> const& elementaryMesh;
    std::vector<Array<T,3> > cellOrigins;
};


template<typename T, template<typename U> class Descriptor>
class RandomPositionCellParticlesForGrowth3D : public BoxProcessingFunctional3D
{
public:
    RandomPositionCellParticlesForGrowth3D (
            TriangularSurfaceMesh<T> const& elementaryMesh_, T hematocrit_, T & ratio_):
                elementaryMesh(elementaryMesh_), hematocrit(hematocrit_), ratio(ratio_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual RandomPositionCellParticlesForGrowth3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    TriangularSurfaceMesh<T> const& elementaryMesh;
    T hematocrit;
    T & ratio;
};



template<typename T, template<typename U> class Descriptor>
void randomPositionCellFieldsForGrowth3D(std::vector<CellField3D<T, Descriptor>* > & cellFields);


template<typename T, template<typename U> class Descriptor>
class RandomPositionCellFieldsForGrowth3D : public BoxProcessingFunctional3D
{
public:
    RandomPositionCellFieldsForGrowth3D (std::vector<CellField3D<T, Descriptor>* > & cellFields_):
                cellFields(cellFields_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual RandomPositionCellFieldsForGrowth3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<CellField3D<T, Descriptor>* > & cellFields;
};


template<typename T, template<typename U> class Descriptor>
class OrderedPositionCellField3D : public BoxProcessingFunctional3D
{
public:
    OrderedPositionCellField3D (std::vector<CellField3D<T, Descriptor>* > & cellFields_):
                cellFields(cellFields_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual OrderedPositionCellField3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<CellField3D<T, Descriptor>* > & cellFields;
};



template<typename T, template<typename U> class Descriptor>
class PrintParticles : public BoxProcessingFunctional3D
{
public:
    PrintParticles() { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual PrintParticles<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
};





#include "initializationCellField3D.hh"
#endif


