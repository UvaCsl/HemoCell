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
// Mesh is supposed to start at (0,0,0)
template<typename T>
void meshRandomRotation (TriangularSurfaceMesh<T> * mesh, Array<T,3> randomAngles=Array<T,3>(-1,-1,-1) ) {
    if (randomAngles[0]<0) {
        randomAngles[0] = guessRandomNumber()* 2 * pi;
        randomAngles[1] = guessRandomNumber()* pi;
        randomAngles[2] = guessRandomNumber()* 2 * pi;
    }

    const T pi = 4.*atan(1.);
    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    mesh->translate(-1.0 * meshCenter);
    mesh->rotate(randomAngles[0], randomAngles[1], randomAngles[2]);

    mesh->computeBoundingBox (xRange, yRange, zRange);
    mesh->translate(Array<T,3>(-xRange[0], -yRange[0], -zRange[0]));
}

// Rotates the mesh in a given angle around its centerpoint.
// Mesh is supposed to start at (0,0,0), however, this function will not change the center position!
template<typename T>
void meshRotation (TriangularSurfaceMesh<T> * mesh, Array<T,3> rotationAngles) {
    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    mesh->translate(-1.0 * meshCenter);
    mesh->rotateXYZ(rotationAngles[0], rotationAngles[1], rotationAngles[2]);
    mesh->translate(meshCenter);
}

// For positioning that uses meshes starting at (0,0,0)
template<typename T>
void meshPositionToOrigin (TriangularSurfaceMesh<T> * mesh) {
    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);    
    mesh->translate(Array<T,3>(-xRange[0], -yRange[0], -zRange[0]));
}

template<typename T, template<typename U> class Descriptor>
void positionCellInParticleField(ParticleField3D<T,Descriptor>& particleField, BlockLattice3D<T,Descriptor>& fluid,
                                            TriangularSurfaceMesh<T> * mesh, Array<T,3> startingPoint, plint cellId) {
    plint nVertices=mesh->getNumVertices();
    Dot3D relativeDisplacement = computeRelativeDisplacement(fluid, particleField);
    //Dot3D fluidLocationDot3D = fluid.getLocation();
    plint iX, iY, iZ;
    plint maxNx = fluid.getNx(), maxNy = fluid.getNy(), maxNz = fluid.getNz();
    for (plint iVertex=0; iVertex < nVertices; ++iVertex) {
        Array<T,3> vertex = startingPoint + mesh->getVertex(iVertex);
        particleField.computeGridPosition (vertex-0.5, iX, iY, iZ);
        Dot3D fluidDomainCell = Dot3D(iX, iY, iZ) - relativeDisplacement;

        bool isInsideFluidDomain = (fluidDomainCell.x >= 0 and fluidDomainCell.y >= 0 and fluidDomainCell.z >= 0) and
        		(fluidDomainCell.x < maxNx and fluidDomainCell.y < maxNy and fluidDomainCell.z < maxNz);

        // Test if particle is inside and in a boundary -> dont add this particle
        if (isInsideFluidDomain and fluid.get(fluidDomainCell.x, fluidDomainCell.y, fluidDomainCell.z).getDynamics().isBoundary())
        	break;
        
        // If all neighbours are boundaries or denied cells
        bool neighboringBoundariesAnywhere = false;  

        // Deny particles that are in the outer most layer, aka. the "shear layer"
        int denyLayerSize = 1; // Size of the outer shear layer to deny particles from (= create a starting CFL). This should scale with dx and be <= 1um.    
        for (int px = -denyLayerSize; px <= denyLayerSize; ++px) {  for (int py = -denyLayerSize; py <= denyLayerSize; ++py) { for (int pz = -denyLayerSize; pz <= denyLayerSize; ++pz) {
                    bool isInsideDomain = (fluidDomainCell.x+px >= 0 and fluidDomainCell.y+py >= 0 and fluidDomainCell.z+pz >= 0) and
                        (fluidDomainCell.x+px < maxNx and fluidDomainCell.y+py < maxNy and fluidDomainCell.z+pz < maxNz);
                    if(isInsideDomain) {
                        neighboringBoundariesAnywhere = neighboringBoundariesAnywhere or fluid.get(fluidDomainCell.x+px, fluidDomainCell.y+py, fluidDomainCell.z+pz).getDynamics().isBoundary();
                    }
                }  
            }  
        }
        
        // This cell does not satisfy all requirements.
        if(neighboringBoundariesAnywhere)
            break; 

        // Finally, if all checks are passed, add the particle.
        particleField.addParticle(particleField.getBoundingBox(), new SurfaceParticle3D<T,Descriptor>(vertex, cellId, iVertex) );

    }
}


template<typename T, template<typename U> class Descriptor>
plint deleteIncompleteCells(ParticleField3D<T,Descriptor>& particleField, BlockLattice3D<T,Descriptor>& fluid, Box3D domain, plint nVertices) {
    // DELETE CELLS THAT ARE NOT WHOLE
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::map<plint, plint> cellIdToNumberOfVertices;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        SurfaceParticle3D<T,Descriptor>* particle =
            dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (particles[iP]);
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
    OrderedPositionCellField3D (std::vector<CellField3D<T, Descriptor>* > & cellFields_, Dot3D latticeSize_):
                cellFields(cellFields_), latticeSize(latticeSize_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual OrderedPositionCellField3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<CellField3D<T, Descriptor>* > & cellFields;
    Dot3D latticeSize;
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


