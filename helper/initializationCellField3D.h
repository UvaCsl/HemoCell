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






#include "initializationCellField3D.hh"
#endif


