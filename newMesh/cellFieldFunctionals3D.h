#ifndef CELL_FIELD_FUNCTIONALS_3D_H
#define CELL_FIELD_FUNCTIONALS_3D_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedBoundaryMethod3D.h"
#include "immersedCellParticle3D.h"
#include "reductionParticle3D.h"
#include "cell3D.h"
#include "meshMetrics.h"
#include "diagonalize.hpp"
#include <stdlib.h>     /* srand, rand */


using namespace std;
using namespace plb;


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
class ForceToFluid3D : public BoxProcessingFunctional3D
{
public:
    ForceToFluid3D (plint ibmKernel_=4 );
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ForceToFluid3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint ibmKernel;
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
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint ibmKernel;
};

template<typename T, template<typename U> class Descriptor>
class ViscousPositionUpdate3D : public BoxProcessingFunctional3D
{
public:
    ViscousPositionUpdate3D (T ratio_=1.0);
    ViscousPositionUpdate3D (ViscousPositionUpdate3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ViscousPositionUpdate3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T ratio;
};


template<typename T, template<typename U> class Descriptor>
class ComputeCellForce3D : public BoxProcessingFunctional3D
{
public:
    ComputeCellForce3D (ConstitutiveModel<T,Descriptor>* cellModel_, std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D_) ;
    ~ComputeCellForce3D() {
//        std::cout <<" ~ComputeCellForce3D() " << global::mpi().getRank() << std::endl;
    }; // { delete cellModel; } ;
    ComputeCellForce3D(ComputeCellForce3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeCellForce3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    ConstitutiveModel<T,Descriptor>* cellModel;
    std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D;
};


template<typename T, template<typename U> class Descriptor>
class FillCellMap : public BoxProcessingFunctional3D
{
public:
    FillCellMap (TriangularSurfaceMesh<T> & mesh_, std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D_);
    virtual ~FillCellMap();
    FillCellMap(FillCellMap<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FillCellMap<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    TriangularSurfaceMesh<T> & mesh;
	std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D;
};


/* ******** CountGlobalNumberOfCells *********************************** */
template<typename T, template<typename U> class Descriptor>
class CountGlobalNumberOfCells : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CountGlobalNumberOfCells(plint numVertices_);
    /// Argument: Particle-field.
    virtual CountGlobalNumberOfCells<T,Descriptor>* clone() const { return new CountGlobalNumberOfCells<T,Descriptor>(*this); };
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; };
    virtual pluint getValue();
    virtual BlockDomain::DomainT appliesTo() const;
protected:
    plint numVertices;
    plint qId;
};


/* ******** ComputeRequiredQuantities *********************************** */
template<typename T, template<typename U> class Descriptor>
class ComputeRequiredQuantities : public BoxProcessingFunctional3D
{
public:
    ComputeRequiredQuantities (std::vector<plint> ccrRequirements_, std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D_);
    virtual ~ComputeRequiredQuantities();
    ComputeRequiredQuantities(ComputeRequiredQuantities<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeRequiredQuantities<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    std::vector<plint> ccrRequirements;
    std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D;
};


template<typename T, template<typename U> class Descriptor>
class SyncCellQuantities : public BoxProcessingFunctional3D
{
public:
    SyncCellQuantities (std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D_);
    virtual ~SyncCellQuantities() { };
    SyncCellQuantities(SyncCellQuantities<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual SyncCellQuantities<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D;
};

/* ================================================================================ */
/* ******** computeEllipsoidFit *********************************** */
/* ================================================================================ */
template<typename T>
void computeEllipsoidFit (std::vector<T> & cellInertia,
                          std::vector<T> & cellsEllipsoidFitAngles,
                          std::vector<T> & cellsEllipsoidFitSemiAxes,
                          T & difference,
                          T cellVolume);


#include "cellFieldFunctionals3D.hh"
#endif  // CELL_FIELD_FUNCTIONALS_3D_H

