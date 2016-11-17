#ifndef CELL_FIELD_FUNCTIONALS_3D_H
#define CELL_FIELD_FUNCTIONALS_3D_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedBoundaryMethod3D.h"
#include "surfaceParticle3D.h"
#include "cellParticle3D.h"
#include "cell3D.h"
#include "meshMetrics.h"
#include "diagonalize.hpp"
#include <stdlib.h>     /* srand, rand */


using namespace std;
using namespace plb;



template<typename T, template<typename U> class Descriptor>
class ChangeParticleUpdateScheme : public BoxProcessingFunctional3D
{
public:
    ChangeParticleUpdateScheme (plint scheme_, T cellTimeStep = 1.0) : scheme(scheme_), dt(cellTimeStep) { } ;
    ~ChangeParticleUpdateScheme() { };
    ChangeParticleUpdateScheme(ChangeParticleUpdateScheme<T,Descriptor> const& rhs) :
        scheme(rhs.scheme), dt(rhs.dt) { } ;
    /// Arguments: [0] Particle-field
    virtual ChangeParticleUpdateScheme<T,Descriptor>* clone() const { return new ChangeParticleUpdateScheme<T,Descriptor>(*this); };
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; } ;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::allVariables; } ;
    void getModificationPattern(std::vector<bool>& isWritten) const { isWritten[0] = true; };
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
        PLB_PRECONDITION( blocks.size()==1 );
        ParticleField3D<T,Descriptor>& particleField =
            *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
        std::vector<Particle3D<T,Descriptor>*> particles;
        particleField.findParticles(domain, particles); // Gets particle only from the bulk
        for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
            //castParticleToICP3D(particles[iParticle])->get_scheme() = scheme;
            //castParticleToICP3D(particles[iParticle])->get_dt() = dt;
            SurfaceParticle3D<T,Descriptor>* particle = dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (particles[iParticle]);
            particle->get_scheme() = scheme;
            particle->get_dt() = dt;
        }
    };
private:
    T dt;
    plint scheme;
};



template<typename T, template<typename U> class Descriptor>
class SyncParticleFieldEnvelope : public BoxProcessingFunctional3D
{
public:
    // Does nothing but signals Palabos, that the field has been modified
    SyncParticleFieldEnvelope ();
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    SyncParticleFieldEnvelope(SyncParticleFieldEnvelope<T,Descriptor> const& rhs);
    virtual SyncParticleFieldEnvelope<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, template<typename U> class Descriptor>
class SyncParticleFieldBulk : public BoxProcessingFunctional3D
{
public:
    // Does nothing but signals Palabos, that the field has been modified
    SyncParticleFieldBulk ();
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    SyncParticleFieldBulk(SyncParticleFieldBulk<T,Descriptor> const& rhs);
    virtual SyncParticleFieldBulk<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};



template<typename T, template<typename U> class Descriptor>
class ForceToFluid3D : public BoxProcessingFunctional3D
{
public:
    // kernelSize gives the envelope size in each direction -> kernelSize=1 will yield a 3x3x3 kernel
    ForceToFluid3D (plint kernelSize = 1, plint ibmKernel_= 2 );
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ForceToFluid3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint ibmKernel;
    plint kernelSize;
};


template<typename T, template<typename U> class Descriptor>
class FluidVelocityToImmersedCell3D : public BoxProcessingFunctional3D
{
public:
    // kernelSize gives the envelope size in each direction -> kernelSize=1 will yield a 3x3x3 kernel
    FluidVelocityToImmersedCell3D (plint kernelSize = 1, plint ibmKernel_= 2 );
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FluidVelocityToImmersedCell3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint ibmKernel;
    plint kernelSize;
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

/* ******** GetGlobalMaxForce *********************************** */
template<typename T, template<typename U> class Descriptor>
class GetGlobalMaxForce : public PlainReductiveBoxProcessingFunctional3D
{
public:
    GetGlobalMaxForce();
    /// Argument: Particle-field.
    virtual GetGlobalMaxForce<T,Descriptor>* clone() const { return new GetGlobalMaxForce<T,Descriptor>(*this); };
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; };
    virtual T getValue();
    virtual BlockDomain::DomainT appliesTo() const;
protected:
    plint qId;
};

/* ******** ComputeRequiredQuantities *********************************** */
template<typename T, template<typename U> class Descriptor>
class DeleteIncompleteCells : public BoxProcessingFunctional3D
{
public:
    DeleteIncompleteCells (std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D_, plint numberOfVertices_=0);
    virtual ~DeleteIncompleteCells() { };
    DeleteIncompleteCells(DeleteIncompleteCells<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual DeleteIncompleteCells<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    std::map<plint, Cell3D<T,Descriptor>* > & cellIdToCell3D;
    plint numberOfVertices;
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

