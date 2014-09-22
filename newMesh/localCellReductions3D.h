#ifndef COLLECTIVE_CELL_REDUCTIONS_H
#define COLLECTIVE_CELL_REDUCTIONS_H

#include "diagonalize.hpp"
#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellReductionTypes.h"
#include "cellField3D.hh"
#include "cell3D.h"


/*
 *    LocalCellReductions3D is the Functional Reductive Box that does the actual job.
 *
 */

template<typename T, template<typename U> class Descriptor>
class LocalCellReductions3D : public BoxProcessingFunctional3D, BlockStatisticsForCellQuantityHolder<T>
{
public:
    LocalCellReductions3D(TriangleBoundary3D<T> const& triangleBoundary_,
            CellField3D<T,Descriptor> & chq_,
            plint numVerticesPerCell_,
            std::vector<plint> const& subscribedQuantities_);
    LocalCellReductions3D(LocalCellReductions3D<T,Descriptor> const& rhs) :
        triangleBoundary(rhs.triangleBoundary),
        chq(rhs.chq),
        numVerticesPerCell(rhs.numVerticesPerCell),
        subscribedQuantities(rhs.subscribedQuantities) { };
    virtual ~LocalCellReductions3D() {};
    virtual LocalCellReductions3D<T,Descriptor>* clone() const { return new LocalCellReductions3D<T,Descriptor>(*this); }
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing; // Particle field.
        modified[1] = modif::allVariables; // Reduction particles;
    }
    BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
private:
    void subscribeParticles(std::vector<Particle3D<T,Descriptor>*> const& particles);
    TriangleBoundary3D<T> const& triangleBoundary;
    CellField3D<T,Descriptor> & chq;
    plint numVerticesPerCell;
    std::vector<plint> subscribedQuantities;
    std::map<plint, plint> particlesPerCellId;
    std::vector<plint> cellIds;
    plint nCellIds;
private:
    T computeQuantity1D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
    Array<T,3> computeQuantity3D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
    std::vector<T> computeQuantityND (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
};



template<typename T, template<typename U> class Descriptor>
class SyncReductionParticles3D : public BoxProcessingFunctional3D
{
public:
    SyncReductionParticles3D (CellField3D<T,Descriptor> & chq_);
    virtual ~SyncReductionParticles3D() { };
    SyncReductionParticles3D(SyncReductionParticles3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual SyncReductionParticles3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    CellField3D<T,Descriptor> & chq;
    plint numVerticesPerCell;
};


/* ================================================================================ */
/* ******** computeEllipsoidFit *********************************** */
/* ================================================================================ */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeEllipsoidFit (std::vector<T> & cellInertia,
                          std::vector<T> & cellsEllipsoidFitAngles,
                          std::vector<T> & cellsEllipsoidFitSemiAxes,
                          T & difference,
                          T cellVolume);


#include "localCellReductions3D.hh"

#endif  // COLLECTIVE_CELL_REDUCTIONS_H
