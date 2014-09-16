#ifndef COLLECTIVE_CELL_REDUCTIONS_H
#define COLLECTIVE_CELL_REDUCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellReductionTypes.h"
#include "cellField3D.hh"
#include "diagonalize.hpp"


/*
 *    LocalCellReductions3D is the Functional Reductive Box that does the actual job.
 *
 */

template<typename T, template<typename U> class Descriptor>
class LocalCellReductions3D : public BoxProcessingFunctional3D
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
    // quantityBins contain the "whichSum"/"whichAverage"/etc
    std::map<plint, std::map<plint, plint >  > quantityBins1D; // quantityBins1D[CCR_EDGE_DISTANCE_MEAN][cellId]
    std::map<plint, std::map<plint, Array<plint,3> >  > quantityBins3D; // quantityBins3D[CCR_VELOCITY_MEAN][cellId]
    std::map<plint, std::map<plint, std::vector<plint> >  > quantityBinsND; // quantityBinsND[CCR_INERTIA][cellId]

    T computeQuantity1D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
    Array<T,3> computeQuantity3D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
    std::vector<T> computeQuantityND (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
private:
    vector<T> sumV, averageV, maxV;
    vector<plint> averageQV;
    plint subscribeSum() { sumV.push_back(T()); return sumV.size() - 1; } ;
    plint subscribeAverage() { averageV.push_back(T()); averageQV.push_back(0); return averageV.size()  - 1; } ;
    plint subscribeMax() { maxV.push_back( -std::numeric_limits<double>::max() ); return maxV.size() - 1; } ;

    void gatherSum(plint qBin, T value) { sumV[qBin] += value; } ;
    void gatherAverage(plint qBin, T value) { averageV[qBin] += value; averageQV[qBin]+=1; } ;
    void gatherMax(plint qBin, T value) { maxV[qBin] = max(maxV[qBin], value); } ;

    T getSum(plint qBin) { return sumV[qBin]; } ;
    T getAverage(plint qBin) { return averageV[qBin]*1.0/averageQV[qBin]; } ;
    T getMax(plint qBin) { return maxV[qBin]; } ;

private: /* ReductiveBoxFunctions */
    plint subscribeReduction1D(plint reductionType);
    Array<plint,3> subscribeReduction3D(plint reductionType);
    std::vector<plint> subscribeReductionND(plint reductionType, plint dimensions);

    void gatherReduction1D(plint reductionType, plint whatQ, T value);
    void gatherReduction3D(plint reductionType, Array<plint,3> whatQ, Array<T,3> value);
    void gatherReductionND(plint reductionType, std::vector<plint> whatQ, std::vector<T> value);

    T getReduction1D(plint reductionType, plint whatQ);
    Array<T,3> getReduction3D(plint reductionType, Array<plint,3> whatQ);
    std::vector<T> getReductionND(plint reductionType, std::vector<plint> whatQ);
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
