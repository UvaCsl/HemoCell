#ifndef COLLECTIVE_CELL_REDUCTIONS_H
#define COLLECTIVE_CELL_REDUCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"


/* IDs for the reduction types,
 * they are created based on the following rules:
 * ==============================
 * Last digit, Dimension:
 *      1D : 1
 *      2D : 2
 *      3D : 3
 *      ND : 4,5,6,7,8,9 // Not above 10
 * Second to last digit, type of reduction:
 *      Sum : 0
 *      Mean: 1
 *      Min : 2
 *      Max : 3
 *      STD : 4 // Still to be implemented
 * ID for the quantity of interest, starting from 1. Can be grouped together like:
 *      Volume         : 1
 *      Angle          : 2
 *      Area           : 3
 *      Edge Distance  : 4
 *      Edge Tile Span : 5
 *      Position       : 6 // Periodic boundary position
 *      Velocity       : 7
 *      Inertia        : 8
 *      Energy         : 9
 *      Position       : 0
 */
#define CCR_NO_PBC_POSITION_MEAN    13 // 3d // BEWARE OF 0 IN FRONT! GOES TO OCT
#define CCR_NO_PBC_POSITION_MIN     23 // 3d // BEWARE OF 0 IN FRONT! GOES TO OCT
#define CCR_NO_PBC_POSITION_MAX     33 // 3d // BEWARE OF 0 IN FRONT! GOES TO OCT
#define CCR_VOLUME                 101 // 1d
#define CCR_ANGLE_MEAN             211 // 1d
#define CCR_ANGLE_MIN              221 // 1d
#define CCR_ANGLE_MAX              231 // 1d
#define CCR_SURFACE                301 // 1d
#define CCR_TRIANGLE_AREA_MEAN     311 // 1d, Better use 301 and divide by the number of triangles
#define CCR_TRIANGLE_AREA_MIN      321 // 1d
#define CCR_TRIANGLE_AREA_MAX      331 // 1d
#define CCR_EDGE_DISTANCE_MEAN     411 // 1d
#define CCR_EDGE_DISTANCE_MIN      421 // 1d
#define CCR_EDGE_DISTANCE_MAX      431 // 1d
#define CCR_TILE_SPAN_MEAN         511 // 1d
#define CCR_TILE_SPAN_MIN          521 // 1d
#define CCR_TILE_SPAN_MAX          531 // 1d
#define CCR_POSITION_MEAN          613 // 3d
#define CCR_POSITION_MIN           623 // 3d
#define CCR_POSITION_MAX           633 // 3d
#define CCR_VELOCITY_MEAN          713 // 3d
#define CCR_VELOCITY_MIN           723 // 3d
#define CCR_VELOCITY_MAX           733 // 3d
#define CCR_INERTIA                809 // 9d, Not working
#define CCR_ENERGY                 901 // 1d
// #define CCR_MAX               41 // 9d


/*
 *    CollectiveCellReductionBox3D is the Functional Reductive Box that does the actual job.
 *
 */
template<typename T, template<typename U> class Descriptor>
class CollectiveCellReductionBox3D : public BoxProcessingFunctional3D
{
public:
    CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint numVerticesPerCell_,
            std::vector<plint> const& subscribedQuantities_);
    CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint numVerticesPerCell_,
            plint subscribedQuantities_);
    CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint numVerticesPerCell_);
    /// Argument: Particle-field.
    virtual ~CollectiveCellReductionBox3D() {};
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);

    virtual CollectiveCellReductionBox3D<T,Descriptor>* clone() const { return new CollectiveCellReductionBox3D<T,Descriptor>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing; // Particle field.
        modified[1] = modif::dynamicVariables; // Reduction particles;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = false;
        isWritten[1] = true;
    }
    BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    void postProcess(std::vector<MultiBlock3D*> blocks);

private:
    void subscribeParticles(std::vector<Particle3D<T,Descriptor>*> const& particles);

    TriangleBoundary3D<T> const& triangleBoundary;
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

    std::map<plint, std::map<plint, T > * > quantitiesFromReductionParticles1D; // carryOnQuantities1D[cellId][CCR_EDGE_DISTANCE_STD]= MEAN_EDGE_DISTANCE
    std::map<plint, std::map<plint, Array<T,3> > *  > quantitiesFromReductionParticles3D;
    std::map<plint, std::map<plint, std::vector<T> > * > quantitiesFromReductionParticlesND;


    T computeQuantity1D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
    Array<T,3> computeQuantity3D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
    std::vector<T> computeQuantityND (plint q, ImmersedCellParticle3D<T,Descriptor>* particle);
private:
    vector<T> sumV, averageV, maxV;
    vector<plint> averageQV;
    plint subscribeSum() { sumV.push_back(T()); return sumV.size(); } ;
    plint subscribeAverage() { averageV.push_back(T()); averageQV.push_back(0); return averageV.size(); } ;
    plint subscribeMax() { maxV.push_back( -std::numeric_limits<double>::max() ); return maxV.size(); } ;

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





#include "collectiveCellReductions3D.hh"

#endif  // COLLECTIVE_CELL_REDUCTIONS_H
