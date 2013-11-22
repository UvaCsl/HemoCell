#ifndef COLLECTIVE_CELL_REDUCTIONS_HH
#define COLLECTIVE_CELL_REDUCTIONS_HH

#include "collectiveCellReductions.h"


template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
            std::vector<plint> subscribedQuantities_) :
            triangleBoundary(triangleBoundary_),
            maxNumberOfCells(maxNumberOfCells_), numVerticesPerCell(numVerticesPerCell_), numTrianglesPerCell(numTrianglesPerCell_),
            subscribedQuantities(subscribedQuantities_)
{
    // Subscribe Quantities to whatQuantity
    for (plint id = 0; id < subscribedQuantities.size(); ++id) { // whatQuantity1D[CCR_EDGE_DISTANCE_MEAN][cellId]
        plint ccrId=subscribedQuantities[id];
        plint dim = ccrId%10; // Find dimension: last digit
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        for (plint cellId = 0; cellId < maxNumberOfCells; ++cellId) {
            if (1 == dim) { whatQuantity1D[ccrId][cellId] = subscribeReduction1D(reductionType); }
            if (3 == dim) { whatQuantity3D[ccrId][cellId] = subscribeReduction3D(reductionType); }
            else          { whatQuantityXD[ccrId][cellId] = subscribeReductionXD(reductionType, dim); }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
//    this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], edgeAngle); // every edgeAngle is evaluated 2 times


}


/* Helper functions for the handling of the BlockStatistics */
/*
 * Subscribe Reductions
 */
template<typename T, template<typename U> class Descriptor>
plint CollectiveCellReductionBox3D<T,Descriptor>::subscribeReduction1D(plint reductionType) {
    if (0 == reductionType)      { return this->getStatistics().subscribeSum(); }
    else if (1 == reductionType) { return this->getStatistics().subscribeAverage(); }
    else if (2 == reductionType) { return this->getStatistics().subscribeMax(); }
    else if (3 == reductionType) { return this->getStatistics().subscribeMax(); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { return this->getStatistics().subscribeAverage(); } // Std is essentially an average
    else { return -1; }
}

/*
 * 1D Gather Reductions
 */
template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::gatherReduction1D(plint reductionType, plint whatQ, T value) {
    if (0 == reductionType)      { this->getStatistics().gatherSum(whatQ, value); }
    else if (1 == reductionType) { this->getStatistics().gatherAverage(whatQ, value); }
    else if (2 == reductionType) { this->getStatistics().gatherMax(whatQ, value); }
    else if (3 == reductionType) { this->getStatistics().gatherMax(whatQ, -value); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { this->getStatistics().gatherAverage(whatQ, value); } // Std is essentially an average
}

/*
 * 1D Get Reductions
 */
template<typename T, template<typename U> class Descriptor>
T CollectiveCellReductionBox3D<T,Descriptor>::getReduction1D(plint reductionType, plint whatQ) {
    if (0 == reductionType)      { return this->getStatistics().getSum(whatQ); }
    else if (1 == reductionType) { return this->getStatistics().getAverage(whatQ); }
    else if (2 == reductionType) { return this->getStatistics().getMax(whatQ); }
    else if (3 == reductionType) { return -this->getStatistics().getMax(whatQ); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { return this->getStatistics().getAverage(whatQ); } // Std is essentially an average
}



/*
 * 3D and XD Reduction operations
 */
// subscribeReduction
template<typename T, template<typename U> class Descriptor>
Array<plint,3> CollectiveCellReductionBox3D<T,Descriptor>::subscribeReduction3D(plint reductionType) {
    plint x = subscribeReduction1D(reductionType);
    plint y = subscribeReduction1D(reductionType);
    plint z = subscribeReduction1D(reductionType);
    return Array<plint,3>(x, y, z);
}

template<typename T, template<typename U> class Descriptor>
std::vector<plint> CollectiveCellReductionBox3D<T,Descriptor>::subscribeReductionXD(plint reductionType, plint dimensions) {
    std::vector<plint> ret;
    for (i = 0; i < dimensions; ++i) { ret.push_back( subscribeReduction1D(reductionType) ); }
    return ret;
}

// gatherReduction
template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::gatherReduction3D(plint reductionType, Array<plint,3> whatQ, Array<T,3> value) {
    for (int i = 0; i < 3; ++i) { gatherReduction1D(reductionType, whatQ[i], value[i]); }
}

template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::gatherReductionXD(plint reductionType, std::vector<plint> whatQ, std::vector<T> value) {
    for (int i = 0; i < whatQ.size(); ++i) { gatherReduction1D(reductionType, whatQ[i], value[i]); }
}

// getReduction
template<typename T, template<typename U> class Descriptor>
Array<T,3> CollectiveCellReductionBox3D<T,Descriptor>::getReduction3D(plint reductionType, Array<plint,3> whatQ) {
    T x = getReduction1D(reductionType, whatQ[0]);
    T y = getReduction1D(reductionType, whatQ[1]);
    T z = getReduction1D(reductionType, whatQ[2]);
    return Array<T,3>(x, y, z);
}

template<typename T, template<typename U> class Descriptor>
std::vector<T> CollectiveCellReductionBox3D<T,Descriptor>::getReductionXD(plint reductionType, std::vector<plint> whatQ) {
    std::vector<T> ret;
    for (i = 0; i < dimensions; ++i) { ret.push_back( subscribeReduction1D(reductionType, whatQ[i]) ); }
    return ret;
}





#endif  // COLLECTIVE_CELL_REDUCTIONS_HH
