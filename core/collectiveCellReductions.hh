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
    // Subscribe Quantities to quantityBins
    for (plint id = 0; id < subscribedQuantities.size(); ++id) { // quantityBins1D[CCR_EDGE_DISTANCE_MEAN][cellId]
        plint ccrId=subscribedQuantities[id];
        plint dim = ccrId%10; // Find dimension: last digit
        plint reductionType = (ccrId%100)/10; // Find reduction type (min,max,etc): second to last digit
        for (plint cellId = 0; cellId < maxNumberOfCells; ++cellId) {
            if (1 == dim) { quantityBins1D[ccrId][cellId] = subscribeReduction1D(reductionType); }
            if (3 == dim) { quantityBins3D[ccrId][cellId] = subscribeReduction3D(reductionType); }
            else          { quantityBinsND[ccrId][cellId] = subscribeReductionND(reductionType, dim); }
        }
    }
}

// * ID for the quantity of interest, starting from 1. Can be grouped together like:
// *      Volume         : 1 // 1D
// *      Angle          : 2 // 1D
// *      Area           : 3 // 1D
// *      Edge Distance  : 4 // 1D
// *      Edge Tile Span : 5 // 1D
// *      Position       : 6 // 3D
// *      Velocity       : 7 // 3D
// *      Inertia        : 8 // ND
// *      Energy         : 9 // 1D

template<typename T, template<typename U> class Descriptor>
T CollectiveCellReductionBox3D<T,Descriptor>::computeQuantity1D (plint q, ImmersedCellParticle3D<T,Descriptor>* particle) {
/*
 *    Calculates Volume, Angle, Area, Edge Distance and Edge Tile Span for each particle.
 *    Input:
 *      q is the ID for the quantity of interest
 *      particle corresponds to the particle of interest.
 */
    T quantity1D;
    plint iVertex = particle->getTag();
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
    if (q==2) {
    // Calculate ANGLE
        T edgeAngle = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            edgeAngle += calculateSignedAngle(triangleMesh, iVertex, neighbors[iB]);
        }
        quantity1D = edgeAngle = edgeAngle*1.0/neighbors.size();
    } else if (q==3) {
    // Calculate AREA
        triangleSurface = triangleMesh.computeVertexArea(iVertex);
    } else if (q==4) {
    // Calculate EDGE DISTANCE
        T edgeDistance = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            quantity1D = edgeDistance +=  triangleMesh.computeEdgeLength(iVertex, neighbors[iB]) ;
        }
        quantity1D = edgeDistance = edgeDistance*1.0/neighbors.size();
    } else if (q==5) {
    // Calculate EDGE TILE SPAN
        T edgeTileSpan = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            edgeTileSpan +=  triangleMesh.computeEdgeTileSpan(iVertex, neighbors[iB]) ;
        }
        quantity1D = edgeTileSpan = edgeTileSpan*1.0/neighbors.size();
    } else if (q==9) {
    // Return Energy of particle
        quantity1D = particle->get_Energy();
    } else if (q==1) {
        // Calculate VOLUME
            for (pluint iB = 0; iB < neighbors.size(); ++iB) {
                plint vId = triangleMesh.getVertexId(neighbors[iB],0);
                ImmersedCellParticle3D<T,Descriptor>* neigboringParticle =
                        dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D.find(vId)->second);
                Array<T,3> v0 = neigboringParticle->get_pbcPosition();
                vId = triangleMesh.getVertexId(neighbors[iB],1);
                neigboringParticle =
                        dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D.find(vId)->second);
                Array<T,3> v1 = neigboringParticle->get_pbcPosition();
                vId = triangleMesh.getVertexId(neighbors[iB],2);
                neigboringParticle =
                        dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (iVertexToParticle3D.find(vId)->second);
                Array<T,3> v2 = neigboringParticle->get_pbcPosition();
                /* Calculating the volume contibution of a face based on the formula:
                 * V[j] = 1.0/6.0 * (X3[j] cross X2[j])*X1[j]  */
                Array<T,3> tmp;
                crossProduct(v1, v2, tmp);
                quantity1D  =  VectorTemplate<T,Descriptor>::scalarProduct(v0,tmp); // * (1.0/6.0)
        }
    }
    return quantity1D;
}


template<typename T, template<typename U> class Descriptor>
Array<T,3> CollectiveCellReductionBox3D<T,Descriptor>::computeQuantityND (plint q, ImmersedCellParticle3D<T,Descriptor>* particle) {
    /*
     *    Calculates Position and Velocity for each particle.
     */
    Array<T,3> quantity3D;
    plint iVertex = particle->getTag();
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
    if (q==6) {
    // POSITION
        quantity3D = particle->getPosition();
    } else if (q==7) {
    // VELOCITY
        quantity3D = particle->get_v();
    }
    return quantity3D;
}


template<typename T, template<typename U> class Descriptor>
std::vector<T> CollectiveCellReductionBox3D<T,Descriptor>::computeQuantityND (plint q, ImmersedCellParticle3D<T,Descriptor>* particle) {
    /*
     *    Calculates Position and Velocity for each particle.
     */
    std::vector<T> quantityND;
    plint iVertex = particle->getTag();
    plit cellId = particle->get_cellId();
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
    if (q==8) {
    // INTERTIA
        T rx, ry, rz;
        T Ixx=0, Ixy=0, Ixz=0;
        T Iyx=0, Iyy=0, Iyz=0;
        T Izx=0, Izy=0, Izz=0;
        Array<T,3> r0 = *(carryOnQuantities3D[CCR_INERTIA][cellId]); // Get the Cell Center;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            T Aj = triangleMesh.computeTriangleArea(neighbors[iB]);
            Array<T,3> nj = triangleMesh.computeTriangleNormal(neighbors[iB]);
            Array<T,3> rj = (triangleMesh.getVertex(neighbors[iB],0) +
                            triangleMesh.getVertex(neighbors[iB],1) + triangleMesh.getVertex(neighbors[iB],2))/3.0;
            rj = rj - r0; // Subtract the Cell Center;
            T dV = 1.0/5.0 * Aj * dot(nj, rj);
            rx = rj[0]; ry = rj[1]; rz = rj[2];
            Ixx += (rz*rz+ry*ry)*dV;
            Iyy += (rz*rz+rx*rx)*dV;
            Izz += (rx*rx+ry*ry)*dV;
            Ixy += -rx*ry*dV;
            Iyx = Ixy;
            Ixz += -rx*rz*dV;
            Izx += Ixz;
            Iyz += -ry*rz*dV;
            Izy = Iyz;
        }
        quantityND.clear();
        quantityND.push_back(Ixx/3.0); // [0] every element is evaluated 3 times
        quantityND.push_back(Ixy/3.0); // [1] every element is evaluated 3 times
        quantityND.push_back(Ixz/3.0); // [2] every element is evaluated 3 times
        quantityND.push_back(Iyx/3.0); // [3] every element is evaluated 3 times
        quantityND.push_back(Iyy/3.0); // [4] every element is evaluated 3 times
        quantityND.push_back(Iyz/3.0); // [5] every element is evaluated 3 times
        quantityND.push_back(Izx/3.0); // [6] every element is evaluated 3 times
        quantityND.push_back(Izy/3.0); // [7] every element is evaluated 3 times
        quantityND.push_back(Izz/3.0); // [8] every element is evaluated 3 times
    }
    return quantityND;
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
    plint ccrId, dim, reductionType, cellId;
    std::map<plint, T > values1D;
    std::map<plint, Array<T,3> > values3D;
    std::map<plint, std::vector<T> > valuesND;

    for (pluint iA = 0; iA < particles.size(); ++iA) {
        values1D.clear(); values3D.clear(); valuesND.clear();
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iA]);
        cellId = particle->get_cellId();

        for (pluint i1d = 0; i1d < quantitiesToReduce.size(); ++i1d) {
            ccrId=quantitiesToReduce[id];
            dim = ccrId%10; // Dimension: last digit
            reductionType = (ccrId%100)/10; // Reduction type (min,max,etc): second to last digit
            quantity = ccrId/100; // Quantity of interest (area, velocity,etc): first digit
            if (dim == 1) {
                T value;
                if (values1D.count(quantity) > 0) {value = values1D[quantity];}
                else { values1D[quantity] = value = computeQuantity1D(quantity, particle); };
                plint qBin = quantityBins1D[ccrId][cellId];
                gatherReduction1D(reductionType, qBin, value);
            } else if (dim == 3) {
                Array<T,3> value;
                if (values3D.count(quantity) > 0) {value = values3D[quantity];}
                else { values3D[quantity] = value = computeQuantity3D(quantity, particle); };
                Array<plint,3> qBin = quantityBins3D[ccrId][cellId];
                gatherReduction3D(reductionType, qBin, value);
            } else {
                std::vector<T> value;
                if (valuesND.count(quantity) > 0) {value = valuesND[quantity];}
                else { valuesND[quantity] = value = computeQuantityND(quantity, particle); };
                std::vector<plint> qBin = quantityBinsND[ccrId][cellId];
                gatherReductionND(reductionType, qBin, value);
            }
        }
    }

}





/* **********************************************************/
/* Helper functions for the handling of the BlockStatistics */
/* **********************************************************/


/*******************************************************************
 *                      1D Subscribe Reductions                    *
 *******************************************************************/
template<typename T, template<typename U> class Descriptor>
plint CollectiveCellReductionBox3D<T,Descriptor>::subscribeReduction1D(plint reductionType) {
    if (0 == reductionType)      { return this->getStatistics().subscribeSum(); }
    else if (1 == reductionType) { return this->getStatistics().subscribeAverage(); }
    else if (2 == reductionType) { return this->getStatistics().subscribeMax(); }
    else if (3 == reductionType) { return this->getStatistics().subscribeMax(); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { return this->getStatistics().subscribeAverage(); } // Std is essentially an average
    else { return -1; }
}

/*******************************************************************
 *                         1D Gather Reductions                    *
 *******************************************************************/

template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::gatherReduction1D(plint reductionType, plint qBin, T value) {
    if (0 == reductionType)      { this->getStatistics().gatherSum(qBin, value); }
    else if (1 == reductionType) { this->getStatistics().gatherAverage(qBin, value); }
    else if (2 == reductionType) { this->getStatistics().gatherMax(qBin, value); }
    else if (3 == reductionType) { this->getStatistics().gatherMax(qBin, -value); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { this->getStatistics().gatherAverage(qBin, value); } // Std is essentially an average
}


/*******************************************************************
 *                          1D Get Reductions                      *
 *******************************************************************/

template<typename T, template<typename U> class Descriptor>
T CollectiveCellReductionBox3D<T,Descriptor>::getReduction1D(plint reductionType, plint qBin) {
    if (0 == reductionType)      { return this->getStatistics().getSum(qBin); }
    else if (1 == reductionType) { return this->getStatistics().getAverage(qBin); }
    else if (2 == reductionType) { return this->getStatistics().getMax(qBin); }
    else if (3 == reductionType) { return -this->getStatistics().getMax(qBin); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { return this->getStatistics().getAverage(qBin); } // Std is essentially an average
}



/*
 * 3D and ND Reduction operations
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
std::vector<plint> CollectiveCellReductionBox3D<T,Descriptor>::subscribeReductionND(plint reductionType, plint dimensions) {
    std::vector<plint> ret;
    for (i = 0; i < dimensions; ++i) { ret.push_back( subscribeReduction1D(reductionType) ); }
    return ret;
}

// gatherReduction
template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::gatherReduction3D(plint reductionType, Array<plint,3> qBin, Array<T,3> value) {
    for (int i = 0; i < 3; ++i) { gatherReduction1D(reductionType, qBin[i], value[i]); }
}

template<typename T, template<typename U> class Descriptor>
void CollectiveCellReductionBox3D<T,Descriptor>::gatherReductionND(plint reductionType, std::vector<plint> qBin, std::vector<T> value) {
    for (int i = 0; i < qBin.size(); ++i) { gatherReduction1D(reductionType, qBin[i], value[i]); }
}

// getReduction
template<typename T, template<typename U> class Descriptor>
Array<T,3> CollectiveCellReductionBox3D<T,Descriptor>::getReduction3D(plint reductionType, Array<plint,3> qBin) {
    T x = getReduction1D(reductionType, qBin[0]);
    T y = getReduction1D(reductionType, qBin[1]);
    T z = getReduction1D(reductionType, qBin[2]);
    return Array<T,3>(x, y, z);
}

template<typename T, template<typename U> class Descriptor>
std::vector<T> CollectiveCellReductionBox3D<T,Descriptor>::getReductionND(plint reductionType, std::vector<plint> qBin) {
    std::vector<T> ret;
    for (i = 0; i < dimensions; ++i) { ret.push_back( subscribeReduction1D(reductionType, qBin[i]) ); }
    return ret;
}





#endif  // COLLECTIVE_CELL_REDUCTIONS_HH
