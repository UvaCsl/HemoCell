#ifndef CELL_IN_SHEAR_FLOW_3D_HH
#define CELL_IN_SHEAR_FLOW_3D_HH

#include "cellInShearFlow3D.h"

namespace plb {

/* ******** InertiaTensorCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
InertiaTensorCellReduceFunctional3D<T,Descriptor>::InertiaTensorCellReduceFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_,
        std::vector< Array<T,3> >& cellCenters_)
    : triangleBoundary(triangleBoundary_)
{
    numberOfCells = cellIds_.size();
    for (pluint i=0; i< (pluint) (numberOfCells); ++i) {
        quantityIds[ cellIds_[i] ].resetToZero();
        cellCenters[ cellIds_[i] ] = cellCenters_[i];
        Array<plint,9> inertiaTensorId; inertiaTensorId.resetToZero();
        for (pluint j=0; j< (pluint) 9; ++j) {
            inertiaTensorId[j] = this->getStatistics().subscribeSum();
        }
        quantityIds[ cellIds_[i] ] = inertiaTensorId;
    }
}


template<typename T, template<typename U> class Descriptor>
InertiaTensorCellReduceFunctional3D<T,Descriptor>::InertiaTensorCellReduceFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
    : triangleBoundary(triangleBoundary_)
{
    numberOfCells = cellIds_.size();
    for (pluint i=0; i< (pluint) (numberOfCells); ++i) {
        quantityIds[ cellIds_[i] ].resetToZero();
        cellCenters[ cellIds_[i] ] = Array<T,3>(0.,0.,0.);
        Array<plint,9> inertiaTensorId; inertiaTensorId.resetToZero();
        for (pluint j=0; j< (pluint) 9; ++j) {
            inertiaTensorId[j] = this->getStatistics().subscribeSum();
        }
        quantityIds[ cellIds_[i] ] = inertiaTensorId;
    }
}


template<typename T, template<typename U> class Descriptor>
void InertiaTensorCellReduceFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    calculateQuantity(triangleMesh, particles, quantityIds);
}


template<typename T, template<typename U> class Descriptor>
void InertiaTensorCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
        std::vector<Particle3D<T,Descriptor>*> & particles, std::map<plint, Array<plint,9> > & quantityIds_)
{
    T rx, ry, rz;
    T Ixx, Ixy, Ixz;
    T Iyx, Iyy, Iyz;
    T Izx, Izy, Izz;
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        plint cellId = particle->get_cellId();
        Array<T,3> r0 = cellCenters[cellId];
        std::vector<plint> neighbors = triangleMesh.getNeighborTriangleIds(iVertex);
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            T Aj = triangleMesh.computeTriangleArea(neighbors[iB]);
            Array<T,3> nj = triangleMesh.computeTriangleNormal(neighbors[iB]);
            Array<T,3> rj = (triangleMesh.getVertex(neighbors[iB],0) +
                            triangleMesh.getVertex(neighbors[iB],1) + triangleMesh.getVertex(neighbors[iB],2))/3.0;
            rj = rj - r0;
            T dV = 1.0/5.0 * Aj * dot(nj, rj);
            rx = rj[0]; ry = rj[1]; rz = rj[2];
            Ixx = (rz*rz+ry*ry)*dV;
            Iyy = (rz*rz+rx*rx)*dV;
            Izz = (rx*rx+ry*ry)*dV;
            Iyx = Ixy = -rx*ry*dV;
            Izx = Ixz = -rx*rz*dV;
            Izy = Iyz = -ry*rz*dV;


            this->getStatistics().gatherSum(quantityIds_[cellId][0], Ixx/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][1], Ixy/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][2], Ixz/3.0); // every element is evaluated 3 times

            this->getStatistics().gatherSum(quantityIds_[cellId][3], Iyx/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][4], Iyy/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][5], Iyz/3.0); // every element is evaluated 3 times

            this->getStatistics().gatherSum(quantityIds_[cellId][6], Izx/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][7], Izy/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][8], Izz/3.0); // every element is evaluated 3 times
        }
    }
}


template<typename T, template<typename U> class Descriptor>
void InertiaTensorCellReduceFunctional3D<T,Descriptor>::getCellQuantityArray(std::vector< std::vector<T> > & cellQuantity, std::vector<plint> cellIds_) const {
    cellQuantity.clear();
    for (pluint i=0; i< (pluint) (cellIds_.size()); ++i) {
        std::vector<T> tmp(9);
        Array<plint,9> inertiaTensorId = this->quantityIds.find(cellIds_[i])->second;
        for (pluint j=0; j < 9; ++j) {
            tmp[j] = this->getStatistics().getSum( inertiaTensorId[j] );
        }
        cellQuantity.push_back(tmp);
    }
}


/* ******** computeCellInertia *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeCellInertia (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,3> >& cellCenters, std::vector< std::vector<T> >& cellInertia) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    InertiaTensorCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds, cellCenters);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellInertia, cellIds);
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeCellInertia (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,9> >& cellInertia) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    InertiaTensorCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellInertia, cellIds);
}


/* ******** computeCellInertia *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeEllipsoidFit (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
                std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume,
                std::vector< std::vector<T> > & cellsEllipsoidFitAngles,
                std::vector< std::vector<T> > & cellsEllipsoidFitSemiAxes) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    std::vector< std::vector<T> > cellInertia;
    InertiaTensorCellReduceFunctional3D<T,Descriptor> inertiaTensorFunctional(Cells, cellIds, cellCenters);
    applyProcessingFunctional(inertiaTensorFunctional, particles.getBoundingBox(), particleArg);
    inertiaTensorFunctional.getCellQuantityArray(cellInertia, cellIds);

    for (pluint i = 0; i < cellIds.size(); ++i) {
        vector<T> semiAxes, ellipsoidAngles;
        getLambdasAndAngles(cellInertia[i], semiAxes, ellipsoidAngles);

        T f1 = semiAxes[0] * 5 / cellsVolume[i];
        T f2 = semiAxes[1] * 5 / cellsVolume[i];
        T f3 = semiAxes[2] * 5 / cellsVolume[i];
        semiAxes[0] = (sqrt( (f1 + f2 - f3)/2 ));
        semiAxes[1] = (sqrt( (f2 + f3 - f1)/2 ));
        semiAxes[2] = (sqrt( (f3 + f1 - f2)/2 ));
        std::sort(semiAxes.begin(), semiAxes.end());
        cellsEllipsoidFitSemiAxes.push_back(semiAxes);
        cellsEllipsoidFitAngles.push_back(ellipsoidAngles);
    }
}


}
#endif // CELL_IN_SHEAR_FLOW_3D_HH
