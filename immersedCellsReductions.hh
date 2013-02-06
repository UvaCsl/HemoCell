#ifndef IMMERSED_CELLS_REDUCTIONS_HH
#define IMMERSED_CELLS_REDUCTIONS_HH

#include "immersedCellsReductions.h"

#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

/* ******** CellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
CellReduceFunctional3D<T,Descriptor>::CellReduceFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_, bool findMax_)
    : triangleBoundary(triangleBoundary_), findMax(findMax_)
{
    numberOfCells = 0;
    for (pluint var = 0; var < cellIds_.size(); ++var)
        if (cellIds_[var] > numberOfCells)
            numberOfCells = cellIds_[var];
    if (findMax) {
        for (pluint i=0; i< (pluint) numberOfCells+1; ++i)
            quantityIds.push_back(this->getStatistics().subscribeMax());
    } else {
        for (pluint i=0; i< (pluint) numberOfCells+1; ++i)
            quantityIds.push_back(this->getStatistics().subscribeSum());
    }
}

template<typename T, template<typename U> class Descriptor>
void CellReduceFunctional3D<T,Descriptor>::processGenericBlocks (
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
void CellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{ }

template<typename T, template<typename U> class Descriptor>
void CellReduceFunctional3D<T,Descriptor>::getCellQuantityArray(std::vector<T>& cellQuantity, std::vector<plint> cellIds) const {
    if (findMax) {
        for (pluint i = 0; i < (pluint) quantityIds.size(); ++i)
            cellQuantity.push_back(this->getStatistics().getMax(quantityIds[i]));
    } else {
        for (pluint i = 0; i < (pluint) quantityIds.size(); ++i)
            cellQuantity.push_back(this->getStatistics().getSum(quantityIds[i]));
    }
}

template<typename T, template<typename U> class Descriptor>
CellReduceFunctional3D<T,Descriptor>* CellReduceFunctional3D<T,Descriptor>::clone() const {
    return new CellReduceFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CellReduceFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}



/* ******** countCellVolume ** countCellSurface *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVolume (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVolumes) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    VolumeCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellVolumes, cellIds);
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellSurface (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellSurfaces) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    SurfaceCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellSurfaces, cellIds);
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanTriangleArea (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanTriangleArea) //Perhaps add TAGS
{
    // We assume homogeneous Cells
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    SurfaceCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellMeanTriangleArea, cellIds);
    std::vector<T> cellNumTriangles;
    NumTrianglesCellReduceFunctional3D<T,Descriptor> nfunctional(Cells, cellIds);
    applyProcessingFunctional(nfunctional, domain, particleArg);
    nfunctional.getCellQuantityArray(cellNumTriangles, cellIds);
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
            cellMeanTriangleArea[iA] /= cellNumTriangles[iA];
    }
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanEdgeDistance (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanEdgeDistance) //Perhaps add TAGS
{
    // We assume homogeneous Cells
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    EdgeDistanceCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellMeanEdgeDistance, cellIds);
    std::vector<T> cellNumVertices;
    NumVerticesCellReduceFunctional3D<T,Descriptor> nfunctional(Cells, cellIds);
    applyProcessingFunctional(nfunctional, domain, particleArg);
    nfunctional.getCellQuantityArray(cellNumVertices, cellIds);
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        cellMeanEdgeDistance[iA] /= cellNumVertices[iA];
    }
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanTileSpan (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanTileSpan) //Perhaps add TAGS
{
    // We assume homogeneous Cells
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    TileSpanCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellMeanTileSpan, cellIds);
    std::vector<T> cellNumVertices;
    NumVerticesCellReduceFunctional3D<T,Descriptor> nfunctional(Cells, cellIds);
    applyProcessingFunctional(nfunctional, domain, particleArg);
    nfunctional.getCellQuantityArray(cellNumVertices, cellIds);
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
    	cellMeanTileSpan[iA] /= cellNumVertices[iA];
    }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMaxEdgeDistance (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMaxEdgeDistance) //Perhaps add TAGS
{
    // We assume homogeneous Cells
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    MaxEdgeDistanceCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellMaxEdgeDistance, cellIds);
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanAngle (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanAngle) //Perhaps add TAGS
{
    // We assume homogeneous Cells
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    AngleCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellMeanAngle, cellIds);
    std::vector<T> cellNumVertices;
    NumVerticesCellReduceFunctional3D<T,Descriptor> nfunctional(Cells, cellIds);
    applyProcessingFunctional(nfunctional, domain, particleArg);
    nfunctional.getCellQuantityArray(cellNumVertices, cellIds);
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        cellMeanAngle[iA] /= cellNumVertices[iA];
    }
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVertices (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVertices) //Perhaps add TAGS
{
    // We assume homogeneous Cells
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    AngleCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellVertices, cellIds);
}



/* ******** VolumeCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void VolumeCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        std::vector<plint> neighbors = triangleMesh.getNeighborTriangleIds(iVertex);
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            Array<T,3> v0 = triangleMesh.getVertex(neighbors[iB],0);
            Array<T,3> v1 = triangleMesh.getVertex(neighbors[iB],1);
            Array<T,3> v2 = triangleMesh.getVertex(neighbors[iB],2);
//       /* ************* Other Calculation ********* */
//          Array<T,3> areaTimesNormal = triangleMesh.computeTriangleNormal(neighbors[iB], true);
//          T triangleVolumeT6 = VectorTemplate<T,Descriptor>::scalarProduct(areaTimesNormal, ((v0+v1+v2)/3.0)) ;
//          this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], triangleVolumeT6/6.0/3.0); // every volume is evaluated 3 times
//       /* ********************************************* */
            /* Calculating the volume contibution of a face based on the formula:
             * V[j] = 1.0/6.0 * (X3[j] cross X2[j])*X1[j]  */
            Array<T,3> tmp;
            crossProduct(v1, v2, tmp);
            T triangleVolumeT6 =  VectorTemplate<T,Descriptor>::scalarProduct(v0,tmp); // * (1.0/6.0)
            this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], triangleVolumeT6/6.0/3.0); // every volume is evaluated 3 times
        }
    }
}


/* ******** SurfaceCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void SurfaceCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        T triangleSurface = triangleMesh.computeVertexArea(iVertex); // .computeVertexArea returns the amount per Vertex
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], triangleSurface);
    }
}


/* ******** EdgeDistanceCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void EdgeDistanceCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
        T edgeDistance = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            edgeDistance +=  triangleMesh.computeEdgeLength(iVertex, neighbors[iB]) ;
        }
        edgeDistance = edgeDistance*1.0/neighbors.size();
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], edgeDistance);
    }
}


/* ******** TileSpanCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void TileSpanCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
        T edgeDistance = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            edgeDistance +=  triangleMesh.computeEdgeTileSpan(iVertex, neighbors[iB]) ;
        }
        edgeDistance = edgeDistance*1.0/neighbors.size();
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], edgeDistance);
    }
}


/* ******** MaxEdgeDistanceCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void MaxEdgeDistanceCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
        T edgeDistance = 0.0;
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            edgeDistance =  max(edgeDistance, triangleMesh.computeEdgeLength(iVertex, neighbors[iB])) ;
        }
        this->getStatistics().gatherMax(quantityIds_[particle->get_cellId()], edgeDistance);
    }
}


/* ******** AnglesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void AngleCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        std::vector<plint> neighborTriangles = triangleMesh.getNeighborTriangleIds(iVertex);
        T edgeAngle = 0.0;
        pluint iEdgeAngle=0;
        std::vector<plint> neighbors = triangleMesh.getNeighborVertexIds(iVertex);
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            std::vector<plint>  adjacentTriangles = triangleMesh.getAdjacentTriangleIds(iVertex, neighbors[iB]) ;
            Array<T,3> V1 = triangleMesh.computeTriangleNormal(adjacentTriangles[0]);
            Array<T,3> V2 = triangleMesh.computeTriangleNormal(adjacentTriangles[1]);
            edgeAngle += angleBetweenVectors(V1, V2) + pi ; // - pi + 2*pi
            iEdgeAngle++;
        }
        edgeAngle /= 1.0*iEdgeAngle;
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], edgeAngle); // every edgeAngle is evaluated 2 times
    }
}


/* ******** NumTrianglesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void NumTrianglesCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        T NumTriangles = triangleMesh.getNeighborTriangleIds(iVertex).size()*1.0/3.0; // .computeVertexArea returns the amount per Vertex
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], NumTriangles);
    }
}


/* ******** NumEdgesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void NumEdgesCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        T NumEdges = triangleMesh.getNeighborVertexIds(iVertex).size();
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], NumEdges*0.5);
    }
}


/* ******** NumVerticesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
void NumVerticesCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh, std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_)
{
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        this->getStatistics().gatherSum(quantityIds_[particle->get_cellId()], 1.0);
    }
}


#endif  // IMMERSED_CELLS_REDUCTIONS_HH

