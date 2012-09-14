#ifndef IMMERSED_CELLS_FUNCTIONAL_3D_HH
#define IMMERSED_CELLS_FUNCTIONAL_3D_HH

#include "immersedCellsFunctional3D.h"

#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"


/* ******** countCellVolume *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVolume (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVolumes) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    ComputeCellVolumeParticlesFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellVolumeArray(cellVolumes, cellIds);
}

/* ******** ComputeCellVolumeParticlesFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::ComputeCellVolumeParticlesFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
    : triangleBoundary(triangleBoundary_)
{
    numberOfCells = 0;
    for (pluint var = 0; var < cellIds_.size(); ++var)
        if (cellIds_[var] > numberOfCells)
            numberOfCells = cellIds_[var];

    for (pluint i=0; i< (pluint) numberOfCells+1; ++i)
        this->getStatistics().subscribeSum() ;
//    pcout << "Done subscribing" << std::endl;
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    std::map<plint, T> tagsVolume;
    typename std::map<plint, T>::iterator tagsVolumeIterator;

//lmount: Find particles within this domain
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::map<plint, plint> tagsNr;
    typename std::map<plint, plint>::iterator tagsNrIterator;
//lmount
    for (pluint iA = 0; iA < particles.size(); ++iA) {
            pcout << "iA: " << iA << "id " << (dynamic_cast<ImmersedCellParticle3D<T,DESCRIPTOR>*> (particles[iA]))->get_cellId() << std::endl;
	}


    plint meshID = 1, cellId;
    triangleBoundary.pushSelect(0,meshID);
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    T triangleVolumeT6;
    Array<T,3> areaTimesNormal;
    Array<T,3> v0, v1, v2, tmp;
    for (plint iTriangle=0; iTriangle<triangleMesh.getNumTriangles(); ++iTriangle) {
        v0 = triangleMesh.getVertex(iTriangle, 0);
        v1 = triangleMesh.getVertex(iTriangle, 1);
        v2 = triangleMesh.getVertex(iTriangle, 2);
        /* Calculating the volume contibution of a face based on the formula:
         * V[j] = 1.0/6.0 * (X3[j] cross X2[j])*X1[j]  */
        crossProduct(v1, v2, tmp);
        triangleVolumeT6 =  VectorTemplate<T,Descriptor>::scalarProduct(v0,tmp); // * (1.0/6.0)
        /* ************* Other Calculation ********* */
//        areaTimesNormal = triangleMesh.computeTriangleNormal(iTriangle, true);
//        triangleVolumeT6 = 2.0 * VectorTemplate<T,Descriptor>::scalarProduct(areaTimesNormal, ((v0+v1+v2)/3.0)) ;
//        triangleVolume = 1.0/3.0 * VectorTemplate<T,Descriptor>::scalarProduct(areaTimesNormal, ((v0+v1+v2)/3.0)) ;
        /* ********************************************* */

        // Update
        plint ctrlId = triangleMesh.getVertexId(iTriangle, 0);
        plint ctrlCellId = (dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[ctrlId]))->get_cellId();
//        pcout << iTriangle << ", id = " << ctrlId << ", cell iD = " << ctrlCellId << " ";
//        for (plint iA = 1; iA < 3; ++iA) {
//        	plint id = triangleMesh.getVertexId(iTriangle, iA);
//        	plint tmpCellId = (dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[id]))->get_cellId();
//        	pcout << ", id = " << id << ", cell iD = " << tmpCellId << " ";
//        	PLB_ASSERT(tmpCellId==ctrlCellId);
//        }
//        pcout << std::endl;
        cellId = triangleMesh.getVertexId(iTriangle, 0);
        cellId = (dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[cellId]))->get_cellId();
        tagsVolume[cellId] += triangleVolumeT6;
        tagsNr[cellId] += 1;
    }
    triangleBoundary.popSelect();

    pcout << "triangleMesh.getNumTriangles " << triangleMesh.getNumTriangles() << std::endl;
    pcout << "triangleMesh.getNumVertices " << triangleMesh.getNumVertices() << std::endl;
    pcout << "NumberOfCells " << triangleMesh.getNumTriangles() << std::endl;
    for ( tagsNrIterator=tagsNr.begin() ; tagsNrIterator != tagsNr.end(); tagsNrIterator++ )
        pcout << "cellId " << (*tagsNrIterator).first << " num " << (*tagsNrIterator).second <<std::endl;
    for ( tagsVolumeIterator=tagsVolume.begin() ; tagsVolumeIterator != tagsVolume.end(); tagsVolumeIterator++ ) {
        this->getStatistics().gatherSum((*tagsVolumeIterator).first, (*tagsVolumeIterator).second/6.0);
    }

}


template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::getCellVolumeArray(std::vector<T>& cellVolumes, std::vector<plint> cellIds) const {
    for (pluint i = 0; i < (pluint) cellIds.size(); ++i)
        cellVolumes.push_back(this->getStatistics().getSum(cellIds[i]));
}


template<typename T, template<typename U> class Descriptor>
ComputeCellVolumeParticlesFunctional3D<T,Descriptor>* ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::clone() const {
    return new ComputeCellVolumeParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}



//template<typename T, template<typename U> class Descriptor>
//CountCellVolumeFunctional3D<T,Descriptor>::CountCellVolumeFunctional3D(plint tag_)
//    : cellVolumeId(this->getStatistics().subscribeIntSum()), tag(tag_)
//{ }
//
//template<typename T, template<typename U> class Descriptor>
//void CountCellVolumeFunctional3D<T,Descriptor>::processGenericBlocks (
//        Box3D domain, std::vector<AtomicBlock3D*> blocks )
//{
//    PLB_PRECONDITION( blocks.size()==1 );
//    ParticleField3D<T,Descriptor>& particleField
//        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
//    std::vector<Particle3D<T,Descriptor>*> particles;
//    particleField.findParticles(domain, particles);
//    T cellVolume = 0;
//    for (pluint iA = 0; iA < particles.size(); ++iA) {
//        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
//        ImmersedCellParticle3D<T,Descriptor>* particle =
//            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
//
//        if (particle->get_cellId() == tag) {
////         if (particle->getTag() == tag) {
//        	cellVolume += 0.;
//        }
//    }
//    this->getStatistics().gatherIntSum(cellVolumeId, cellVolume);
//}
//
//template<typename T, template<typename U> class Descriptor>
//CountCellVolumeFunctional3D<T,Descriptor>* CountCellVolumeFunctional3D<T,Descriptor>::clone() const {
//    return new CountCellVolumeFunctional3D<T,Descriptor>(*this);
//}
//
//template<typename T, template<typename U> class Descriptor>
//void CountCellVolumeFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
//    modified[0] = modif::nothing;
//}
//
//template<typename T, template<typename U> class Descriptor>
//T CountCellVolumeFunctional3D<T,Descriptor>::getVollumeCell() const {
//    return this->getStatistics().getSum(cellVolumeId);
//}
//
//template< typename T, template<typename U> class Descriptor,
//          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
//T countCellVolume (
//                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, plint tag)
//{
//    std::vector<MultiBlock3D*> particleArg;
//    particleArg.push_back(&particles);
//
//    CountCellVolumeFunctional3D<T,Descriptor> functional(tag);
//    applyProcessingFunctional(functional, domain, particleArg);
//    return functional.getVollumeCell();
//}
//
//
//





#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
