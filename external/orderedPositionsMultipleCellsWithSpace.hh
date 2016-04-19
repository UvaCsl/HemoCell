#ifndef ORDERED_POSISIONS_OF_MULTIPLE_CELLS_WITH_SPACE_HH
#define ORDERED_POSISIONS_OF_MULTIPLE_CELLS_WITH_SPACE_HH
#include "orderedPositionsMultipleCellsWithSpace.h"
#include <algorithm>    // std::random_shuffle


template<typename T, template<typename U> class Descriptor>
void OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    pluint numberOfCellFields=cellFields.size();
    PLB_PRECONDITION( blocks.size()==(numberOfCellFields+1) );
    int mpiSize = global::mpi().getSize();;
    int mpiRank = global::mpi().getRank();
    srand (mpiSize * mpiRank);
    T ratio;
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);

    std::vector<ParticleField3D<T,Descriptor>* > particleFields(numberOfCellFields);
    std::vector<T> volumes(numberOfCellFields);
    std::vector<T> volumeFractions(numberOfCellFields);
    std::vector<TriangularSurfaceMesh<T>* > meshes(numberOfCellFields);
    std::vector<ElementsOfTriangularSurfaceMesh<T> > emptyEoTSM(numberOfCellFields);
    std::vector<Box3D> relativeDomains(numberOfCellFields);
	std::vector<plint> Np(numberOfCellFields);

    T totalVolumeFraction=0;

    T Vdomain = domain.getNx() * domain.getNy() * domain.getNz();

    Dot3D fLocation(fluid.getLocation());
    Box3D realDomain(
                        domain.x0 + fLocation.x, domain.x1 + fLocation.x,
                        domain.y0 + fLocation.y, domain.y1 + fLocation.y,
                        domain.z0 + fLocation.z, domain.z1 + fLocation.z );
    Array<T,2> xRange, yRange, zRange;

    for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {
        T vf = cellFields[iCF]->getVolumeFraction();
        if (vf > 1.0) { vf /= 100; }
        volumeFractions[iCF] = vf;

        TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(cellFields[iCF]->getMesh(), emptyEoTSM[iCF]);
        mesh->computeBoundingBox (xRange, yRange, zRange);
        mesh->translate(Array<T,3>(-xRange[0], -yRange[0], -zRange[0]));
        meshes[iCF] = mesh;
        volumes[iCF] = MeshMetrics<T>(*mesh).getVolume();

        Np[iCF] = (1 + 0.05*(guessRandomNumber() -0.5)) * (Vdomain/volumes[iCF] * volumeFractions[iCF]) ;
        totalVolumeFraction += volumes[iCF];
        particleFields[iCF] = ( dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[iCF+1]) );
        particleFields[iCF]->removeParticles(particleFields[iCF]->getBoundingBox());
        Dot3D relativeDisplacement = computeRelativeDisplacement(fluid, *(particleFields[iCF]));
        relativeDomains[iCF] = ( Box3D(
                            domain.x0 + relativeDisplacement.x, domain.x1 + relativeDisplacement.x,
                            domain.y0 + relativeDisplacement.y, domain.y1 + relativeDisplacement.y,
                            domain.z0 + relativeDisplacement.z, domain.z1 + relativeDisplacement.z ) );
    }

    std::vector<std::vector<Array<T,3> > > positions;
    std::vector<std::vector<plint> > cellIds;
    std::vector<std::vector<Array<T,3> > > randomAngles;
    getOrderedPositionsMultipleCellsVector(realDomain, meshes, Np, positions, cellIds, randomAngles);

    for (pluint iCF = 0; iCF < positions.size(); ++iCF)
    {
        cout << "(OrderedPositionMultipleCellField3DWithSpace) ";
    	for (pluint c = 0; c < positions[iCF].size(); ++c)
    	{
//    	    ElementsOfTriangularSurfaceMesh<T> emptyEoTSMCopy;
//    	    TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*meshes[iCF], emptyEoTSMCopy);
//    	    meshRandomRotation(meshCopy, randomAngles[iCF][c]);
			positionCellInParticleField(*(particleFields[iCF]), fluid, 
			        meshes[iCF], positions[iCF][c]-0.5, cellIds[iCF][c]);
//			delete meshCopy;
    	}

        // Deleting all LSPs from Box3D ToDelete
        for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {
           particleFields[iCF]->removeParticles(ToDelete);
        }

	    // DELETE CELLS THAT ARE NOT WHOLE
        plint nVertices=meshes[iCF]->getNumVertices();
        cout << mpiRank;
        plint cellsDeleted = deleteIncompleteCells(*(particleFields[iCF]), fluid, relativeDomains[iCF], nVertices); 
        std::vector<Particle3D<T,Descriptor>*> particles;
	particleFields[iCF]->findParticles(particleFields[iCF]->getBoundingBox(),   particles);
	//cout    << ", number of particles/nVertices " << particles.size()*1.0/nVertices
	//        << " (deleted:" << cellsDeleted << ") for particleId:" << iCF << std::endl;
	delete meshes[iCF];
   }
}


template<typename T, template<typename U> class Descriptor>
OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>* OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>::clone() const {
    return new OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        modified[1+iField] = modif::dynamicVariables; // Particle fields.
    }
}

template<typename T, template<typename U> class Descriptor>
void OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        isWritten[1+iField] = true; // Particle fields.
    }

}


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}




// Function
template<typename T, template<typename U> class Descriptor>
void orderedPositionMultipleCellField3DWithSpace(std::vector<CellField3D<T, Descriptor>* > & cellFields, Box3D ToDelete) {
//void orderedPositionMultipleCellField3DWithSpace(CellField3D<T, Descriptor>  & cellField, Box3D ToDelete) {
    int mpiRank = global::mpi().getRank();

    global::timer("CellInit").start();
    std::vector<MultiBlock3D*> fluidAndParticleFieldsArg;
    fluidAndParticleFieldsArg.push_back( &(cellFields[0]->getFluidField3D()) );
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back( &(cellFields[icf]->getParticleField3D()) );
    }

    std::vector<MultiBlock3D*> particleFieldsArg;
    particleFieldsArg.push_back( &(cellFields[0]->getParticleField3D()) );

    applyProcessingFunctional (
        new OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>(cellFields, ToDelete),
        cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg );
    pcout << "Ready to start.." << std::endl;

        
	    // DELETE CELLS THAT ARE NOT WHOLE
        plint nVertices=cellFields[0]->getMesh().getNumVertices();
        cout << mpiRank;
        // FOLLOWING SHOULD BE INSIDE THE OrderedPositionMultipleCellField3DWithSpace FUNCTIONAL (AND THEY ARE)
        // plint cellsDeleted = deleteIncompleteCells((cellFields[icf]->getParticleField3D()), cellFields[icf]->getFluidField3D(), relativeDomains[icf], nVertices); 
        // std::vector<Particle3D<T,Descriptor>*> particles;
	// cellFields[icf]->getParticleField3D().findParticles(cellFields[icf]->getParticleField3D().getBoundingBox(),   particles);
	// cout    << ", number of particles/nVertices " << particles.size()*1.0/nVertices
	        // << " (deleted:" << cellsDeleted << ") for particleId:" << icf << std::endl;
	// delete meshes[icf];

        cellFields[0]->advanceParticles();
        cellFields[0]->synchronizeCellQuantities();



    global::timer("CellInit").stop();
}

#endif // ORDERED_POSISIONS_OF_MULTIPLE_CELLS_HH

