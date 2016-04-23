#ifndef RANDOM_POSISIONS_OF_MULTIPLE_CELLS_HH
#define RANDOM_POSISIONS_OF_MULTIPLE_CELLS_HH

#include "randomPositionsMultipleCells.h"


template<typename T>
void getRandomPositionsMultipleCellsVector(Box3D realDomain,
                                            std::vector<TriangularSurfaceMesh<T>* > & meshes,
                                            std::vector<plint> & Np,
                                            std::vector<std::vector<Array<T,3> > > & positions,
                                            std::vector<std::vector<plint> > & cellIds,
                                            std::vector<std::vector<Array<T,3> > > & randomAngles,
                                            T packingDensity, plint maxPackIter)
{
    PLB_PRECONDITION( meshes.size()==Np.size() );


    vector<int> domainSize;
    domainSize.push_back((int)realDomain.getNx());
    domainSize.push_back((int)realDomain.getNy());
    domainSize.push_back((int)realDomain.getNz());

    vector<vector3> diameters;
    vector<int> nPartsPerComponent;

    nPartsPerComponent.clear(); nPartsPerComponent.resize(Np.size());
    diameters.clear(); diameters.resize(Np.size());

    for(int i = 0; i < Np.size(); i++)
    {
        T dx, dy, dz;
        Array<T,2> xRange, yRange, zRange;
        meshes[i]->computeBoundingBox (xRange, yRange, zRange);
        dx = (xRange[1] - xRange[0]);
        dy = (yRange[1] - yRange[0]);
        dz = (zRange[1] - zRange[0]);

        diameters[i] = vector3(dx,dy,dz);
        nPartsPerComponent[i] = Np[i];
    }

    pcout << "Executing packing dynamics..." << std::endl;

    Packing pack;

    pack.initSuspension(nPartsPerComponent, diameters, domainSize, packingDensity, maxPackIter);
    pack.execute();

    vector<vector<vector3> > packPositions;
    vector<vector<vector3> > packAngles;

    pack.getOutput(packPositions, packAngles);
    
    //Debug
    //pack.testOutput(); pack.savePov("ellipsoids.pov");

    pcout << "Packing Done." << std::endl;


    // Copy results of packing to appropriate arrays for ficsion
    positions.clear();	positions.resize(Np.size());
    randomAngles.clear(); randomAngles.resize(Np.size());
    cellIds.clear();	cellIds.resize(Np.size());

    plint ni=0;
    int mpiRank = global::mpi().getRank();
    Array<T,3> rdomain(realDomain.x0, realDomain.y0, realDomain.z0);

    for (pluint i = 0; i < Np.size(); ++i)
    {
        positions[i].resize(Np[i]);
        cellIds[i].resize(Np[i]);
        randomAngles[i].resize(Np[i]);

        for (pluint j = 0; j < Np[i]; ++j)
        {
            // Store mesh positions and rotations
            //randomAngles[i][j] = Array<T, 3>(1.0,0.0,0.0);
            randomAngles[i][j] = Array<T, 3>(packAngles[i][j][0], packAngles[i][j][1], packAngles[i][j][2]);
            positions[i][j] = rdomain + Array<T,3>(packPositions[i][j][0], packPositions[i][j][1], packPositions[i][j][2]);
            cellIds[i][j] = ni + 1000*mpiRank;

            // Rotate mesh

            ni++;
        }
    }

}



/* ******** OrderedPositionMultipleCellField3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void RandomPositionMultipleCellField3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
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

    getRandomPositionsMultipleCellsVector(realDomain, meshes, Np, positions, cellIds, randomAngles, packingDensity, maxPackIter);

    for (pluint iCF = 0; iCF < positions.size(); ++iCF)
    {
        cout << "(RandomPositionMultipleCellField3D) ";
        for (pluint c = 0; c < positions[iCF].size(); ++c)
        {
            ElementsOfTriangularSurfaceMesh<T> emptyEoTSMCopy;
    	    TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*meshes[iCF], emptyEoTSMCopy);
//    	    meshRandomRotation(meshCopy, randomAngles[iCF][c]);
            meshRotation (meshCopy, randomAngles[iCF][c]);
            //positionCellInParticleField(*(particleFields[iCF]), fluid,
            //                            meshes[iCF], positions[iCF][c]-0.5, cellIds[iCF][c]);
            positionCellInParticleField(*(particleFields[iCF]), fluid,
                                         meshCopy, positions[iCF][c]-0.5, cellIds[iCF][c]);
			delete meshCopy;
        }

        // DELETE CELLS THAT ARE NOT WHOLE
        plint nVertices=meshes[iCF]->getNumVertices();
        cout << mpiRank;
        plint cellsDeleted = deleteIncompleteCells(*(particleFields[iCF]), fluid, relativeDomains[iCF], nVertices);
        std::vector<Particle3D<T,Descriptor>*> particles;
        particleFields[iCF]->findParticles(particleFields[iCF]->getBoundingBox(),   particles);
        cout    << ", number of particles/nVertices " << particles.size()*1.0/nVertices
        << " (deleted:" << cellsDeleted << ") for particleId:" << iCF << std::endl;
        delete meshes[iCF];
    }
}


template<typename T, template<typename U> class Descriptor>
RandomPositionMultipleCellField3D<T,Descriptor>* RandomPositionMultipleCellField3D<T,Descriptor>::clone() const {
    return new RandomPositionMultipleCellField3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionMultipleCellField3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        modified[1+iField] = modif::dynamicVariables; // Particle fields.
    }
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionMultipleCellField3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        isWritten[1+iField] = true; // Particle fields.
    }

}


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT RandomPositionMultipleCellField3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor>
void randomPositionMultipleCellField3D(std::vector<CellField3D<T, Descriptor> *> &cellFields, T packingDensity, plint maxPackIter = 25000) {
    global::timer("CellInit").start();
    std::vector<MultiBlock3D *> fluidAndParticleFieldsArg;

    fluidAndParticleFieldsArg.push_back(&(cellFields[0]->getFluidField3D()));

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back(&(cellFields[icf]->getParticleField3D()));
    }

    //std::vector<MultiBlock3D *> particleFieldsArg;
    //particleFieldsArg.push_back(&(cellFields[0]->getParticleField3D()));

    applyProcessingFunctional(
            new RandomPositionMultipleCellField3D<T, Descriptor>(cellFields, packingDensity, maxPackIter),
            cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg);

    pcout << "Ready to start.." << std::endl;

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        cellFields[icf]->advanceParticles();
        cellFields[icf]->synchronizeCellQuantities();
    }

    global::timer("CellInit").stop();
}

#endif