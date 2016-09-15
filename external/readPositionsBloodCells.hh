#ifndef READ_POSISIONS_OF_MULTIPLE_CELLS_HH
#define READ_POSISIONS_OF_MULTIPLE_CELLS_HH

#include "readPositionsBloodCells.h"


template<typename T>
void getReadPositionsBloodCellsVector(Box3D realDomain,
                                            std::vector<TriangularSurfaceMesh<T>* > & meshes,
                                            std::vector<plint> & Np,
                                            std::vector<std::vector<Array<T,3> > > & positions,
                                            std::vector<std::vector<plint> > & cellIds,
                                            std::vector<std::vector<Array<T,3> > > & randomAngles,
                                            const char* positionsFileName)
{
    PLB_PRECONDITION( meshes.size() == 2 );


    if(global::mpi().getSize()>1)
    {
        pcout << "WARNING! (readPositionsBloodCels) You should run this type of initialisation in a single process!" << endl
             << "WARNING! If you need multiprocessor initialisation use (randomPositionsMultipleCells)!"
             << endl << "WARNING! You should cancel this run (unless you are really sure you want this), as only the master will be initialised this way!" << endl;
    }

    pcout << "(readPositionsBloodCels) Reading particle positions..." << std::endl;


    vector<vector3> packPositions[2];
    vector<vector3> packAngles[2];


    // Reading data from file

    plb_ifstream fIn(positionsFileName);

    //uint nCell[2];
    Np.resize(2);
    fIn >> Np[0] >> Np[1];

    pcout << "(readPositionsBloodCels) Particle count (RBCs, PLTs): " << Np[0] << ", " << Np[1] << endl;

    // TODO: proper try-catch
    for(pluint j = 0; j < 2; j++) {

        packPositions[j].resize(Np[j]); packAngles[j].resize(Np[j]);

        for (pluint i = 0; i < Np[j]; i++) {
            fIn >> packPositions[j][i][0] >> packPositions[j][i][1] >> packPositions[j][i][2] >> packAngles[j][i][0]
                >> packAngles[j][i][1] >> packAngles[j][i][2];
            packAngles[j][i] *= pi/180.0; // Deg to Rad
            packAngles[j][i] *= -1.0;  // Right- to left-handed coordinate system

        }
    }
    //

    pcout << "(readPositionsBloodCels) Reading done." << std::endl;

    pcout << "(readPositionsBloodCels) Domain: " << (int)realDomain.getNx() << " " << (int)realDomain.getNy() << " " << (int)realDomain.getNz() << endl;

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



/* ******** ReadPositionMultipleCellField3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void ReadPositionsBloodCellField3D<T,Descriptor>::processGenericBlocks (
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
        mesh->translate(Array<T,3>(-(xRange[0]+xRange[1])/2.0, -(yRange[0]+yRange[1])/2.0, -(zRange[0]+zRange[1])/2.0));
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

    getReadPositionsBloodCellsVector(realDomain, meshes, Np, positions, cellIds, randomAngles, positionsFileName);

    // Change positions to match dx (it is in um originally)
    T posRatio = 1e-6/dx;

    for (pluint iCF = 0; iCF < positions.size(); ++iCF)
    {
        cout << "(ReadPositionsBloodCellField3D) ";
        for (pluint c = 0; c < positions[iCF].size(); ++c)
        {
            ElementsOfTriangularSurfaceMesh<T> emptyEoTSMCopy;
    	    TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*meshes[iCF], emptyEoTSMCopy);
//    	    meshRandomRotation(meshCopy, randomAngles[iCF][c]);
            meshRotation (meshCopy, randomAngles[iCF][c]);
            //positionCellInParticleField(*(particleFields[iCF]), fluid,
            //                            meshes[iCF], positions[iCF][c]-0.5, cellIds[iCF][c]);
            positionCellInParticleField(*(particleFields[iCF]), fluid,
                                         meshCopy, positions[iCF][c]*posRatio, cellIds[iCF][c]);
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
ReadPositionsBloodCellField3D<T,Descriptor>* ReadPositionsBloodCellField3D<T,Descriptor>::clone() const {
    return new ReadPositionsBloodCellField3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ReadPositionsBloodCellField3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        modified[1+iField] = modif::dynamicVariables; // Particle fields.
    }
}

template<typename T, template<typename U> class Descriptor>
void ReadPositionsBloodCellField3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        isWritten[1+iField] = true; // Particle fields.
    }

}


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ReadPositionsBloodCellField3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor>
void readPositionsBloodCellField3D(std::vector<CellField3D<T, Descriptor> *> &cellFields, T dx, const char* positionsFileName) {
    global::timer("CellInit").start();
    std::vector<MultiBlock3D *> fluidAndParticleFieldsArg;

    fluidAndParticleFieldsArg.push_back(&(cellFields[0]->getFluidField3D()));

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back(&(cellFields[icf]->getParticleField3D()));
    }

    //std::vector<MultiBlock3D *> particleFieldsArg;
    //particleFieldsArg.push_back(&(cellFields[0]->getParticleField3D()));

    applyProcessingFunctional(
            new ReadPositionsBloodCellField3D<T, Descriptor>(cellFields, dx, positionsFileName),
            cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg);

    pcout << "Ready to start.." << std::endl;

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        cellFields[icf]->advanceParticles();
        cellFields[icf]->synchronizeCellQuantities();
    }

    global::timer("CellInit").stop();
}

#endif