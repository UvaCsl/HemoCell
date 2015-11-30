#ifndef INITIALIZATION_CELL_FIELD_3D_HH
#define INITIALIZATION_CELL_FIELD_3D_HH
#include "initializationCellField3D.h"
#include "ParticleHdf5IO.h"



template<typename T, template<typename U> class Descriptor>
void randomPositionCellFieldsForGrowth3D(std::vector<CellField3D<T, Descriptor>* > & cellFields) {
    global::timer("CellInit").start();
    std::vector<MultiBlock3D*> fluidAndParticleFieldsArg;
    fluidAndParticleFieldsArg.push_back( &(cellFields[0]->getFluidField3D()) );
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back( &(cellFields[icf]->getParticleField3D()) );
    }

    std::vector<MultiBlock3D*> particleFieldsArg;
    particleFieldsArg.push_back( &(cellFields[0]->getParticleField3D()) );

    applyProcessingFunctional (
        new RandomPositionCellFieldsForGrowth3D<T,Descriptor>(cellFields),
        cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg );
    pcout << "Ready to start.." << std::endl;
//    applyProcessingFunctional (
//        new PrintParticles<T,Descriptor>(),
//        cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg );

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        cellFields[icf]->advanceParticles();
        cellFields[icf]->synchronizeCellQuantities();
    }

    global::timer("CellInit").stop();
}

/* ******** PositionCellParticles3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void PositionCellParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    /* ######## Center mesh ######## */
    ElementsOfTriangularSurfaceMesh<T> emptyEoTSM;
    TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(elementaryMesh, emptyEoTSM);
    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    T dx = xRange[1] - xRange[0], dy = yRange[1] - yRange[0], dz = zRange[1] - zRange[0];
    T maxSide = max( max(dx, dy), dz);
    // Bring mesh to the center
    Array<T,3> mvp(dx, dy, dz);
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    mesh->translate(-1.0 * meshCenter);
    /* ############################# */

    plint nVertices = mesh->getNumVertices();
    for (pluint iCO=0; iCO < cellOrigins.size(); ++iCO) {
        Array<T,3> & cellOrigin = cellOrigins[iCO];
        plint cellId = iCO;
        for (plint iVertex=0; iVertex < nVertices; ++iVertex) {
            Array<T,3> vertex = cellOrigin + mesh->getVertex(iVertex);
            particleField.addParticle(domain, new SurfaceParticle3D<T,Descriptor>(vertex, cellId, iVertex) );
        }
    }
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);

   pcout << "(PositionCellParticles3D) pid " << global::mpi().getRank()
        << " particles.size " << particles.size()
        << " nVertices " << nVertices
        << " cellOrigin.size " << cellOrigins.size() << std::endl;
   delete mesh;
}

template<typename T, template<typename U> class Descriptor>
PositionCellParticles3D<T,Descriptor>* PositionCellParticles3D<T,Descriptor>::clone() const {
    return new PositionCellParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void PositionCellParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::allVariables; // Particle field.
}

template<typename T, template<typename U> class Descriptor>
void PositionCellParticles3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT PositionCellParticles3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** RandomPositionCellParticlesForGrowth3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellParticlesForGrowth3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    int mpiSize = global::mpi().getSize();;
    int mpiRank = global::mpi().getRank();

    if (hematocrit > 1.0) { hematocrit /= 100; }
    PLB_PRECONDITION( hematocrit < 1.0 );

    srand (mpiSize * mpiRank);
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);
//    fluid.get(1, 1, 1).getDynamics().isBoundary();
    ElementsOfTriangularSurfaceMesh<T> emptyEoTSM;
    TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(elementaryMesh, emptyEoTSM);

    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    T dx = xRange[1] - xRange[0], dy = yRange[1] - yRange[0], dz = zRange[1] - zRange[0];
    T maxSide = max( max(dx, dy), dz);
    // Bring mesh to the center
    Array<T,3> mvp(dx, dy, dz);
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    mesh->translate(-1.0 * meshCenter);
    MeshMetrics<T> mm(*mesh);
    T volume = mm.getVolume();

    plint DeltaX = domain.getNx();
    plint DeltaY = domain.getNy();
    plint DeltaZ = domain.getNz();

    plint step = plint( pow(volume/hematocrit, 1.0/3.0) );
    T Vdomain = DeltaX * DeltaY * DeltaZ;
    T Vsubdomain = plint((DeltaX-0.5)/step) * plint((DeltaY-0.5)/step) * plint((DeltaZ-0.5)/step) * step*step*step ;
    T rescaledHematocrit = hematocrit * Vdomain / Vsubdomain;
    T prob = step * step * step * 1.0 / volume * rescaledHematocrit;
    while (prob > 1.0) {
        step -= 1.0;
        Vsubdomain = plint(DeltaX/step) * plint(DeltaY/step) * plint(DeltaZ/step) * step*step*step ;
        rescaledHematocrit = hematocrit * Vdomain / Vsubdomain;
        prob = step * step * step * 1.0 / volume * rescaledHematocrit;
    }

    T scale = step*1.0/maxSide;
    if (scale > 1.0) { scale = 1.0; };
    mesh->scale(scale);
    ratio = scale;

    // Access the position of the atomic-block inside the multi-block.
    Dot3D relativePosition = particleField.getLocation();
    Array<T,3> relativeCoordinate(relativePosition.x, relativePosition.y, relativePosition.z);
    Dot3D fluidLocationDot3D = fluid.getLocation();
    Array<T,3> fluidLocation(fluidLocationDot3D.x, fluidLocationDot3D.y, fluidLocationDot3D.z);

//    cout << global::mpi().getRank()
//             << "domain "
//                 << domain.x0 << " " << domain.x1 << ", "
//                 << domain.y0 << " " << domain.y1 << ", "
//                 << domain.z0 << " " << domain.z1 << ", "
//             << "fluidLocation "
//                 << fluidLocationDot3D.x << ", "
//                 << fluidLocationDot3D.y << ", "
//                 << fluidLocationDot3D.z << ", "
//             << "particlePosition "
//                 << relativePosition.x << ", "
//                 << relativePosition.y << ", "
//                 << relativePosition.z << ", "
//             << std::endl;


    plint nVertices = mesh->getNumVertices();
    plint cellId = 0;
    // Loop through the domain and place cell depending on the
    for (plint iX=domain.x0; iX<=domain.x1-step; iX+=step) {
        for (plint iY=domain.y0; iY<=domain.y1-step; iY+=step) {
            for (plint iZ=domain.z0; iZ<=domain.z1-step; iZ+=step) {
                T rn = guessRandomNumber();
                if (rn <= prob) {
                    meshRandomRotation(mesh);
                    for (plint iVertex=0; iVertex < nVertices; ++iVertex) {
                        Array<T,3> vertex = Array<T,3>(iX*1.0, iY*1.0, iZ*1.0) + relativeCoordinate + mesh->getVertex(iVertex);
                        Dot3D cellInFluidDomain = Dot3D(vertex[0], vertex[1], vertex[2]) - fluidLocationDot3D;
                        if (fluid.get(cellInFluidDomain.x, cellInFluidDomain.y, cellInFluidDomain.z).getDynamics().isBoundary()) {
//                            cout << global::mpi().getRank() << "vertex (" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")" <<std::endl;
                            break;
                        } else {
                            particleField.addParticle(domain, new SurfaceParticle3D<T,Descriptor>(vertex, cellId + 1000*mpiRank, iVertex) );
                        }

                    }
                    cellId+=1;
                }
            }
        }
    }

    // DELETE CELLS THAT ARE NOT WHOLE
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::map<plint, plint> cellIdToNumberOfVertices;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        SurfaceParticle3D<T,Descriptor>* particle =
            dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (particles[iP]);
        cellIdToNumberOfVertices[particle->get_cellId()]++;
    }
    plint cellsDeleted=0;
    typename std::map<plint, plint >::iterator itrt;
    for (itrt  = cellIdToNumberOfVertices.begin(); itrt != cellIdToNumberOfVertices.end(); ++itrt) {
        if (itrt->second < nVertices) {
            cellsDeleted++;
            particleField.removeParticles(domain, itrt->first);
        }
    }
    particleField.findParticles(domain, particles);
    pcout << "(RandomPositionCellParticlesForGrowth3D) "
            << mpiRank << "Number of particles/nVertices " << particles.size()*1.0/nVertices
            << " scale " << scale << " step " << step
            << " prob " << prob
            << " h " << particles.size()*1.0/nVertices * volume * 1.0/ (DeltaX*DeltaY*DeltaZ) << " (deleted:" << cellsDeleted << ")" << std::endl;
    delete mesh;
}

template<typename T, template<typename U> class Descriptor>
RandomPositionCellParticlesForGrowth3D<T,Descriptor>* RandomPositionCellParticlesForGrowth3D<T,Descriptor>::clone() const {
    return new RandomPositionCellParticlesForGrowth3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellParticlesForGrowth3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing; // Fluid field.
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellParticlesForGrowth3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
    isWritten[1] = false;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT RandomPositionCellParticlesForGrowth3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}







/* ******** RandomPositionCellFieldsForGrowth3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellFieldsForGrowth3D<T,Descriptor>::processGenericBlocks (
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
    T totalVolumeFraction=0;


    Dot3D fLocation(fluid.getLocation());
    Box3D realDomain(
                        domain.x0 + fLocation.x, domain.x1 + fLocation.x,
                        domain.y0 + fLocation.y, domain.y1 + fLocation.y,
                        domain.z0 + fLocation.z, domain.z1 + fLocation.z );
//    bool particleIsInBulk = particleField.isContained(particle->getPosition(), realDomain);

    for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {
        T vf = cellFields[iCF]->getVolumeFraction();
        if (vf > 1.0) { vf /= 100; }
        volumeFractions[iCF] = vf;

        TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(cellFields[iCF]->getMesh(), emptyEoTSM[iCF]);
        meshes[iCF] = mesh;
        volumes[iCF] = MeshMetrics<T>(*mesh).getVolume();

        totalVolumeFraction += volumes[iCF];
        particleFields[iCF] = ( dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[iCF+1]) );
        particleFields[iCF]->removeParticles(particleFields[iCF]->getBoundingBox());

        Dot3D relativeDisplacement = computeRelativeDisplacement(fluid, *(particleFields[iCF]));
        relativeDomains[iCF] = ( Box3D(
                            domain.x0 + relativeDisplacement.x, domain.x1 + relativeDisplacement.x,
                            domain.y0 + relativeDisplacement.y, domain.y1 + relativeDisplacement.y,
                            domain.z0 + relativeDisplacement.z, domain.z1 + relativeDisplacement.z ) );
    }

    plint DeltaX = realDomain.getNx();
    plint DeltaY = realDomain.getNy();
    plint DeltaZ = realDomain.getNz();


    Array<T,2> xRange, yRange, zRange;
    meshes[0]->computeBoundingBox (xRange, yRange, zRange);
    T dx = xRange[1] - xRange[0], dy = yRange[1] - yRange[0], dz = zRange[1] - zRange[0];
    T maxSide = max( max(dx, dy), dz);
    // Bring mesh to the center
    Array<T,3> mvp(dx, dy, dz);
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    meshes[0]->translate(-1.0 * meshCenter);
    MeshMetrics<T> mm(*(meshes[0]));
    T volume = mm.getVolume();

    plint step = plint( pow(volume/volumeFractions[0], 1.0/3.0) );
    T Vdomain = DeltaX * DeltaY * DeltaZ;
    T Vsubdomain = plint((DeltaX-0.5)/step) * plint((DeltaY-0.5)/step) * plint((DeltaZ-0.5)/step) * step*step*step ;
    T rescaledHematocrit = volumeFractions[0] * Vdomain / Vsubdomain;
    T prob = step * step * step * 1.0 / volume * rescaledHematocrit;
    while (prob > 1.0) {
        step -= 1.0;
        Vsubdomain = plint(DeltaX/step) * plint(DeltaY/step) * plint(DeltaZ/step) * step*step*step ;
        rescaledHematocrit = volumeFractions[0] * Vdomain / Vsubdomain;
        prob = step * step * step * 1.0 / volume * rescaledHematocrit;
    }

    T scale = step*1.0/maxSide;
    if (scale > 1.0) { scale = 1.0; };
    meshes[0]->scale(scale);
    ratio = scale;

    // Access the position of the atomic-block inside the multi-block.
    plint nVertices = meshes[0]->getNumVertices();
    plint cellId = 0;

    // Loop through the domain and place cell depending on the
    std::vector<Array<T,3> > freeStartingPoints;
    T smalleps = 0.01;
    T stepX = step * (1 + fmod(DeltaX-smalleps, step)*1.0/DeltaX );
    T stepY = step * (1 + fmod(DeltaY-smalleps, step)*1.0/DeltaY );
    T stepZ = step * (1 + fmod(DeltaZ-smalleps, step)*1.0/DeltaZ );
    for (T iX=realDomain.x0; iX<realDomain.x1-smalleps; iX+=stepX) {
        for (T iY=realDomain.y0; iY<realDomain.y1-smalleps; iY+=stepY) {
            for (T iZ=realDomain.z0; iZ<realDomain.z1-smalleps; iZ+=stepZ) {
                T rn = guessRandomNumber();
                Array<T,3> startingPoint(iX*1.0, iY*1.0, iZ*1.0);
                if (rn <= prob) {
                    meshRandomRotation(meshes[0]);
                    positionCellInParticleField(*(particleFields[0]), fluid, meshes[0], startingPoint, cellId + 1000*mpiRank);
                    cellId+=1;
                } else {
                    freeStartingPoints.push_back(startingPoint);
                }
            }
        }
    }

    // DELETE CELLS THAT ARE NOT WHOLE
    plint cellsDeleted = deleteIncompleteCells(*(particleFields[0]), fluid, relativeDomains[0], nVertices);
    // Position the deleted cells to the free positions
    plint nCellsToReposition = cellsDeleted<freeStartingPoints.size()?cellsDeleted:freeStartingPoints.size();
    for (int iC = 0; iC < nCellsToReposition; ++iC) {
        meshRandomRotation(meshes[0]);
        positionCellInParticleField(*(particleFields[0]), fluid, meshes[0], freeStartingPoints[iC], cellId + 1000*mpiRank);
        cellId+=1;
    }

    cellsDeleted = deleteIncompleteCells(*(particleFields[0]), fluid, relativeDomains[0], nVertices);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleFields[0]->findParticles(particleFields[0]->getBoundingBox(),   particles);
    pcout << "(RandomPositionCellFieldsForGrowth3D) "
            << mpiRank << "Number of particles/nVertices " << particles.size()*1.0/nVertices
            << " scale " << scale << " step " << step
            << " prob " << prob
            << " h " << particles.size()*1.0/nVertices * volume * 1.0/ (DeltaX*DeltaY*DeltaZ) << " (deleted:" << cellsDeleted << ")" << std::endl;
    delete meshes[0];
}


template<typename T, template<typename U> class Descriptor>
RandomPositionCellFieldsForGrowth3D<T,Descriptor>* RandomPositionCellFieldsForGrowth3D<T,Descriptor>::clone() const {
    return new RandomPositionCellFieldsForGrowth3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellFieldsForGrowth3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        modified[1+iField] = modif::dynamicVariables; // Particle fields.
    }
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellFieldsForGrowth3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        isWritten[1+iField] = true; // Particle fields.
    }

}


template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT RandomPositionCellFieldsForGrowth3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T>
void getOrderedPositionsVector(Box3D realDomain, Array<T,3> step,
            std::vector<Array<T,3> > & positions, std::vector<plint> & cellIds) {
    plint maxNumberPerSide=1290; // max long int == 2147483647^(1/3). Yields 1cm, (1200 * 8) or 2.5 mm (1200 * 2)
    set<plint> uniqueCellIds;
    for (T iX=realDomain.x0-2.5*step[0]; iX<=realDomain.x1+2*step[0]; iX+=step[0]) {
        for (T iY=realDomain.y0-2.5*step[1]; iY<=realDomain.y1+2*step[1]; iY+=step[1]) {
            for (T iZ=realDomain.z0-2.5*step[2]; iZ<=realDomain.z1+2*step[2]; iZ+=step[2]) {
                int nx = floor(iX/step[0]);
                int ny = floor(iY/step[1]);
                int nz = floor(iZ/step[2]);
                if (nx>=0 and ny>=0 and nz>=0) {
                    T x = nx*step[0];
                    T y = ny*step[1];
                    T z = nz*step[2];
                    plint cellId =  nx + ny*maxNumberPerSide + nz*maxNumberPerSide*maxNumberPerSide;
                    if (uniqueCellIds.count(cellId) == 0) {
                    	positions.push_back( Array<T,3>(x,y,z) );
                    	cellIds.push_back( cellId );
                    	uniqueCellIds.insert(cellId);
                    }
                }
            }
        }
    }

}


template<typename T, template<typename U> class Descriptor>
void orderedPositionCellField3D(std::vector<CellField3D<T, Descriptor>* > & cellFields, Dot3D latticeSize) {
    global::timer("CellInit").start();
    std::vector<MultiBlock3D*> fluidAndParticleFieldsArg;
    fluidAndParticleFieldsArg.push_back( &(cellFields[0]->getFluidField3D()) );
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back( &(cellFields[icf]->getParticleField3D()) );
    }

    std::vector<MultiBlock3D*> particleFieldsArg;
    particleFieldsArg.push_back( &(cellFields[0]->getParticleField3D()) );

    applyProcessingFunctional (
        new OrderedPositionCellField3D<T,Descriptor>(cellFields, latticeSize),
        cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg );
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        cellFields[icf]->createCellMap();
//        writeCellField3D_HDF5(*cellFields[icf], 1., 1., 1);
        cellFields[icf]->deleteIncompleteCells();
//        cellFields[icf]->createCellMap();
//        writeCellField3D_HDF5(*cellFields[icf], 1., 1., 2);
    }
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
//        cout << global::mpi().getRank() << " N_cells = " << cellFields[icf]->getNumberOfCells_Global() << std::endl;
        cellFields[icf]->advanceParticles();
        cellFields[icf]->synchronizeCellQuantities();
    }
    pcout << "Ready to start.." << std::endl;
    global::timer("CellInit").stop();
}



/* ******** OrderedPositionCellField3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void OrderedPositionCellField3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    pluint numberOfCellFields=cellFields.size();
    PLB_PRECONDITION( blocks.size()==(numberOfCellFields+1) );
    int mpiSize = global::mpi().getSize();;
    int mpiRank = global::mpi().getRank();
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>* particleField = ( dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]) );

    T phi = cellFields[0]->getVolumeFraction();
    if (phi > 1.0) { phi /= 100; }

    ElementsOfTriangularSurfaceMesh<T> emptyEoTSM;
    TriangularSurfaceMesh<T>* mesh = copyTriangularSurfaceMesh(cellFields[0]->getMesh(), emptyEoTSM);;
    T volume = MeshMetrics<T>(*mesh).getVolume();


    Dot3D fLocation(fluid.getLocation());
    Box3D realDomain(
                        domain.x0 + fLocation.x, domain.x1 + fLocation.x,
                        domain.y0 + fLocation.y, domain.y1 + fLocation.y,
                        domain.z0 + fLocation.z, domain.z1 + fLocation.z );

    plint DeltaX = realDomain.getNx();
    plint DeltaY = realDomain.getNy();
    plint DeltaZ = realDomain.getNz();


    Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    T dx = xRange[1] - xRange[0], dy = yRange[1] - yRange[0], dz = zRange[1] - zRange[0];
    // Bring mesh to the center
    Array<T,3> meshCenter = Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * 0.5;
    Array<T,3> mvp(dx, dy, dz);
    mesh->translate((0.5*mvp)-meshCenter);

    T phi_max = volume/(1.0*dx*dy*dz);
//    pcout <<"phi " << phi_max << " " << phi << std::endl;
    T stepFactor = pow(phi_max*1.0/phi, 1.0/3.0);
    PLB_ASSERT(stepFactor > 1.0);
    if (stepFactor < 1.) {
    	stepFactor = 1.05;
    }
    Array<T,3> step = stepFactor * mvp;

    T Vdomain = latticeSize.x * latticeSize.y * latticeSize.z;
	T Vsubdomain = plint((latticeSize.x-0.5)/step[0]) * plint((latticeSize.y-0.5)/step[1]) * plint((latticeSize.z-0.5)/step[2]) * step[0]*step[1]*step[2];
	phi = phi * Vdomain / Vsubdomain;
	stepFactor = pow(phi_max*1.0/phi, 1.0/3.0);
	step = stepFactor * mvp;


    Array<T,3> offset(-0.49, -0.49, -0.49);
    T smalleps = 0.01;
    step[0] = step[0] * (1 + fmod(latticeSize.x-smalleps-offset[0], step[0])*1.0/latticeSize.x );
    step[1] = step[1] * (1 + fmod(latticeSize.y-smalleps-offset[1], step[1])*1.0/latticeSize.y );
    step[2] = step[2] * (1 + fmod(latticeSize.z-smalleps-offset[2], step[2])*1.0/latticeSize.z );

    std::vector<Array<T,3> > positions;
    std::vector<plint > cellIds;
    getOrderedPositionsVector(realDomain, step, positions, cellIds);
//    cout << "(OrderedPositionCellField3D) " << positions.size() << " particles" << std::endl;
    // Loop through the domain and place cell depending on the
    for (pluint iC = 0; iC < positions.size(); ++iC) {
//        cout << "(OrderedPositionCellField3D) cellId: " << cellIds[iC] << " r:" << mpiRank
//                << " (" << positions[iC][0] << ", "
//                 << "" << positions[iC][1] << ", "
//                 << "" << positions[iC][2] << "), "
//                 << std::endl;
        positionCellInParticleField(*particleField, fluid, mesh, positions[iC] + offset, cellIds[iC]);
    }
    pcout << "(OrderedPositionCellField3D) Steps: "
            << "(" << step[0] << ", "
             << "" << step[1] << ", "
             << "" << step[2] << "), "
             << ", Cells " << positions.size()
             << ", Volume " << volume
             << ", stepFactor " << stepFactor
             << ", mvp (" << mvp[0] << ", "
              << "" << mvp[1] << ", "
              << "" << mvp[2] << "), "
            << std::endl;

//    plint nVertices = mesh->getNumVertices();
//    std::vector<Particle3D<T,Descriptor>*> particles;
//    particleField->findParticles(particleField->getBoundingBox(),   particles);
//    cout << "(OrderedPositionCellField3D) "
//            << mpiRank << "Number of particles/nVertices " << particles.size()*1.0/nVertices
//            << " h " << particles.size()*1.0/nVertices * volume * 1.0/ (DeltaX*DeltaY*DeltaZ) << " (deleted:" << 0 << ")" << std::endl;
    delete mesh;
}


template<typename T, template<typename U> class Descriptor>
OrderedPositionCellField3D<T,Descriptor>* OrderedPositionCellField3D<T,Descriptor>::clone() const {
    return new OrderedPositionCellField3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void OrderedPositionCellField3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        modified[1+iField] = modif::dynamicVariables; // Particle fields.
    }
}

template<typename T, template<typename U> class Descriptor>
void OrderedPositionCellField3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (int iField = 0; iField < cellFields.size(); ++iField) {
        isWritten[1+iField] = true; // Particle fields.
    }

}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT OrderedPositionCellField3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}





/* ******** PrintParticles *********************************** */

template<typename T, template<typename U> class Descriptor>
void PrintParticles<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    int mpiSize = global::mpi().getSize();;
    int mpiRank = global::mpi().getRank();

    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);
    cout << "Yes, now we print!" << std::endl;
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(particleField.getBoundingBox(), particles);
    for (pluint iP = 0; iP < particles.size(); ++iP)
        cout << "(PRINT PARTICLES)" << mpiRank << " ("
             << particles[iP]->getPosition()[0] << ", "
             << particles[iP]->getPosition()[1] << ", "
             << particles[iP]->getPosition()[2] << ") "
             << std::endl;


}


template<typename T, template<typename U> class Descriptor>
PrintParticles<T,Descriptor>* PrintParticles<T,Descriptor>::clone() const {
    return new PrintParticles<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void PrintParticles<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
}

template<typename T, template<typename U> class Descriptor>
void PrintParticles<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT PrintParticles<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


#endif


