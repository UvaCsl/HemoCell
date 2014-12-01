#ifndef INITIALIZATION_CELL_FIELD_3D_HH
#define INITIALIZATION_CELL_FIELD_3D_HH
#include "initializationCellField3D.h"



template<typename T, template<typename U> class Descriptor>
void randomPositionCellFieldsForGrowth3D(std::vector<CellField3D<T, Descriptor>* > & cellFields) {
    global::timer("CellInit").start();
    std::vector<MultiBlock3D*> fluidAndParticleFieldsArg;
    fluidAndParticleFieldsArg.push_back( &(cellFields[0]->getFluidField3D()) );
    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back( &(cellFields[icf]->getParticleField3D()) );
    }

    applyProcessingFunctional (
        new RandomPositionCellFieldsForGrowth3D<T,Descriptor>(cellFields),
        cellFields[0]->getFluidField3D().getBoundingBox(), fluidAndParticleFieldsArg );

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
            particleField.addParticle(domain, new ImmersedCellParticle3D<T,Descriptor>(vertex, cellId, iVertex) );
        }
    }
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);

    cout << "(PositionCellParticles3D) pid " << global::mpi().getRank()
        << " particles.size " << particles.size()
        << " nVertices " << nVertices
        << " cellOrigin.size " << cellOrigins.size() << std::endl;
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

    plint DeltaX = domain.x1-domain.x0 + 1;
    plint DeltaY = domain.y1-domain.y0 + 1;
    plint DeltaZ = domain.z1-domain.z0 + 1;

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
                            particleField.addParticle(domain, new ImmersedCellParticle3D<T,Descriptor>(vertex, cellId + 1000*mpiRank, iVertex) );
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
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iP]);
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
    int mpiSize = global::mpi().getSize();;
    int mpiRank = global::mpi().getRank();

    T hematocrit = cellFields[0]->getVolumeFraction();
    TriangularSurfaceMesh<T> & elementaryMesh = cellFields[0]->getMesh();
    T ratio;
    if (hematocrit > 1.0) { hematocrit /= 100; }
    PLB_PRECONDITION( hematocrit < 1.0 );

    srand (mpiSize * mpiRank);
    PLB_PRECONDITION( blocks.size()==(cellFields.size()+1) );
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);
    particleField.removeParticles(domain);

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

    Dot3D pfLocation(fluid.getLocation());
    Box3D realDomain(
                        domain.x0 + pfLocation.x, domain.x1 + pfLocation.x,
                        domain.y0 + pfLocation.y, domain.y1 + pfLocation.y,
                        domain.z0 + pfLocation.z, domain.z1 + pfLocation.z );
//    bool particleIsInBulk = particleField.isContained(particle->getPosition(), realDomain);

    plint DeltaX = domain.x1-domain.x0 + 1;
    plint DeltaY = domain.y1-domain.y0 + 1;
    plint DeltaZ = domain.z1-domain.z0 + 1;

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
    Dot3D particlePositionDot3D = particleField.getLocation();
    Array<T,3> particlePosition(particlePositionDot3D.x, particlePositionDot3D.y, particlePositionDot3D.z);
    Dot3D fluidLocationDot3D = fluid.getLocation();
    Array<T,3> fluidLocation(fluidLocationDot3D.x, fluidLocationDot3D.y, fluidLocationDot3D.z);

    plint nVertices = mesh->getNumVertices();
    plint cellId = 0;
    cout << global::mpi().getRank()
            << "domain "
                << domain.x0 << " " << domain.x1 << ", "
                << domain.y0 << " " << domain.y1 << ", "
                << domain.z0 << " " << domain.z1 << ", "
            << "fluidLocation "
                << fluidLocationDot3D.x << ", "
                << fluidLocationDot3D.y << ", "
                << fluidLocationDot3D.z << ", "
            << "particlePosition "
                << particlePositionDot3D.x << ", "
                << particlePositionDot3D.y << ", "
                << particlePositionDot3D.z << ", "
            << std::endl;


    Box3D newDomain(domain.x0-100, domain.x1+100, domain.y0-100, domain.y1+100, domain.z0-100, domain.z1+100);
    // Loop through the domain and place cell depending on the
    for (plint iX=newDomain.x0; iX<=newDomain.x1-step; iX+=step) {
        for (plint iY=newDomain.y0; iY<=newDomain.y1-step; iY+=step) {
            for (plint iZ=newDomain.z0; iZ<=newDomain.z1-step; iZ+=step) {
                T rn = guessRandomNumber();
                if (rn <= prob) {
                    meshRandomRotation(mesh);
                    Array<T,3> center = Array<T,3>(iX*1.0, iY*1.0, iZ*1.0) + fluidLocation;
                    for (plint iVertex=0; iVertex < nVertices; ++iVertex) {
                        Array<T,3> vertex = center + mesh->getVertex(iVertex);
                        Dot3D cellInFluidDomain = Dot3D(vertex[0], vertex[1], vertex[2]) - fluidLocationDot3D;
//                        if (fluid.get(cellInFluidDomain.x, cellInFluidDomain.y, cellInFluidDomain.z).getDynamics().isBoundary()) {
//                            cout << global::mpi().getRank() << "vertex (" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")" <<std::endl;
//                            break;
//                        } else {
                            particleField.addParticle(newDomain, new ImmersedCellParticle3D<T,Descriptor>(vertex, cellId + 1000*mpiRank, iVertex) );
//                        }

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
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iP]);
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
    cout << "(RandomPositionCellParticlesForGrowth3D) "
            << mpiRank << "Number of particles/nVertices " << particles.size()*1.0/nVertices
            << " scale " << scale << " step " << step
            << " prob " << prob
            << " h " << particles.size()*1.0/nVertices * volume * 1.0/ (DeltaX*DeltaY*DeltaZ) << " (deleted:" << cellsDeleted << ")" << std::endl;
}


template<typename T, template<typename U> class Descriptor>
RandomPositionCellFieldsForGrowth3D<T,Descriptor>* RandomPositionCellFieldsForGrowth3D<T,Descriptor>::clone() const {
    return new RandomPositionCellFieldsForGrowth3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellFieldsForGrowth3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing; // Fluid field.
}

template<typename T, template<typename U> class Descriptor>
void RandomPositionCellFieldsForGrowth3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
    isWritten[1] = false;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT RandomPositionCellFieldsForGrowth3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



#endif


