#ifndef CELLFIELD_3D_HH
#define CELLFIELD_3D_HH
#include "cellField3D.h"
#include "initializationCellField3D.h"


template<typename T, template<typename U> class Descriptor>
CellField3D<T, Descriptor>::CellField3D(
    MultiBlockLattice3D<T, Descriptor> & lattice_, 
	  TriangularSurfaceMesh<T> & elementaryMesh_, 
    T hematocrit_,
	  ShellModel3D<T> * cellModel_, 
    std::string identifier_) 
    :
		lattice(lattice_), elementaryMesh(elementaryMesh_),
		hematocrit(hematocrit_), cellModel(cellModel_), identifier(identifier_)
{
    plint maxEdgeLengthLU = ceil(cellModel->getMaximumEdgeExtensionLengthLU());
    plint maxCellDiameterLU = ceil(cellModel->getMaxCellDiameterLU());
    pluint particleEnvelopeWidth = maxEdgeLengthLU;
    pluint reductionParticleEnvelopeWidth = maxCellDiameterLU;

    pcout << "(CellField3D) " << identifier << "-> particle envelope [lu]: " << particleEnvelopeWidth << std::endl;
    pcout << "(CellField3D) " << identifier << "-> reduction particle envelope width [lu]: " << reductionParticleEnvelopeWidth << std::endl;

    MultiBlockManagement3D const& latticeManagement(lattice.getMultiBlockManagement());
	MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            particleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );

    immersedParticles = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(
            particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
    immersedParticles->periodicity().toggleAll(true);
    immersedParticles->toggleInternalStatistics(false);

    MultiBlockManagement3D reductionParticleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            reductionParticleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );

    reductionParticles = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(reductionParticleManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );
    reductionParticles->periodicity().toggleAll(true);
    reductionParticles->toggleInternalStatistics(false);

    ccrRequirements.insert(cellModel->getSyncRequirements());
    ccrRequirements.insert(CCR_NO_PBC_POSITION_MEAN);
    particleArg.clear(); particleLatticeArg.clear(); particleReductionParticleArg.clear();
    particleArg.push_back(immersedParticles);
    particleLatticeArg.push_back(immersedParticles);
    particleLatticeArg.push_back(&lattice);
    particleReductionParticleArg.push_back(immersedParticles);
    particleReductionParticleArg.push_back(reductionParticles);
    reductionParticleArg.push_back(reductionParticles);
	/* Default values*/
	coupleWithIBM = true;

}

template<typename T, template<typename U> class Descriptor>
CellField3D<T, Descriptor>::~CellField3D() {
    typename std::map<plint, Cell3D<T,Descriptor>* >::iterator iter;
    for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
        delete (iter->second);
    }
    cellIdToCell3D.clear();

    delete cellModel;
	delete immersedParticles;
    delete reductionParticles;
}



template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::initialize(std::vector<Array<T,3> > & centers) {
    global::timer("CellInit").start();
    applyProcessingFunctional (
        new PositionCellParticles3D<T,Descriptor>(elementaryMesh, centers),
        immersedParticles->getBoundingBox(), particleArg );
    advanceParticles();
    synchronizeCellQuantities();
    global::timer("CellInit").stop();
}


template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::grow(plint growIterations) {
    global::timer("CellInit").start();
    T ratio;
    coupleWithIBM = false;
    T k_int = 0.25, DeltaX=1.0, R=3.75, k=1.5;
    PowerLawForce<T> PLF(k_int, DeltaX, R, k);
    plint nIter = (1-ratio)*growIterations;

    applyProcessingFunctional (
        new RandomPositionCellParticlesForGrowth3D<T,Descriptor>(elementaryMesh, hematocrit, ratio),
        lattice.getBoundingBox(), particleLatticeArg );
    advanceParticles();
    synchronizeCellQuantities();
    for (plint i = 0; i < growIterations; ++i) {
        T iRatio = 1.0; 
        iRatio = ratio + i*1.0 / (growIterations*0.5) ;
        if (iRatio > 1.0) { iRatio = 1.0; }
        if (i%100 == 0) {
            pcout << "growth iter:" << i<< ", " <<  i*100.0/growIterations << "%" <<std::endl;
//            writeCellField3D_HDF5(*this, 1.0, 1.0, i, "init_");
        }
        // #1# Calculate forces
        applyConstitutiveModel();
        applyCellCellForce(PLF, R); // #1b# Repulsive force
        // #2# Force to velocity
        computeVelocity(ratio);
        // #3# Update position of particles
        advanceParticles();
        synchronizeCellQuantities();
    }
    coupleWithIBM = true;
    global::timer("CellInit").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::createCellMap() {
    global::timer("Quantities").start();
    applyProcessingFunctional (
        new FillCellMap<T,Descriptor> (elementaryMesh, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    global::timer("Quantities").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::advanceParticles() {
    global::timer("IBM").start();
    applyProcessingFunctional ( // advance particles in time according to velocity
        new AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>,
        immersedParticles->getBoundingBox(), particleArg );
    global::timer("IBM").stop();
    createCellMap();

}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::spreadForceIBM() {
    global::timer("IBM").start();
    if (coupleWithIBM != 0) { // Force from the Cell dynamics to the Fluid
        applyProcessingFunctional ( // compute force applied on the fluid by the particles
                new ForceToFluid3D<T,Descriptor> (),
                immersedParticles->getBoundingBox(), particleLatticeArg );
    }
    global::timer("IBM").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::setFluidExternalForce(Array<T,3> force) {
    global::timer("IBM").start();
        setExternalVector( lattice, lattice.getBoundingBox(),
                       Descriptor<T>::ExternalField::forceBeginsAt, 
                       force);
    global::timer("IBM").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::setFluidExternalForce(T forceScalar) {
    setFluidExternalForce(Array<T,3>(forceScalar,0.0,0.0) );
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::interpolateVelocityIBM() {
    global::timer("IBM").start();
    if (coupleWithIBM != 0) { // Force from the Cell dynamics to the Fluid
        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedCell3D<T,Descriptor>(),
            immersedParticles->getBoundingBox(), particleLatticeArg);
    }
    global::timer("IBM").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::computeVelocity(T ratio) {
    interpolateVelocityIBM();
//    if (coupleWithIBM != 0) { // Force from the Cell dynamics to the Fluid
//        interpolateVelocityIBM();
//    } else {
//        global::timer("IBM").start();
//        applyProcessingFunctional ( // copy fluid velocity on particles
//            new ViscousPositionUpdate3D<T,Descriptor>(ratio),
//            immersedParticles->getBoundingBox(), particleLatticeArg);
//        global::timer("IBM").stop();
//    }
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::applyConstitutiveModel() {
    global::timer("Model").start();
    // #1# Membrane Model + Stretching
    applyProcessingFunctional (
                    new ComputeCellForce3D<T,Descriptor> (cellModel, cellIdToCell3D),
                    immersedParticles->getBoundingBox(), particleArg );
    global::timer("Model").stop();
}


template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::applyCellCellForce(CellCellForce3D<T> & calcForce, T cutoffRadius) {
    global::timer("CellCellForce").start();
    applyProcessingFunctional (
       new ComputeCellCellForces3D<T,Descriptor> (calcForce, cutoffRadius),
       immersedParticles->getBoundingBox(), particleArg );
    global::timer("CellCellForce").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::applyWallCellForce(CellCellForce3D<T> & calcForce, T cutoffRadius, MultiParticleField3D<DenseParticleField3D<T,Descriptor> > *wallParticles) {
    std::vector<MultiBlock3D*> wallParticleArg;
    wallParticleArg.push_back(wallParticles);
    wallParticleArg.push_back(immersedParticles);
    global::timer("WallCellForce").start();
    applyProcessingFunctional (
       new ComputeWallCellForces3D<T,Descriptor> (calcForce, cutoffRadius),
       immersedParticles->getBoundingBox(), wallParticleArg );
    global::timer("WallCellForce").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::applyDifferentCellForce(CellCellForce3D<T> & calcForce, T cutoffRadius,
        MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * otherCellParticles) {

    std::vector<MultiBlock3D*> differentCellCellParticleArg;
    differentCellCellParticleArg.push_back(otherCellParticles);
    differentCellCellParticleArg.push_back(immersedParticles);
    global::timer("CellCellForce").start();
    applyProcessingFunctional (
       new ComputeDifferentCellForces3D<T,Descriptor> (calcForce, cutoffRadius),
       immersedParticles->getBoundingBox(), differentCellCellParticleArg );
    global::timer("CellCellForce").stop();
}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::deleteIncompleteCells() {
	createCellMap();
    SyncRequirements ccReq;
    ccReq.insert(CCR_NO_PBC_POSITION_MEAN);
    applyProcessingFunctional (
        new ComputeRequiredQuantities<T,Descriptor> (ccReq.getSyncRequirements(), cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleReductionParticleArg );

    std::vector<MultiBlock3D*> reductionParticleParticleArg;
    reductionParticleParticleArg.push_back(reductionParticles);
    reductionParticleParticleArg.push_back(immersedParticles);
	plint Nv = getMesh().getNumVertices() ;
    applyProcessingFunctional (
        new DeleteIncompleteCells<T,Descriptor> (cellIdToCell3D, Nv ),
        reductionParticles->getBoundingBox(), reductionParticleParticleArg );

}



template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::synchronizeSyncRequirements_Local(SyncRequirements ccReq) {
    ccReq.insert(CCR_NO_PBC_POSITION_MEAN);
    SyncRequirements ccrIndependent, ccrDependent;
    separateDependencies(ccReq, ccrIndependent, ccrDependent);

    // First do the independent reductions (volume surface)...
    global::timer("Model").start();
    applyProcessingFunctional (
        new ComputeRequiredQuantities<T,Descriptor> (ccrIndependent.getSyncRequirements(), cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleReductionParticleArg );
    global::timer("Model").stop();

    global::timer("Quantities").start();
    applyProcessingFunctional (
        new SyncCellQuantities<T,Descriptor> (cellIdToCell3D),
        reductionParticles->getBoundingBox(), reductionParticleArg );
    global::timer("Quantities").stop();

    // ... and then the ones that depend on something else (Inertia etc)
    if (ccrDependent.size() > 0) {
        global::timer("Model").start();
        applyProcessingFunctional (
            new ComputeRequiredQuantities<T,Descriptor> (ccrDependent.getSyncRequirements(), cellIdToCell3D),
            immersedParticles->getBoundingBox(), particleReductionParticleArg );
        global::timer("Model").stop();

        global::timer("Quantities").start();
        applyProcessingFunctional (
            new SyncCellQuantities<T,Descriptor> (cellIdToCell3D),
            reductionParticles->getBoundingBox(), reductionParticleArg );
        global::timer("Quantities").stop();
    }

}

template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::synchronizeCellQuantities_Local(SyncRequirements ccReq) {
    ccReq.insert(ccrRequirements);
    synchronizeSyncRequirements_Local(ccReq);
}


#endif


