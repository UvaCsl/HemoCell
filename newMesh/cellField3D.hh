#ifndef CELLFIELD_3D_HH
#define CELLFIELD_3D_HH
#include "cellField3D.h"

template<typename T, template<typename U> class Descriptor>
CellField3D<T, Descriptor>::CellField3D(MultiBlockLattice3D<T, Descriptor> & lattice_, 
	TriangularSurfaceMesh<T> & elementaryMesh_, T hematocrit_,
	ConstitutiveModel<T, Descriptor> * cellModel_, std::string identifier_) :
		lattice(lattice_), elementaryMesh(elementaryMesh_),
		hematocrit(hematocrit_), cellModel(cellModel_), identifier(identifier_)
{
    plint maxEdgeLengthLU = ceil(cellModel->getMaximumEdgeExtensionLengthLU());
    plint maxRadiusLU = ceil(cellModel->getCellRadiusLU());
    pluint particleEnvelopeWidth = maxEdgeLengthLU;
    pluint reductionParticleEnvelopeWidth = 5*maxRadiusLU;
    pcout << "particleEnvelopeWidth " << particleEnvelopeWidth << std::endl;
    pcout << "reductionParticleEnvelopeWidth " << reductionParticleEnvelopeWidth << std::endl;
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
    particleArg.clear(); particleLatticeArg.clear(); particleReductioParticleArg.clear();
    particleArg.push_back(immersedParticles);
    particleLatticeArg.push_back(immersedParticles);
    particleLatticeArg.push_back(&lattice);
    particleReductioParticleArg.push_back(immersedParticles);
    particleReductioParticleArg.push_back(reductionParticles);
    reductionParticleArg.push_back(reductionParticles);
	/* Default values*/
	ibmKernel = 2;
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
void CellField3D<T, Descriptor>::initialize() {
    global::timer("CellInit").start();
    std::vector<Array<T,3> > cellOrigins;
    applyProcessingFunctional (
//        new RandomPositionCellParticlesForGrowth3D<T,Descriptor>(elementaryMesh, .4),
        new PositionCellParticles3D<T,Descriptor>(elementaryMesh, cellOrigins),
        lattice.getBoundingBox(), particleLatticeArg );
    applyProcessingFunctional (
        new FillCellMap<T,Descriptor> (elementaryMesh, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    advanceParticles();
    synchronizeCellQuantities();
    global::timer("CellInit").stop();
}


template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::grow() {
    global::timer("CellInit").start();
    T ratio;
    applyProcessingFunctional (
        new RandomPositionCellParticlesForGrowth3D<T,Descriptor>(elementaryMesh, hematocrit, ratio),
        lattice.getBoundingBox(), particleLatticeArg );
    applyProcessingFunctional (
        new FillCellMap<T,Descriptor> (elementaryMesh, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    writeCellField3D_HDF5(*this, 1.0, 1.0, 0);
    synchronizeCellQuantities();

    T k_int = 0.00025, DeltaX=1.0, R=0.75, k=1.5;
    PowerLawForce<T> PLF(k_int, DeltaX, R, k);
    plint nIter = (1-ratio)*10000 + 1000;
    for (plint i = 0; i < nIter; ++i) {
        T iRatio = ratio + i*1.0 / 10000.0 ;
        if (iRatio > 1.0) { iRatio = 1.0; }
        // #1# Calculate forces
        applyProcessingFunctional ( // #1a# Membrane Model
                        new ComputeCellForce3D<T,Descriptor> (cellModel, cellIdToCell3D),
                        immersedParticles->getBoundingBox(), particleArg );
        applyProcessingFunctional (  // #1b# Repulsive force
           new ComputeCellCellForces3D<T,Descriptor> (PLF, R),
           immersedParticles->getBoundingBox(), particleArg );

        // #2# Force to velocity
        applyProcessingFunctional ( // copy fluid velocity on particles
            new ViscousPositionUpdate3D<T,Descriptor>(iRatio),
            immersedParticles->getBoundingBox(), particleLatticeArg);
        // #3# Update position of particles
        applyProcessingFunctional ( // advance particles in time according to velocity
            new AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>,
            immersedParticles->getBoundingBox(), particleArg );
        applyProcessingFunctional (
            new FillCellMap<T,Descriptor> (elementaryMesh, cellIdToCell3D),
            immersedParticles->getBoundingBox(), particleArg );
        synchronizeCellQuantities();
        if (i%100 == 0) {
            pcout << "growth iter:" << i<<std::endl;
            writeCellField3D_HDF5(*this, 1.0, 1.0, i);
        }
    }

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
                new ForceToFluid3D<T,Descriptor> (ibmKernel),
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
	        new FluidVelocityToImmersedCell3D<T,Descriptor>(ibmKernel),
    	    immersedParticles->getBoundingBox(), particleLatticeArg);
    }
    global::timer("IBM").stop();
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
void CellField3D<T, Descriptor>::synchronizeCellQuantities_Local(SyncRequirements ccReq) {
    ccReq.insert(ccrRequirements);
    SyncRequirements ccrIndependent, ccrDependent;
    separateDependencies(ccReq, ccrIndependent, ccrDependent);

    // First do the independent reductions (volume surface)...
    global::timer("Model").start();
    applyProcessingFunctional (
        new ComputeRequiredQuantities<T,Descriptor> (ccrIndependent.getSyncRequirements(), cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleReductioParticleArg );
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
            immersedParticles->getBoundingBox(), particleReductioParticleArg );
        global::timer("Model").stop();

        global::timer("Quantities").start();
        applyProcessingFunctional (
            new SyncCellQuantities<T,Descriptor> (cellIdToCell3D),
            reductionParticles->getBoundingBox(), reductionParticleArg );
        global::timer("Quantities").stop();
    }

}


#endif


