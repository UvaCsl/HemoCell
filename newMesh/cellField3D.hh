#ifndef CELLFIELD_3D_HH
#define CELLFIELD_3D_HH
#include "cellField3D.h"

template<typename T, template<typename U> class Descriptor>
CellField3D<T, Descriptor>::CellField3D(MultiBlockLattice3D<T, Descriptor> & lattice_, 
	TriangularSurfaceMesh<T> & elementaryMesh_,  pluint numberOfCells_,
	ConstitutiveModel<T, Descriptor> * cellModel_) : 
		lattice(lattice_), elementaryMesh(elementaryMesh_), cellModel(cellModel_)
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
    global::timer("Quantities").start();
	std::vector<Array<T,3> > cellOrigins;
    applyProcessingFunctional (
        new RandomPositionCellParticlesForGrowth3D<T,Descriptor>(elementaryMesh, .4),
        lattice.getBoundingBox(), particleLatticeArg );
    applyProcessingFunctional (
        new FillCellMap<T,Descriptor> (elementaryMesh, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    global::timer("Quantities").stop();
    advanceParticles();
//    synchronizeCellQuantities();
}


template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::advanceParticles() {
    global::timer("Quantities").start();
    applyProcessingFunctional ( // advance particles in time according to velocity
        new AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>,
        immersedParticles->getBoundingBox(), particleArg );

    applyProcessingFunctional (
        new FillCellMap<T,Descriptor> (elementaryMesh, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    global::timer("Quantities").stop();

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


