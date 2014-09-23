#ifndef CELLFIELD_3D_HH
#define CELLFIELD_3D_HH
#include "CellField3D.h"

template<typename T, template<typename U> class Descriptor>
CellField3D<T, Descriptor>::CellField3D(MultiBlockLattice3D<T, Descriptor> & lattice_, 
	TriangularSurfaceMesh<T> const& elementaryMesh_,  pluint numberOfCells_, 
	ConstitutiveModel<T, Descriptor> * cellModel_) : 
		lattice(lattice_), elementaryMesh(elementaryMesh_) 
{
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

    particleArg.push_back(immersedParticles);
  	particleLatticeArg.push_back(immersedParticles);
	particleLatticeArg.push_back(&lattice);
	/* Default values*/
	ibmKernel = 2;
	coupleWithIBM = true;

}

template<typename T, template<typename U> class Descriptor>
CellField3D<T, Descriptor>::~CellField3D() {
	delete [] cellModel; 
	delete [] immersedParticles;
}




template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::advanceParticles() {
    global::timer("Quantities").start();
    applyProcessingFunctional ( // advance particles in time according to velocity
        new AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>,
        immersedParticles->getBoundingBox(), particleArg );

    applyProcessingFunctional (
        new FillCellMap<T,Descriptor> (mesh, cellIdToCell3D),
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
    global::timer("IBM").start();
        setExternalVector( lattice, lattice.getBoundingBox(),
                       Descriptor<T>::ExternalField::forceBeginsAt, 
                       Array<T,3>(forceScalar,0.0,0.0) );
    global::timer("IBM").stop();
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
                    new ComputeImmersedElasticForce3D<T,Descriptor> (Cells, cellModel->clone(), RBCField),
                    immersedParticles->getBoundingBox(), particleArg );
    global::timer("Model").stop();
}


template<typename T, template<typename U> class Descriptor>
void CellField3D<T, Descriptor>::synchronizeCellQuantities_Local() {
    global::timer("Model").start();
    applyProcessingFunctional (
        new ComputeRequiredQuantities<T,Descriptor> (ccrRequirements, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    global::timer("Model").stop();

    global::timer("Quantities").start();
    applyProcessingFunctional (
        new SyncCellQuantities<T,Descriptor> (mesh, cellIdToCell3D),
        immersedParticles->getBoundingBox(), particleArg );
    global::timer("Quantities").stop();

}


#endif