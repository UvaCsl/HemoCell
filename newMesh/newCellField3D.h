#ifndef CELLFIELD_3D_H
#define CELLFIELD_3D_H


particleEnvelopeWidth = 2;

template<typename T, template<typename U> class Descriptor>
class CellField3D
{
public:
	CellField3D(MultiBlockLattice3D<T, Descriptor> const& lattice_, TriangularSurfaceMesh<T> const& elementaryMesh_, 
			plint npar_, ConstitutiveModel<T, Descriptor> * cellModel_) : 
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
    	immersedParticles.periodicity().toggleAll(true);
	    immersedParticles.toggleInternalStatistics(false);

	    particleArg.push_back(&immersedParticles);
  	  	particleLatticeArg.push_back(&immersedParticles);
    	particleLatticeArg.push_back(&lattice);
    	/* Default values*/
    	ibmKernel = 2;
    	coupleWithIBM = true;

	}
	~CellField3D() { delete [] cellModel; };
public:
	void setIBMKernel(plint ibmKernel_) { ibmKernel = ibmKernel_; }
	void setIBMCoupling(bool coupleWithIBM_) { coupleWithIBM = coupleWithIBM_; }
    Cell3D<T,Descriptor> operator[](plint cellId) { return CellIdToCell3D[cellId]; }

public:
	void advanceParticles() {
	    applyProcessingFunctional ( // advance particles in time according to velocity
	        new AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>,
	        immersedParticles.getBoundingBox(), particleArg );
	}

	void spreadForceIBM() {
        global::timer("IBM").start();
        if (coupleWithIBM != 0) { // Force from the Cell dynamics to the Fluid
            applyProcessingFunctional ( // compute force applied on the fluid by the particles
                    new ForceToFluid3D<T,Descriptor> (ibmKernel),
                    immersedParticles.getBoundingBox(), particleLatticeArg );
        }
        global::timer("IBM").stop();
	}

	void interpolateVelocityIBM() {
        global::timer("IBM").start();
        if (coupleWithIBM != 0) { // Force from the Cell dynamics to the Fluid
	        applyProcessingFunctional ( // copy fluid velocity on particles
    	        new FluidVelocityToImmersedCell3D<T,Descriptor>(ibmKernel),
        	    immersedParticles.getBoundingBox(), particleLatticeArg);
	    }
        global::timer("IBM").stop();
	}
	void applyConstitutiveModel() {
        global::timer("Model").start();
        // #1# Membrane Model + Stretching
        applyProcessingFunctional (
                        new ComputeImmersedElasticForce3D<T,Descriptor> (Cells, cellModel->clone(), RBCField),
                        immersedParticles.getBoundingBox(), particleArg );
        global::timer("Model").stop();
	}


private:
	MultiBlockLattice3D<T, Descriptor> const& lattice
	TriangularSurfaceMesh<T> & elementaryMesh;
	ConstitutiveModel<T, Descriptor> * cellModel;
	plint ibmKernel;
	bool coupleWithIBM;

	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * immersedParticles;
    std::vector<MultiBlock3D*> particleArg;
    std::vector<MultiBlock3D*> particleLatticeArg;

    std::map<plint, Cell3D<T,Descriptor> > CellIdToCell3D;

	/* data */
};


#include "CellField3D.hh"
#endif