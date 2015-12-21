#ifndef FICSION_HELPER_FUNCTIONALS_H
#define FICSION_HELPER_FUNCTIONALS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "surfaceParticle3D.h"


//    fluid.get(1, 1, 1).getDynamics().isBoundary();
template<typename T, template<typename U> class Descriptor>
class PositionBoundaryParticles : public BoxProcessingFunctional3D
{
public:
    PositionBoundaryParticles ();
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual PositionBoundaryParticles<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
};



template<typename T, template<typename U> class Descriptor>
MultiParticleField3D<DenseParticleField3D<T,Descriptor> > *
		createBoundaryParticleField3D(MultiBlockLattice3D<T, Descriptor> & lattice) {

	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * boundaryParticleField3D
				= new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(lattice);
	boundaryParticleField3D->periodicity().toggleAll(true);
	boundaryParticleField3D->toggleInternalStatistics(false);

	std::vector<MultiBlock3D*> particleLatticeArg;
    particleLatticeArg.push_back(boundaryParticleField3D);
    particleLatticeArg.push_back(&lattice);


    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new PositionBoundaryParticles<T,Descriptor> (),
            boundaryParticleField3D->getBoundingBox(), particleLatticeArg );

	return boundaryParticleField3D;
}


template<typename T, template<typename U> class Descriptor>
MultiParticleField3D<DenseParticleField3D<T,Descriptor> > *
		createBoundaryParticleField3D(MultiBlockLattice3D<T, Descriptor> & lattice, Box3D domain) {

	MultiParticleField3D<DenseParticleField3D<T,Descriptor> > * boundaryParticleField3D
				= new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(lattice);
	boundaryParticleField3D->periodicity().toggleAll(true);
	boundaryParticleField3D->toggleInternalStatistics(false);

	std::vector<MultiBlock3D*> particleLatticeArg;
    particleLatticeArg.push_back(boundaryParticleField3D);
    particleLatticeArg.push_back(&lattice);


    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new PositionBoundaryParticles<T,Descriptor> (),
            domain, particleLatticeArg );

	return boundaryParticleField3D;
}




#include "helperFunctionals.hh"
#endif  // CELL_3D_HH

