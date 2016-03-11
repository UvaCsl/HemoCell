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

//    fluid.get(1, 1, 1).getDynamics().isBoundary();
template<typename T, template<typename U> class Descriptor>
class DeleteParticles3D : public BoxProcessingFunctional3D
{
public:
	DeleteParticles3D ();
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual DeleteParticles3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
};


template<typename T, template<typename U> class Descriptor>
void addParticlesToBoundaryParticleField3D(
				MultiBlockLattice3D<T, Descriptor> & lattice,
				MultiParticleField3D<LightParticleField3D<T,Descriptor> > & boundaryParticleField3D,
				Box3D domain) {
	std::vector<MultiBlock3D*> particleLatticeArg;
    particleLatticeArg.push_back(&boundaryParticleField3D);
    particleLatticeArg.push_back(&lattice);

    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new PositionBoundaryParticles<T,Descriptor> (),
            domain, particleLatticeArg );
}


template<typename T, template<typename U> class Descriptor>
void removeParticlesFromBoundaryParticleField3D(
				MultiBlockLattice3D<T, Descriptor> & lattice,
				MultiParticleField3D<LightParticleField3D<T,Descriptor> > & boundaryParticleField3D,
				Box3D domain) {
	std::vector<MultiBlock3D*> particleLatticeArg;
    particleLatticeArg.push_back(&boundaryParticleField3D);
    particleLatticeArg.push_back(&lattice);

    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new DeleteParticles3D<T,Descriptor> (),
            domain, particleLatticeArg );
}



template<typename T, template<typename U> class Descriptor>
MultiParticleField3D<LightParticleField3D<T,Descriptor> > *
		createBoundaryParticleField3D(MultiBlockLattice3D<T, Descriptor> & lattice, Box3D domain) {

	MultiParticleField3D<LightParticleField3D<T,Descriptor> > * boundaryParticleField3D
				= new MultiParticleField3D<LightParticleField3D<T,Descriptor> >(lattice);
	boundaryParticleField3D->periodicity().toggleAll(true);
	boundaryParticleField3D->toggleInternalStatistics(false);
    addParticlesToBoundaryParticleField3D(lattice, *boundaryParticleField3D, domain);
	return boundaryParticleField3D;
}


template<typename T, template<typename U> class Descriptor>
MultiParticleField3D<LightParticleField3D<T,Descriptor> > *
		createBoundaryParticleField3D(MultiBlockLattice3D<T, Descriptor> & lattice) {
	return createBoundaryParticleField3D(lattice, lattice.getBoundingBox() );
}



#include "helperFunctionals.hh"
#endif  // CELL_3D_HH

