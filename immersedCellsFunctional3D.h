#ifndef IMMERSED_CELLS_FUNCTIONAL_3D_H
#define IMMERSED_CELLS_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "immersedWallParticle3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "shellModel3D.h"
#include <map>

/* ******** AdvanceCellsFunctional3D *********************************** */

/// Execute the iteration step during which particles advance.
template<typename T, template<typename U> class Descriptor>
class AdvanceCellsFunctional3D : public BoxProcessingFunctional3D
{
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
	AdvanceCellsFunctional3D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual AdvanceCellsFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T cutOffValue;
};

/* ******** FluidVelocityToImmersedCells3D *********************************** */

template<typename T, template<typename U> class Descriptor>
class FluidVelocityToImmersedCells3D : public BoxProcessingFunctional3D
{
public:
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FluidVelocityToImmersedCells3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};


#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
