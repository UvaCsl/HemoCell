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


/// Count the number of particles, no matter which kind, found inside the domain.
template<typename T, template<typename U> class Descriptor>
class CountCellVolumeFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CountCellVolumeFunctional3D(plint tag_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountCellVolumeFunctional3D<T,Descriptor>* clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint cellVolumeId, tag;
};


#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
