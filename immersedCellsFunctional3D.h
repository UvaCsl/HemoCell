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

/* ******** countCellVolume *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVolume (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVolumes);

/* ******** ComputeCellVolumeParticlesFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class ComputeCellVolumeParticlesFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    ComputeCellVolumeParticlesFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeCellVolumeParticlesFunctional3D<T,Descriptor>* clone() const;
    void  getCellVolumeArray(std::vector<T>& cellVolumes, std::vector<plint> cellIds) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint numberOfCells;
    TriangleBoundary3D<T> const& triangleBoundary;
};


//
///// Count the number of particles, no matter which cellId, found inside the domain.
//template<typename T, template<typename U> class Descriptor>
//class CountCellVolumeFunctional3D : public PlainReductiveBoxProcessingFunctional3D
//{
//public:
//    CountCellVolumeFunctional3D(plint tag_);
//    /// Argument: Particle-field.
//    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
//    virtual CountCellVolumeFunctional3D<T,Descriptor>* clone() const;
//    T getVollumeCell() const;
//    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
//private:
//    plint cellVolumeId, tag;
//};
//


#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
