#ifndef IMMERSED_CELLS_FUNCTIONAL_3D_HH
#define IMMERSED_CELLS_FUNCTIONAL_3D_HH

#include "immersedCellsFunctional3D.h"

#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"


/* ******** CountCellVolumeFunctional3D *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T>&  countCellVolume (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, plint numberOfCells )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    ComputeCellVolumeParticlesFunctional3D<T,Descriptor> functional(Cells, numberOfCells);
    applyProcessingFunctional(functional, domain, particleArg);
    return functional.getCellVolumeArray();
}



//template<typename T, template<typename U> class Descriptor>
//CountCellVolumeFunctional3D<T,Descriptor>::CountCellVolumeFunctional3D(plint tag_)
//    : cellVolumeId(this->getStatistics().subscribeIntSum()), tag(tag_)
//{ }
//
//template<typename T, template<typename U> class Descriptor>
//void CountCellVolumeFunctional3D<T,Descriptor>::processGenericBlocks (
//        Box3D domain, std::vector<AtomicBlock3D*> blocks )
//{
//    PLB_PRECONDITION( blocks.size()==1 );
//    ParticleField3D<T,Descriptor>& particleField
//        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
//    std::vector<Particle3D<T,Descriptor>*> particles;
//    particleField.findParticles(domain, particles);
//    T cellVolume = 0;
//    for (pluint iA = 0; iA < particles.size(); ++iA) {
//        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
//        ImmersedWallParticle3D<T,Descriptor>* particle =
//            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
//
//        if (particle->get_cellId() == tag) {
////         if (particle->getTag() == tag) {
//        	cellVolume += 0.;
//        }
//    }
//    this->getStatistics().gatherIntSum(cellVolumeId, cellVolume);
//}
//
//template<typename T, template<typename U> class Descriptor>
//CountCellVolumeFunctional3D<T,Descriptor>* CountCellVolumeFunctional3D<T,Descriptor>::clone() const {
//    return new CountCellVolumeFunctional3D<T,Descriptor>(*this);
//}
//
//template<typename T, template<typename U> class Descriptor>
//void CountCellVolumeFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
//    modified[0] = modif::nothing;
//}
//
//template<typename T, template<typename U> class Descriptor>
//T CountCellVolumeFunctional3D<T,Descriptor>::getVollumeCell() const {
//    return this->getStatistics().getSum(cellVolumeId);
//}
//
//template< typename T, template<typename U> class Descriptor,
//          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
//T countCellVolume (
//                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, plint tag)
//{
//    std::vector<MultiBlock3D*> particleArg;
//    particleArg.push_back(&particles);
//
//    CountCellVolumeFunctional3D<T,Descriptor> functional(tag);
//    applyProcessingFunctional(functional, domain, particleArg);
//    return functional.getVollumeCell();
//}
//
//
//

#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
