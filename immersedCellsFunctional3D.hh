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

/* ******** AdvanceCellsFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
AdvanceCellsFunctional3D<T,Descriptor>::AdvanceCellsFunctional3D (
        T cutOffValue_ )
  : cutOffValue(cutOffValue_)
{ }

template<typename T, template<typename U> class Descriptor>
void AdvanceCellsFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    particleField.advanceParticles(domain, cutOffValue);
}

template<typename T, template<typename U> class Descriptor>
AdvanceCellsFunctional3D<T,Descriptor>* AdvanceCellsFunctional3D<T,Descriptor>::clone() const {
    return new AdvanceCellsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AdvanceCellsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;  // Important: access envelope as well,
                                          // because particles are streamed from the
                                          // envelope into the bulk.
}

template<typename T, template<typename U> class Descriptor>
void AdvanceCellsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** FluidVelocityToImmersedCells3D *********************************** */


template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedCells3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    std::vector<Cell<T,Descriptor>*> cells(8);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iParticle];
        ImmersedWallParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedWallParticle3D<T,Descriptor>*> (nonTypedParticle);
        PLB_ASSERT( particle );
        Array<T,3> position(particle->getPosition());
        linearInterpolationCoefficients(fluid, position, cellPos, weights);

        // Use copy constructor in order to initialize dynamics object.
        Cell<T,Descriptor>* cellOnVertex;
        for (plint iCell=0; iCell<8; ++iCell) {
            cells[iCell] = &fluid.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z);
        }
        cellOnVertex = new Cell<T,Descriptor>(*cells[0]);
        for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
            (*cellOnVertex)[iPop] =
                weights[0]*(*cells[0])[iPop] + weights[1]*(*cells[1])[iPop] + weights[2]*(*cells[2])[iPop] +
                weights[3]*(*cells[3])[iPop] + weights[4]*(*cells[4])[iPop] + weights[5]*(*cells[5])[iPop] +
                weights[6]*(*cells[6])[iPop] + weights[7]*(*cells[7])[iPop];
        }
        cellOnVertex->computeVelocity(particle->get_v());
        delete cellOnVertex;
    }
}

template<typename T, template<typename U> class Descriptor>
FluidVelocityToImmersedCells3D<T,Descriptor>* FluidVelocityToImmersedCells3D<T,Descriptor>::clone() const {
    return new FluidVelocityToImmersedCells3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FluidVelocityToImmersedCells3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing; // Particle field.
}




#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
