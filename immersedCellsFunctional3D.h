#ifndef IMMERSED_CELLS_FUNCTIONAL_3D_H
#define IMMERSED_CELLS_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "immersedCellParticle3D.h"
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
    std::vector<plint> volumeIds;
};


/// Remove all particles of a certain tag from a given domain.
template<typename T, template<typename U> class Descriptor>
class TranslateTaggedParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    TranslateTaggedParticlesFunctional3D(Array<T,3> const& translation_, plint tag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    /// Argument: Particle-field.
    virtual TranslateTaggedParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint tag;
	Array<T,3> const& translation;
};





/* ******** TranslateTaggedParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
TranslateTaggedParticlesFunctional3D<T,Descriptor>::
    TranslateTaggedParticlesFunctional3D(Array<T,3> const& translation_, plint tag_) : tag(tag_), translation(translation_)
{ }

template<typename T, template<typename U> class Descriptor>
void TranslateTaggedParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<plint> cellIds;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iP];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);

        if (particle->get_cellId() == tag) {
            particle->getPosition() = particle->getPosition() + translation;
        }
    }

//     pcout << cellIds.size() << std::endl;
//     particleField.removeParticles(domain,tag);
}

template<typename T, template<typename U> class Descriptor>
TranslateTaggedParticlesFunctional3D<T,Descriptor>* TranslateTaggedParticlesFunctional3D<T,Descriptor>::clone() const {
    return new TranslateTaggedParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void TranslateTaggedParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
