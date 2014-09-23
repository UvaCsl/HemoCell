#ifndef CELL_FIELD_FUNCTIONALS_3D_HH
#define CELL_FIELD_FUNCTIONALS_3D_HH
#include "cellFieldFunctionals3D.h"




/* ******** FillCellMap *********************************** */

template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>::FillCellMap (
	TriangularSurfaceMesh<T>& mesh_,
	std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D_)
: mesh(mesh_), cellIdToCell3D(cellIdToCell3D_)
{ }

template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>::~FillCellMap()
{
}

template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>::FillCellMap (
            FillCellMap<T,Descriptor> const& rhs)
    : cellIdToCell3D(rhs.cellIdToCell3D), mesh(rhs.mesh)
{ }


template<typename T, template<typename U> class Descriptor>
void FillCellMap<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    cellIdToCell3D.clear();

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(particleField.getBoundingBox(), found); // Gets the whole domain.

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        ImmersedCellParticle3D<T,Descriptor>* particle = 
        		dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (found[iParticle]);

        plint cellId = particle->get_cellId();
        if (cellIdToCell3D.count(cellId) == 0) {
        	cellIdToCell3D[cellId] = Cell3D<T,Descriptor>(mesh, cellId);
        }
        cellIdToCell3D[cellId].push_back(particle);
        // if (tagToParticle3D.count(vertexId) == 0) {
        //     tagToParticle3D[vertexId] = nonTypedParticle;
        // } else if (contained(nonTypedParticle->getPosition(), domain)) {
        //     tagToParticle3D[vertexId] = nonTypedParticle;
        // }
    }

    typename std::map<plint, T >::iterator iter;
    for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
    	(iter->second).close();
    }
}

template<typename T, template<typename U> class Descriptor>
FillCellMap<T,Descriptor>*
    FillCellMap<T,Descriptor>::clone() const
{
    return new FillCellMap<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FillCellMap<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT FillCellMap<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk; // This changes in processGenericBlocks with getBoundingBox();
}

template<typename T, template<typename U> class Descriptor>
void FillCellMap<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Particle field.
}








#endif  // CELL_FIELD_FUNCTIONALS_3D_HH

