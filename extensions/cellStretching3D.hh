#if 0
#ifndef CELL_STRETCHING_3D_HH
#define CELL_STRETCHING_3D_HH

#include "cellStretching3D.h"

/* ================================================================================ */
/* ******** applyForceToCells *********************************** */
/* ================================================================================ */

template<typename T, template<typename U> class Descriptor>
void applyForceToCells(HemoCellField & cellField,
        std::vector<plint> const& cellIds, std::vector<std::vector<plint> > iVertices,
        std::vector<Array<T,3> > const& forces)
{
    applyProcessingFunctional (
        new ApplyForce3D<T,Descriptor>(cellField, cellIds, iVertices, forces),
        cellField.getBoundingBox(), cellField.getParticleArg() );
}


template<typename T, template<typename U> class Descriptor>
void applyForceToCells(HemoCellField & cellField,
        std::vector<plint> const& cellIds,
        std::vector<Array<T,3> > const& forces)
{
    applyProcessingFunctional (
        new ApplyForce3D<T,Descriptor>(cellField, cellIds, forces),
        cellField.getBoundingBox(), cellField.getParticleArg() );
}

/* ================================================================================ */
/* ******** ApplyForce3D *********************************** */
/* ================================================================================ */

template<typename T, template<typename U> class Descriptor>
ApplyForce3D<T,Descriptor>::ApplyForce3D (ApplyForce3D<T,Descriptor> const& rhs) :
    cellField(rhs.cellField), cellIds(rhs.cellIds), iVertices(rhs.iVertices), forces(rhs.forces)
{ }


template<typename T, template<typename U> class Descriptor>
ApplyForce3D<T,Descriptor>::ApplyForce3D (HemoCellField & cellField_,
      std::vector<plint> const& cellIds_, std::vector<std::vector<plint> > iVertices_,
      std::vector<Array<T,3> > const& forces_)
      :   cellField(cellField_), cellIds(cellIds_), iVertices(iVertices_), forces(forces_)
{
    PLB_ASSERT(  (cellIds.size() == iVertices.size()) && (cellIds.size() == forces.size()));
};

template<typename T, template<typename U> class Descriptor>
ApplyForce3D<T,Descriptor>::ApplyForce3D (HemoCellField & cellField_,
      std::vector<plint> const& cellIds_,
      std::vector<Array<T,3> > const& forces_)
      :   cellField(cellField_), cellIds(cellIds_), forces(forces_)
{
    pluint numVertices = cellField.getMesh().getNumVertices();
    for (pluint iCell = 0; iCell < cellIds.size(); ++iCell) {
        vector<plint> vertexVector;
        for (pluint iV = 0; iV < numVertices; ++iV) {
            vertexVector.push_back(iV);
        }
        iVertices.push_back(vertexVector);
    }
    PLB_ASSERT(  (cellIds.size() == iVertices.size()) && (cellIds.size() == forces.size()));
};


template<typename T, template<typename U> class Descriptor>
void ApplyForce3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );

    for (pluint i = 0; i < cellIds.size(); ++i) {
        plint cellId = cellIds[i];
        std::vector<plint> const& iiVertex= iVertices[i];
        Array<T,3> const& force = forces[i];
        if (cellField.count(cellId) > 0) {
            Cell3D<T,Descriptor> * cell = cellField[cellId];
            pluint nVertex = iiVertex.size();
            for (pluint iV = 0; iV < nVertex; ++iV) {
                plint iVertex = iiVertex[iV];
                if (cell->isInBulk(iVertex)) {
                    castParticleToICP3D( cell->getParticle3D(iVertex) )->get_force() += (force*1.0/(1.0*nVertex));
                }
            }
        }
    }

}

template<typename T, template<typename U> class Descriptor>
ApplyForce3D<T,Descriptor>*
    ApplyForce3D<T,Descriptor>::clone() const
{
    return new ApplyForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ApplyForce3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ApplyForce3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ApplyForce3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::allVariables; // Particle field.
}





#endif  // CELL_STRETCHING_FORCES_3D_HG

#endif