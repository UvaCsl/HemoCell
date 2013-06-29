#ifndef CELL_IN_SHEAR_FLOW_3D_H
#define CELL_IN_SHEAR_FLOW_3D_H


#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedCellsReductions.h"
#include <map>
#include <algorithm>
#include "external/diagonalize.cpp"

namespace plb {

/* ******** InertiaTensorCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class InertiaTensorCellReduceFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    InertiaTensorCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_);
    InertiaTensorCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_,
                                        std::vector< Array<T,3> >& cellCenters_);
    /// Argument: Particle-field.
    virtual InertiaTensorCellReduceFunctional3D<T,Descriptor>* clone() const { return new InertiaTensorCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual void getCellQuantityArray(std::vector< std::vector<T> > & cellQuantity, std::vector<plint> cellIds_) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::map<plint, Array<plint,9> > & quantityIds_);
protected:
    plint numberOfCells;
    std::map<plint, Array<plint,9> > quantityIds;
    std::map<plint, Array<T,3> > cellCenters;
    TriangleBoundary3D<T> const& triangleBoundary;
    bool findMax;
};


/* ******** computeCellInertia *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeCellInertia (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,3> >& cellCenters, std::vector< Array<T,9> >& cellInertia);

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeCellInertia (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,9> >& cellInertia);
}

#include "cellInShearFlow3D.hh"

#endif // CELL_IN_SHEAR_FLOW_3D_H
