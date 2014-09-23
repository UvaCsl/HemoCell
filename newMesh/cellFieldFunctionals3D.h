#ifndef CELL_FIELD_FUNCTIONALS_3D_H
#define CELL_FIELD_FUNCTIONALS_3D_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "immersedCellParticle3D.h"
#include "cell3D.h"

template<typename T, template<typename U> class Descriptor>
class FillCellMap : public BoxProcessingFunctional3D
{
public:
    FillCellMap (TriangularSurfaceMesh<T>& mesh_, std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D_);
    virtual ~FillCellMap();
    FillCellMap(FillCellMap<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FillCellMap<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
	std::map<plint, Cell3D<T,Descriptor> > & cellIdToCell3D;
	TriangularSurfaceMesh<T>& mesh;
};


#include "cellFieldFunctionals3D.hh"
#endif  // CELL_FIELD_FUNCTIONALS_3D_H

