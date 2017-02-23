#if 0
#ifndef ORDERED_POSISIONS_OF_MULTIPLE_CELLS_H
#define ORDERED_POSISIONS_OF_MULTIPLE_CELLS_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "initializationCellField3D.h"
#include "initializationCellField3D.hh"
#include "binpacking/3dbpp.cpp"
#include <vector>

using namespace std;
using namespace plb;

template<typename T, template<typename U> class Descriptor>
void orderedPositionMultipleCellField3D(std::vector<HemoCellField* > & cellFields);


template<typename T>
void getOrderedPositionsMultipleCellsVector(Box3D realDomain, 
    std::vector<TriangularSurfaceMesh<T>* > & meshes,
    std::vector<plint> & Np,
    std::vector<std::vector<Array<T,3> > > & positions, 
    std::vector<std::vector<plint> > & cellIds) ;


template<typename T, template<typename U> class Descriptor>
class OrderedPositionMultipleCellField3D : public BoxProcessingFunctional3D
{
public:
    OrderedPositionMultipleCellField3D (std::vector<HemoCellField* > & cellFields_):
                cellFields(cellFields_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual OrderedPositionMultipleHemoCellField* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<HemoCellField* > & cellFields;
};




#include "orderedPositionsMultipleCells.hh"
#endif // ORDERED_POSISIONS_OF_MULTIPLE_CELLS_H
#endif

