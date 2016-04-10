#ifndef ORDERED_POSISIONS_OF_MULTIPLE_CELLS_WITH_SPACE_H
#define ORDERED_POSISIONS_OF_MULTIPLE_CELLS_WITH_SPACE_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "initializationCellField3D.h"
#include "initializationCellField3D.hh"
// #include "3dbpp.cpp"
#include <vector>
#include "orderedPositionsMultipleCells.h"

using namespace std;
using namespace plb;

template<typename T, template<typename U> class Descriptor>
void orderedPositionMultipleCellField3DWithSpace(std::vector<CellField3D<T, Descriptor>* > & cellFields, Box3D ToDelete);
//void orderedPositionMultipleCellField3DWithSpace(CellField3D<T, Descriptor> & cellFields, Box3D ToDelete);



template<typename T, template<typename U> class Descriptor>
class OrderedPositionMultipleCellField3DWithSpace : public BoxProcessingFunctional3D
{
public:
    OrderedPositionMultipleCellField3DWithSpace (std::vector<CellField3D<T, Descriptor>* > & cellFields_, Box3D ToDelete_):
   //  OrderedPositionMultipleCellField3DWithSpace (CellField3D<T, Descriptor> & cellField_, Box3D ToDelete_):
                cellFields(cellFields_), ToDelete(ToDelete_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual OrderedPositionMultipleCellField3DWithSpace<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<CellField3D<T, Descriptor>* > & cellFields;
   // CellField3D<T, Descriptor> & cellField;
    Box3D ToDelete;
};




#include "orderedPositionsMultipleCellsWithSpace.hh"
#endif // ORDERED_POSISIONS_OF_MULTIPLE_CELLS_WITH_SPACE_H


