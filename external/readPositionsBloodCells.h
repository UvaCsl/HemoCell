#ifndef READ_POSISIONS_OF_BLOOD_CELLS_H
#define READ_POSISIONS_OF_BLOOD_CELLS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "initializationCellField3D.h"
#include "initializationCellField3D.hh"
#include <vector>

using namespace std;
using namespace plb;

template<typename T, template<typename U> class Descriptor>
void readPositionsBloodCellField3D(std::vector<CellField3D<T, Descriptor>* > & cellFields, T dx, const char* positionsFileName);

template<typename T>
void getReadPositionsBloodCellsVector(Box3D realDomain,
                                           std::vector<TriangularSurfaceMesh<T>* > & meshes,
                                           std::vector<plint> & Np,
                                           std::vector<std::vector<Array<T,3> > > & positions,
                                           std::vector<std::vector<plint> > & cellIds,
                                           std::vector<std::vector<Array<T,3> > > & randomAngles,
                                           const char* positionsFileName);

template<typename T, template<typename U> class Descriptor>
class ReadPositionsBloodCellField3D : public BoxProcessingFunctional3D
{
public:
    ReadPositionsBloodCellField3D (std::vector<CellField3D<T, Descriptor>* > & cellFields_, T dx_, const char* positionsFileName_):
            cellFields(cellFields_), dx(dx_), positionsFileName(positionsFileName_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ReadPositionsBloodCellField3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<CellField3D<T, Descriptor>* > & cellFields;
    const char* positionsFileName;
    T dx;
};


#include "readPositionsBloodCells.hh"
#endif