#ifndef READ_POSISIONS_OF_BLOOD_CELLS_H
#define READ_POSISIONS_OF_BLOOD_CELLS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include "initializationCellField3D.h"
#include "initializationCellField3D.hh"

using namespace std;
using namespace plb;

void readPositionsBloodCellField3D(CellFields3D & cellFields, double dx, const char* positionsFileName);

template<typename T>
void getReadPositionsBloodCellsVector(Box3D realDomain,
                                           std::vector<TriangularSurfaceMesh<T>* > & meshes,
                                           std::vector<plint> & Np,
                                           std::vector<std::vector<Array<T,3> > > & positions,
                                           std::vector<std::vector<plint> > & cellIds,
                                           std::vector<std::vector<Array<T,3> > > & randomAngles,
                                           const char* positionsFileName);

class ReadPositionsBloodCellField3D : public BoxProcessingFunctional3D
{
public:
    ReadPositionsBloodCellField3D (CellFields3D & cellFields_, double dx_, const char* positionsFileName_)
            : cellFields(cellFields_) {dx = dx_; positionsFileName=positionsFileName_;}
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ReadPositionsBloodCellField3D* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    const char* positionsFileName;
    double dx;
    CellFields3D & cellFields;
};


#include "readPositionsBloodCells.hh"
#endif
