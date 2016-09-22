#ifndef RANDOM_POSISIONS_OF_MULTIPLE_CELLS_H
#define RANDOM_POSISIONS_OF_MULTIPLE_CELLS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "initializationCellField3D.h"
#include "initializationCellField3D.hh"
#include "dynpacking/packEllipsoids.cpp"
#include <vector>

using namespace std;
using namespace plb;

template<typename T, template<typename U> class Descriptor>
void randomPositionMultipleCellField3D(std::vector<CellField3D<T, Descriptor>* > & cellFields, T packingDensity, T dx, plint maxPackIter = 25000);

template<typename T>
void getRandomPositionsMultipleCellsVector(Box3D realDomain,
                                           std::vector<TriangularSurfaceMesh<T>* > & meshes,
                                           std::vector<plint> & Np,
                                           std::vector<std::vector<Array<T,3> > > & positions,
                                           std::vector<std::vector<plint> > & cellIds,
                                           std::vector<std::vector<Array<T,3> > > & randomAngles,
                                           T packingDensity);

template<typename T, template<typename U> class Descriptor>
class RandomPositionMultipleCellField3D : public BoxProcessingFunctional3D
{
public:
    RandomPositionMultipleCellField3D (std::vector<CellField3D<T, Descriptor>* > & cellFields_, T packingDensity_, T dx_, plint maxPackIter_):
            cellFields(cellFields_), packingDensity(packingDensity_), dx(dx_), maxPackIter(maxPackIter_) { }
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual RandomPositionMultipleCellField3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<CellField3D<T, Descriptor>* > & cellFields;
    T packingDensity;
    plint maxPackIter;
    T dx;
};


#include "randomPositionsMultipleCells.hh"
#endif