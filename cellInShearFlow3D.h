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


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
class SingleCellInShearFlow
{
public:
    SingleCellInShearFlow(bool store_=false);
    SingleCellInShearFlow(plint iteration, TriangleBoundary3D<T> Cells,
            MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
            std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume, bool store_=false);
    SingleCellInShearFlow(TriangleBoundary3D<T> Cells,
            MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
            std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume, bool store_=false);
    ~SingleCellInShearFlow() {};
    void addData(plint iteration, TriangleBoundary3D<T> Cells,
            MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
            std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume);
    void write(plb_ofstream & shearResultFile);
    void writeHeader(plb_ofstream & shearResultFile);
    void write() { this->write(*shearResultFile); } ;
    std::vector<T> const& getTaylorDeformationIndex();
    std::vector<T> const& getDeformationIndex() { return deformationIndex; } ;
    std::vector<T> const& getIterations() { return iterations; } ;
    std::vector< Array<T,3> > const& getAngles() { return angles; } ;
    std::vector< vector<T> > const& getInertiaTensor() { return inertiaTensor; } ;
    std::vector<T> const& getSymmetryDeviation() { return symmetryDeviation; } ;
    void setShearResultFile(plb_ofstream & shearResultFile_) { shearResultFile = &shearResultFile_; } ;
    void set_dx(T dx_) { dx = dx_; } ;
    void set_dt(T dt_) { dt = dt_; } ;
    void set_dm(T dm_) { dm = dm_; } ;
    void set_dxdtdm(T dx_, T dt_, T dm_) { dx = dx_; dt = dt_; dm = dm_;} ;
private:
    std::vector<T> deformationIndex;
    std::vector<T> iterations;
    std::vector< Array<T,3> > angles;
    std::vector< Array<T,3> > diameters;
    std::vector< vector<T> > inertiaTensor;
    std::vector<T> symmetryDeviation;
    bool store;
    T maxDiameter;
    plb_ofstream * shearResultFile;
    T dx, dt, dm;
};


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


/* ******** computeCellInertia *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeEllipsoidFit (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
                std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume,
                std::vector< std::vector<T> > & cellsEllipsoidFitAngles,
                std::vector< std::vector<T> > & cellsEllipsoidFitSemiAxes,
                std::vector< std::vector<T> > & cellInertia,
                std::vector<T> & difference);

}

#include "cellInShearFlow3D.hh"

#endif // CELL_IN_SHEAR_FLOW_3D_H
