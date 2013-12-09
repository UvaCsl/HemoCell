#ifndef CELL_IN_SHEAR_FLOW_3D_H
#define CELL_IN_SHEAR_FLOW_3D_H


#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedCellsReductions.h"
#include "immersedCellParticleVtk3D.hh"
#include "cellStretchingForces3D.h"
#include <map>
#include <algorithm>
#include "external/diagonalize.cpp"

namespace plb {


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
class SingleCellInShearFlow
{
public:
    SingleCellInShearFlow(TriangleBoundary3D<T> const& Cells,
            MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
            std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume,
            plint numParticlesPerSide_, plint flowType_, T dx_, T dt_, T dNewton_,
            std::map<plint, Particle3D<T,DESCRIPTOR>*> * tagToParticle3D_,
            bool store_=false);
    ~SingleCellInShearFlow() {};
    void updateQuantities(plint iteration, std::vector<plint> cellIds,
            std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume);
    void writeHeader(bool writeOutput=true);
    void write(bool writeOutput=true);
    std::vector<T> const& getIterations() { return iterations; } ;
    std::vector< Array<T,3> > const& getTumblingAngles() { return tumblingAngles; } ;
    std::vector< Array<T,3> > const& getTankTreadingAngles() { return tankTreadingAngles; } ;
    std::vector< Array<T,3> > const& getDiameters() { return diameters; } ;
    std::vector< vector<T> > const& getInertiaTensor() { return inertiaTensor; } ;
    std::vector<T> const& getSymmetryDeviation() { return symmetryDeviation; } ;
    std::vector<T> const& getDeformationIndex() { return deformationIndex; } ;
    std::vector<T> const& getTaylorDeformationIndex();
    void set_dx(T dx_) { dx = dx_; } ;
    void set_dt(T dt_) { dt = dt_; } ;
    void set_dm(T dm_) { dm = dm_; } ;
    void set_dxdtdm(T dx_, T dt_, T dm_) { dx = dx_; dt = dt_; dm = dm_;} ;
private:
    TriangleBoundary3D<T> const& Cells;
    MultiParticleField3D<ParticleFieldT<T,Descriptor> > * particles;
    std::vector<T> iterations;
    std::vector< Array<T,3> > tumblingAngles;
    std::vector< Array<T,3> > diameters;
    std::vector< vector<T> > inertiaTensor;
    std::vector<T> symmetryDeviation;
    std::vector<T> deformationIndex;
    T maxDiameter;
    plb_ofstream shearResultFile;
    plint numParticlesPerSide;
    plint flowType;
    T dx, dt, dNewton, dm;
    std::map<plint, Particle3D<T,DESCRIPTOR>*> * tagToParticle3D;
    bool store;

    std::vector< Array<T,3> > tankTreadingAngles; // This quantity is temporary (is not stored)
    std::vector<T> tankTreadingLengths; // This quantity is temporary (is not stored)

    std::vector<plint> outerLeftTags, outerRightTags;
    std::vector<plint> outerUpTags, outerDownTags;
    std::vector<plint> outerFrontTags, outerBackTags;
    std::vector<std::vector<plint>*> lateralCellParticleTags;

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
