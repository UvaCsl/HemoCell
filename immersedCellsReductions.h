#ifndef IMMERSED_CELLS_REDUCTIONS_H
#define IMMERSED_CELLS_REDUCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "shellModel3D.h"
#include <map>

/* ******** CellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class CellReduceFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_, bool findMax_=false);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CellReduceFunctional3D<T,Descriptor>* clone() const;
    virtual void  getCellQuantityArray(std::vector<T>& cellQuantity, std::vector<plint> cellIds) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
protected:
    plint numberOfCells;
    std::vector<plint> quantityIds;
    TriangleBoundary3D<T> const& triangleBoundary;
    bool findMax;
};


/* ******** countCellVolume ** countCellSurface ********************************* */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVolume (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVolumes);

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellSurface (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellSurfaces) ;//Perhaps add TAGS

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanTriangleArea (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanTriangleArea) ;

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanEdgeDistance (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanEdgeDistance) ;

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanTileSpan (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanTileSpan) ;

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellCenters(TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,3> > & cellsCenter, std::vector<T> & cellNumVertices); //Perhaps add TAGS

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVelocity(TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,3> > & cellsVelocity, std::vector<T> & cellNumVertices); //Perhaps add TAGS

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMaxEdgeDistance (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMaxEdgeDistance) ;

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellMeanAngle (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellMeanAngle) ;

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVertices (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVertices);


/* ******** CenterCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class CenterCellReduceFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CenterCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_);
    /// Argument: Particle-field.
    virtual CenterCellReduceFunctional3D<T,Descriptor>* clone() const { return new CenterCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual void getCellQuantityArray(std::vector< Array<T,3> > & cellsCenter, std::vector<plint> cellIds) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
protected:
    plint numberOfCells;
    std::vector<plint> quantityIds;
    TriangleBoundary3D<T> const& triangleBoundary;
    bool findMax;
};


/* ******** VelocityCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class VelocityCellReduceFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    VelocityCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_);
    /// Argument: Particle-field.
    virtual VelocityCellReduceFunctional3D<T,Descriptor>* clone() const { return new VelocityCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual void getCellQuantityArray(std::vector< Array<T,3> > & cellsVelocity, std::vector<plint> cellIds) const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
protected:
    plint numberOfCells;
    std::vector<plint> quantityIds;
    TriangleBoundary3D<T> const& triangleBoundary;
    bool findMax;
};

/* ******** VolumeCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class VolumeCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    VolumeCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual VolumeCellReduceFunctional3D<T,Descriptor>* clone() const { return new VolumeCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** SurfaceCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class SurfaceCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    SurfaceCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual SurfaceCellReduceFunctional3D<T,Descriptor>* clone() const { return new SurfaceCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** EdgeDistanceCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class EdgeDistanceCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    EdgeDistanceCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual EdgeDistanceCellReduceFunctional3D<T,Descriptor>* clone() const { return new EdgeDistanceCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** TileSpanCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class TileSpanCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
	TileSpanCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual TileSpanCellReduceFunctional3D<T,Descriptor>* clone() const { return new TileSpanCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};

/* ******** EdgeDistanceCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class MaxEdgeDistanceCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    MaxEdgeDistanceCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_, true) { }; // true refers to the max Reduction
    virtual MaxEdgeDistanceCellReduceFunctional3D<T,Descriptor>* clone() const { return new MaxEdgeDistanceCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** AngleCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class AngleCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    AngleCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual AngleCellReduceFunctional3D<T,Descriptor>* clone() const { return new AngleCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** NumTrianglesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class NumTrianglesCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    NumTrianglesCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual NumTrianglesCellReduceFunctional3D<T,Descriptor>* clone() const { return new NumTrianglesCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** NumEdgesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class NumEdgesCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    NumEdgesCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual NumEdgesCellReduceFunctional3D<T,Descriptor>* clone() const { return new NumEdgesCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


/* ******** NumVerticesCellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class NumVerticesCellReduceFunctional3D : public CellReduceFunctional3D<T, Descriptor> {
public:
    NumVerticesCellReduceFunctional3D(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
     : CellReduceFunctional3D<T, Descriptor>(triangleBoundary_, cellIds_) { };
    virtual NumVerticesCellReduceFunctional3D<T,Descriptor>* clone() const { return new NumVerticesCellReduceFunctional3D<T,Descriptor>(*this); };
    virtual void calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
            std::vector<Particle3D<T,Descriptor>*> & particles, std::vector<plint> & quantityIds_);
private:
};


#include "immersedCellsReductions.hh"

#endif  // IMMERSED_CELLS_REDUCTIONS_H
