#ifndef COLLECTIVE_CELL_REDUCTIONS_H
#define COLLECTIVE_CELL_REDUCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"

#define CCR_VOLUME               0 // 1d
#define CCR_SURFACE              1 // 1d
#define CCR_ANGLE_MEAN           2 // 1d
#define CCR_ANGLE_MIN            3 // 1d
#define CCR_ANGLE_MAX            4 // 1d
#define CCR_TRIANGLE_AREA_MEAN   5 // 1d
#define CCR_TRIANGLE_AREA_MIN    6 // 1d
#define CCR_TRIANGLE_AREA_MAX    7 // 1d
#define CCR_EDGE_DISTANCE_MEAN   8 // 1d
#define CCR_EDGE_DISTANCE_MIN    9 // 1d
#define CCR_EDGE_DISTANCE_MAX    10 // 1d
#define CCR_TILE_SPAN_MEAN       11 // 1d
#define CCR_TILE_SPAN_MIN        12 // 1d
#define CCR_TILE_SPAN_MAX        13 // 1d
#define CCR_POSITION_MEAN        14 // 3d
#define CCR_POSITION_MIN         17 // 3d
#define CCR_POSITION_MAX         20 // 3d
#define CCR_VELOCITY_MEAN        23 // 3d
#define CCR_VELOCITY_MIN         26 // 3d
#define CCR_VELOCITY_MAX         29 // 3d
#define CCR_INERTIA              32 // 9d
// #define CCR_MAX               41 // 9d



/* ******** CellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class CellReductorWrapper
{
    public:
        CellReductorWrapper(std::vector<MultiBlock3D*> & particleArg_,
                            TriangleBoundary3D<T> const& triangleBoundary_,
                            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
                            std::vector<plint> & quantitiesToReduce_);
        CellReductorWrapper(std::vector<MultiBlock3D*> & particleArg_,
                            TriangleBoundary3D<T> const& triangleBoundary_,
                            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
                            plint quantitiesToReduce_=-1);

        void reduce(plint quantitiesToReduce_=-1);
        void reduce(std::vector<plint> & quantitiesToReduce);
    private:
        TriangleBoundary3D<T> const& triangleBoundary;
        plint maxNumberOfCells;
        plint numVerticesPerCell;
        plint numTrianglesPerCell;
        std::vector<plint> * quantitiesToReduce;

        std::vector<plint> cellIDs;
        std::map<plint, std::map<plint, T >  > quantities1D; // quantities1D[cellId][CCR_EDGE_DISTANCE_MEAN]
        std::map<plint, std::map<plint, Array<T,3> >  > quantities3D; // quantities3D[cellId][CCR_VELOCITY_MEAN]
        std::map<plint, std::map<plint, std::vector<T> >  > quantitiesND; // quantitiesND[cellId][CCR_INERTIA]
};



/* ******** CellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class CollectiveCellReductions : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CollectiveCellReductions(TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_, plint numberOfCells_, bool findMax_=false);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CellReduceFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

    void calculateQuantity();
    void getCellIDs(std::vector<plint> & cellIDs);
private:
    std::vector<plint> quantityIds;
    TriangleBoundary3D<T> const& triangleBoundary;
    plint maxNumberOfCells;
    plint numVerticesPerCell;
    plint numTrianglesPerCell;
    std::vector<plint> * quantitiesToReduce;

};





#include "collectiveCellReductions.hh"

#endif  // COLLECTIVE_CELL_REDUCTIONS_H
