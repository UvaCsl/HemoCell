#ifndef COLLECTIVE_CELL_REDUCTIONS_H
#define COLLECTIVE_CELL_REDUCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"


/* IDs for the reduction types,
 * they are created based on the following rules:
 * ==============================
 * Last digit, Dimension:
 *      1D : 1
 *      2D : 2
 *      3D : 3
 *      XD : 4,5,6,7,8,9
 * Second to last digit, type of reduction:
 *      Sum : 0
 *      Mean: 1
 *      Min : 2
 *      Max : 3
 *      STD : 4 // Still to be implemented
 * ID for the quantity of interest, starting from 1. Can be grouped together like:
 *      Angle          : 2
 *      Area           : 3
 *      Edge Distance  : 4
 *      Edge Tile Span : 5
 *      Position       : 6
 *      Velocity       : 7
 */
#define CCR_VOLUME                 101 // 1d
#define CCR_ANGLE_MEAN             211 // 1d
#define CCR_ANGLE_MIN              221 // 1d
#define CCR_ANGLE_MAX              231 // 1d
#define CCR_SURFACE                301 // 1d
#define CCR_TRIANGLE_AREA_MEAN     311 // 1d, Better use 301 and divide by the number of triangles
#define CCR_TRIANGLE_AREA_MIN      321 // 1d
#define CCR_TRIANGLE_AREA_MAX      331 // 1d
#define CCR_EDGE_DISTANCE_MEAN     411 // 1d
#define CCR_EDGE_DISTANCE_MIN      421 // 1d
#define CCR_EDGE_DISTANCE_MAX      431 // 1d
#define CCR_TILE_SPAN_MEAN         511 // 1d
#define CCR_TILE_SPAN_MIN          521 // 1d
#define CCR_TILE_SPAN_MAX          531 // 1d
#define CCR_POSITION_MEAN          613 // 3d
#define CCR_POSITION_MIN           623 // 3d
#define CCR_POSITION_MAX           633 // 3d
#define CCR_VELOCITY_MEAN          713 // 3d
#define CCR_VELOCITY_MIN           723 // 3d
#define CCR_VELOCITY_MAX           733 // 3d
#define CCR_INERTIA                809 // 9d, Not working
// #define CCR_MAX               41 // 9d




/* ******** CellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class CellReductorWrapper
{
    public:
        CellReductorWrapper(std::vector<MultiBlock3D*> & particleArg_,
                            TriangleBoundary3D<T> const& triangleBoundary_,
                            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
                            std::vector<plint> quantitiesToReduce_);
        CellReductorWrapper(std::vector<MultiBlock3D*> & particleArg_,
                            TriangleBoundary3D<T> const& triangleBoundary_,
                            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
                            plint quantitiesToReduce_=-1);
        virtual ~CellReductorWrapper() { };

        void reduce(plint quantitiesToReduce_=-1);
        void reduce(std::vector<plint> & quantitiesToReduce);
        void setCarryOnQuantities1D(std::map<plint, std::map<plint, T > * > const& carryOnQuantities1D);
        void setCarryOnQuantities3D(std::map<plint, std::map<plint, Array<T,3> > *  > const& carryOnQuantities3D);
        // carryOnQuantities3D[CCR_INERTIA][cellId] = Positions
        void setCarryOnQuantitiesXD(std::map<plint, std::map<plint, std::vector<T> > * > const& carryOnQuantitiesXD);
    private:
        CollectiveCellReductions<T,Descriptor> CCR;
        TriangleBoundary3D<T> const& triangleBoundary;
        plint maxNumberOfCells;
        plint numVerticesPerCell;
        plint numTrianglesPerCell;
        std::vector<plint> quantitiesToReduce;

        std::vector<plint> cellIDs;
        std::map<plint, std::map<plint, T > * > quantities1D; // quantities1D[CCR_EDGE_DISTANCE_MEAN][cellId]
        std::map<plint, std::map<plint, Array<T,3> > * > quantities3D; // quantities3D[CCR_VELOCITY_MEAN][cellId]
        std::map<plint, std::map<plint, std::vector<T> > * > quantitiesXD; // quantitiesXD[CCR_INERTIA][cellId]
};

/*
 *    CollectiveCellReductionBox3D is the Functional Reductive Box that does the actual job.
 *
 */


/* ******** CellReduceFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
class CollectiveCellReductionBox3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
            std::vector<plint> subscribedQuantities_);
    CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_,
            plint subscribedQuantities_);
    CollectiveCellReductionBox3D(TriangleBoundary3D<T> const& triangleBoundary_,
            plint maxNumberOfCells_, plint numVerticesPerCell_, plint numTrianglesPerCell_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CellReduceFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

    void reduceQuantity();
    void reduceQuantity(plint quantityId);
    void reduceQuantity(std::vector<plint> quantitiesToReduce);
    void getCellIDs(std::vector<plint> & cellIDs);
    /* This method is overridden multiple times */
    void getCellQuantity(std::map<plint, T > & quantity1D,
                         std::map<plint, std::map<plint, Array<T,3> >  > & quantity3D,
                         std::map<plint, std::map<plint, std::vector<T> >  > & quantitiesXD);
private:
    std::map<plint, std::map<plint, plint >  > whatQuantity1D; // whatQuantity1D[CCR_EDGE_DISTANCE_MEAN][cellId]
    std::map<plint, std::map<plint, Array<plint,3> >  > whatQuantity3D; // whatQuantity3D[CCR_VELOCITY_MEAN][cellId]
    std::map<plint, std::map<plint, std::vector<plint> >  > whatQuantityXD; // whatQuantityXD[CCR_INERTIA][cellId]

    TriangleBoundary3D<T> const& triangleBoundary;
    std::vector<plint> quantitiesToReduce;
    std::set<plint> subscribedQuantities;
    std::vector<plint> & cellIDsInsideTheDomain;
    plint maxNumberOfCells, numVerticesPerCell, numTrianglesPerCell;
public:
    void getCellQuantity(std::map<plint, std::map<plint, Array<T,3> > * > & quantity3D);
    void getCellQuantity(std::map<plint, std::map<plint, std::vector<T> > * > & quantitiesXD);
    void getCellQuantity(std::map<plint, T > & quantity1D,
                         std::map<plint, std::map<plint, Array<T,3> > * > & quantity3D);
private: /* ReductiveBoxFunctions */
    plint subscribeReduction1D(plint reductionType);
    Array<plint,3> subscribeReduction3D(plint reductionType);
    std::vector<plint> subscribeReductionXD(plint reductionType, plint dimensions);

    void gatherReduction1D(plint reductionType, plint whatQ, T value);
    void gatherReduction3D(plint reductionType, Array<plint,3> whatQ, Array<T,3> value);
    void gatherReductionXD(plint reductionType, std::vector<plint> whatQ, std::vector<T> value);

    T getReduction1D(plint reductionType, plint whatQ);
    Array<T,3> getReduction3D(plint reductionType, Array<plint,3> whatQ);
    std::vector<T> getReductionXD(plint reductionType, std::vector<plint> whatQ);



};





#include "CollectiveCellReductions.hh"

#endif  // COLLECTIVE_CELL_REDUCTIONS_H
