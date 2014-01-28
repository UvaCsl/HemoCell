#ifndef CELL_REDUCTION_TYPES_H
#define CELL_REDUCTION_TYPES_H

#include "palabos3D.h"
#include "palabos3D.hh"


/* IDs for the reduction types,
 * they are created based on the following rules:
 * ==============================
 * Last digit, Dimension:
 *      1D : 1
 *      2D : 2
 *      3D : 3
 *      ND : 4,5,6,7,8,9 // Not above 10
 * Second to last digit, type of reduction:
 *      Sum : 0
 *      Mean: 1
 *      Min : 2
 *      Max : 3
 *      STD : 4 // Still to be implemented
 * ID for the quantity of interest, starting from 1. Can be grouped together like:
 *      Volume         : 1
 *      Angle          : 2
 *      Area           : 3
 *      Edge Distance  : 4
 *      Edge Tile Span : 5
 *      Position       : 6 // Periodic boundary position
 *      Velocity       : 7
 *      Inertia        : 8
 *      Energy         : 9
 *      Position       : 0
 */
#define CCR_NO_PBC_POSITION_MEAN    13 // 3d // BEWARE OF 0 IN FRONT! GOES TO OCT
#define CCR_NO_PBC_POSITION_MIN     23 // 3d // BEWARE OF 0 IN FRONT! GOES TO OCT
#define CCR_NO_PBC_POSITION_MAX     33 // 3d // BEWARE OF 0 IN FRONT! GOES TO OCT
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
#define CCR_ENERGY                 901 // 1d
// #define CCR_MAX               41 // 9d


#endif
