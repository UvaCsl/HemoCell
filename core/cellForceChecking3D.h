#ifndef CELL_FORCE_CHECKING_3D_H
#define CELL_FORCE_CHECKING_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <cmath>
#include <map>
#include "computeCellForces3D.h"

using namespace plb;
using namespace std;

typedef double T;

void testInPlane();
void testBending() ;

void calculateInPlaneForce(Array<T,3> x[3], Array<T,3> fx[3], T eqLengthRatio, T eqLength, T k_inPlane);

void writeInPlaneDataVTK(
        Array<T,3> const& x1, Array<T,3> const& x2,
        Array<T,3> const& x3,
        Array<T,3> const& fx1, Array<T,3> const& fx2,
        Array<T,3> const& fx3,
        std::string fName);

void writeBendingDataVTK(
        Array<T,3> const& x1, Array<T,3> const& x2,
        Array<T,3> const& x3, Array<T,3> const& x4,
        Array<T,3> const& fx1, Array<T,3> const& fx2,
        Array<T,3> const& fx3, Array<T,3> const& fx4,
        std::string fName);



#include "cellForceChecking3D.hh"
#endif  // CELL_FORCE_CHECKING_3D_H
