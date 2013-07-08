#ifndef MESH_GENERATING_FUNCTIONS_H
#define MESH_GENERATING_FUNCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"

namespace plb {

template<typename T>
TriangleSet<T> constructSphereIcosahedron(Array<T,3> const& center, T radius, plint minNumOfTriangles);

template<typename T>
Array<T,3> spherePointToRBCPoint(const Array<T,3> point, T R=1.0);


template<typename T>
TriangleSet<T> constructRBC(Array<T,3> const& center, T radius, plint minNumOfTriangles, std::vector<T> const& eulerAngles);


template<typename T>
TriangleSet<T> constructRBCFromSphere(Array<T,3> const& center, T radius, plint minNumOfTriangles,
        std::vector<T> const& eulerAngles);


template<typename T>
TriangleSet<T> constructCell(Array<T,3> const& center, T radius, std::string cellFilename, std::vector<T> const& eulerAngles);



} // namespace plb

#include "meshGeneratingFunctions.hh"

#endif  // MESH_GENERATING_FUNCTIONS_H
