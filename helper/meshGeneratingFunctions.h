#ifndef MESH_GENERATING_FUNCTIONS_H
#define MESH_GENERATING_FUNCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"

namespace plb {

template<typename T>
TriangleSet<T> constructCell(Array<T,3> const& center, T radius, std::string cellFilename, Array<T,3> const& eulerAngles);


template<typename T>
TriangleSet<T> constructSphereIcosahedron(Array<T,3> const& center, T radius, plint minNumOfTriangles);


template<typename T>
Array<T,3> mapMeshAsRBC(const Array<T,3> point, const Array<T,3> center, T R) ;


template<typename T>
Array<T,3> spherePointToRBCPoint(const Array<T,3> point, T R=1.0);


template<typename T>
TriangleSet<T> constructRBC(Array<T,3> const& center, T radius, plint minNumOfTriangles, Array<T,3> const& eulerAngles);


// initialSphereShape: [0] Octahedron (PLB Sphere ) [1] Icosahedron
template<typename T>
TriangleSet<T> constructRBCFromSphere(Array<T,3> const& center, T radius, plint minNumOfTriangles,
        Array<T,3> const& eulerAngles, pluint initialSphereShape=0);

template<typename T>
TriangleBoundary3D<T> constructMeshElement(plint shape, T radius, plint cellNumTriangles, T dx, std::string& cellPath, Array<T,3> const& eulerAngles) {
    Array<T,3> center(2*radius, 2*radius, 2*radius);
    std::vector<TriangleSet<T> > allTriangles;
    TriangleSet<T> wholeTriangleSet;
    if (shape == 0) {
        wholeTriangleSet.append(constructSphereIcosahedron<T>(center, radius, cellNumTriangles));
    } else if (shape == 1) {
        wholeTriangleSet.append(constructRBCFromSphere<T>(center, radius, cellNumTriangles, eulerAngles, 1));
    } else if (shape == 2) {
        wholeTriangleSet.append(constructCell<T>(center, radius, cellPath, eulerAngles));
    } else if (shape == 3) {
        wholeTriangleSet.append(constructRBC<T>(center, radius, cellNumTriangles, eulerAngles));
    } else if (shape == 4) {
        wholeTriangleSet.append(constructRBCFromSphere<T>(center, radius, cellNumTriangles, eulerAngles, 0));
    } else if (shape == 5) {
        wholeTriangleSet.append(constructSphere<T>(center, radius, cellNumTriangles));
    }
    /* Unknown code, was there from the beginning, let it be */
    DEFscaledMesh<T> defMesh (wholeTriangleSet);
    defMesh.setDx(dx);
    TriangleBoundary3D<T> Cells(defMesh);
    Cells.getMesh().inflate();
    return Cells;
}


} // namespace plb

#include "meshGeneratingFunctions.hh"

#endif  // MESH_GENERATING_FUNCTIONS_H
