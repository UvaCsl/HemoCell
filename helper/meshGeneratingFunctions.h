#ifndef MESH_GENERATING_FUNCTIONS_H
#define MESH_GENERATING_FUNCTIONS_H

#include "hemocell_internal.h"
#include "config.h"

namespace plb {

template<typename T>
TriangleSet<T> constructCell(hemo::Array<T,3> const& center, T radius, std::string cellFilename, hemo::Array<T,3> const& eulerAngles);


template<typename T>
TriangleSet<T> constructSphereIcosahedron(hemo::Array<T,3> const& center, T radius, plint minNumOfTriangles);


template<typename T>
hemo::Array<T,3> mapMeshAsRBC(const hemo::Array<T,3> point, const hemo::Array<T,3> center, T R) ;


template<typename T>
hemo::Array<T,3> spherePointToRBCPoint(const hemo::Array<T,3> point, T R=1.0);


template<typename T>
TriangleSet<T> constructRBC(hemo::Array<T,3> const& center, T radius, plint minNumOfTriangles, hemo::Array<T,3> const& eulerAngles);


// initialSphereShape: [0] Octahedron (PLB Sphere ) [1] Icosahedron
template<typename T>
TriangleSet<T> constructRBCFromSphere(hemo::Array<T,3> const& center, T radius, plint minNumOfTriangles,
        hemo::Array<T,3> const& eulerAngles, pluint initialSphereShape=0);

template<typename T>
TriangleSet<T> constructEllipsoidFromSphere(hemo::Array<T,3> const& center, T radius, T aspectRatio, plint minNumOfTriangles,
        hemo::Array<T,3> const& eulerAngles, pluint initialSphereShape);


template<typename T>
TriangleBoundary3D<T> constructMeshElement(plint shape, T radius, plint cellNumTriangles, T dx, std::string cellPath, hemo::Array<T,3> const& eulerAngles, T aspectRatio=0.3) {
    hemo::Array<T,3> center(0.0, 0.0, 0.0);
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
    } else if (shape == 6) {
        wholeTriangleSet.append(constructEllipsoidFromSphere<T>(center, radius, aspectRatio, cellNumTriangles, eulerAngles, 0));
    }
    /* Unknown code, was there from the beginning, let it be */
    DEFscaledMesh<T> defMesh (wholeTriangleSet);
    defMesh.setDx(dx);
    TriangleBoundary3D<T> Cells(defMesh);
    Cells.getMesh().inflate();
    return Cells;
}


} // namespace plb

TriangularSurfaceMesh<double> * constructStringMeshFromConfig(Config & materialCfg);

#include "meshGeneratingFunctions.hh"

#endif  // MESH_GENERATING_FUNCTIONS_H
