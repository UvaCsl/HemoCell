/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MESH_GENERATING_FUNCTIONS_H
#define MESH_GENERATING_FUNCTIONS_H

#include "constant_defaults.h"
#include "config.h"

#include "offLattice/triangleSet.hh"
#include "offLattice/triangleToDef.hh"
#include "offLattice/triangleSetGenerator.hh"
#include "offLattice/triangularSurfaceMesh.hh"
#include "offLattice/triangleBoundary3D.hh"

namespace plb {

template<typename T>
TriangleSet<T> constructCell(plb::Array<T,3> const& center, T radius, std::string cellFilename, plb::Array<T,3> const& eulerAngles);


template<typename T>
TriangleSet<T> constructSphereIcosahedron(plb::Array<T,3> const& center, T radius, plint minNumOfTriangles);


template<typename T>
plb::Array<T,3> mapMeshAsRBC(const plb::Array<T,3> point, const plb::Array<T,3> center, T R) ;


template<typename T>
plb::Array<T,3> spherePointToRBCPoint(const plb::Array<T,3> point, T R=1.0);


template<typename T>
TriangleSet<T> constructRBC(plb::Array<T,3> const& center, T radius, plint minNumOfTriangles, plb::Array<T,3> const& eulerAngles);


// initialSphereShape: [0] Octahedron (PLB Sphere ) [1] Icosahedron
template<typename T>
TriangleSet<T> constructRBCFromSphere(plb::Array<T,3> const& center, T radius, plint minNumOfTriangles,
        plb::Array<T,3> const& eulerAngles, pluint initialSphereShape=0);

template<typename T>
TriangleSet<T> constructEllipsoidFromSphere(plb::Array<T,3> const& center, T radius, T aspectRatio, plint minNumOfTriangles,
        plb::Array<T,3> const& eulerAngles, pluint initialSphereShape);


template<typename T>
TriangleBoundary3D<T> constructMeshElement(plint shape, T radius, plint cellNumTriangles, T dx, std::string cellPath, plb::Array<T,3> const& eulerAngles, T aspectRatio=0.3) {
    plb::Array<T,3> center({0.0, 0.0, 0.0});
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

plb::TriangularSurfaceMesh<T> * constructStringMeshFromConfig(hemo::Config & materialCfg);

#include "meshGeneratingFunctions.hh"

#endif  // MESH_GENERATING_FUNCTIONS_H
