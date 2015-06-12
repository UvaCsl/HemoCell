/* This file is part of the bloodyCells library.
 * Copyright (C) 2009, 2010 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IMMERSED_CELLS_3D_HH
#define IMMERSED_CELLS_3D_HH

#include "immersedCells3D.h"

using namespace plb;
using namespace std;


#ifndef PI__
#define PI__
const T pi = 4.*atan(1.);
#endif  // PI__



plint xDirection  = 0;
plint yDirection  = 1;
plint zDirection  = 2;
plint margin          = 3;  // Extra margin of allocated cells around the obstacle,
                             //   for moving walls.
//const plint extendedEnvelopeWidth = 2;  // Because Guo needs 2-cell neighbor access.
//const plint particleEnvelopeWidth = 6;




template<typename T, template<typename U> class Descriptor>
void deleteCell(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                   const Box3D &outlet, std::vector<plint> &numParts,
                   std::vector<plint> &cellIds,
                   std::vector<Array<T,3> > &centers, std::vector<T> &radii )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        // count all particles of a certain tag in a buffer zone
        plint numPartsPerTag = countParticles(particleField,outlet,cellIds[iA]);
        // if all the particle of a certain tag are in a buffer zone
        if (numPartsPerTag > 0) {
            // then delete these particles in all the buffer zones
            plint before = countParticles(particleField,outlet,cellIds[iA]);

            applyProcessingFunctional (
                new AbsorbTaggedParticlesFunctional3D<T,Descriptor>(cellIds[iA]),
                particleField.getBoundingBox(), particleArg );

            plint after = countParticles(particleField,outlet,cellIds[iA]);

            pcout << "erased particles = " << before << ", " << after << std::endl;

//             delete Cells[iA];                     // delete the iA-th pointer mesh
//             numParts.erase(numParts.begin()+iA);        // erase the iA-th number of particles
//             cellIds.erase(cellIds.begin()+iA);                // erase the iA-th tag
//             Cells.erase(Cells.begin()+iA);  // erase the iA-th mesh
//             centers.erase(centers.begin()+iA);          // delete the iA-th center
//             radii.erase(radii.begin()+iA);              // delete the iA-th radius
//             --iA;

        }
    }
//    if (erased) {
//        pcout << "Particles absorbed : Number of particles per tag." << std::endl;
//        for (pluint iA = 0; iA < centers.size(); ++iA) {
//            pcout << cellIds[iA] << " , " <<  countParticles(particleField, particleField.getBoundingBox(), cellIds[iA]) << std::endl;
//        }
//    }
}

template<typename T, template<typename U> class Descriptor>
bool generateCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                      const Box3D &inlet, std::vector<plint> &cellIds, TriangleBoundary3D<T> &Cells,
                      plint numPartsPerCell, plint numOfCellsPerInlet, plint &slice )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        applyProcessingFunctional (
            new CreateTaggedImmersedCellParticle3D<T,Descriptor>(Cells, cellIds[iA], numPartsPerCell),
            particleField.getBoundingBox(), particleArg );
    }
    return true;
}

template<typename T>
TriangleBoundary3D<T> createCompleteMesh(
    const std::vector<Array<T,3> > &centers, const std::vector<T> &radii, Array<T,3> const& eulerAngles,
    std::vector<plint> &cellIds, IncomprFlowParam<T> const& parameters,
    plint shape, std::string cellPath, plint &cellNumTriangles, plint &numPartsPerCell)
{
//	shape 0:Sphere, 1:RBC
    PLB_ASSERT(centers.size() == radii.size());

    // creation of the whole mesh including the one that is not used by particles
    std::vector<TriangleSet<T> > allTriangles;
    TriangleSet<T> wholeTriangleSet; 
    for (pluint iA = 0; iA < centers.size(); ++iA) {
        Array<T,3> center(centers[iA]);
        T radius = radii[iA];
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
        cellIds.push_back(iA);
    }
//    wholeTriangleSet.merge(allTriangles);
//    wholeTriangleSet.writeAsciiSTL("/tmp/test.stl");
//    Dot3D location(centers[0][0]-radii[0],centers[0][1]-radii[0],centers[0][2]-radii[0]);
//    DEFscaledMesh<T> defMesh (
//            wholeTriangleSet, parameters.getResolution(), xDirection, margin, location );
    DEFscaledMesh<T> defMesh (wholeTriangleSet);
    defMesh.setDx(parameters.getDeltaX());
    TriangleBoundary3D<T> Cells(defMesh);
    Cells.getMesh().inflate();
    pcout << "Original sphere at location [" << centers[0][0] << ","
             << centers[0][1] << "," << centers[0][2] << "] " << std::endl;
//    pcout << "defMesh.getDx = " <<  defMesh.getDx() << std::endl;

    numPartsPerCell = Cells.getMesh().getNumVertices() / centers.size();
    plint modulo = Cells.getMesh().getNumVertices() % centers.size();
//     pcout << "num parts per triangles = " << numPartsPerCell << ", modulo = " << modulo << std::endl;
    pcout << "num Vertices = " << Cells.getMesh().getNumVertices() ;
    pcout << " (" << Cells.getMesh().getNumVertices()/centers.size() << ")" ;
    pcout << ", num Triangles = " << Cells.getMesh().getNumTriangles();
    pcout << " (" << Cells.getMesh().getNumTriangles()/centers.size() << ")" ;
    pcout << ", centers = " << centers.size() ;
    pcout << ", modulo = " << modulo << std::endl;

    PLB_ASSERT(modulo == 0);
    Cells.cloneVertexSet(0);
    
    return Cells;
}


#endif  // IMMERSED_CELLS_3D_HH


