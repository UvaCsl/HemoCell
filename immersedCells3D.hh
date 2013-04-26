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
#include "immersedCellsFunctional3D.h"
#include "immersedCellsFunctional3D.hh"
#include <algorithm>
#include <string>


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
void createImmersedCellParticles (
        MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
        TriangleBoundary3D<T>& boundary, plint tag, plint numPartsPerCell )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);
    applyProcessingFunctional (
        new CreateTaggedImmersedCellParticle3D<T,Descriptor>(boundary,tag,numPartsPerCell),
        particleField.getBoundingBox(), particleArg );
}

template<typename T, template<typename U> class Descriptor>
void translateCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                   const Box3D &outlet, std::vector<plint> &numParts,
                   std::vector<plint> &cellIds,
                   std::vector<Array<T,3> > &centers, std::vector<T> &radii, Array<T,3> const& translation )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);

    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        // count all particles of a certain tag in a buffer zone
        plint numPartsPerTag = countParticles(particleField, outlet, cellIds[iA]);
        // if all the particle of a certain tag are in a buffer zone
        if (numPartsPerTag > 0) {
            // then delete these particles in all the buffer zones
            plint before = countParticles(particleField,outlet,cellIds[iA]);

            applyProcessingFunctional ( // copy fluid velocity on particles
                            new TranslateTaggedParticlesFunctional3D<T,DESCRIPTOR>(translation, cellIds[iA]),
                            particleField.getBoundingBox(), particleArg);
//            applyProcessingFunctional (
//                new AbsorbTaggedParticlesFunctional3D<T,Descriptor>(translation, cellIds[iA]),
//                    outlet, particleArg );
            plint after = countParticles(particleField,outlet,cellIds[iA]);
            pcout << "translated particles = " << before << ", " << after << std::endl;
        }
    }
//    applyProcessingFunctional ( // update mesh position
//        new CopyParticleToVertex3D<T,DESCRIPTOR>(Cells.getMesh()),
//        particleField.getBoundingBox(), particleArg);
//    if (erased) {
//        pcout << "Particles absorbed : Number of particles per tag." << std::endl;
//        for (pluint iA = 0; iA < centers.size(); ++iA) {
//            pcout << cellIds[iA] << " , " <<  countParticles(particleField, particleField.getBoundingBox(), cellIds[iA]) << std::endl;
//        }
//    }
}


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
    bool created = false;
    plint numPartsPerTag = 0;
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        // count all particles of a certain tag in a buffer zone
        numPartsPerTag += countParticles(particleField,inlet,cellIds[iA]);
    }

    std::vector<plint> newcellIds;
    for (plint iA = slice*numOfCellsPerInlet; iA < (slice+1)*numOfCellsPerInlet; ++iA) {
        newcellIds.push_back(cellIds[iA]);
    }
    
    if (numPartsPerTag == 0) {
        createCells(Cells, newcellIds, numPartsPerCell, particleField);
        created = true;
        ++slice;
    }

//    if (created) {
//        pcout << "Particles created : Number of particles per tag." << std::endl;
//        for (pluint iA = 0; iA < newcellIds.size(); ++iA) {
//            pcout << newcellIds[iA] + 1<< " : " <<  countParticles(particleField, particleField.getBoundingBox(), newcellIds[iA]) << std::endl;
//        }
//    }
    return created;
}
template<typename T, template<typename U> class Descriptor>
void createCells(TriangleBoundary3D<T> &Cells,
                    const std::vector<plint> &cellIds,
                    plint numPartsPerCell,
                    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > &immersedParticles)
{
    for (pluint iA = 0; iA < cellIds.size(); ++iA) {
        createImmersedCellParticles(immersedParticles, Cells, cellIds[iA], numPartsPerCell);
    }
}
//
//template<typename T>
//TriangleSet<T> constructRBC(Array<T,3> const& center, T radius, plint minNumOfTriangles) {
//    TriangleSet<T> RBC("./lib/RBC.stl");
//    RBC.scale(radius/3.93);
//    RBC.rotate(pi/2.0, pi/2.0, 0.);
//    RBC.translate(center);
//    return RBC;
//}


template<typename T>
TriangleSet<T> constructRBC(Array<T,3> const& center, T radius, plint minNumOfTriangles) {
    return constructCell(center, radius, "./lib/RBC.stl");
}

template<typename T>
TriangleSet<T> constructCell(Array<T,3> const& center, T radius, std::string cellFilename) {
//    Cuboid<T> boundingCuboid;
    TriangleSet<T> Cell(cellFilename);
    Cuboid<T> cb = Cell.getBoundingCuboid();
    Array<T,3> dr = (cb.upperRightCorner - cb.lowerLeftCorner);
    T scaleFactor = std::max(dr[0],std::max(dr[1],dr[2]));
    Cell.scale(radius*2.0/scaleFactor);
    Cell.rotate(pi/2.0, pi/2.0, 0.);
    Cell.translate(center);
    return Cell;
}

template<typename T>
TriangleBoundary3D<T> createCompleteMesh(
    const std::vector<Array<T,3> > &centers, const std::vector<T> &radii,
    std::vector<plint> &cellIds, plint &numPartsPerCell, IncomprFlowParam<T> const& parameters,
    plint shape, std::string cellPath, plint minNumOfTriangles=100)
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
        	allTriangles.push_back(constructSphere<T>(center, radius, minNumOfTriangles));
        }
        else if (shape == 1) {
        	allTriangles.push_back(constructRBC<T>(center, radius, minNumOfTriangles));
        }
        else if (shape == 2) {
            allTriangles.push_back(constructCell<T>(center, radius, cellPath));
        }
        cellIds.push_back(iA);
    }
    wholeTriangleSet.merge(allTriangles);

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
    pcout << "num Vertices = " << Cells.getMesh().getNumVertices() << ", centers = " << centers.size() << std::endl;
    pcout << "num Triangles = " << Cells.getMesh().getNumTriangles() << ", centers = " << centers.size() << std::endl;
    pcout << "num Vertices per Cell= " << numPartsPerCell << ", modulo = " << modulo << std::endl;
    pcout << "num Triangles per Cell= " << Cells.getMesh().getNumTriangles()/centers.size() << ", modulo = " << modulo << std::endl;

    PLB_ASSERT(modulo == 0);

    Cells.cloneVertexSet(0);
    
    return Cells;
}


#endif  // IMMERSED_CELLS_3D_HH


