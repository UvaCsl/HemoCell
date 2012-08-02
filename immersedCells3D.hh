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

#include "core/globalDefs.h"
#include "immersedCells3D.h"
using namespace plb;
using namespace std;

plint xDirection  = 0;
plint yDirection  = 1;
plint zDirection  = 2;
plint borderWidth     = 1;  // Because Guo acts in a one-cell layer.
// Requirement: margin>=borderWidth.
plint margin          = 3;  // Extra margin of allocated cells around the obstacle,
                             //   for moving walls.
plint extraLayer      = 0;  // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
const plint blockSize = 0; // Zero means: no sparse representation.
const plint extendedEnvelopeWidth = 2;  // Because Guo needs 2-cell neighbor access.
const plint particleEnvelopeWidth = 6;



template<typename T, template<typename U> class Descriptor>
void createImmersedWallParticles (
        MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
        TriangleBoundary3D<T>& boundary, plint tag, plint numPartsPerBloodCell )
{
    boundary.pushSelect(0,1);
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);
    applyProcessingFunctional (
        new CreateTaggedImmersedWallParticle3D<T,Descriptor>(boundary,tag,numPartsPerBloodCell),
        particleField.getBoundingBox(), particleArg );
    boundary.popSelect();
}

template<typename T, template<typename U> class Descriptor>
void deleteBloodCell(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                   const Box3D &outlet, std::vector<plint> &numParts,
                   std::vector<plint> &tags, TriangleBoundary3D<T> &bloodCells,
                   std::vector<Array<T,3> > &centers, std::vector<plint > &radii )
{
    bool erased = false;
    for (pluint iA = 0; iA < tags.size(); ++iA) {
        // count all particles of a certain tag in a buffer zone
        plint numPartsPerTag = countParticles(particleField,outlet,tags[iA]);
        // if all the particle of a certain tag are in a buffer zone
        if (numPartsPerTag == numParts[iA] && numPartsPerTag > 0) {
            // then delete these particles in all the buffer zones
            plint before = countParticles(particleField,outlet,tags[iA]);
            std::vector<MultiBlock3D*> particleArg;
            particleArg.push_back(&particleField);

            applyProcessingFunctional (
                new AbsorbTaggedParticlesFunctional3D<T,Descriptor>(tags[iA]),
                    outlet, particleArg );

            plint after = countParticles(particleField,outlet,tags[iA]);

            pcout << "erased particles = " << before << ", " << after << std::endl;
            
//             delete bloodCells[iA];                     // delete the iA-th pointer mesh
//             numParts.erase(numParts.begin()+iA);        // erase the iA-th number of particles
//             tags.erase(tags.begin()+iA);                // erase the iA-th tag
//             bloodCells.erase(bloodCells.begin()+iA);  // erase the iA-th mesh
//             centers.erase(centers.begin()+iA);          // delete the iA-th center
//             radii.erase(radii.begin()+iA);              // delete the iA-th radius
//             --iA;

            erased = true;
        }
    }

    if (erased) {
        pcout << "Particles absorbed : Number of particles per tag." << std::endl;
        for (pluint iA = 0; iA < centers.size(); ++iA) {
            pcout << tags[iA] << " , " <<  countParticles(particleField, particleField.getBoundingBox(), tags[iA]) << std::endl;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
bool generateBloodCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                      const Box3D &inlet, std::vector<plint> &tags, TriangleBoundary3D<T> &bloodCells,
                      plint numPartsPerBloodCell, plint numOfBloodCellsPerInlet, plint &slice )
{
    bool created = false;
    plint numPartsPerTag = 0;
    for (pluint iA = 0; iA < tags.size(); ++iA) {
        // count all particles of a certain tag in a buffer zone
        numPartsPerTag += countParticles(particleField,inlet,tags[iA]);
    }

    std::vector<plint> newTags;
    for (plint iA = slice*numOfBloodCellsPerInlet; iA < (slice+1)*numOfBloodCellsPerInlet; ++iA) {
        newTags.push_back(tags[iA]);
    }
    
    if (numPartsPerTag == 0) {
        createBloodCells(bloodCells, newTags, numPartsPerBloodCell, particleField);
        created = true;
        ++slice;
    }

    if (created) {
        pcout << "Particles created : Number of particles per tag." << std::endl;
        for (pluint iA = 0; iA < newTags.size(); ++iA) {
            pcout << newTags[iA] << " , " <<  countParticles(particleField, particleField.getBoundingBox(), newTags[iA]) << std::endl;
        }
    }
    return created;
}
template<typename T, template<typename U> class Descriptor>
void createBloodCells(TriangleBoundary3D<T> &bloodCells,
                    const std::vector<plint> &tags,
                    plint numPartsPerBloodCell,
                    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > &immersedParticles)
{
    for (pluint iA = 0; iA < tags.size(); ++iA) {
        createImmersedWallParticles(immersedParticles, bloodCells, tags[iA], numPartsPerBloodCell);
    }
}

template<typename T>
TriangleBoundary3D<T> createCompleteMesh(
    const std::vector<Array<T,3> > &centers, const std::vector<plint> &radii,
    std::vector<plint> &tags, plint &numPartsPerBloodCell)
{
    PLB_ASSERT(centers.size() == radii.size());

    // creation of the whole mesh including the one that is not used by particles
    std::vector<TriangleSet<T> > allTriangles;
    TriangleSet<T> wholeTriangleSet; 
    for (pluint iA = 0; iA < centers.size(); ++iA) {
        Array<T,3> center(centers[iA]);
        plint radius = radii[iA];
        allTriangles.push_back(constructSphere<T>(center, radius, 100));
        tags.push_back(iA);
    }
    wholeTriangleSet.merge(allTriangles);

    Dot3D location(centers[0][0]-radii[0],centers[0][1]-radii[0],centers[0][2]-radii[0]);
    DEFscaledMesh<T> defMesh (
            wholeTriangleSet, 10, xDirection, margin, location );
//         pcout << "Original sphere at location [" << location.x << ","
//             << location.y << "," << location.z << "] with radius " << radius << std::endl;

    TriangleBoundary3D<T> bloodCells(defMesh);
    bloodCells.getMesh().inflate();

    numPartsPerBloodCell = bloodCells.getMesh().getNumVertices() / centers.size();
    plint modulo = bloodCells.getMesh().getNumVertices() % centers.size();
//     pcout << "num parts per triangles = " << numPartsPerBloodCell << ", modulo = " << modulo << std::endl;
    pcout << "num vertices = " << bloodCells.getMesh().getNumVertices() << ", centers = " << centers.size() << std::endl;
    pcout << "num parts per triangles = " << numPartsPerBloodCell << ", modulo = " << modulo << std::endl;
    PLB_ASSERT(modulo == 0);

    bloodCells.cloneVertexSet(0);
    
    return bloodCells;
}

//template<typename T>
//TriangleBoundary3D<T> createCompleteMeshRBCs(
//	    const std::vector<Array<T,3> > &centers, const std::vector<plint> &radii,
//	std::vector<plint> &tags, plint &numPartsPerBloodCell)
//{
//	PLB_ASSERT(centers.size() == radii.size());
//
//	// creation of the whole mesh including the one that is not used by particles
//	std::vector<TriangleSet<T> > allTriangles;
//	TriangleSet<T> wholeTriangleSet;
//	for (pluint iA = 0; iA < centers.size(); ++iA) {
//		Array<T,3> center(centers[iA]);
//		plint radius = radii[iA];
//		constructRBC<T> RBC(center, radius, 100);
//		allTriangles.push_back(RBC);
//		tags.push_back(iA);
//	}
//	wholeTriangleSet.merge(allTriangles);
//
//	Dot3D location(centers[0][0]-radii[0],centers[0][1]-radii[0],centers[0][2]-radii[0]);
//	DEFscaledMesh<T> defMesh (
//			wholeTriangleSet, 10, xDirection, margin, location );
////         pcout << "Original sphere at location [" << location.x << ","
////             << location.y << "," << location.z << "] with radius " << radius << std::endl;
//
//	TriangleBoundary3D<T> bloodCells(defMesh);
//	bloodCells.getMesh().inflate();
//
//	numPartsPerBloodCell = bloodCells.getMesh().getNumVertices() / centers.size();
//	plint modulo = bloodCells.getMesh().getNumVertices() % centers.size();
////     pcout << "num parts per triangles = " << numPartsPerBloodCell << ", modulo = " << modulo << std::endl;
//	pcout << "num vertices = " << bloodCells.getMesh().getNumVertices() << ", centers = " << centers.size() << std::endl;
//	pcout << "num parts per triangles = " << numPartsPerBloodCell << ", modulo = " << modulo << std::endl;
//	PLB_ASSERT(modulo == 0);
//
//	bloodCells.cloneVertexSet(0);
//
//	return bloodCells;
//}

#endif  // IMMERSED_CELLS_3D_HH


