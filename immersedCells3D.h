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

#ifndef IMMERSED_CELLS_3D_H
#define IMMERSED_CELLS_3D_H

#include "core/globalDefs.h"
#include "core/array.h"

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedWallParticle3D.h"
#include "immersedWallParticle3D.hh"
#include "immersedWallParticleFunctional3D.h"
#include "immersedWallParticleFunctional3D.hh"
#include "immersedWallParticleVtk3D.h"
#include "immersedWallParticleVtk3D.hh"
#include "shellModel3D.h"
#include "shellModel3D.hh"
#include "immersedCells3D.h"
#include "immersedCells3D.hh"
// #include "lib/immersedRBCs.hh"
#include <vector>
using namespace plb;
using namespace std;


template<typename T, template<typename U> class Descriptor>
void createImmersedWallParticles (
        MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
        TriangleBoundary3D<T>& boundary, plint tag, plint numPartsPerBloodCell );



template<typename T, template<typename U> class Descriptor>
void deleteBloodCell(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                   const Box3D &outlet, std::vector<plint> &numParts,
                   std::vector<plint> &tags, TriangleBoundary3D<T> &bloodCells,
                   std::vector<Array<T,3> > &centers, std::vector<plint > &radii );



template<typename T, template<typename U> class Descriptor>
bool generateBloodCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                      const Box3D &inlet, std::vector<plint> &tags, TriangleBoundary3D<T> &bloodCells,
                      plint numPartsPerBloodCell, plint numOfBloodCellsPerInlet, plint &slice );

template<typename T, template<typename U> class Descriptor>
void createBloodCells(TriangleBoundary3D<T> &bloodCells,
                    const std::vector<plint> &tags,
                    plint numPartsPerBloodCell,
                    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > &immersedParticles);

template<typename T>
TriangleBoundary3D<T> createCompleteMesh(
    const std::vector<Array<T,3> > &centers, const std::vector<plint> &radii,
    std::vector<plint> &tags, plint &numPartsPerBloodCell) ;

template<typename T>
TriangleBoundary3D<T> createCompleteMeshRBCs	(
    const std::vector<Array<T,3> > &centers, const std::vector<plint> &radii,
    std::vector<plint> &tags, plint &numPartsPerBloodCell) ;

#endif  // IMMERSED_CELLS_3D_H


