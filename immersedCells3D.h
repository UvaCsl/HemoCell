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

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedCellParticle3D.hh"
#include "immersedCellParticleFunctional3D.h"
#include "immersedCellParticleFunctional3D.hh"
#include "immersedCellParticleVtk3D.h"
#include "immersedCellParticleVtk3D.hh"
#include "shellModel3D.h"
#include "shellModel3D.hh"
#include "immersedCells3D.h"
#include "immersedCells3D.hh"
// #include "lib/immersedRBCs.hh"
#include <vector>
using namespace plb;
using namespace std;


template<typename T, template<typename U> class Descriptor>
void createImmersedCellParticles (
        MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
        TriangleBoundary3D<T>& boundary, plint tag, plint numPartsPerCell );


template<typename T, template<typename U> class Descriptor>
void translateCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                   const Box3D &outlet, std::vector<plint> &numParts,
                   std::vector<plint> &cellIds, TriangleBoundary3D<T> &Cells,
                   std::vector<Array<T,3> > &centers, std::vector<T> &radii, Array<T,3> const& translation );


template<typename T, template<typename U> class Descriptor>
void deleteCell(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                   const Box3D &outlet, std::vector<plint> &numParts,
                   std::vector<plint> &cellIds, TriangleBoundary3D<T> &Cells,
                   std::vector<Array<T,3> > &centers, std::vector<T > &radii );


template<typename T, template<typename U> class Descriptor>
bool generateCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
                      const Box3D &inlet, std::vector<plint> &cellIds, TriangleBoundary3D<T> &Cells,
                      plint numPartsPerCell, plint numOfCellsPerInlet, plint &slice );


template<typename T, template<typename U> class Descriptor>
void createCells(TriangleBoundary3D<T> &Cells,
                    const std::vector<plint> &cellIds,
                    plint numPartsPerCell,
                    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > &immersedParticles);


template<typename T>
TriangleSet<T> constructRBC(Array<T,3> const& center, T radius, plint minNumOfTriangles) ;


template<typename T>
TriangleSet<T> constructCell(Array<T,3> const& center, T radius, std::string cellFilename);


template<typename T>
TriangleBoundary3D<T> createCompleteMesh(
    const std::vector<Array<T,3> > &centers, const std::vector<T> &radii,
    std::vector<plint> &cellIds, IncomprFlowParam<T> const& parameters,
    plint shape, std::string cellPath, plint &cellNumTriangles, plint &numPartsPerCell);


#endif  // IMMERSED_CELLS_3D_H


