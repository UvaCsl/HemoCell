/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011 FlowKit Sarl
 * Avenue de Chailly 23
 * 1012 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IMMERSED_PARTICLE_VTK_3D_H
#define IMMERSED_PARTICLE_VTK_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <string>
#include <map>
#include "surfaceParticle3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void writeImmersedSurfaceVTK( TriangleBoundary3D<T> const& boundary,
                      MultiParticleField3D<LightParticleField3D<T,Descriptor> >& particles,
                      std::string const& fName);

template<typename T, template<typename U> class Descriptor>
void vtkForImmersedVertices(std::vector<Particle3D<T,Descriptor>*> const& particles,
                    TriangleBoundary3D<T> const& boundary,
                    std::string fName,
                    plint nx, plint ny, plint nz);

template<typename T>
void writeImmersedPointsVTK(TriangleBoundary3D<T> const& boundary, std::vector<plint> const& vertices, T const& dx,
                    std::string fName);

template<typename T>
void writeImmersedPointsVTK(std::vector<Array<T,3> > const& positions, std::vector<plint> const& tags, T const& dx,
                    std::string fName);
}  // namespace plb

#include "immersedCellParticleVtk3D.hh"

#endif  // IMMERSED_PARTICLE_VTK_3D_H
