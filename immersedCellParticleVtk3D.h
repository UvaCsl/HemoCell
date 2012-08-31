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

#include "core/globalDefs.h"
#include "particles/multiParticleField3D.h"
#include "offLattice/triangleBoundary3D.h"
#include <vector>
#include <string>

namespace plb {

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void writeImmersedSurfaceVTK( TriangleBoundary3D<T> const& boundary,
                      MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles,
                      std::vector<std::string> const& scalars,
                      std::vector<std::string> const& vectors,
                      std::string const& fName, bool dynamicMesh, plint tag );


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void writeImmersedSurfaceVTK( TriangleBoundary3D<T> const& boundary,
                      MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles,
                      std::vector<std::string> const& scalars,
                      std::vector<std::string> const& vectors,
                      std::string const& fName, bool dynamicMesh, plint tag,
                      std::vector<T> const& scalarFactor, std::vector<T> const& vectorFactor );


template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
void vtkForImmersedVertices(std::vector<Particle3D<T,Descriptor>*> const& particles,
                    TriangleBoundary3D<T> const& boundary,
                    std::vector<std::string> const& scalars,
                    std::vector<std::string> const& vectors,
                    std::string fName, bool dynamicMesh, plint tag,
                    std::vector<T> const& scalarFactor, std::vector<T> const& vectorFactor );

}  // namespace plb

#endif  // IMMERSED_PARTICLE_VTK_3D_H