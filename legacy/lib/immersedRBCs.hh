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

/* Main author: Lampros Mountrakis */

#ifndef IMMERSED_RBCS_3D_HH
#define IMMERSED_RBCS_3D_HH

#include "core/globalDefs.h"
#include "offLattice/triangleSetGenerator.h"
#include "offLattice/triangleSetGenerator.hh"
#include "offLattice/triangleSet.h"
#include "offLattice/triangleSet.hh"

namespace plb {

TriangleSet<T> * rbcTriangleSet = 0;

template<typename T>
TriangleSet<T> constructRBC(Array<T,3> const& center, T radius, plint minNumOfTriangles)
{
    if (0 == rbcTriangleSet) {
        rbcTriangleSet = new TriangleSet<T>("lib/RBC.stl");
    }
    TriangleSet<T> RBC(*rbcTriangleSet);
    RBC.scale(radius);
    RBC.translate(center);
    return RBC;
};

} // namespace plb

#endif  // IMMERSED_RBCS_3D_HH


