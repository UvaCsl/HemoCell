/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
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

#ifndef IMMERSEDBOUNDARYMETHOD_3D_H
#define IMMERSEDBOUNDARYMETHOD_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedBoundaryMethod3D.hh"
#include <vector>

namespace plb {

template<typename T>
T phi2 (T x) ;

template<typename T>
T phi3 (T x) ;

template<typename T>
T phi4 (T x) ;

template<typename T>
T phi4c (T x) ;

template<typename T>
void interpolationCoefficients (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights,
        plint ibmKernel=2);

template<typename T>
void interpolationCoefficientsPhi2 (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights );

template<typename T>
void interpolationCoefficientsPhi3 (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights );

template<typename T>
void interpolationCoefficientsPhi4 (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights );

template<typename T>
void interpolationCoefficientsPhi4c (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights );


}
#endif  // IMMERSEDBOUNDARYMETHOD_3D_H

