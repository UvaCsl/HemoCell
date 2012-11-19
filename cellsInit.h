/* This file is part of the ficsion framework.

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

#ifndef CELLS_INIT_H
#define CELLS_INIT_H


#include "ficsionInit.hh"
#include "immersedCells3D.h"
#include "immersedCells3D.hh"
#include "immersedCellsFunctional3D.h"
#include "immersedCellsFunctional3D.hh"

void positionCells(plint shape, T radius, plint npar, IncomprFlowParam<T> const& parameters,
        std::vector<Array<T,3> > & centers, std::vector<T> & radii) ;




#endif  // CELLS_INIT_H


