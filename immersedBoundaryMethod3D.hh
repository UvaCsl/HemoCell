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

#ifndef IMMERSEDBOUNDARYMETHOD_3D_HH
#define IMMERSEDBOUNDARYMETHOD_3D_HH

#include "core/globalDefs.h"
#include "core/util.h"
#include "immersedBoundaryMethod3D.h"
#include <vector>

/* ******** Function linearInterpolationCoefficients27PointStensil ********************* */
namespace plb {

template<typename T>
void linearInterpolationCoefficients8PointStensil (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights )
{
    cellPos.resize(8);
    weights.resize(8);
    plint i = 0;
    for (int dx = 0; dx < 2; ++dx) {
        for (int dy = 0; dy < 2; ++dy) {
            for (int dz = 0; dz < 2; ++dz) {
                cellPos[i] = Dot3D( (plint) position[0] + dx, (plint) position[1] + dy, (plint) position[2] + dz);
                T u = fabs(position[0] - (T)cellPos[i].x);
                T v = fabs(position[1] - (T)cellPos[i].y);
                T w = fabs(position[2] - (T)cellPos[i].z);
                if ( (u > 1.0) || (v > 1.0) || (w > 1.0)) {
                    weights[i] = 0.0;
                } else {
                    weights[i] = (1.-u) * (1.-v) * (1.-w);
                }
//                pcout << "i" << i<< " dx " << dx << " " <<
//                        "dy " << dy << " " <<
//                        "dz " << dz << " " <<
//                        "weight " << weights[i] << "\n";
            i+=1;
            }
        }
    }
    // Convert cell position to local coordinates.
    for (plint iPos=0; iPos<8; ++iPos) {
        cellPos[iPos] -= block.getLocation();
    }
}

}
#endif  // IMMERSEDBOUNDARYMETHOD_3D_HH
