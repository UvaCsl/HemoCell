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

#ifndef CELLS_INIT_HH
#define CELLS_INIT_HH


#include "cellsInit.hh"


void positionCells(plint shape, T radius, plint npar, IncomprFlowParam<T> const& parameters,
        std::vector<Array<T,3> > & centers, std::vector<T> & radii) {

    std::vector<T> posY, posZ;
    T diameter = 2*radius;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    for (T r =1 + diameter; r < ny-diameter-1; r+=1.05*diameter) {
        posZ.push_back(r);
    }
    if (shape==1) { // Y and Z for RBC are not symmetric, pack it more dense
        T dY =    0.265106361*radius;
        for (T r = 1 + 2*dY; r < ny-2*dY - 1; r+=2*1.05*2*dY) {
            posY.push_back(r);
        }
    } else {
        posY = posZ;
    }
    plint n=0;
    for (T iN = 3+diameter; (iN < nx-(diameter) - 1); iN+=1.05*diameter) { // create as many particles as possible
        for (pluint iA = 0; iA < posY.size(); ++iA) {
            for (pluint iB = 0; iB < posZ.size() && (n < npar); ++iB) {
                centers.push_back(Array<T,3>(iN,posY[iA],posZ[iB]));
                radii.push_back(radius);
                n++;
            }
        }
    }
    if (nz == 0) {}
}





#endif  // CELLS_INIT_HH


