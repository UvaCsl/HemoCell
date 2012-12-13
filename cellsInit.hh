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


#include "cellsInit.h"


void positionCells(plint shape, T radius, plint & npar, IncomprFlowParam<T> const& parameters,
        std::vector<Array<T,3> > & centers, std::vector<T> & radii, plint flowType) {

    std::vector<T> posX, posY, posZ;
    const plint nx = parameters.getNx() ;
    const plint ny = parameters.getNy()  ;
    const plint nz = parameters.getNz()  ;
    const T dX = 2.1 * radius ;
    const T dY = 2.1 * radius * ( (shape==1) ? 0.265106361 : 1 );
    const T dZ = 2.1 * radius;

    plint NdX = (nx-1)*1.0/dX;
    plint NdY = (ny-1)*1.0/dY;
    plint NdZ = (nz-1)*1.0/dZ;
    npar = npar<(NdX*NdY*NdZ)?npar:(NdX*NdY*NdZ);
    plint slices = npar/(NdY*NdZ);

    for (plint i = 0; i < slices; ++i) {
        for (plint iy = 0; iy < NdY; ++iy) {
            for (plint iz = 0; iz < NdZ; ++iz) {
                posX.push_back((i+0.5)*dX);
                posY.push_back((iy+0.5)*dY);
                posZ.push_back((iz+0.5)*dZ);
            }
        }
    }
    plint lastSlice = npar%(NdY*NdZ);
    plint rows = lastSlice/NdZ;
    for (plint iy = 1; iy <= rows; ++iy) {
        for (plint iz = 0; iz < NdZ; ++iz) {
            posX.push_back((slices+0.5)*dX);
            posY.push_back(ny * iy*1.0/(rows+1 + 1.0));
            posZ.push_back((iz+0.5)*dZ);
        }
    }

    plint mods = lastSlice%NdZ;
    for (plint iz = 1; iz <= mods; ++iz) {
        posX.push_back((slices+0.5)*dX);
        posY.push_back(ny * (rows+1)*1.0/(rows+1 + 1.0));
        posZ.push_back(nz * iz*1.0/(mods + 1.0));
    }

    T addToX = 0.0;
    if (flowType == 1) {
        addToX = (NdX - slices) * dX * 0.5;
    }


    for (pluint iA = 0; iA < posX.size(); ++iA) {
        centers.push_back(Array<T,3>(posX[iA]+addToX,posY[iA],posZ[iA]));
        radii.push_back(radius);
    }
}





#endif  // CELLS_INIT_HH


