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
T phi2 (T x) {
    x = fabs(x);
    if (x <= 1.0) {
        return (1.0 - x);
    } else {
        return 0;
    }
}

template<typename T>
T phi3 (T x) {
    x = fabs(x);
    if (x <= 0.5) {
        return (1 + sqrt(1 + 3*x*x))/3.0;
    } else if (x<= 1.5) {
        return (5 - 3*x - sqrt(-2 + 6*x - 3*x*x) );
    } else {
        return 0;
    }
}

template<typename T>
T phi4 (T x) {
    x = fabs(x);
    if (x <= 1.0) {
        return 1.0/8.0 * (3 - 2*x + sqrt(1 + 4*x - 4*x*x));
    } else if (x<= 2.0) {
        return 1.0/8.0 * (5 - 2*x + sqrt(-7 + 12*x - 4*x*x));
    } else {
        return 0;
    }
}

template<typename T>
void interpolationCoefficientsPhi2 (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights )
{
    cellPos.clear();
    weights.clear();
    plint i = 0;
    plint x0=-1, x1=2;
    Box3D boundingBox(block.getBoundingBox());
    for (int dx = x0; dx < x1; ++dx) {
        for (int dy = x0; dy < x1; ++dy) {
            for (int dz = x0; dz < x1; ++dz) {
                Dot3D cellPosition(Dot3D( (plint) position[0] + dx, (plint) position[1] + dy, (plint) position[2] + dz));
                Dot3D cellPositionInDomain = cellPosition - block.getLocation(); // Convert cell position to local coordinates.
                if (contained(cellPositionInDomain,boundingBox)) {
                    T phi[3];
                    phi[0] = (position[0] - (T)cellPosition.x);
                    phi[1] = (position[1] - (T)cellPosition.y);
                    phi[2] = (position[2] - (T)cellPosition.z);
                    T weight = phi2(phi[0]) * phi2(phi[1]) * phi2(phi[2]);
                    if (weight>0) {
                        weights.push_back(weight);
                        cellPos.push_back(cellPositionInDomain);
                        i+=1;
                    }
                }
            }
        }
    }
}


template<typename T>
void interpolationCoefficientsPhi3 (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights )
{
    cellPos.clear();
    weights.clear();
    plint i = 0;
    plint x0=-2, x1=3;
    Box3D boundingBox(block.getBoundingBox());
//    pcout << "location( " << block.getLocation().x <<
//            block.getLocation().y <<
//            block.getLocation().z <<
//            "), bB =(" << boundingBox.x0 <<
//            ", " << boundingBox.x1 <<
//            "), (" << boundingBox.y0 <<
//            ", " << boundingBox.y1 <<
//            "), " << boundingBox.z0 <<
//            "," << boundingBox.z1 <<
//             "), p.x = (" << position[0] <<
//             ", " << position[1] <<
//             ", " << position[2] << ")" << std::endl;
    for (int dx = x0; dx < x1; ++dx) {
        for (int dy = x0; dy < x1; ++dy) {
            for (int dz = x0; dz < x1; ++dz) {
                Dot3D cellPosition(Dot3D( (plint) position[0] + dx, (plint) position[1] + dy, (plint) position[2] + dz));
                Dot3D cellPositionInDomain = cellPosition - block.getLocation(); // Convert cell position to local coordinates.
                if (contained(cellPositionInDomain,boundingBox)) {
                    T phi[3];
                    phi[0] = (position[0] - (T)cellPosition.x);
                    phi[1] = (position[1] - (T)cellPosition.y);
                    phi[2] = (position[2] - (T)cellPosition.z);
                    T weight = phi3(phi[0]) * phi3(phi[1]) * phi3(phi[2]);
                    if (weight>0) {
                        weights.push_back(weight);
                        cellPos.push_back(cellPositionInDomain);
                        i+=1;
                    }
                }
            }
        }
    }
}


template<typename T>
void interpolationCoefficientsPhi4 (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights )
{
    cellPos.clear();
    weights.clear();
    plint i = 0;
    plint x0=-2, x1=3;
    Box3D boundingBox(block.getBoundingBox());
//    pcout << "location( " << block.getLocation().x <<
//            block.getLocation().y <<
//            block.getLocation().z <<
//            "), bB =(" << boundingBox.x0 <<
//            ", " << boundingBox.x1 <<
//            "), (" << boundingBox.y0 <<
//            ", " << boundingBox.y1 <<
//            "), " << boundingBox.z0 <<
//            "," << boundingBox.z1 <<
//             "), p.x = (" << position[0] <<
//             ", " << position[1] <<
//             ", " << position[2] << ")" << std::endl;
    for (int dx = x0; dx < x1; ++dx) {
        for (int dy = x0; dy < x1; ++dy) {
            for (int dz = x0; dz < x1; ++dz) {
                Dot3D cellPosition(Dot3D( (plint) position[0] + dx, (plint) position[1] + dy, (plint) position[2] + dz));
                Dot3D cellPositionInDomain = cellPosition - block.getLocation(); // Convert cell position to local coordinates.
                if (contained(cellPositionInDomain,boundingBox)) {
                    T phi[3];
                    phi[0] = (position[0] - (T)cellPosition.x);
                    phi[1] = (position[1] - (T)cellPosition.y);
                    phi[2] = (position[2] - (T)cellPosition.z);
                    T weight = phi4(phi[0]) * phi4(phi[1]) * phi4(phi[2]);
                    if (weight>0) {
                        weights.push_back(weight);
                        cellPos.push_back(cellPositionInDomain);
                        i+=1;
                    }
                }
            }
        }
    }
}


}
#endif  // IMMERSEDBOUNDARYMETHOD_3D_HH
