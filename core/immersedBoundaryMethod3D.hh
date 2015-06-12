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

#ifndef PI__
#define PI__
const double pi = 4.*atan(1.);
#endif  // PI__

#include "immersedBoundaryMethod3D.h"

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
        return 1.0/3.0 * (1 + sqrt(1 - 3*x*x));
    } else if (x<= 1.5) {
        return 1.0/6.0 * (5 - 3*x - sqrt(-2 + 6*x - 3*x*x) );
    } else {
        return 0;
    }
}

template<typename T>
T phi4 (T x) {
    x = fabs(x);
    if (x <= 1.0) {
        return 0.125 * (3 - 2*x + sqrt(1 + 4*x - 4*x*x));
    } else if (x<= 2.0) {
        return 0.125 * (5 - 2*x - sqrt(-7 + 12*x - 4*x*x));
    } else {
        return 0;
    }
}

template<typename T>
T phi4c (T x) {
    x = fabs(x);
    if (x<=2)
        return 0.25*(1 + cos(pi*x*0.5));
    else
        return 0.0;
}

template<typename T>
void interpolationCoefficients (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights,
        plint ibmKernel) {
    /*
     * ibmKernel == 2, Phi2 ! Default
     * ibmKernel == 3, Phi3
     * ibmKernel == 4, Phi4
     * ibmKernel == 5, Phi4c
     */
    if (ibmKernel == 2) {
        interpolationCoefficientsPhi2(block, position, cellPos, weights);
    } else if (ibmKernel == 3) {
        interpolationCoefficientsPhi3(block, position, cellPos, weights);
    } else if (ibmKernel == 4) {
        interpolationCoefficientsPhi4(block, position, cellPos, weights);
    } else if (ibmKernel == 5) {
        interpolationCoefficientsPhi4c(block, position, cellPos, weights);
    } else {
        interpolationCoefficientsPhi4(block, position, cellPos, weights);
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

template<typename T>
void interpolationCoefficientsPhi4c (
        AtomicBlock3D const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights )
{
    cellPos.clear();
    weights.clear();
    plint i = 0;
    plint x0=-2, x1=3;
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
                    T weight = phi4c(phi[0]) * phi4c(phi[1]) * phi4c(phi[2]);
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

/*
 * In case one of the interpolating boundary nodes is a boundary,
 * the force is spread to the rest of the nodes
 * */
template<typename T, template<typename U> class Descriptor>
void curateInterpolationCoefficients (BlockLattice3D<T,Descriptor>& fluid,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights) {
    Cell<T,Descriptor>* cell;
    T percentageInBoundaries=1;
    std::vector<T> newWeights = weights;
    std::vector<Dot3D> newCellPos = cellPos;
    weights.clear(); cellPos.clear();
    for (pluint iCell = 0; iCell < newWeights.size() ; ++iCell) {
        Dot3D cellPosition = newCellPos[iCell];
        cell = &fluid.get(cellPosition.x, cellPosition.y, cellPosition.z);
        if (not cell->getDynamics().isBoundary()) {
            percentageInBoundaries -= newWeights[iCell];
            weights.push_back(newWeights[iCell]);
            cellPos.push_back(newCellPos[iCell]);
        }
    }
    plint nRemaining = weights.size();
    for (plint iCell = 0; iCell < nRemaining; ++iCell) {
        weights[iCell] += percentageInBoundaries*1.0/nRemaining;
    }
}

}
#endif  // IMMERSEDBOUNDARYMETHOD_3D_HH


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
