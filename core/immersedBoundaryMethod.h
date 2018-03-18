/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef IMMERSEDBOUNDARYMETHOD_H
#define IMMERSEDBOUNDARYMETHOD_H

#include "hemocell_internal.h"
#include <vector>

namespace plb {

/// Decide if a Lagrangian point is contained in 3D box, boundaries exclusive
inline bool contained_sane(hemo::Array<plint,3> const& x, Box3D const& box) {
    return x[0]>=box.x0 && x[0]<=box.x1 &&
           x[1]>=box.y0 && x[1]<=box.y1 &&
           x[2]>=box.z0 && x[2]<=box.z1;
}
inline double phi2 (double x) {
    x = fabs(x);
    x = 1.0 - x;
    return max(x,0.0);
}

template<typename T>
T phi3 (T x) ;

template<typename T>
T phi4 (T x) ;

template<typename T>
T phi4c (T x) ;

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficients (
        BlockLattice3D<T,Descriptor> const& block, hemo::Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficientsPhi1 (
        BlockLattice3D<T,Descriptor> const& block, hemo::Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

inline void interpolationCoefficientsPhi2 (
        BlockLattice3D<double,DESCRIPTOR> & block, HemoCellParticle * particle)
{
    //Clean current
    particle->kernelWeights.clear();
    particle->kernelLocations.clear();
    
    // Fixed kernel size
    const plint x0=-1, x1=2; //const for nice loop unrolling

    //Coordinates are relative
    const Dot3D tmpDot = block.getLocation(); 
    const hemo::Array<plint,3> relLoc = {tmpDot.x, tmpDot.y, tmpDot.z};

    //Get position, relative
    const hemo::Array<double,3> position_tmp = particle->position;
    const hemo::Array<double,3> position = {position_tmp[0] -relLoc[0], position_tmp[1]-relLoc[1],position_tmp[2]-relLoc[2]};

    //Get our reference node (0,0)
    const hemo::Array<plint,3> center({plint(position[0] + 0.5), plint(position[1] + 0.5), plint(position[2] + 0.5)}); 
    
    //Boundingbox of lattice
    const Box3D boundingBox = block.getBoundingBox();
    
    //Prealloc is better than JItalloc
    hemo::Array<plint,3> posInBlock;
    double phi[3];
    double weight;
    double total_weight = 0;
    
    for (int dx = x0; dx < x1; ++dx) {
        for (int dy = x0; dy < x1; ++dy) {
            for (int dz = x0; dz < x1; ++dz) {
                posInBlock = {center[0] + dx, center[1] + dy, center[2] + dz};
                
                //Sanity checks, skip if outside domain or boundary
                if (!contained_sane(posInBlock,boundingBox)) {
                  continue;
                }

                phi[0] = (position[0] - posInBlock[0]); //Get absolute distance
                phi[1] = (position[1] - posInBlock[1]);
                phi[2] = (position[2] - posInBlock[2]);
                weight = phi2(phi[0]) * phi2(phi[1]) * phi2(phi[2]);
                
                if (weight  == 0.0){
                  continue;
                }
                if (block.get(posInBlock[0],posInBlock[1],posInBlock[2]).getDynamics().isBoundary()) {
                  continue;
                }              
                
                total_weight+=weight;

                particle->kernelWeights.push_back(weight);
                particle->kernelLocations.push_back(&block.get(posInBlock[0],posInBlock[1],posInBlock[2]));
            }
        }
    }
    const double weight_coeff = 1.0 / total_weight;
    for(double & weight_ : particle->kernelWeights) { //Normalize weight to 1
      weight_ *= weight_coeff;
    }
}

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficientsPhi3 (
    
        BlockLattice3D<T,Descriptor> const& block, hemo::Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

void interpolationCoefficientsPhi4 (
        BlockLattice3D<double,DESCRIPTOR> const& block, hemo::Array<double,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<double>& weights);

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficientsPhi4c (
        BlockLattice3D<T,Descriptor> const& block, hemo::Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

/*
 * In case one of the interpolating boundary nodes is a boundary,
 * the force is spread to the rest of the nodes
 * */
template<typename T, template<typename U> class Descriptor>
void curateInterpolationCoefficients (BlockLattice3D<T,Descriptor>& fluid,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights) ;


}
#endif  // IMMERSEDBOUNDARYMETHOD_3D_H

