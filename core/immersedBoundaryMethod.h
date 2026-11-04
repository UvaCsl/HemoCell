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

#include <cmath>
#include <vector>

namespace hemo {

/// Decide if a Lagrangian point is contained in 3D box, boundaries exclusive
inline bool contained_sane(hemo::Array<plint, 3> const & x, Box3D const & box) {
  return x[0] >= box.x0 && x[0] <= box.x1 && x[1] >= box.y0 && x[1] <= box.y1 && x[2] >= box.z0 && x[2] <= box.z1;
}
inline T phi2(T x) {
  x = fabs(x);
  x = 1.0 - x;
  return max(x, (T)0.0);
}

template <typename T>
T phi3(T x);

template <typename T>
T phi4(T x);

template <typename T>
T phi4c(T x);

template <typename T, template <typename U> class Descriptor>
void interpolationCoefficients(BlockLattice3D<T, Descriptor> const & block, hemo::Array<T, 3> const & position,
                               std::vector<Dot3D> & cellPos, std::vector<T> & weights);

template <typename T, template <typename U> class Descriptor>
void interpolationCoefficientsPhi1(BlockLattice3D<T, Descriptor> const & block, hemo::Array<T, 3> const & position,
                                   std::vector<Dot3D> & cellPos, std::vector<T> & weights);

/**
 * @brief Compute the IBM kernel of a particle with the 2-point hat function
 * phi2 (trilinear weights): `kernelLocations`/`kernelWeights` (and, with
 * INTERIOR_VISCOSITY, `kernelCoordinates`).
 *
 * phi2 is zero for distances >= 1, so along each axis only the nodes
 * floor(p) and floor(p)+1 can contribute (at most 2x2x2 nodes). Nodes outside
 * the block, with zero weight, or with boundary dynamics are skipped; the
 * remaining weights are normalised to sum to one. Nodes are visited in
 * ascending x, then y, then z order.
 *
 * @param block    fluid block of the particle (block-local node coordinates)
 * @param particle particle with absolute lattice position `sv.position`
 */
inline void interpolationCoefficientsPhi2(BlockLattice3D<T, DESCRIPTOR> & block, HemoCellParticle & particle) {
  particle.kernelWeights.clear();
  particle.kernelWeights.reserve(8);
  particle.kernelLocations.clear();
#ifdef INTERIOR_VISCOSITY
  particle.kernelCoordinates.clear();
  particle.kernelCoordinates.reserve(8);
#endif

  const Dot3D location = block.getLocation();
  const Box3D boundingBox = block.getBoundingBox();
  const plint lower[3] = {boundingBox.x0, boundingBox.y0, boundingBox.z0};
  const plint upper[3] = {boundingBox.x1, boundingBox.y1, boundingBox.z1};

  // Position relative to the block.
  const T position[3] = {particle.sv.position[0] - location.x, particle.sv.position[1] - location.y,
                         particle.sv.position[2] - location.z};

  // Candidate nodes and their one-dimensional weights per axis.
  plint node[3][2];
  T weight1D[3][2];
  bool inside[3][2];
  for (int d = 0; d < 3; ++d) {
    const plint base = plint(std::floor(position[d]));
    for (int i = 0; i < 2; ++i) {
      node[d][i] = base + i;
      weight1D[d][i] = phi2(position[d] - node[d][i]);
      inside[d][i] = node[d][i] >= lower[d] && node[d][i] <= upper[d];
    }
  }

  T total_weight = 0;
  for (int i = 0; i < 2; ++i) {
    if (!inside[0][i]) { continue; }
    for (int j = 0; j < 2; ++j) {
      if (!inside[1][j]) { continue; }
      for (int k = 0; k < 2; ++k) {
        if (!inside[2][k]) { continue; }
        const T weight = weight1D[0][i] * weight1D[1][j] * weight1D[2][k];
        if (weight == 0.0) { continue; }

        plb::Cell<T, DESCRIPTOR> & cell = block.get(node[0][i], node[1][j], node[2][k]);
        if (cell.getDynamics().isBoundary()) { continue; }

        total_weight += weight;
        particle.kernelWeights.push_back(weight);
        particle.kernelLocations.push_back(&cell);
#ifdef INTERIOR_VISCOSITY
        particle.kernelCoordinates.push_back({node[0][i], node[1][j], node[2][k]});
#endif
      }
    }
  }
  const T weight_coeff = 1.0 / total_weight;
  for (T & weight_ : particle.kernelWeights) { //Normalize weight to 1
    weight_ *= weight_coeff;
  }
}

template <typename T, template <typename U> class Descriptor>
void interpolationCoefficientsPhi3(

    BlockLattice3D<T, Descriptor> const & block, hemo::Array<T, 3> const & position, std::vector<Dot3D> & cellPos,
    std::vector<T> & weights);

void interpolationCoefficientsPhi4(BlockLattice3D<double, DESCRIPTOR> const & block,
                                   hemo::Array<double, 3> const & position, std::vector<Dot3D> & cellPos,
                                   std::vector<double> & weights);

template <typename T, template <typename U> class Descriptor>
void interpolationCoefficientsPhi4c(BlockLattice3D<T, Descriptor> const & block, hemo::Array<T, 3> const & position,
                                    std::vector<Dot3D> & cellPos, std::vector<T> & weights);

/*
 * In case one of the interpolating boundary nodes is a boundary,
 * the force is spread to the rest of the nodes
 * */
template <typename T, template <typename U> class Descriptor>
void curateInterpolationCoefficients(BlockLattice3D<T, Descriptor> & fluid, std::vector<Dot3D> & cellPos,
                                     std::vector<T> & weights);


}
#endif // IMMERSEDBOUNDARYMETHOD_3D_H
