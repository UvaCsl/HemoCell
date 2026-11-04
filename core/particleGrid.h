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
#ifndef HEMO_PARTICLEGRID_H
#define HEMO_PARTICLEGRID_H

#include "hemoCellParticle.h"

#include <vector>

namespace hemo {

/**
 * @brief Cell list that bins the particles of one atomic block onto its
 * lattice nodes.
 *
 * Used for short-range neighbour searches (particle-particle repulsion,
 * particle-wall repulsion, solidification). Each particle is assigned to the
 * node `int(position - location + 0.5)` per axis, i.e. the same truncating
 * rounding that HemoCell used before this class existed. Particles whose node
 * falls outside `[0,nx) x [0,ny) x [0,nz)` are not binned.
 *
 * Storage is compressed (CSR-like): `nodeStart` holds `nx*ny*nz+1` offsets
 * into `particleIndices`, so the particles of node `n` are
 * `particleIndices[nodeStart[n] .. nodeStart[n+1])`. Unlike the previous
 * fixed-capacity grid there is no limit on the number of particles per node,
 * and no capacity checks are needed while filling it. Within a node, particle
 * indices are stored in ascending order.
 *
 * The stored indices refer to the particle vector passed to build(). Any
 * insertion, removal or reordering of that vector, and any change of particle
 * positions, invalidates the grid; the owner must rebuild it before use.
 *
 * Not thread safe; one instance per atomic block.
 */
class ParticleGrid {
public:
  /**
   * @brief (Re)build the grid.
   * @param[in] particles  Particles of the atomic block (absolute lattice coordinates).
   * @param[in] location   Absolute lattice position of node (0,0,0) of the block.
   * @param[in] nx,ny,nz   Number of lattice nodes of the block along each axis.
   *
   * Cost: O(nx*ny*nz + particles.size()). Buffers are reused between builds.
   */
  void build(const std::vector<HemoCellParticle> & particles, const plb::Dot3D & location,
             plint nx, plint ny, plint nz);

  /// Linear index of block-local node (x,y,z); the caller guarantees it is inside the grid.
  inline unsigned int index(plint x, plint y, plint z) const {
    return z + nz*(y + ny*x);
  }
  /// First particle index of node @p node (see class description).
  inline const unsigned int * begin(unsigned int node) const {
    return particleIndices.data() + nodeStart[node];
  }
  /// One past the last particle index of node @p node.
  inline const unsigned int * end(unsigned int node) const {
    return particleIndices.data() + nodeStart[node+1];
  }
  /// Number of particles binned in node @p node.
  inline unsigned int count(unsigned int node) const {
    return nodeStart[node+1] - nodeStart[node];
  }
  inline plint getNx() const { return nx; }
  inline plint getNy() const { return ny; }
  inline plint getNz() const { return nz; }

private:
  plint nx = 0, ny = 0, nz = 0;
  std::vector<unsigned int> nodeStart;        ///< nx*ny*nz+1 offsets into particleIndices
  std::vector<unsigned int> particleIndices;  ///< particle indices, grouped per node
  std::vector<unsigned int> nodeOfParticle;   ///< scratch: node of each particle (or sentinel)
};

/**
 * @brief Accumulate the pairwise inter-cell repulsion force into
 * `sv.force_repulsion` of all binned particles.
 *
 * For every pair of particles of *different* cells that lie in the same or in
 * neighbouring grid nodes and are closer than @p r_cutoff, a force of
 * magnitude `r_const * r_cutoff / distance` is applied along the connecting
 * line, with equal and opposite contributions to both particles (each pair is
 * visited exactly once). Pairs further apart than one grid node are not found,
 * so @p r_cutoff should not exceed one lattice spacing.
 *
 * Nodes on the upper x, y and z faces of the grid only act as neighbours of
 * interior nodes; pairs entirely within those outermost envelope layers are
 * skipped, as in the original implementation.
 *
 * The caller resets `force_repulsion` beforehand. Coincident particles
 * (distance 0) yield a non-finite force; this is not checked.
 *
 * @param[in,out] particles  Particle vector the grid was built from.
 * @param[in] grid           Grid built from @p particles.
 * @param[in] r_const        Force constant [lattice force units].
 * @param[in] r_cutoff       Cut-off distance [lattice units].
 */
void computePairRepulsion(std::vector<HemoCellParticle> & particles, const ParticleGrid & grid,
                          T r_const, T r_cutoff);

}
#endif // HEMO_PARTICLEGRID_H
