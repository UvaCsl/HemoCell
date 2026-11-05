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
#include "particleGrid.h"

#include <cmath>

namespace hemo {

void ParticleGrid::build(const std::vector<HemoCellParticle> & particles, const plb::Dot3D & location, plint nx_,
                         plint ny_, plint nz_) {
  nx = nx_;
  ny = ny_;
  nz = nz_;
  const unsigned int nNodes = nx * ny * nz;
  const unsigned int nParticles = particles.size();
  const unsigned int notBinned = nNodes; // sentinel node index

  nodeStart.assign(nNodes + 1, 0);
  nodeOfParticle.resize(nParticles);

  // Pass 1: find the node of every particle and count particles per node
  // (the count of node n is stored in nodeStart[n]).
  for (unsigned int i = 0; i < nParticles; i++) {
    const hemo::Array<T, 3> & pos = particles[i].sv.position;
    const plint x = nearestNode(pos[0] - location.x);
    const plint y = nearestNode(pos[1] - location.y);
    const plint z = nearestNode(pos[2] - location.z);
    if ((x >= 0) && (x < nx) && (y >= 0) && (y < ny) && (z >= 0) && (z < nz)) {
      const unsigned int node = index(x, y, z);
      nodeOfParticle[i] = node;
      nodeStart[node]++;
    } else {
      nodeOfParticle[i] = notBinned;
    }
  }

  // Inclusive prefix sum: nodeStart[n] becomes the end offset of node n.
  for (unsigned int n = 1; n < nNodes; n++) { nodeStart[n] += nodeStart[n - 1]; }
  nodeStart[nNodes] = (nNodes > 0) ? nodeStart[nNodes - 1] : 0;
  particleIndices.resize(nodeStart[nNodes]);

  // Pass 2: scatter the particle indices. Walking backwards and decrementing
  // the end offsets leaves nodeStart[n] at the start offset of node n and keeps
  // the indices within a node in ascending order.
  for (unsigned int i = nParticles; i-- > 0;) {
    const unsigned int node = nodeOfParticle[i];
    if (node != notBinned) { particleIndices[--nodeStart[node]] = i; }
  }
}

namespace {

inline void repulseParticlePair(HemoCellParticle & lParticle, HemoCellParticle & nParticle, const T r_const,
                                const T r_cutoff) {
  if (lParticle.sv.cellId == nParticle.sv.cellId) { return; }
  const hemo::Array<T, 3> dv = lParticle.sv.position - nParticle.sv.position;
  const T distance = std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
  // Coincident particles have no direction to repel along (and 1/distance would
  // be infinite); they are left to the other forces.
  if (distance < r_cutoff && distance > 0) {
    const hemo::Array<T, 3> rfm = r_const * (1 / (distance / r_cutoff)) * (dv / distance);
    lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;
    nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - rfm;
  }
}

// All pairs between two *different* nodes.
inline void repulseNodes(std::vector<HemoCellParticle> & particles, const ParticleGrid & grid,
                         const unsigned int l_index, const unsigned int n_index, const T r_const, const T r_cutoff) {
  const unsigned int * const l_end = grid.end(l_index);
  const unsigned int * const n_begin = grid.begin(n_index);
  const unsigned int * const n_end = grid.end(n_index);
  for (const unsigned int * i = grid.begin(l_index); i != l_end; ++i) {
    HemoCellParticle & lParticle = particles[*i];
    for (const unsigned int * j = n_begin; j != n_end; ++j) {
      repulseParticlePair(lParticle, particles[*j], r_const, r_cutoff);
    }
  }
}

// All unordered pairs within one node, each visited once.
inline void repulseWithinNode(std::vector<HemoCellParticle> & particles, const ParticleGrid & grid,
                              const unsigned int index, const T r_const, const T r_cutoff) {
  const unsigned int * const end = grid.end(index);
  for (const unsigned int * i = grid.begin(index); i != end; ++i) {
    HemoCellParticle & lParticle = particles[*i];
    for (const unsigned int * j = i + 1; j != end; ++j) {
      repulseParticlePair(lParticle, particles[*j], r_const, r_cutoff);
    }
  }
}

} // namespace

void computePairRepulsion(std::vector<HemoCellParticle> & particles, const ParticleGrid & grid, const T r_const,
                          const T r_cutoff) {
  const plint nx = grid.getNx();
  const plint ny = grid.getNy();
  const plint nz = grid.getNz();

  // Half-shell neighbour search: each node is paired with itself and with 13
  // of its 26 neighbours, so every pair of nodes is visited once.
  for (plint x = 0; x < nx - 1; x++) {
    for (plint y = 0; y < ny - 1; y++) {
      for (plint z = 0; z < nz - 1; z++) {
        const unsigned int l_index = grid.index(x, y, z);
        if (grid.count(l_index) == 0) { continue; } // no pairs start here

        // (0,0,0)
        repulseWithinNode(particles, grid, l_index, r_const, r_cutoff);
        // (0,0,1) (0,1,0) (0,1,1)
        repulseNodes(particles, grid, l_index, grid.index(x, y, z + 1), r_const, r_cutoff);
        repulseNodes(particles, grid, l_index, grid.index(x, y + 1, z), r_const, r_cutoff);
        repulseNodes(particles, grid, l_index, grid.index(x, y + 1, z + 1), r_const, r_cutoff);
        // (1,-1..1,-1..1)
        for (plint yy = y - 1; yy <= y + 1; yy++) {
          if (yy < 0) { continue; }
          for (plint zz = z - 1; zz <= z + 1; zz++) {
            if (zz < 0) { continue; }
            repulseNodes(particles, grid, l_index, grid.index(x + 1, yy, zz), r_const, r_cutoff);
          }
        }
        // (0,1,-1)
        if (z > 0) { repulseNodes(particles, grid, l_index, grid.index(x, y + 1, z - 1), r_const, r_cutoff); }
      }
    }
  }
}

}
