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
#include "hemoCellStretch.h"

using namespace std;
using namespace plb;
namespace hemo {


HemoCellStretch::FindForcedLsps * HemoCellStretch::FindForcedLsps::clone() const {
  return new HemoCellStretch::FindForcedLsps(*this);
}

void HemoCellStretch::FindForcedLsps::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  vector<HemoCellParticle *> found;
  // Every block that holds the cell gives the same vertices: take the first.
  if (!stretch.lower_lsps.empty()) { return; }
  HemoCellParticleField * pf = dynamic_cast<HemoCellParticleField *>(blocks[0]);
  const map<plint, vector<int>> & ppc = pf->get_particles_per_cell();
  if (ppc.find(0) == ppc.end()) { return; }

  const vector<int> & p_indices = ppc.at(0);
  for (int p_index : p_indices) {
    if (p_index == -1) { fatal("(HemoCellStretch) The stretched cell is incomplete on this block."); }
    found.push_back(&pf->particles[p_index]);
  }
  //sort found on first dimension
  //Use simple sort, dont want to overload < operator of particle
  HemoCellParticle * tmp;
  for (unsigned int i = 0; i < found.size() - 1; i++) {
    for (unsigned int j = 1; j < found.size() - i; j++) {
      if (found[j - 1]->sv.position[0] > found[j]->sv.position[0]) {
        tmp = found[j - 1];
        found[j - 1] = found[j];
        found[j] = tmp;
      }
    }
  }
  for (unsigned int i = 0; i < stretch.n_forced_lsps; i++) {
    stretch.lower_lsps.push_back(found[i]->sv.vertexId);
    stretch.upper_lsps.push_back(found[found.size() - 1 - i]->sv.vertexId);
  }
}

HemoCellStretch::ForceForcedLsps * HemoCellStretch::ForceForcedLsps::clone() const {
  return new HemoCellStretch::ForceForcedLsps(*this);
}

void HemoCellStretch::ForceForcedLsps::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  const map<plint, std::vector<int>> & ppc = dynamic_cast<HemoCellParticleField *>(blocks[0])->get_particles_per_cell();
  vector<HemoCellParticle> * particles = &dynamic_cast<HemoCellParticleField *>(blocks[0])->particles;

  hemo::Array<T, 3> ex_force = {stretch.external_force * stretch.scale, 0., 0.};
  for (unsigned int vi : stretch.lower_lsps) {
    if (ppc.find(0) == ppc.end()) { continue; }
    if (ppc.at(0)[vi] < 0) { continue; }
    (*particles)[ppc.at(0)[vi]].sv.force -= ex_force;
  }
  for (unsigned int vi : stretch.upper_lsps) {
    if (ppc.find(0) == ppc.end()) { continue; }
    if (ppc.at(0)[vi] < 0) { continue; }
    (*particles)[ppc.at(0)[vi]].sv.force += ex_force;
  }
}

HemoCellStretch::HemoCellStretch(HemoCellField & cellfield_, unsigned int n_forced_lsps_, T external_force_)
    : cellfield(cellfield_), n_forced_lsps(n_forced_lsps_), external_force(external_force_ / n_forced_lsps_) {
  //Sanity Checks
  if (cellfield.cellFields.number_of_cells != 1) {
    fatal("(HemoCellStretch) Needs exactly 1 cell, there are ", cellfield.cellFields.number_of_cells, ".");
  }
  if (n_forced_lsps == 0) { fatal("(HemoCellStretch) The number of pulled vertices must be positive."); }
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(cellfield.getParticleField3D());
  cellfield.cellFields.syncEnvelopes();
  applyProcessingFunctional(new FindForcedLsps(*this), cellfield.getParticleField3D()->getBoundingBox(), wrapper);
}

void HemoCellStretch::applyForce() {
  if (cellfield.timescale != 1) {
    fatal("(HemoCellStretch) Needs a material update interval of 1, it is ", cellfield.timescale, ".");
  }
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(cellfield.getParticleField3D());
  applyProcessingFunctional(new ForceForcedLsps(*this), cellfield.getParticleField3D()->getBoundingBox(), wrapper);
}

}