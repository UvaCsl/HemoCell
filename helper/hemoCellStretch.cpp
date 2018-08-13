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
namespace hemo {
  

HemoCellStretch::FindForcedLsps * HemoCellStretch::FindForcedLsps::clone() const { return new HemoCellStretch::FindForcedLsps(*this);}

void HemoCellStretch::FindForcedLsps::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  vector<HemoCellParticle*> found;
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  const map<int,vector<int>> & ppc = pf->get_particles_per_cell();
  
  const vector<int> & p_indices = ppc.at(0);
  for (int p_index : p_indices) {
    if (p_index == -1) {
      cout << "Error -1 found in cell, exiting" << endl;
      exit(1);
    }
    found.push_back(&pf->particles[p_index]);
  }
  //sort found on first dimension
  //Use simple sort, dont want to overload < operator of particle
  HemoCellParticle * tmp;
  for (unsigned int i = 0 ; i <  found.size() - 1 ; i++) {
    for (unsigned int j = 1 ; j < found.size() - i ; j++) {
      if (found[j-1]->sv.position[0] > found[j]->sv.position[0]) {
        tmp = found[j-1];
        found[j-1] = found[j];
        found[j] = tmp;
      }
    }
  }
  for (unsigned int i = 0 ; i < n_forced_lsps; i ++) {
    lower_lsps.push_back(found[i]->sv.vertexId);
    upper_lsps.push_back(found[found.size()-1-i]->sv.vertexId);
  }
}

HemoCellStretch::ForceForcedLsps * HemoCellStretch::ForceForcedLsps::clone() const { return new HemoCellStretch::ForceForcedLsps(*this);}

void HemoCellStretch::ForceForcedLsps::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  const map<int,std::vector<int>> & ppc = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->get_particles_per_cell();
  vector<HemoCellParticle> * particles = &dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->particles;

  hemo::Array<T,3> ex_force = {external_force*scale,0.,0.};
  for (unsigned int vi : lower_lsps) {
    if (ppc.find(0) == ppc.end()) { continue; }
    if (ppc.at(0)[vi] < 0) { continue; }
    (*particles)[ppc.at(0)[vi]].sv.force -= ex_force;
  }
  for (unsigned int vi : upper_lsps) {
    if (ppc.find(0) == ppc.end()) { continue; }
    if (ppc.at(0)[vi] < 0) { continue; }
    (*particles)[ppc.at(0)[vi]].sv.force += ex_force;
  }
}

HemoCellStretch::HemoCellStretch(HemoCellField & cellfield_, unsigned int n_forced_lsps_, T external_force_) 
                : cellfield(cellfield_)
{
  //Sanity Checks
  if (cellfield.cellFields.number_of_cells != 1) {
    pcout << "(HemoCellStretch) Refusing to run with more or less than 1 cell" << endl;
    exit(1);
  }
  n_forced_lsps = n_forced_lsps_;
  external_force = external_force_/n_forced_lsps;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfield.getParticleField3D());
  cellfield.cellFields.syncEnvelopes();
  applyProcessingFunctional(new FindForcedLsps(),cellfield.getParticleField3D()->getBoundingBox(),wrapper);
}

vector<plint> HemoCellStretch::lower_lsps = vector<plint>();
vector<plint> HemoCellStretch::upper_lsps = vector<plint>();
unsigned int HemoCellStretch::n_forced_lsps = 0;
T HemoCellStretch::external_force = 0.0;
T HemoCellStretch::scale = 1.0;

void HemoCellStretch::applyForce() {
  if (cellfield.timescale != 1) {
    pcout << "Refusing to stretch with particle update timestep larger than 1" << endl;
    exit(1);
  }
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfield.getParticleField3D());
  applyProcessingFunctional(new ForceForcedLsps(),cellfield.getParticleField3D()->getBoundingBox(),wrapper);
}

}