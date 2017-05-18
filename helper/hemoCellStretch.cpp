#ifndef HEMOCELLSTRETCH_CPP
#define HEMOCELLSTRETCH_CPP

#include "hemoCellStretch.h"
HemoCellStretch::FindForcedLsps * HemoCellStretch::FindForcedLsps::clone() const { return new HemoCellStretch::FindForcedLsps(*this);}

void HemoCellStretch::FindForcedLsps::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  vector<HemoCellParticle*> found;
  dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->findParticles(domain,found);
  
  //sort found on first dimension
  //Use simple sort, dont want to overload < operator of particle
  HemoCellParticle * tmp;
  for (unsigned int i = 0 ; i <  found.size() - 1 ; i++) {
    for (unsigned int j = 1 ; j < found.size() - i ; j++) {
      if (found[j-1]->position[0] > found[j]->position[0]) {
        tmp = found[j-1];
        found[j-1] = found[j];
        found[j] = tmp;
      }
    }
  }
  for (unsigned int i = 0 ; i < n_forced_lsps; i ++) {
    lower_lsps.push_back(found[i]->vertexId);
    upper_lsps.push_back(found[found.size()-1-i]->vertexId);
  }
}

HemoCellStretch::ForceForcedLsps * HemoCellStretch::ForceForcedLsps::clone() const { return new HemoCellStretch::ForceForcedLsps(*this);}

void HemoCellStretch::ForceForcedLsps::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  map<int,std::vector<HemoCellParticle*>> * ppc = &dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0])->particles_per_cell;
  Array<double,3> ex_force = {external_force,0.,0.};
  for (unsigned int vi : lower_lsps) {
    (*ppc)[0][vi]->force -= ex_force;
  }
  for (unsigned int vi : upper_lsps) {
    (*ppc)[0][vi]->force += ex_force;
  }
}

HemoCellStretch::HemoCellStretch(HemoCellField & cellfield_, unsigned int n_forced_lsps_, double external_force_) 
                : cellfield(cellfield_)
{
  n_forced_lsps = n_forced_lsps_;
  external_force = external_force_/n_forced_lsps;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfield.getParticleField3D());
  applyProcessingFunctional(new FindForcedLsps(),cellfield.getParticleField3D()->getBoundingBox(),wrapper);
}

vector<plint> HemoCellStretch::lower_lsps = vector<plint>();
vector<plint> HemoCellStretch::upper_lsps = vector<plint>();
unsigned int HemoCellStretch::n_forced_lsps = 0;
double HemoCellStretch::external_force = 0.0;

void HemoCellStretch::applyForce() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfield.getParticleField3D());
  applyProcessingFunctional(new ForceForcedLsps(),cellfield.getParticleField3D()->getBoundingBox(),wrapper);
}

#endif
