#ifndef LOADBALANCER_CPP
#define LOADBALANCER_CPP

#include "loadBalancer.h"

LoadBalancer::LoadBalancer(HemoCell & hemocell_) : hemocell(hemocell_) {

}


void LoadBalancer::GatherNumberOfLSPs::processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  vector<HemoCellParticle *> found;
  pf->findParticles(pf->localDomain,found);
  gatherValues[pf->atomicBlockId] = found.size();
}

double LoadBalancer::calculateFractionalLoadImbalance() {
  //We need the total number of atomic blocks for the Gather Functional
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell.cellfields->immersedParticles);

  map<int,int> gatherValues;
  applyProcessingFunctional(new GatherNumberOfLSPs(gatherValues),hemocell.cellfields->immersedParticles->getBoundingBox(), wrapper);

  int numAtomicBlock = hemocell.lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
  HemoCellGatheringFunctional<int>::gather(gatherValues,numAtomicBlock);

  for (auto const & entry : gatherValues) {
    pcout << "Atomic block " << entry.first << " has " << entry.second << " lsps." << endl;
  }

  return 0.;
}


void LoadBalancer::doLoadBalance() {

}

//Necessary C++ crap
LoadBalancer::GatherNumberOfLSPs * LoadBalancer::GatherNumberOfLSPs::clone() const { return new LoadBalancer::GatherNumberOfLSPs(*this); }
#endif
