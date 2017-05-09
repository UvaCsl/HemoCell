#ifndef LOADBALANCER_CPP
#define LOADBALANCER_CPP

#include "loadBalancer.h"

LoadBalancer::LoadBalancer(HemoCell & hemocell_) : hemocell(hemocell_) {

}

void LoadBalancer::GatherTimeOfAtomicBlocks::processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  BlockLattice3D<double,DESCRIPTOR> * ff = dynamic_cast<BlockLattice3D<double,DESCRIPTOR>*>(blocks[1]);

  gatherValues[pf->atomicBlockId].particle_time = pf->timer.getTime();
  gatherValues[pf->atomicBlockId].fluid_time = ff->timer.getTime();
  gatherValues[pf->atomicBlockId].mpi_proc = global::mpi().getRank();
  
  vector<HemoCellParticle *> found;
  pf->findParticles(pf->localDomain,found);
  gatherValues[pf->atomicBlockId].n_lsp = found.size();
  

 
  gatherValues[pf->atomicBlockId].n_neighbours = pf->neighbours.size();
  for (unsigned int i = 0 ; i < pf->neighbours.size(); i++) {
    gatherValues[pf->atomicBlockId].neighbours[i] = pf->neighbours[i];
  }
    
  pf->timer.reset();
  ff->timer.reset();
}

double LoadBalancer::calculateFractionalLoadImbalance() {
  //set FLI_iscalled
  FLI_iscalled = true;
  //We need the total number of atomic blocks for the Gather Functional
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell.cellfields->immersedParticles);
  int numAtomicBlock = hemocell.lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();

  //Do a more intricate functional (The above one could be part of this one, but
  //it is just for example
  wrapper.push_back(hemocell.lattice); //push back fluid field as well
  map<int,TOAB_t> gatherValues;
  applyProcessingFunctional(new GatherTimeOfAtomicBlocks(gatherValues),hemocell.cellfields->immersedParticles->getBoundingBox(), wrapper);
  HemoCellGatheringFunctional<TOAB_t>::gather(gatherValues,numAtomicBlock);

  for (auto const & entry : gatherValues) {
    pcout << "Atomic block " << entry.first << " is on proc " << entry.second.mpi_proc << " spent " << entry.second.particle_time << " for the particle field, spent " << entry.second.fluid_time << " for the fluid field" << endl;
  }

  this->gatherValues = gatherValues;

  return 0.;
}


void LoadBalancer::doLoadBalance() {
  if(!FLI_iscalled) {
    pcerr << "You did not calculate the fractional load imbalance before trying to balance, this is required! Exiting...";
    exit(0);
  }
  
  //Unfortunately, due to inconsequent atomic block naming, it is easier to start everything on the root process
  //Variable naming according to Parmetis Manual
  vector<idx_t> vtxdist(global::mpi().getSize()+1);
  for (unsigned int i = 1 ; i < vtxdist.size(); i++){
    vtxdist[i] = gatherValues.size(); 
  }
  
  if (global::mpi().getRank()) {
    goto parmetis; // Do not do anything if not ROOT
  }
  
    
  
  parmetis: //LABEL PARMETIS
  return;
 /*idx_t  
  int ret = ParMETIS_V3_PartGeomKway(
       idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
       idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
       idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
       idx_t *edgecut, idx_t *part, MPI_Comm *comm);*/
}

//Necessary C++ crap
LoadBalancer::GatherTimeOfAtomicBlocks * LoadBalancer::GatherTimeOfAtomicBlocks::clone() const { return new LoadBalancer::GatherTimeOfAtomicBlocks(*this); }
#endif
