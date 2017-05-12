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
  
  gatherValues[pf->atomicBlockId].location[0] = pf->getLocation().x;
  gatherValues[pf->atomicBlockId].location[1] = pf->getLocation().y;
  gatherValues[pf->atomicBlockId].location[2] = pf->getLocation().z;

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
  
  //Variable naming according to Parmetis Manual
  idx_t wgtflag = 2;
  idx_t numflag = 0;
  idx_t ndims = 3;
  idx_t ncon = 1;
  idx_t nparts = global::mpi().getSize();
  idx_t options[3] = {1,PARMETIS_DBGLVL_TIME|PARMETIS_DBGLVL_INFO|PARMETIS_DBGLVL_PROGRESS,0};
  idx_t edgecut = 0;
  MPI_Comm mc = MPI_COMM_WORLD;
  unsigned int remainder = gatherValues.size()%nparts;
  unsigned int partsize = gatherValues.size()/nparts;
  unsigned int rank = global::mpi().getRank();
  
  vector<idx_t> vtxdist(nparts+1);
  for (unsigned int i = 1 ; i < vtxdist.size(); i++){
    vtxdist[i] = vtxdist[i-1] + partsize;
    if (i <= remainder) {
      vtxdist[i]++;
    }
  }
  
  unsigned int nv = vtxdist[rank+1] - vtxdist[rank];
  unsigned int ofs= vtxdist[rank];
  vector<idx_t> part(nv);

  
  vector<idx_t> xadj(nv+1);
  xadj[0] = 0;
  for (unsigned int i = 0 ; i + 1 < xadj.size() ; i ++) {
    xadj[i+1] = xadj[i] + gatherValues[ofs+i].n_neighbours;
  }

  vector<idx_t> adjncy(xadj.back());
  unsigned int entry = 0;
  for (unsigned int i = 0 ; i < nv ; i ++) {
    for (int j = 0 ; j < gatherValues[ofs+i].n_neighbours; j++) {
      adjncy[entry] = gatherValues[ofs+i].neighbours[j];
      entry++;
    }
  }

  vector<real_t> xyz(nv*ndims);
  for (unsigned int i = 0 ; i < xyz.size()/ndims ; i++) {
    for (int j = 0 ; j < ndims ; j++ ) {
      xyz[i*ndims + j] = gatherValues[ofs+i].location[j];
    }
  }

  vector<idx_t> vwgt(nv);
  for (unsigned int i = 0 ; i < vwgt.size() ; i++) {
    vwgt[i] = gatherValues[ofs+i].n_lsp;
  }
  
  vector<real_t> tpwghts(ncon*nparts,1.0/(nparts*ncon));
  vector<real_t> ubvec(ncon,1.05);

  ParMETIS_V3_PartGeomKway(&vtxdist[0], &xadj[0], &adjncy[0], &vwgt[0], NULL, &wgtflag, &numflag,  &ndims, &xyz[0], 
                           &ncon, &nparts, &tpwghts[0], &ubvec[0], &options[0],&edgecut, &part[0], &mc);
  
  map<int,plint> newProc; //Gather the results to all mpi processes, we can use the gathering functional for that as well!
  for (unsigned int i = 0 ; i < part.size() ; i++){
    newProc[ofs+i] = part[i];
  }
  int numAtomicBlock = hemocell.lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
  HemoCellGatheringFunctional<plint>::gather(newProc,numAtomicBlock);
  
  for (auto & pair : newProc) {
    pcout << "Atomic block " << pair.first << " is assigned to processor " << pair.second << endl;
  }
  
  hemocell.saveCheckPoint(); // Save Checkpoint
  pcout << "(LoadBalancer) Recreating Fluid field with new Distribution of Atomic Blocks" << endl;

  map<plint,plint> nTA; //Conversion is necessary for next function
  for (auto & pair : newProc) { nTA[pair.first] = pair.second; }
  ThreadAttribution* newThreadAttribution = new ExplicitThreadAttribution(nTA);        
  MultiBlockLattice3D<double,DESCRIPTOR> * newlattice = new
                MultiBlockLattice3D<double,DESCRIPTOR>(MultiBlockManagement3D (
                hemocell.lattice->getSparseBlockStructure(),
                newThreadAttribution,
                //hemocell.lattice->getMultiBlockManagement().getThreadAttribution().clone(),
                hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth(),
                hemocell.lattice->getMultiBlockManagement().getRefinementLevel() ),
                hemocell.lattice->getBlockCommunicator().clone(),
                hemocell.lattice->getCombinedStatistics().clone(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<double,DESCRIPTOR>(),
                hemocell.lattice->getBackgroundDynamics().clone() );
  newlattice->periodicity().toggle(0, hemocell.lattice->periodicity().get(0));
  newlattice->periodicity().toggle(1, hemocell.lattice->periodicity().get(1));
  newlattice->periodicity().toggle(2, hemocell.lattice->periodicity().get(2));
  newlattice->toggleInternalStatistics(hemocell.lattice->isInternalStatisticsOn());

  delete hemocell.lattice;
  hemocell.lattice = newlattice;
  hemocell.cellfields->lattice = newlattice;
  
  delete hemocell.cellfields->immersedParticles;
  hemocell.cellfields->createParticleField();
  
  hemocell.loadCheckPoint();
  pcout << "(LoadBalancer) Continuing simulation with rebalanced application" << endl;

  return;
}

//Necessary C++ crap
LoadBalancer::GatherTimeOfAtomicBlocks * LoadBalancer::GatherTimeOfAtomicBlocks::clone() const { return new LoadBalancer::GatherTimeOfAtomicBlocks(*this); }
#endif
