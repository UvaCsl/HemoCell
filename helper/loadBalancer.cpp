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
#include "loadBalancer.h"
LoadBalancer::LoadBalancer(HemoCell & hemocell_) : hemocell(hemocell_), original_block_structure(hemocell_.lattice->getSparseBlockStructure().clone()),original_thread_attribution(hemocell_.lattice->getMultiBlockManagement().getThreadAttribution().clone()) { 

}

#ifdef HEMO_PARMETIS
#include <parmetis.h>

void LoadBalancer::reloadCheckpoint() {
  //Firstly reload the config
  delete hemocell.documentXML;
  hemocell.documentXML = new XMLreader(global::directories().getOutputDir() + "checkpoint.xml");
  hemocell.cellfields->syncEnvelopes();
  hemocell.loadCheckPoint();
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
  //Error Checking
  std::vector<plint> const& blocks = hemocell.cellfields->immersedParticles->getLocalInfo().getBlocks();
  for (pluint iBlock=0; iBlock<blocks.size(); ++iBlock) {
    if (hemocell.cellfields->immersedParticles->getComponent(blocks[iBlock]).neighbours.size() > HEMOCELL_MAX_NEIGHBOURS) {
      cerr << "(HemoCell) ERROR: More neighbours of atomic block exist than allowed: " << hemocell.cellfields->immersedParticles->getComponent(blocks[iBlock]).neighbours.size() << " instead of " << HEMOCELL_MAX_NEIGHBOURS << "!. Check the atomic block structure!\n";
      exit(0);
    }
  }
    
  //set FLI_iscalled
  FLI_iscalled = true;
  int size = global::mpi().getSize();
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
  
  vector<double> times(size);
  vector<double> lsps(size);
  
  /*for (auto const & entry : gatherValues) {
    pcout << "Atomic block " << entry.first << " is on proc " << entry.second.mpi_proc << " spent " << entry.second.particle_time << " for the particle field, spent " << entry.second.fluid_time << " for the fluid field" << endl;
  }*/

  this->gatherValues = gatherValues;
  
  for (auto const & entry : gatherValues) {
    times[entry.second.mpi_proc] += entry.second.particle_time + entry.second.fluid_time;
    lsps[entry.second.mpi_proc] += entry.second.n_lsp;
  }
  
  // lsps.
  double sum = 0;
  double average = 0;
  double max = 0;
  double fli = 0;
  
  // this one is for time
  double sum2 = 0;
  double average2 = 0;
  double max2 = 0;
  double fli2 = 0;
  
  for (unsigned int i = 0 ; i < times.size(); i++){
      sum = sum + times[i];
      sum2 = sum2 + lsps[i];
  }
  
  average = sum / size ;
  max = *std::max_element(times.begin(),times.end());

  fli = (max/average)-1;
  
  average2 = sum2 / size ;
  max2 = *std::max_element(lsps.begin(),lsps.end());

  fli2 = (max2/average2)-1;
  
  pcout << "fli (time):  " << fli << " fli (lsp): "<< fli2 << std::endl;
  
  return fli2;
}

void LoadBalancer::doLoadBalance() {
  if(!FLI_iscalled) {
    pcerr << "You did not calculate the fractional load imbalance before trying to balance, this is required! Exiting...";
    exit(0);
  }
  hemocell.saveCheckPoint(); // Save Checkpoint

  if (original_block_stored) {
    plint envelopeWidth = hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth();
    plint refinementLevel = hemocell.lattice->getMultiBlockManagement().getRefinementLevel();

    Dynamics<double,DESCRIPTOR> * dynamics = hemocell.lattice->getBackgroundDynamics().clone();
    bool perX = hemocell.lattice->periodicity().get(0);
    bool perY = hemocell.lattice->periodicity().get(1);
    bool perZ = hemocell.lattice->periodicity().get(2);
    bool internalStat = hemocell.lattice->isInternalStatisticsOn();

    delete hemocell.lattice;

    MultiBlockLattice3D<double,DESCRIPTOR> * newlattice = new
                MultiBlockLattice3D<double,DESCRIPTOR>(MultiBlockManagement3D (
                *original_block_structure->clone(),
                original_thread_attribution->clone(),
                envelopeWidth,
                refinementLevel ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),                
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<double,DESCRIPTOR>(),
                dynamics );
    
    
  newlattice->periodicity().toggle(0, perX);
  newlattice->periodicity().toggle(1, perY);
  newlattice->periodicity().toggle(2, perZ);
  newlattice->toggleInternalStatistics(internalStat);

  hemocell.lattice = newlattice;
  hemocell.cellfields->lattice = newlattice;
  
  delete hemocell.cellfields->immersedParticles;
  hemocell.cellfields->createParticleField();
  
  reloadCheckpoint();
  
  pcout << "(LoadBalancer) Re-Calculating FLI of original block structure" << endl;
  calculateFractionalLoadImbalance();
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
  
  /*for (auto & pair : newProc) {
    pcout << "Atomic block " << pair.first << " is assigned to processor " << pair.second << endl;
  }*/

  
  pcout << "(LoadBalancer) Recreating Fluid field with new Distribution of Atomic Blocks" << endl;

  map<plint,plint> nTA; //Conversion is necessary for next function
  for (auto & pair : newProc) { 
      nTA[pair.first] = pair.second; 
  }
  ExplicitThreadAttribution* newThreadAttribution = new ExplicitThreadAttribution(nTA);
  delete original_thread_attribution;
  original_thread_attribution = newThreadAttribution->clone();
  
  plint envelopeWidth = hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth();
  plint refinementLevel = hemocell.lattice->getMultiBlockManagement().getRefinementLevel();

  Dynamics<double,DESCRIPTOR> * dynamics = hemocell.lattice->getBackgroundDynamics().clone();
  bool perX = hemocell.lattice->periodicity().get(0);
  bool perY = hemocell.lattice->periodicity().get(1);
  bool perZ = hemocell.lattice->periodicity().get(2);
  bool internalStat = hemocell.lattice->isInternalStatisticsOn();

  delete hemocell.lattice;

  MultiBlockLattice3D<double,DESCRIPTOR> * newlattice = new
            MultiBlockLattice3D<double,DESCRIPTOR>(MultiBlockManagement3D (
            *original_block_structure->clone(),
            newThreadAttribution,
            envelopeWidth,
            refinementLevel ),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),                
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double,DESCRIPTOR>(),
            dynamics );
  
  newlattice->periodicity().toggle(0, perX);
  newlattice->periodicity().toggle(1, perY);
  newlattice->periodicity().toggle(2, perZ);
  newlattice->toggleInternalStatistics(internalStat);
 
  hemocell.lattice = newlattice;
  hemocell.cellfields->lattice = newlattice;
  
  delete hemocell.cellfields->immersedParticles;
  hemocell.cellfields->createParticleField(original_block_structure->clone(),newThreadAttribution->clone());

  reloadCheckpoint();
  pcout << "(LoadBalancer) Continuing simulation with balanced application" << endl;

  return;
}

//Necessary C++ crap
LoadBalancer::GatherTimeOfAtomicBlocks * LoadBalancer::GatherTimeOfAtomicBlocks::clone() const { return new LoadBalancer::GatherTimeOfAtomicBlocks(*this); }

void LoadBalancer::restructureBlocks(bool checkpoint_available) {

  ThreadAttribution * oldThreads = original_thread_attribution->clone();
  SparseBlockStructure3D * old_structure = original_block_structure->clone();
  vector<plint> localBlocks = old_structure->getLocalBlocks(*oldThreads);
  
  vector<Box3D> boxes(localBlocks.size());
  for(unsigned int bid = 0 ; bid < localBlocks.size() ; bid ++) {
    old_structure->getBulk(localBlocks[bid],boxes[bid]);
  }
  
  //Simple algorithm, check in three directions (the other three are from neigh -> root)
  //If we find a matching face, extend root block, delete neighbour block, and restart until we find no merges
  restart:
  for (unsigned int root = 0 ; root < boxes.size() ; root++) {
    for (unsigned int neigh = 0 ; neigh < boxes.size() ; neigh++) {
      if (root == neigh) continue;
      if (boxes[root].x0 == boxes[neigh].x1 + 1&& 
          boxes[root].y0 == boxes[neigh].y0 && boxes[root].y1 == boxes[neigh].y1 &&
          boxes[root].z0 == boxes[neigh].z0 && boxes[root].z1 == boxes[neigh].z1)
      {
        boxes[root].x0 = boxes[neigh].x0;
        boxes[neigh] = boxes.back();
        boxes.pop_back();
        goto restart;
      }
      if (boxes[root].y0 == boxes[neigh].y1 + 1 && 
          boxes[root].x0 == boxes[neigh].x0 && boxes[root].x1 == boxes[neigh].x1 &&
          boxes[root].z0 == boxes[neigh].z0 && boxes[root].z1 == boxes[neigh].z1)
      {
        boxes[root].y0 = boxes[neigh].y0;
        boxes[neigh] = boxes.back();
        boxes.pop_back();
        goto restart;
      }
      if (boxes[root].z0 == boxes[neigh].z1 + 1 && 
          boxes[root].x0 == boxes[neigh].x0 && boxes[root].x1 == boxes[neigh].x1 &&
          boxes[root].y0 == boxes[neigh].y0 && boxes[root].y1 == boxes[neigh].y1)
      {
        boxes[root].z0 = boxes[neigh].z0;
        boxes[neigh] = boxes.back();
        boxes.pop_back();
        goto restart;
      }
    }
  }
  
  //Communicate new blocks to all procs
  map<int,Box3D_simple> newBlocks;
  for (unsigned int box = 0 ; box < boxes.size() ; box++) {
    newBlocks[localBlocks[box]] = boxes[box];
  }
  
  int numAtomicBlock = old_structure->getNumBlocks();
  HemoCellGatheringFunctional<Box3D_simple>::gather(newBlocks,numAtomicBlock);
  
  pcout << "(LoadBalancer) (Restructure) went from " << numAtomicBlock << " to " << newBlocks.size() << " atomic blocks (whole domain)" << endl;
  
  //The new structure we can fill, Initialize such because palabos
  SparseBlockStructure3D * new_structure = old_structure->clone();
  delete old_structure;
  for (const auto & pair : new_structure->getBulks()) {
    new_structure->removeBlock(pair.first); //Lol
  }
  
  map<plint,plint> nTA; //Conversion is necessary for new threadattribution
  plint blockId = 0;
  for (const auto & pair : newBlocks ) {
    nTA[blockId] = oldThreads->getMpiProcess(pair.first);
    new_structure->addBlock(pair.second.getBox3D(),blockId);
    blockId++;
  
  }
  ExplicitThreadAttribution* newThreadAttribution = new ExplicitThreadAttribution(nTA);        

  plint envelopeWidth = hemocell.lattice->getMultiBlockManagement().getEnvelopeWidth();
  plint refinementLevel = hemocell.lattice->getMultiBlockManagement().getRefinementLevel();
  Dynamics<double,DESCRIPTOR> * dynamics = hemocell.lattice->getBackgroundDynamics().clone();
  bool perX = hemocell.lattice->periodicity().get(0);
  bool perY = hemocell.lattice->periodicity().get(1);
  bool perZ = hemocell.lattice->periodicity().get(2);
  bool internalStat = hemocell.lattice->isInternalStatisticsOn();
  
  delete hemocell.lattice;
  
  MultiBlockLattice3D<double,DESCRIPTOR> * newlattice = new
                MultiBlockLattice3D<double,DESCRIPTOR>(MultiBlockManagement3D (
                *new_structure->clone(),
                newThreadAttribution,
                envelopeWidth,
                refinementLevel),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),                
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<double,DESCRIPTOR>(),
                dynamics );
  newlattice->periodicity().toggle(0, perX);
  newlattice->periodicity().toggle(1, perY);
  newlattice->periodicity().toggle(2, perZ);
  newlattice->toggleInternalStatistics(internalStat);
  
  hemocell.lattice = newlattice;
  hemocell.cellfields->lattice = newlattice;
  
  delete hemocell.cellfields->immersedParticles;
  hemocell.cellfields->createParticleField(new_structure->clone(),newThreadAttribution->clone());

  reloadCheckpoint();
  pcout << "(LoadBalancer) (Restructure) Continuing simulation with restructured application" << endl;

  delete newThreadAttribution;
  delete new_structure;
  
  return;
}
#endif
