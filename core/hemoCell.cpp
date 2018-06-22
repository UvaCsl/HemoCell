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
#include "hemocell.h"
#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>

volatile sig_atomic_t interrupted = 0;
void set_interrupt(int signum) {
  interrupted = 1;
}

HemoCell::HemoCell(char * configFileName, int argc, char * argv[])
{
  plbInit(&(argc),&(argv));

  pcout << "(HemoCell) (Config) reading " << configFileName << endl;
  cfg = new Config(configFileName);
  documentXML = new XMLreader(configFileName);
  loadDirectories(configFileName, cfg);
  loadGlobalConfigValues(cfg);
  printHeader();
  
  //Start statistics
  statistics.start();
  
  // start clock for basic performance feedback
  lastOutputAt = 0;
  global::timer("atOutput").start();
  
#ifdef FORCE_LIMIT
    hlogfile << "(HemoCell) WARNING: Force limit active at " << FORCE_LIMIT << " pN. Results can be inaccurate due to force capping." << endl;
#endif
  if (sizeof(T) == sizeof(float)) {
    hlog << "(HemoCell) WARNING: Running with single precision, you might want to switch to double precision" << endl;
  }
  
  ///Set signal handlers to exit gracefully on many signals
  struct sigaction sa;
  memset(&sa, 0, sizeof(struct sigaction));
  sa.sa_handler = set_interrupt;
  sigaction(SIGINT,&sa,0);
  sigaction(SIGTERM,&sa,0);
  sigaction(SIGHUP,&sa,0);
  sigaction(SIGQUIT,&sa,0);
  sigaction(SIGABRT,&sa,0);
  sigaction(SIGUSR1,&sa,0);
  sigaction(SIGUSR2,&sa,0);
 
  
}

void HemoCell::latticeEquilibrium(T rho, hemo::Array<T, 3> vel) {
  hlog << "(HemoCell) (Fluid) Setting Fluid Equilibrium" << endl;
  plb::Array<T,3> vel_plb = {vel[0],vel[1],vel[2]};
  initializeAtEquilibrium(*lattice, (*lattice).getBoundingBox(), rho, vel_plb);
}

void HemoCell::initializeCellfield() {
  cellfields = new HemoCellFields(*lattice,(*cfg)["domain"]["particleEnvelope"].read<int>(),*this);

  //Correct place for init
  loadBalancer = new LoadBalancer(*this);
}

void HemoCell::setOutputs(string name, vector<int> outputs) {
  hlog << "(HemoCell) (CellField) Setting output variables for " << name << " cells" << endl;
  vector<int> outputs_c = outputs;
  (*cellfields)[name]->setOutputVariables(outputs_c);
}

void HemoCell::setFluidOutputs(vector<int> outputs) {
  hlog << "(HemoCell) (Fluid) Setting output variables for fluid field" << endl;
  vector<int> outputs_c = outputs;
  cellfields->desiredFluidOutputVariables = outputs_c;
}

void HemoCell::setSystemPeriodicity(unsigned int axis, bool bePeriodic) {
  if (lattice == 0) {
    pcerr << "(HemoCell) (Periodicity) please create a lattice before trying to set the periodicity" << endl;
    exit(1);    
  }
  if (cellfields->immersedParticles == 0) {
    pcerr << "(HemoCell) (Periodicity) please create a particlefield (hemocell.initializeCellfields()) before trying to set the periodicity" << endl;
    exit(1);   
  }
  lattice->periodicity().toggle(axis,bePeriodic);
  cellfields->immersedParticles->periodicity().toggle(axis, bePeriodic);
  cellfields->InitAfterLoadCheckpoint();
}

void HemoCell::setSystemPeriodicityLimit(unsigned int axis, int limit) {
  hlog << "(HemoCell) (Periodicity) Setting periodicity limit of axis " << axis << " to " << limit << endl;
  cellfields->periodicity_limit[axis] = limit;
  
  //recalculate offsets :
  cellfields->periodicity_limit_offset_y = cellfields->periodicity_limit[0];
  cellfields->periodicity_limit_offset_z = cellfields->periodicity_limit[0]*cellfields->periodicity_limit[1];
}

void HemoCell::loadParticles() {
  hlog << "(HemoCell) (CellField) Loading particle positions "  << endl;
  loadParticlesIsCalled = true;
  readPositionsBloodCellField3D(*cellfields, param::dx, *cfg);
  cellfields->syncEnvelopes();
  cellfields->deleteIncompleteCells(false);
}

void HemoCell::loadCheckPoint() {
  pcout << "(HemoCell) (Saving Functions) Loading Checkpoint"  << endl;
  cellfields->load(documentXML, iter, cfg);
}

void HemoCell::saveCheckPoint() {
  pcout << "(HemoCell) (Saving Functions) Saving Checkpoint at timestep " << iter << endl;
  cellfields->save(documentXML, iter, cfg);
}

void HemoCell::writeOutput() {
  statistics["output"].start();
  // Very naive performance approximation

  pcout << "(HemoCell) (Output) writing output at timestep " << iter << " (" << param::dt * iter<< " s). Approx. performance: " << to_string(stof(statistics.elapsed_string())/(iter-lastOutputAt)) << " s / iteration." << endl;
  lastOutputAt = iter;
  if(repulsionEnabled) {
    statistics["output"]["repulsionForce"].start();
    cellfields->applyRepulsionForce();
    statistics["output"]["repulsionForce"].stop();
  }
  if(boundaryRepulsionEnabled) {
    statistics["output"]["boundaryRepulsionForce"].start();
    cellfields->applyBoundaryRepulsionForce();
    statistics["output"]["boundaryRepulsionForce"].stop();
  }
  statistics["output"]["syncEnvelopes"].start();
  cellfields->syncEnvelopes();
  statistics["output"]["syncEnvelopes"].stop();
  statistics["output"]["deleteCells"].start();
  if (globalConfigValues.cellsDeletedInfo) {
    cellfields->deleteIncompleteCells(true);
  } else {
    cellfields->deleteIncompleteCells(false);   
  }
  statistics["output"]["deleteCells"].stop();

  //Repoint surfaceparticle forces for output
  statistics["output"]["separateForceVectors"].start();
  cellfields->separate_force_vectors();
  statistics["output"]["separateForceVectors"].stop();

  //Recalculate the forces
  statistics["output"]["recalculateConstitutiveModel"].start();
  cellfields->applyConstitutiveModel(true);
  statistics["output"]["recalculateConstitutiveModel"].stop();

  // Creating a new directory per save
  if (global::mpi().isMainProcessor()) {
    string folder = global::directories().getOutputDir() + "/hdf5/" + to_string(iter) ;
    mkpath(folder.c_str(), 0777);
    folder = global::directories().getOutputDir() + "/csv/" + to_string(iter) ;
    mkpath(folder.c_str(), 0777);
  }
  global::mpi().barrier();

  //Write Output
  statistics["output"]["writeOutput"].start();
  statistics["output"]["writeOutput"]["writeCellField"].start();
  writeCellField3D_HDF5(*cellfields,param::dx,param::dt,iter);
  statistics["output"]["writeOutput"]["writeCellField"].stop();
  statistics["output"]["writeOutput"]["writeFluidField"].start();
  writeFluidField_HDF5(*cellfields,param::dx,param::dt,iter);
  statistics["output"]["writeOutput"]["writeFluidField"].stop();
  statistics["output"]["writeOutput"]["writeCellCSVInfo"].start();
  writeCellInfo_CSV(this);
  statistics["output"]["writeOutput"]["writeCellCSVInfo"].stop();
  statistics["output"]["writeOutput"].stop();

  //Repoint surfaceparticle forces for speed
  statistics["output"]["unifyForceVectors"].start();
  cellfields->unify_force_vectors();
  statistics["output"]["unifyForceVectors"].stop();
  statistics["output"]["syncEnvelopes"].start();
  cellfields->syncEnvelopes();
  statistics["output"]["syncEnvelopes"].stop();
  // Continue with performance measurement
  statistics["output"].stop();
}

void HemoCell::checkExitSignals() {
  if (interrupted == 1) {
    cout << endl << "Caught Signal, saving work and quitting!" << endl << std::flush;
    exit(1);
  }
}

void HemoCell::iterate() {
  checkExitSignals();
  statistics["iterate"].start();
  // ### 1 ### Particle Force to Fluid
  if(repulsionEnabled && iter % cellfields->repulsionTimescale == 0) {
    statistics["iterate"]["repulsionForce"].start();
    cellfields->applyRepulsionForce();
    statistics["iterate"]["repulsionForce"].stop();
  }
  if(boundaryRepulsionEnabled && iter % cellfields->boundaryRepulsionTimescale == 0) {
    statistics["iterate"]["boundaryRepulsionForce"].start();
    cellfields->applyBoundaryRepulsionForce();
    statistics["iterate"]["boundaryRepulsionForce"].stop();
  }
  statistics["iterate"]["spreadParticleForce"].start();
  cellfields->spreadParticleForce();
  statistics["iterate"]["spreadParticleForce"].stop();

  // #### 2 #### LBM
  statistics["iterate"]["collideAndStream"].start();
  lattice->collideAndStream();
  statistics["iterate"]["collideAndStream"].stop();

  if(iter %cellfields->particleVelocityUpdateTimescale == 0) {
    // #### 3 #### IBM interpolation
    statistics["iterate"]["interpolateFluidVelocity"].start();
    cellfields->interpolateFluidVelocity();
    statistics["iterate"]["interpolateFluidVelocity"].stop();

    // ### 4 ### sync the particles
    statistics["iterate"]["syncEnvelopes"].start();
    cellfields->syncEnvelopes();
    statistics["iterate"]["syncEnvelopes"].stop();
  }
  
  // ### 5 ###
  statistics["iterate"]["advanceParticles"].start();
  cellfields->advanceParticles();
  statistics["iterate"]["advanceParticles"].stop();

  // ### 6 ###
  statistics["iterate"]["applyConstitutiveModel"].start();
  cellfields->applyConstitutiveModel();    // Calculate Force on Vertices
  statistics["iterate"]["applyConstitutiveModel"].stop();
    
  //We can safely delete non-local cells here, assuming model timestep is divisible by velocity timestep
  if(iter % cellfields->particleVelocityUpdateTimescale == 0) {
    if (globalConfigValues.cellsDeletedInfo) {
      statistics["iterate"]["deleteIncompleteCells"].start();
      cellfields->deleteIncompleteCells(true);
      statistics["iterate"]["deleteIncompleteCells"].stop();
    }
    statistics["iterate"]["deleteNonLocalParticles"].start();
    cellfields->deleteNonLocalParticles(3);
    statistics["iterate"]["deleteNonLocalParticles"].stop();
  }

  statistics["iterate"]["setExternalVector"].start();
  // Reset Forces on the lattice, TODO do own efficient implementation
  setExternalVector(*lattice, (*lattice).getBoundingBox(),
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
  statistics["iterate"]["setExternalVector"].stop();
  iter++;
  statistics["iterate"].stop();
}

T HemoCell::calculateFractionalLoadImbalance() {
  hlog << "(HemoCell) (LoadBalancer) Calculating Fractional Load Imbalance at timestep " << iter << endl;
  return loadBalancer->calculateFractionalLoadImbalance();
}

void HemoCell::setMaterialTimeScaleSeparation(string name, unsigned int separation){
  hlog << "(HemoCell) (Timescale Seperation) Setting seperation of " << name << " to " << separation << " timesteps"<<endl;
  //pcout << "(HemoCell) WARNING if the timescale separation is not dividable by tmeasure, checkpointing is non-deterministic!"<<endl; //not true anymore, with checkpointing remaining force is saved
  (*cellfields)[name]->timescale = separation;
  if (separation%cellfields->particleVelocityUpdateTimescale!=0) {
     pcout << "(HemoCell) Error, Velocity timescale separation cannot divide this material timescale separation, exiting ..." <<endl;
     exit(1);
  }
}

void HemoCell::setParticleVelocityUpdateTimeScaleSeparation(unsigned int separation) {
  hlog << "(HemoCell) (Timescale separation) Setting update separation of all particles to " << separation << " timesteps" << endl;
  hlogfile << "(HemoCell) WARNING this introduces great errors" << endl;

  for (unsigned int i = 0; i < cellfields->size() ; i++) {
    if ((*cellfields)[i]->timescale%separation !=0) {
      pcout << "(HemoCell) Error, Velocity timescale separation cannot divide all material timescale separations, exiting ..." <<endl;
      exit(1);
    }
  }
  cellfields->particleVelocityUpdateTimescale = separation;
}

void HemoCell::setRepulsionTimeScaleSeperation(unsigned int separation){
  hlog << "(HemoCell) (Repulsion Timescale Seperation) Setting seperation to " << separation << " timesteps"<<endl;
  cellfields->repulsionTimescale = separation;
  if (separation%cellfields->particleVelocityUpdateTimescale!=0) {
     pcout << "(HemoCell) Error, Velocity timescale separation cannot divide this repulsion timescale separation, exiting ..." <<endl;
     exit(1);
  }
}

void HemoCell::setMinimumDistanceFromSolid(string name, T distance) {
  hlog << "(HemoCell) (Set Distance) Setting minimum distance from solid to " << distance << " micrometer for " << name << endl; 
  if (loadParticlesIsCalled) {
    pcout << "(HemoCell) (Set Distance) WARNING: this function is called after the particles are loaded, so it probably has no effect" << endl;
  }
  (*cellfields)[name]->minimumDistanceFromSolid = distance;
}

void HemoCell::setRepulsion(T repulsionConstant, T repulsionCutoff) {
  hlog << "(HemoCell) (Repulsion) Setting repulsion constant to " << repulsionConstant << ". repulsionCutoff to" << repulsionCutoff << " Âµm" << endl;
  hlogfile << "(HemoCell) (Repulsion) Enabling repulsion" << endl;
  cellfields->repulsionConstant = repulsionConstant;
  cellfields->repulsionCutoff = repulsionCutoff*(1e-6/param::dx);
  repulsionEnabled = true;
}

void HemoCell::enableBoundaryParticles(T boundaryRepulsionConstant, T boundaryRepulsionCutoff, unsigned int timestep) {
  cellfields->populateBoundaryParticles();
  hlog << "(HemoCell) (Repulsion) Setting boundary repulsion constant to " << boundaryRepulsionConstant << ". boundary repulsionCutoff to" << boundaryRepulsionCutoff << " Âµm" << endl;
  hlogfile << "(HemoCell) (Repulsion) Enabling boundary repulsion" << endl;
  if (timestep%cellfields->particleVelocityUpdateTimescale!=0) {
     pcout << "(HemoCell) Error, Velocity timescale separation cannot divide this repulsion timescale separation, exiting ..." <<endl;
     exit(1);
  }
  cellfields->boundaryRepulsionConstant = boundaryRepulsionConstant;
  cellfields->boundaryRepulsionCutoff = boundaryRepulsionCutoff*(1e-6/param::dx);
  cellfields->boundaryRepulsionTimescale = timestep;
  boundaryRepulsionEnabled = true;
}



#ifdef HEMO_PARMETIS
void HemoCell::doLoadBalance() {
	pcout << "(HemoCell) (LoadBalancer) Balancing Atomic Block over mpi processes" << endl;
  loadBalancer->doLoadBalance();
}
#endif

void HemoCell::doRestructure(bool checkpoint_avail) {
  hlog << "(HemoCell) (LoadBalancer) Restructuring Atomic Blocks on processors" << endl;
  loadBalancer->restructureBlocks(checkpoint_avail);
}
