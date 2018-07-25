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

#include "readPositionsBloodCells.h"
#include "hemoCellFunctional.h"
#include "hemoCellParticle.h"
#include "hemoCellField.h"
#include "ParticleHdf5IO.h"
#include "FluidHdf5IO.h"
#include "writeCellInfoCSV.h"
#include "genericFunctions.h"

#include "core/plbInit.hh"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "dataProcessors/dataInitializerFunctional3D.hh"
#include "dataProcessors/dataInitializerWrapper3D.hh"

using namespace hemo;

volatile sig_atomic_t interrupted = 0;
void set_interrupt(int signum) {
  interrupted = 1;
}

HemoCell::HemoCell(char * configFileName, int argc, char * argv[])
{
  plb::plbInit(&(argc),&(argv));

  if (global.hemoCellInitialized) {
    pcout << "(HemoCell) (Error) Hemocell object already created, refusing to construct another one" << endl;
    exit(1);
  }
  global.hemoCellInitialized = true;
  
  pcout << "(HemoCell) (Config) reading " << configFileName << endl;
  cfg = new Config(configFileName);
  documentXML = new XMLreader(configFileName);
  loadDirectories(configFileName, cfg);
  loadGlobalConfigValues(cfg);
  printHeader();
  
  //Start statistics
  global.statistics.start();
  
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
  plb::initializeAtEquilibrium(*lattice, (*lattice).getBoundingBox(), rho, vel_plb);
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

void HemoCell::setCEPACOutputs(vector<int> outputs) {
  hlog << "(HemoCell) (Fluid) Setting CEPAC output variables for fluid field" << endl;
  vector<int> outputs_c = outputs;
  cellfields->desiredCEPACfieldOutputVariables = outputs_c;
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
  global.statistics["output"].start();
  std::string tpi = ((iter != lastOutputAt) ? Profiler::toString((global.statistics.elapsed()-lastOutput)/(iter-lastOutputAt)):"0.00");
  
  lastOutput = global.statistics.elapsed();
  lastOutputAt  = iter;
  // Very naive performance approximation
  pcout << "(HemoCell) (Output) writing output at timestep " << iter << " (" << param::dt * iter<< " s). Approx. performance: " << tpi << " s / iteration." << endl;
  if(repulsionEnabled) {
    cellfields->applyRepulsionForce();
  }
  if(boundaryRepulsionEnabled) {
    cellfields->applyBoundaryRepulsionForce();
  }
  cellfields->syncEnvelopes();
  if (global.cellsDeletedInfo) {
    cellfields->deleteIncompleteCells(true);
  } else {
    cellfields->deleteIncompleteCells(false);   
  }

  //Repoint surfaceparticle forces for output
  cellfields->separate_force_vectors();

  //Recalculate the forces
  cellfields->applyConstitutiveModel(true);

  // Creating a new directory per save
  if (global::mpi().isMainProcessor()) {
    string folder = global::directories().getOutputDir() + "/hdf5/" + to_string(iter) ;
    mkpath(folder.c_str(), 0777);
    folder = global::directories().getOutputDir() + "/csv/" + to_string(iter) ;
    mkpath(folder.c_str(), 0777);
  }
  global::mpi().barrier();

  //Write Output
  global.statistics.getCurrent()["writeOutput"].start();
  writeCellField3D_HDF5(*cellfields,param::dx,param::dt,iter);
  writeFluidField_HDF5(*cellfields,param::dx,param::dt,iter);
  if (global.enableCEPACfield) {
    writeCEPACField_HDF5(*cellfields,param::dx,param::dt,iter);
  }
  writeCellInfo_CSV(this);
  global.statistics.getCurrent().stop();

  //Repoint surfaceparticle forces for speed
  cellfields->unify_force_vectors();
  cellfields->syncEnvelopes();
  // Continue with performance measurement
  global.statistics.getCurrent().stop();
}

void HemoCell::checkExitSignals() {
  if (interrupted == 1) {
    cout << endl << "Caught Signal, saving work and quitting!" << endl << std::flush;
    exit(1);
  }
}

void HemoCell::iterate() {
  checkExitSignals();
  if (!sanityCheckDone) {
    sanityCheck();
  }
  global.statistics.getCurrent()["iterate"].start();
  // ### 1 ### Particle Force to Fluid
  if(repulsionEnabled && iter % cellfields->repulsionTimescale == 0) {
    cellfields->applyRepulsionForce();
  }
  if(boundaryRepulsionEnabled && iter % cellfields->boundaryRepulsionTimescale == 0) {
    cellfields->applyBoundaryRepulsionForce();
  }
  cellfields->spreadParticleForce();

  // #### 2 #### LBM
  global.statistics.getCurrent()["collideAndStream"].start();
  lattice->collideAndStream();
  global.statistics.getCurrent().stop();
  if (global.enableCEPACfield) {
    global.statistics.getCurrent()["CEPACcollideAndStream"].start();
    cellfields->CEPACfield->collideAndStream();
    global.statistics.getCurrent().stop();
  }

  if(iter %cellfields->particleVelocityUpdateTimescale == 0) {
    // #### 3 #### IBM interpolation
    cellfields->interpolateFluidVelocity();

    // ### 4 ### sync the particles
    cellfields->syncEnvelopes();
  }
  
  if(global.enableSolidifyMechanics && iter%cellfields->solidifyTimescale) {
    cellfields->solidifyCells();
  }
  // ### 5 ###
  cellfields->advanceParticles();

  // ### 6 ###
  cellfields->applyConstitutiveModel();    // Calculate Force on Vertices
    
  //We can safely delete non-local cells here, assuming model timestep is divisible by velocity timestep
  if(iter % cellfields->particleVelocityUpdateTimescale == 0) {
    if (global.cellsDeletedInfo) {
      cellfields->deleteIncompleteCells(true);
    }
    cellfields->deleteNonLocalParticles(3);
  }

  global.statistics.getCurrent()["setExternalVector"].start();
  // Reset Forces on the lattice, TODO do own efficient implementation
  setExternalVector(*lattice, (*lattice).getBoundingBox(),
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
  global.statistics.getCurrent().stop();
  iter++;
  global.statistics.getCurrent().stop();
}

T HemoCell::calculateFractionalLoadImbalance() {
  hlog << "(HemoCell) (LoadBalancer) Calculating Fractional Load Imbalance at timestep " << iter << endl;
  return loadBalancer->calculateFractionalLoadImbalance();
}

void HemoCell::setMaterialTimeScaleSeparation(string name, unsigned int separation){
  hlog << "(HemoCell) (Timescale Seperation) Setting seperation of " << name << " to " << separation << " timesteps"<<endl;
  (*cellfields)[name]->timescale = separation;
}

void HemoCell::setParticleVelocityUpdateTimeScaleSeparation(unsigned int separation) {
  hlog << "(HemoCell) (Timescale separation) Setting update separation of all particles to " << separation << " timesteps" << endl;
  hlogfile << "(HemoCell) WARNING this introduces great errors" << endl;
  cellfields->particleVelocityUpdateTimescale = separation;
}

void HemoCell::setRepulsionTimeScaleSeperation(unsigned int separation){
  hlog << "(HemoCell) (Repulsion Timescale Seperation) Setting seperation to " << separation << " timesteps"<<endl;
  cellfields->repulsionTimescale = separation;

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

void HemoCell::sanityCheck() {
  hlog << "(HemoCell) (SanityCheck) Performing Sanity check on simulation parameters and setup" << endl;
  //Lattice sanity
  if (param::dx != 5e-7) {
    hlog << "(HemoCell) (SanityCheck) WARNING: Fluid dx is not 5e-7 but " << param::dx << " This is unvalidated!" << endl;
  }

  int env_min_width = (12e-6 /param::dx)+1;
  int env_width = cellfields->immersedParticles->getMultiBlockManagement().getEnvelopeWidth();
  if (env_width < env_min_width) {
    hlog << "(HemoCell) (SanityCheck) WARNING: Envelope width is very small: " << env_width << " (" << env_width*param::dx << "µm) Instead of "  << env_min_width << "!" << endl;

  }
  
  //Material sanity
  if (boundaryRepulsionEnabled) {
    if (cellfields->boundaryRepulsionTimescale%cellfields->particleVelocityUpdateTimescale!=0) {
     hlog << "(HemoCell) Error, Particle velocity timescale separation cannot divide this repulsion timescale separation, exiting ..." <<endl;
     exit(1);
    }
  }
  
  if (repulsionEnabled) {
    if (cellfields->repulsionTimescale%cellfields->particleVelocityUpdateTimescale!=0) {
       hlog << "(HemoCell) Error, Velocity timescale separation cannot divide this repulsion timescale separation, exiting ..." <<endl;
       exit(1);
    }
  }
  
  for (unsigned int i = 0; i < cellfields->size() ; i++) {
    if ((*cellfields)[i]->timescale%cellfields->particleVelocityUpdateTimescale !=0) {
      hlog << "(HemoCell) Error, Velocity timescale separation cannot divide all material timescale separations, exiting ..." <<endl;
      exit(1);
    }
  }
  
  //Cellfields Sanity
  //Check number of neighbours
  if (cellfields->max_neighbours > 30) {
    hlog << "(HemoCell) WARNING: The number of atomic neighbours is suspiciously high: " << cellfields->max_neighbours << " Usually it should be < 30 ! Check the atomic block structure!" << endl;
  }
    
  //Parameter Sanity
#ifdef FORCE_LIMIT
  hlog << "(HemoCell) WARNING: Force limit active at " << FORCE_LIMIT << " pN. Results can be inaccurate due to force capping." << endl;
#endif
  if (sizeof(T) == sizeof(float)) {
    hlog << "(HemoCell) WARNING: Running with single precision, you might want to switch to double precision" << endl;
  }

  // Check lattice viscosity [0.01, 0.45]
  if(param::nu_lbm < 0.01 || param::nu_lbm > 0.45) {
        hlog << "(WARNING!!!) lattice viscosity [" << param::nu_lbm << "] is not in the stable range for LBM [0.01, 0.45]!" << std::endl;
  }

  // Check for lattice velocity to ensure low Courant number (LBM is explicit afterall...)
  if(param::u_lbm_max > 0.1) {
    hlog << "(WARNING!!!) lattice velocity [" << param::u_lbm_max << "] is too high [>0.1]!" << std::endl;
  }
     
  sanityCheckDone = true;
}