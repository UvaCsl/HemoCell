#ifndef HEMOCELL_CPP
#define HEMOCELL_CPP

#include "hemocell.h"

HemoCell::HemoCell(char * configFileName, int argc, char * argv[]) {
  plbInit(&(argc),&(argv));
  printHeader();

	//TODO This should be done through hemocell config, not some palabos global
  global::directories().setOutputDir("./tmp/");
  global::directories().setLogOutDir("./log/");
  global::directories().setInputDir("./");
  global::IOpolicy().activateParallelIO(true);
  global::IOpolicy().setStlFilesHaveLowerBound(false);
  mkpath("./tmp/hdf5/", 0777);

  pcout << "(HemoCell) (Config) reading " << configFileName << endl;
  cfg = new Config(configFileName);
  documentXML = new XMLreader(configFileName);

  loadBalancer = new LoadBalancer(*this);

  // start clock for basic performance feedback
  lastOutputAt = 0;
  global::timer("atOutput").start();
  
#ifdef FORCE_LIMIT
  pcout << "(HemoCell) WARNING: Force limit active at " << FORCE_LIMIT << " pN. Results can be inaccurate due to force capping." << endl;
#endif
}

void HemoCell::latticeEquilibrium(double rho, hemo::Array<double, 3> vel) {
  pcout << "(HemoCell) (Fluid) Setting Fluid Equilibrium" << endl;
  initializeAtEquilibrium(*lattice, (*lattice).getBoundingBox(), rho, vel);
}

void HemoCell::initializeCellfield() {
  cellfields = new HemoCellFields(*lattice,(*cfg)["domain"]["particleEnvelope"].read<int>(),*this);
}

void HemoCell::setOutputs(string name, vector<int> outputs) {
  pcout << "(HemoCell) (CellField) Setting output variables for " << name << " cells" << endl;
  vector<int> outputs_c = outputs;
  (*cellfields)[name]->setOutputVariables(outputs_c);
}

void HemoCell::setFluidOutputs(vector<int> outputs) {
  pcout << "(HemoCell) (Fluid) Setting output variables for fluid field" << endl;
  vector<int> outputs_c = outputs;
  cellfields->desiredFluidOutputVariables = outputs_c;
}

void HemoCell::loadParticles() {
  pcout << "(HemoCell) (CellField) Loading particle positions "  << endl;
  loadParticlesIsCalled = true;
  readPositionsBloodCellField3D(*cellfields, param::dx, *cfg);
  cellfields->syncEnvelopes();
  cellfields->deleteIncompleteCells(false);
}

void HemoCell::loadCheckPoint() {
  pcout << "(HemoCell) (Saving Functions) Loading Checkpoint"  << endl;
  cellfields->load(documentXML, iter);
}

void HemoCell::saveCheckPoint() {
  pcout << "(HemoCell) (Saving Functions) Saving Checkpoint at timestep " << iter << endl;
  cellfields->save(documentXML, iter);
}

void HemoCell::writeOutput() {
  // Very naive performance approximation
  double dtSinceLastOutput = global::timer("atOutput").stop();
  double timePerIter = dtSinceLastOutput / (iter - lastOutputAt);
  lastOutputAt = iter;

	pcout << "(HemoCell) (Output) writing output at timestep " << iter << " (" << param::dt * iter<< " s). Approx. performance: " << timePerIter << " s / iteration." << endl;
	//Repoint surfaceparticle forces for output
	cellfields->separate_force_vectors();

  //Recalculate the forces
  cellfields->applyConstitutiveModel(true);

  if(repulsionEnabled) {
    cellfields->applyRepulsionForce(true);
  }

  // Creating a new directory per save
  if (global::mpi().isMainProcessor()) {
    string folder = global::directories().getOutputDir() + "/hdf5/" + zeroPadNumber(iter) ;
    mkpath(folder.c_str(), 0777);
    folder = global::directories().getOutputDir() + "/csv/" + zeroPadNumber(iter) ;
    mkpath(folder.c_str(), 0777);
  }
  global::mpi().barrier();

  //Write Output
  writeCellField3D_HDF5(*cellfields,param::dx,param::dt,iter);
  writeFluidField_HDF5(*cellfields,param::dx,param::dt,iter);
  writeCellInfo_CSV(this);

  //Repoint surfaceparticle forces for speed
  cellfields->unify_force_vectors();

  // Continue with performance measurement
  global::timer("atOutput").restart();
}

void HemoCell::iterate() {
  cellfields->applyConstitutiveModel();    // Calculate Force on Vertices

  // Calculate repulsion Force
  if(repulsionEnabled) {
    cellfields->applyRepulsionForce();
  }

  // ##### Particle Force to Fluid ####
  cellfields->spreadParticleForce();

  // ##### 3 ##### LBM
  lattice->timedCollideAndStream();


  // ##### 4 ##### IBM interpolation
  cellfields->interpolateFluidVelocity();

  //### 6 ### Might be together with interpolation
  cellfields->syncEnvelopes();

  //### 7 ### 
  cellfields->advanceParticles();

  // Reset Forces on the lattice, TODO do own efficient implementation
  setExternalVector(*lattice, (*lattice).getBoundingBox(),
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          hemo::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));

  iter++;
}

double HemoCell::calculateFractionalLoadImbalance() {
	pcout << "(HemoCell) (LoadBalancer) Calculating Fractional Load Imbalance at timestep " << iter << endl;
  return loadBalancer->calculateFractionalLoadImbalance();
}

void HemoCell::setMaterialTimeScaleSeparation(string name, unsigned int separation){
  pcout << "(HemoCell) (Timescale Seperation) Setting seperation of " << name << " to " << separation << " timesteps"<<endl;
  pcout << "(HemoCell) WARNING if the timescale separation is not dividable by tmeasure, checkpointing is non-deterministic!"<<endl;
  (*cellfields)[name]->timescale = separation;
}

void HemoCell::setParticlePositionUpdateTimeScaleSeparation(unsigned int separation) {
  pcout << "(HemoCell) (Timescale separation) Setting update separation of all particles to " << separation << " timesteps" << endl;
  pcout << "(HemoCell) WARNING this introduces great errors, WARNING make sure tmeasure and the material timescale separation are dividable by this value" << endl;
  cellfields->particleUpdateTimescale = separation;
}

void HemoCell::setRepulsionTimeScaleSeperation(unsigned int separation){
  pcout << "(HemoCell) (Repulsion Timescale Seperation) Setting seperation to " << separation << " timesteps"<<endl;
  pcout << "(HemoCell) WARNING if the timescale separation is not dividable by tmeasure, checkpointing is non-deterministic!"<<endl;
  cellfields->repulsionTimescale = separation;
}

void HemoCell::setMinimumDistanceFromSolid(string name, double distance) {
  pcout << "(HemoCell) (Set Distance) Setting minimum distance from solid to " << distance << " micrometer for " << name << endl; 
  if (loadParticlesIsCalled) {
    pcout << "(HemoCell) (Set Distance) WARNING: this function is called after the particles are loaded, so it probably has no effect" << endl;
  }
  (*cellfields)[name]->minimumDistanceFromSolid = distance;
}

void HemoCell::setRepulsion(double repulsionConstant, double repulsionCutoff) {
  pcout << "(HemoCell) (Repulsion) Setting repulsion constant to " << repulsionConstant << ". repulsionCutoff to" << repulsionCutoff << " µm" << endl;
  pcout << "(HemoCell) (Repulsion) Enabling repulsion" << endl;
  cellfields->repulsionConstant = repulsionConstant;
  cellfields->repulsionCutoff = repulsionCutoff*(1e-6/param::dx);
  repulsionEnabled = true;
}


#ifdef HEMO_PARMETIS
void HemoCell::doLoadBalance() {
	pcout << "(HemoCell) (LoadBalancer) Balancing Atomic Block over mpi processes" << endl;
  loadBalancer->doLoadBalance();
}
#endif

#endif
