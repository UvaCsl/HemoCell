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
}

void HemoCell::latticeEquilibrium(double rho, Array<double, 3> vel) {
  pcout << "(HemoCell) (Fluid) Setting Fluid Equilibrium" << endl;
  initializeAtEquilibrium(*lattice, (*lattice).getBoundingBox(), rho, vel);
}

void HemoCell::initializeCellfield() {
  cellfields = new HemoCellFields(*lattice,(*cfg)["domain"]["particleEnvelope"].read<int>());
}

void HemoCell::setOutputs(string name, vector<int> outputs) {
  pcout << "(HemoCell) (CellField) Setting output variables for " << name << " cells" << endl;
  vector<int> outputs_c = outputs;
  (*cellfields)[name]->setOutputVariables(outputs_c);
}

void HemoCell::setFluidOutputs(vector<int> outputs) {
  pcout << "(HemoCell) (Fluid) Setting output variables for Fluid" << endl;
  vector<int> outputs_c = outputs;
  cellfields->desiredFluidOutputVariables = outputs_c;
}

void HemoCell::loadParticles(string filename) {
  pcout << "(HemoCell) (CellField) Loading particle positions from " << filename  << endl;
  readPositionsBloodCellField3D(*cellfields, param::dx, filename.c_str(), *cfg);
  cellfields->syncEnvelopes();
  cellfields->deleteIncompleteCells();
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
	pcout << "(HemoCell) (Output) writing desired output at timestep " << iter << endl;
	//Repoint surfaceparticle forces for output
	cellfields->separate_force_vectors();

	//Recalculate the forces
	cellfields->applyConstitutiveModel();

	//Write Output
	writeCellField3D_HDF5(*cellfields,param::dx,param::dt,iter);
	writeFluidField_HDF5(*cellfields,param::dx,param::dt,iter);

	//Repoint surfaceparticle forces for speed
	cellfields->unify_force_vectors();

}

void HemoCell::iterate() {
	cellfields->applyConstitutiveModel();    // Calculate Force on Vertices

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
		Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));

  iter++;
}

double HemoCell::calculateFractionalLoadImbalance() {
	pcout << "(HemoCell) (LoadBalancer) Calculating Fractional Load Imbalance at timestep " << iter << endl;
  return loadBalancer->calculateFractionalLoadImbalance();
}

void HemoCell::doLoadBalance() {
	pcout << "(HemoCell) (LoadBalancer) Balancing Atomic Block over mpi processes" << endl;
  loadBalancer->doLoadBalance();
}


#endif
