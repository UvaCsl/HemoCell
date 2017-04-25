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
}

void HemoCell::latticeEquilibrium(double rho, Array<double, 3> vel) {
  pcout << "(HemoCell) (Fluid) Setting Fluid Equilibrium" << endl;
  initializeAtEquilibrium(*lattice, (*lattice).getBoundingBox(), rho, vel);
}


#endif
