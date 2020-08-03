#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>
#include "helper/hemocellInit.hh"
#include "helper/leesEdwardsBC.h"
#include "helper/cellInfo.h"

int main (int argc, char * argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1],argc,argv);
  Config * cfg = hemocell.cfg;

  plint nz = 100.0*(1e6*(*cfg)["domain"]["dx"].read<T>());
  plint ny = 100.0*(1e6*(*cfg)["domain"]["dx"].read<T>());
  plint nx = 100.0*(1e6*(*cfg)["domain"]["dx"].read<T>());
  double dt = (*cfg)["domain"]["dt"].read<double>();
  plint extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with with 2.

  pcout << "(CellStretch) Initializing lattice: " << nx <<"x" << ny <<"x" << nz << " [lu]" << std::endl;

  param::lbm_shear_parameters((*cfg),nz);
  param::printParameters();

  hemocell.lattice = new MultiBlockLattice3D<T,DESCRIPTOR>(
    defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
    defaultMultiBlockPolicy3D().getBlockCommunicator(),
	defaultMultiBlockPolicy3D().getCombinedStatistics(),
	defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
    new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

  hemocell.lattice->toggleInternalStatistics(false);
  LeesEdwardsBC<T, DESCRIPTOR> LEbc (*hemocell.lattice, param::shearrate_lbm, dt, &hemocell.LEcurrentDisplacement);
  LEbc.initialize();

  hemocell.lattice->initialize();
  // hemocell.outputInSiUnits = true;

  hemocell.initializeCellfield();
  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", 20);
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(5);
  hemocell.setOutputs("RBC_HO", {OUTPUT_POSITION,OUTPUT_VELOCITY,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA, OUTPUT_FORCE_VISC} );

  hemocell.setFluidOutputs( { OUTPUT_VELOCITY, OUTPUT_BOUNDARY, OUTPUT_SHEAR_STRESS,OUTPUT_STRAIN_RATE, OUTPUT_DENSITY } );

  hemocell.loadParticles();

  if (hemocell.iter == 0)    {
    pcout << "(OneCellShear) fresh start: warming up cell-free fluid domain for " << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt)
    {
      hemocell.lattice->collideAndStream();
    }
  }


  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  while (hemocell.iter < tmax)
  {
    hemocell.iterate();
    LEbc.updateLECurDisplacement(hemocell.iter);

    if (hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();
    }
  }

  return 0;
}
