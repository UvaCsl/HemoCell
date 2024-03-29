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
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include <fenv.h>
#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;


  int nx, ny, nz;
  nx = ny = nz = (*cfg)["domain"]["refDirN"].read<int>() ;
  //nx = 2 * ny ;
  hlog << "(unbounded) (Parameters) calculating flow parameters" << endl;
  param::lbm_pipe_parameters((*cfg),nx);
  param::printParameters();

  hlog << "(unbounded) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>()),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));

  //hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
  //nx+2, ny+2, nz+2, new BGKdynamics<double,DESCRIPTOR>(1.0/param::tau) );

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.lattice->periodicity().toggle(1,true);
  hemocell.lattice->periodicity().toggle(2,true);
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.0,0.0,0.0));

  hlog << getMultiBlockInfo(*hemocell.lattice) << endl;

  //Driving Force
  hlog << "(PipeFlow) (Fluid) Setting up driving Force" << endl;
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
  DESCRIPTOR<T>::ExternalField::forceBeginsAt,
  plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, poiseuilleForce, poiseuilleForce));


  hemocell.lattice->initialize();

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 1); //Micrometer! not LU

 hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  // Turn on periodicity in the all directions
  hemocell.setSystemPeriodicity(0, true);
  hemocell.setSystemPeriodicity(1, true);
  hemocell.setSystemPeriodicity(2, true);

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
    }

  if (hemocell.iter == 0) {
    hlog << "(unbounded) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  hlog << "(unbounded) Starting simulation..." << endl;
  hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC") << endl;
  int ncells = CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
  hlog << " | RBC Volume ratio [x100%]: " << ncells * 77.0 * 100 / (nx * ny * nz) << endl;
  hlog << "(main)   nCells (global) = " << ncells << endl ;

  while (hemocell.iter < tmax ) {
    hemocell.iterate();

    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
		      plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, poiseuilleForce, poiseuilleForce));


       if (hemocell.iter % tmeas == 0) {
      hlog << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
      hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;


      hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
      }
  }

  hlog << "(main) Simulation finished :) " << endl;

  return 0;
}
