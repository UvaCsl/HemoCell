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
#include "cellInfo.h"
#include "fluidInfo.h"
#include "hemocell.h"
#include "particleInfo.h"
#include "pltSimpleModel.h"
#include "rbcHighOrderModel.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

typedef double T;

using namespace hemo;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config *cfg = hemocell.cfg;

  // equal sized-cube
  int nx, ny, nz;
  nx = ny = nz = (*cfg)["domain"]["refDirN"].read<int>();

  hlog << "(unbounded) (Parameters) calculating flow parameters" << endl;
  param::lbm_shear_parameters((*cfg), nx);
  param::printParameters();

  hlog << "(unbounded) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>()),
      defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));


  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
                = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

  hemocell.lattice->toggleInternalStatistics(false);

  // top and bottom sides of the domain
  Box3D top = Box3D(0, nx-1, 0, ny-1, nz-1, nz-1);
  Box3D bottom = Box3D(0, nx-1, 0, ny-1, 0, 0);

  // all directions have periodicity
  hemocell.lattice->periodicity().toggle(0, true);
  hemocell.lattice->periodicity().toggle(1, true);
  hemocell.lattice->periodicity().toggle(2, true);

  // assign velocity at top/bottom
  (*boundaryCondition).setVelocityConditionOnBlockBoundaries ( *hemocell.lattice, top );
  (*boundaryCondition).setVelocityConditionOnBlockBoundaries ( *hemocell.lattice, bottom );

  // define shear rate along z axis
  T vHalf = (nz-1)*param::shearrate_lbm*0.5;
  setBoundaryVelocity(*hemocell.lattice, top, plb::Array<T,3>(-vHalf,0.0,0.0));
  setBoundaryVelocity(*hemocell.lattice, bottom, plb::Array<T,3>(vHalf,0.0,0.0));

  // not sure what this does?
  hemocell.latticeEquilibrium(1., plb::Array<T, 3>(0.0, 0.0, 0.0));

  // ??
  hlog << getMultiBlockInfo(*hemocell.lattice) << endl;

  // assign force vector equal in x, y, z directions
  // FIXME (*hemocell.lattice). to hemocell.lattice->?
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(
                        0.0, 0.0, 0.0));

  // initialise the lattic
  hemocell.lattice->initialize();

  // initialise the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(
      "RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.setParticleVelocityUpdateTimeScaleSeparation(
      (*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE};
  hemocell.setOutputs("RBC", outputs);

  outputs = {OUTPUT_VELOCITY, OUTPUT_DENSITY};
  hemocell.setFluidOutputs(outputs);

  // Turn on periodicity in the all directions
  hemocell.setSystemPeriodicity(0, true);
  hemocell.setSystemPeriodicity(1, true);
  hemocell.setSystemPeriodicity(2, true);

  // loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  if (hemocell.iter == 0) {
    hlog << "(unbounded) fresh start: warming up cell-free fluid domain for "
         << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..."
         << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>();
         ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();

  hlog << "(unbounded) Starting simulation..." << endl;
  int ncells =
      CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
  hlog << " | RBC Volume ratio [x100%]: "
       << ncells * 77.0 * 100 / (nx * ny * nz) << endl;
  hlog << "(main)   nCells (global) = " << ncells << endl;

  while (hemocell.iter < tmax) {
    hemocell.iterate();

    // Set driving force as required after each iteration
//    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
//                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
//                      plb::Array<T, DESCRIPTOR<T>::d>(
//                          poiseuilleForce, 0, 0));
//
    if (hemocell.iter % tmeas == 0) {
      hlog << "(main) Stats. @ " << hemocell.iter << " ("
           << hemocell.iter * param::dt << " s):" << endl;
      hlog << "\t # of cells: "
           << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      hlog << " | # of RBC: "
           << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell,
                                                                   "RBC");
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell);
      double toMpS = param::dx / param::dt;
      hlog << "\t Velocity  -  max.: " << finfo.max * toMpS
           << " m/s, mean: " << finfo.avg * toMpS
           << " m/s, rel. app. viscosity: "
           << (param::u_lbm_max * 0.5) / finfo.avg << endl;

      hemocell.writeOutput();
    }
  }

  hemo::global.statistics.printStatistics();
  hemo::global.statistics.outputStatistics();

  hlog << "(main) Simulation finished :) " << endl;
  return 0;
}
