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
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "helper/hemocellInit.hh"
#include "writeCellInfoCSV.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

typedef double T;
using namespace hemo;
using namespace plb;
using namespace std;

int main(int argc, char * argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  // ----------------- Read in config file & geometry ---------------------------

  std::unique_ptr<MultiScalarField3D<int>> flagMatrix;
  std::unique_ptr<VoxelizedDomain3D<T>> voxelizedDomain;
  getFlagMatrixFromSTL(*cfg, voxelizedDomain, flagMatrix);

  plint nx = (*cfg)["domain"]["refDirN"].read<int>();
  plint ny = 0.5 * nx;
  plint nz = 0.75 * nx;
  param::lbm_shear_parameters((*cfg), nz);
  param::printParameters();

  // ------------------------ Init lattice --------------------------------
  pcout << "(stent_strut) Initializing lattice: " << nx << "x" << ny << "x" << nz << " [lu]" << std::endl;

  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      voxelizedDomain.get()->getMultiBlockManagement(), defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));

  // -------------------------- Define boundary conditions ---------------------

  OnLatticeBoundaryCondition3D<T, DESCRIPTOR> * boundaryCondition = createLocalBoundaryCondition3D<T, DESCRIPTOR>();
  Box3D bb = hemocell.lattice->getBoundingBox();
  Box3D bounceback_box(bb.x0 + 1, bb.x1 - 1, bb.y0 + 1, bb.y1 - 1, bb.z0, bb.z1 - 1);

  defineDynamics(*hemocell.lattice, *flagMatrix.get(), bounceback_box, new BounceBack<T, DESCRIPTOR>(1.), 0);

  // -------------------------- Define shear surface ---------------------
  T surf_velocityLU =
      (*cfg)["domain"]["velocity"].read<T>() *
      ((*cfg)["domain"]["dt"].read<T>() / (*cfg)["domain"]["dx"].read<T>()); // unit conversion from m/s to LUs
  pcout << "(stent_strut) Surface velocity: " << surf_velocityLU << " in LU" << endl;

  Box3D shear_surf(bb.x0, bb.x1, bb.y0, bb.y1, bb.z1, bb.z1);

  boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice, shear_surf);
  setBoundaryVelocity(*hemocell.lattice, shear_surf, plb::Array<T, 3>(surf_velocityLU, 0, 0));


  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0, true);
  hemocell.lattice->periodicity().toggle(1, true); //b
  hemocell.latticeEquilibrium(1., plb::Array<double, 3>(0., 0., 0.));

  hemocell.lattice->initialize();

  // ----------------------- Init cell models --------------------------

  hemocell.initializeCellfield();
  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 0.5); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,      OUTPUT_TRIANGLES,  OUTPUT_FORCE,      OUTPUT_FORCE_VOLUME,
                         OUTPUT_FORCE_BENDING, OUTPUT_FORCE_LINK, OUTPUT_FORCE_AREA, OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,   OUTPUT_DENSITY,     OUTPUT_FORCE,       OUTPUT_BOUNDARY,
             OUTPUT_SHEAR_RATE, OUTPUT_STRAIN_RATE, OUTPUT_SHEAR_STRESS};
  hemocell.setFluidOutputs(outputs);

  // ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    pcout << "(stent_strut) CHECKPOINT found!" << endl;
    hemocell.loadCheckPoint();
  }


  if (hemocell.iter == 0) {
    pcout << "(stent_strut) fresh start: warming up cell-free fluid domain for "
          << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tcsv = (*cfg)["sim"]["tcsv"].read<unsigned int>();


  while (hemocell.iter < tmax) {
    hemocell.iterate();

    if (hemocell.iter % tmeas == 0) {
      hlog << "(main) Stats. @ " << hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
      hlog << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell);
      T toMpS = param::dx / param::dt;
      hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS
           << " m/s, rel. app. viscosity: " << (param::u_lbm_max * 0.5) / finfo.avg << endl;
      ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell);
      T topN = param::df * 1.0e12;
      hlog << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max
           << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

      hemocell.writeOutput();
    }

    if (hemocell.iter % tcsv == 0) {
      hlog << "Saving simple mean cell values to CSV at timestep " << hemocell.iter << endl;
      writeCellInfo_CSV(hemocell);
    }

    if (hemocell.iter % tcheckpoint == 0) { hemocell.saveCheckPoint(); }
  }

  pcout << "(stent_strut) Simulation finished :)" << std::endl;
  return 0;
}
