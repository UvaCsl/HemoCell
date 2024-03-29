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
#include "wbcHighOrderModel.h"
#include "helper/hemocellInit.hh"
#include "helper/cellInfo.h"
#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
      cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
      return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;




// ----------------- Read in config file & calc. LBM parameters ---------------------------
  pcout << "(CellCollision) (Parameters) calculating shear flow parameters" << endl;
  plint nx = 25.0*(1e-6/(*cfg)["domain"]["dx"].read<T>());
  plint ny = nx;
  plint nz = ny*0.6;
  param::lbm_shear_parameters((*cfg),ny);
  param::printParameters();

  // ------------------------ Init lattice --------------------------------

  pcout << "(CellCollision) Initializing lattice: " << nx <<"x" << ny <<"x" << nz << " [lu]" << std::endl;

  plint extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with with 2.

  hemocell.lattice = new MultiBlockLattice3D<T,DESCRIPTOR>(
      defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
      defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

  pcout << "(CellCollision) Re corresponds to u_max = " << (param::re * param::nu_p)/(hemocell.lattice->getBoundingBox().getNy()*param::dx) << " [m/s]" << endl;
  // -------------------------- Define boundary conditions ---------------------

  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
      = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

  hemocell.lattice->toggleInternalStatistics(false);

  iniLatticeSquareCouette(*hemocell.lattice, nx, ny, nz, *boundaryCondition, param::shearrate_lbm);

  hemocell.lattice->initialize();

  // ----------------------- Init cell models --------------------------

  hemocell.initializeCellfield();
  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_INNER_LINKS};
  hemocell.addCellType<WbcHighOrderModel>("ELL", ELLIPSOID_FROM_SPHERE);
  hemocell.setOutputs("ELL", outputs);
  hemocell.addCellType<WbcHighOrderModel>("ELL2", ELLIPSOID_FROM_SPHERE);
  hemocell.setOutputs("ELL2", outputs);

  outputs = {OUTPUT_VELOCITY};
  hemocell.setFluidOutputs(outputs);
  hemocell.outputInSiUnits = true; //HDF5 output in SI units (except location (so fluid location, particle location is still in LU)

// ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    pcout << "(CellCollision) CHECKPOINT found!" << endl;
    hemocell.loadCheckPoint();
  }


  if (hemocell.iter == 0) {
    pcout << "(CellCollision) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  pcout << "(CellCollision) Shear rate: " << (*cfg)["domain"]["shearrate"].read<T>() << " s^-1." << endl;

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();


  while (hemocell.iter < tmax ) {

    hemocell.iterate();

    if (hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  pcout << "(CellCollision) Simulation finished :)" << std::endl;
  return 0;
}
