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
#include "helper/hemocellInit.hh"
#include "helper/cellInfo.h"
#include "gtest/gtest.h"
#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;
using namespace plb;
using namespace std;

// upper and lower bound for the deformation index
const T exp_def_upper = 0.0127919;
const T exp_def_lower = 0.0127918;

TEST(Validation, oneCellShear) {
  char * args[] = {(char *)"test", (char *)"path", NULL};
  char * inp = (char *)"tests/validation/oneCellShear/config.xml";

  hemo::HemoCell hemocell(inp, 0, args, hemo::HemoCell::MPIHandle::External);
  Config * cfg = hemocell.cfg;

  // ----------------- Read in config file & calc. LBM parameters ---------------------------
  pcout << "(OneCellShear) (Parameters) calculating shear flow parameters" << endl;
  plint nz = 10.0 * (1e-6 / (*cfg)["domain"]["dx"].read<T>());
  plint nx = 2 * nz;
  plint ny = 2 * nz;
  param::lbm_shear_parameters((*cfg), ny);
  param::printParameters();
  plint max_iteration = (*cfg)["sim"]["tmax"].read<plint>();   // number of iterations
  plint warmup = (*cfg)["parameters"]["warmup"].read<plint>(); // warm-up iterations for the fluid
  T D0 = 2.0 * (*cfg)["ibm"]["radius"].read<T>() * 1e6;
  T def_idx; // deformation index

  // ------------------------ Init lattice --------------------------------

  pcout << "(CellStretch) Initializing lattice: " << nx << "x" << ny << "x" << nz << " [lu]" << std::endl;

  plint extendedEnvelopeWidth = 2; // Because we might use ibmKernel with with 2.

  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
      defaultMultiBlockPolicy3D().getBlockCommunicator(), defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));

  pcout << "(OneCellShear) Re corresponds to u_max = "
        << (param::re * param::nu_p) / (hemocell.lattice->getBoundingBox().getNy() * param::dx) << " [m/s]" << endl;
  // -------------------------- Define boundary conditions ---------------------

  OnLatticeBoundaryCondition3D<T, DESCRIPTOR> * boundaryCondition = createLocalBoundaryCondition3D<T, DESCRIPTOR>();

  hemocell.lattice->toggleInternalStatistics(false);

  iniLatticeSquareCouette(*hemocell.lattice, nx, ny, nz, *boundaryCondition, param::shearrate_lbm);

  hemocell.lattice->initialize();
  hemocell.outputInSiUnits = true;

  // ---------------------- Initialise particle positions and perform warm-up iterations for the fluid ---------------

  hemocell.initializeCellfield();
  hemocell.addCellType<RbcHighOrderModel>("tests/validation/oneCellShear/RBC", RBC_FROM_SPHERE);

  hemocell.loadParticles();

  if (hemocell.iter == 0) {
    pcout << "(OneCellShear) fresh start: warming up cell-free fluid domain for "
          << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < warmup; ++itrt) { hemocell.lattice->collideAndStream(); }
  }

  pcout << "(OneCellShea) Shear rate: " << (*cfg)["domain"]["shearrate"].read<T>() << " s^-1." << endl;

  // --------------------- Run the simulation --------------------------------------------

  while (hemocell.iter < max_iteration) {
    hemocell.iterate();

    // Fill up the static info structure with desired data
    CellInformationFunctionals::calculateCellVolume(&hemocell);
    CellInformationFunctionals::calculateCellArea(&hemocell);
    CellInformationFunctionals::calculateCellPosition(&hemocell);
    CellInformationFunctionals::calculateCellStretch(&hemocell);
    CellInformationFunctionals::calculateCellBoundingBox(&hemocell);

    T largest_diam = (CellInformationFunctionals::info_per_cell[0].stretch) / (1e-6 / param::dx);
    T rel_D2 = (largest_diam / D0) * (largest_diam / D0);
    def_idx = (rel_D2 - 1.0) / (rel_D2 + 1.0) * 100.0;

    CellInformationFunctionals::clear_list();
  }
  // --------------------- Perform the assertions --------------------------------------

  ASSERT_LT(def_idx, exp_def_upper);
  ASSERT_GT(def_idx, exp_def_lower);
}
