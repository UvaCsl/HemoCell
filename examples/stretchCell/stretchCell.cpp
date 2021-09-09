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
#include "helper/cellInfo.h"
#include "helper/hemoCellStretch.h"
#include "hemocell.h"
#include "rbcHighOrderModel.h"

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

// The used IBM kernel width.
const plb::plint extendedEnvelopeWidth = 2;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  // Initialise HemoCell and parse the parameteres configuration file.
  HemoCell hemocell(argv[1], argc, argv);
  Config &cfg = *hemocell.cfg;

  param::lbm_base_parameters(cfg);
  param::printParameters();

  // Extract the stretch force, both in lattice units and SI.
  auto stretch_force = cfg["parameters"]["stretchForce"].read<T>();
  auto stretch_force_si = stretch_force * 1e-12 / param::df;

  // The resolution of the domain in number of cells to achieve an domain with
  // size 26x13x13 micron (x, y, z).
  plint nz = 13 * (1e-6 / cfg["domain"]["dx"].read<T>());
  plint nx = 2 * nz;
  plint ny = nz;
  hlog << "Lattice: " << nx << "x" << ny << "x" << nz << " [lu]" << std::endl;

  // An default initialisation for the lattice.
  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      defaultMultiBlockPolicy3D().getMultiBlockManagement(
          nx, ny, nz, extendedEnvelopeWidth),
      defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);

  // All boundaries are subjected to zero velocity.
  OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *boundaryCondition =
      createLocalBoundaryCondition3D<T, DESCRIPTOR>();

  boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice);
  setBoundaryVelocity(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      plb::Array<T, 3>(0., 0., 0.));

  // The domain is initialised with zero velocity.
  hemocell.latticeEquilibrium(1., hemo::Array<T, 3>({0., 0., 0.}));
  hemocell.lattice->initialize();

  // Initialise cells particles
  hemocell.initializeCellfield();
  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);

  // Define any output fields of interest to be written to HDF5 output. 
  hemocell.setOutputs("RBC", {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE,
                              OUTPUT_FORCE_VOLUME, OUTPUT_FORCE_BENDING,
                              OUTPUT_FORCE_LINK, OUTPUT_FORCE_AREA,
                              OUTPUT_FORCE_VISC});
  hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_FORCE});

  // HDF5 output in SI units (except location (so fluid location, particle
  // location is still in LU)
  hemocell.outputInSiUnits = true;

  // loading the cellfield
  if (cfg.checkpointed) {
    hemocell.loadCheckPoint();
  } else {
    hemocell.loadParticles();
    hemocell.writeOutput();
  }

  // Setting up the stretching
  auto rbc_field = (*hemocell.cellfields)["RBC"];
  unsigned n_forced_lsps = 1 + 6;
  HemoCellStretch cellStretch(*rbc_field, n_forced_lsps, stretch_force_si);

  hlog << "External stretching force" << stretch_force << " [pN(flb)] ("
       << stretch_force_si << ")" << endl;

  // Get undeformed values
  T volume_lbm = rbc_field->meshmetric->getVolume();
  T surface_lbm = rbc_field->meshmetric->getSurface();
  T volume_eq = volume_lbm / pow(1e-6 / param::dx, 3);
  T surface_eq = surface_lbm / pow(1e-6 / param::dx, 2);

  auto bb = rbc_field->getOriginalBoundingBox();
  hlog << "Original Bounding box:" << endl;
  hlog << "\tx: " << bb[0] << " : " << bb[1] << endl;
  hlog << "\ty: " << bb[2] << " : " << bb[3] << endl;
  hlog << "\tz: " << bb[4] << " : " << bb[5] << endl;

  // Creating output log file
  plb_ofstream fOut;
  auto filename = "stretch-" + std::to_string((int)stretch_force) + ".log";
  if (cfg.checkpointed) {
    fOut.open(filename.c_str(), std::ofstream::app);
  } else {
    // Initialise file with header of data
    fOut.open(filename.c_str());
    fOut << "iteration axial transverse" << std::endl;
  }

  unsigned int tmax = cfg["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = cfg["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = cfg["sim"]["tcheckpoint"].read<unsigned int>();

  while (hemocell.iter < tmax) {
    // NOTE: the `applyForce()` is not done in `hemocell.iterate()` and has to
    // be called manually every iteration to ensure the forces are applied.
    cellStretch.applyForce();
    hemocell.iterate();

    if (hemocell.iter == 1 || hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();

      // Fill up the static info structure with desired data
      CellInformationFunctionals::calculateCellVolume(&hemocell);
      CellInformationFunctionals::calculateCellArea(&hemocell);
      CellInformationFunctionals::calculateCellPosition(&hemocell);
      CellInformationFunctionals::calculateCellStretch(&hemocell);
      CellInformationFunctionals::calculateCellBoundingBox(&hemocell);

      T volume = (CellInformationFunctionals::info_per_cell[0].volume) / pow(1e-6 / param::dx, 3);
      T surface = (CellInformationFunctionals::info_per_cell[0].area) / pow(1e-6 / param::dx, 2);
      hemo::Array<T, 3> position = CellInformationFunctionals::info_per_cell[0].position / (1e-6 / param::dx);
      hemo::Array<T, 6> bbox = CellInformationFunctionals::info_per_cell[0].bbox / (1e-6 / param::dx);
      T largest_diam = (CellInformationFunctionals::info_per_cell[0].stretch) / (1e-6 / param::dx);

      hlog << "\tCell center at: {" << position[0] << ", " << position[1] << ", " << position[2] << "} µm" << endl;
      hlog << "\tDiameters: {" << bbox[1] - bbox[0] << ", " << bbox[3] - bbox[2] << ", " << bbox[5] - bbox[4] << "} µm" << endl;
      hlog << "\tSurface: " << surface << " µm^2" << " (" << surface / surface_eq * 100.0 << "%)" << endl;
      hlog << "\tVolume: " << volume << " µm^3" << " (" << volume / volume_eq * 100.0 << "%)" << endl;
      hlog << "\tLargest diameter: " << largest_diam << " µm." << endl;

      // Track (iteration, axial, transverse) in the output file for
      // possible post-processing.
      fOut << hemocell.iter << " " << bbox[1] - bbox[0] << " "
           << bbox[3] - bbox[2] << endl;

      CellInformationFunctionals::clear_list();
    }

    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  fOut.close();
  hlog << "(CellStretch) Simulation finished :)" << std::endl;
  return 0;
}
