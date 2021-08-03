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
#include "wbcHighOrderModel.h"

#include "lykov.h"
#include "palabos3D.h"
#include "palabos3D.hh"
#include "wedge.h"

using namespace hemo;

/// \brief Possible scenarios for the capillary flow examples.
/// The Lykov variant considiers a split channel with two parallel capillaries,
/// whereas the Wedge variant considers a triangular shaped wedge through which
/// a cell is squeezed.
enum Scenario {
  Lykov, /// The Lykov capillary example.
  Wedge, /// A wedge-like geometry.
};

/// The chosen scenario for the simulation.
const auto scenario = Scenario::Lykov;

/// Driving force for the scenario: Scenario::Wedge.
const double wedge_driving_force = 6.7e5;

/// Extended envelope for the ibm kernel: 4 or 2 might be sufficient. NOTE:
/// also depends on the used resolution dx.
const unsigned extendedEnvelopeWidth = 2;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  // Initialise hemocell and the flow parameters.
  HemoCell hemocell(argv[1], argc, argv);
  Config *cfg = hemocell.cfg;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

  // Define the extend of the simulation domain based on the resolution.
  auto resolution = (*cfg)["domain"]["refDirN"].read<int>();
  auto nx = 0;
  auto ny = 0;
  auto nz = 0;

  switch (scenario) {
  case Scenario::Lykov:
    std::tie(nx, ny, nz) = Lykov::domain_size(resolution);
    break;
  case Scenario::Wedge:
    std::tie(nx, ny, nz) = Wedge::domain_size(resolution);
    break;
  };

  // Initialisation of the fluid domain and block lattice configuration.
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
      defaultMultiBlockPolicy3D().getMultiBlockManagement(
          nx, ny, nz, extendedEnvelopeWidth),
      defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0 / param::tau));

  // Define the geometry layout for each scenario. This enforces solid (i.e.
  // bounce back dynamics) to specific regions of the fluid field.
  int error = 0;
  switch (scenario) {
  case Scenario::Lykov: {
    auto capillary_diameter = (*cfg)["domain"]["capillaryD"].read<double>();
    capillary_diameter /= param::dx;
    error = Lykov::geometry(hemocell.lattice, resolution, capillary_diameter);
    break;
  }
  case Scenario::Wedge:
    error = Wedge::geometry(hemocell.lattice, resolution);
    break;
  };
  if (error != 0)
    return error;

  // Both scenarios consider periodicity along x-direction.
  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggle(0, true);
  hemocell.latticeEquilibrium(1., plb::Array<double, 3>(0., 0., 0.));

  // Define driving force and cell properties
  auto driving_force = plb::Array<double, 3>(0, 0, 0);
  std::string cell_type;
  switch (scenario) {
  case Scenario::Lykov: {
    auto dForce = (*cfg)["domain"]["drivingForce"].read<double>();
    driving_force = Lykov::driving_force(dForce);
    cell_type = "WBC_lykov";
    break;
  }
  case Scenario::Wedge:
    driving_force = Wedge::driving_force(wedge_driving_force);
    cell_type = "WBC_wedge";
    break;
  };

  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
                    driving_force);

  // Initialise the lattice and add all the cells.
  hemocell.lattice->initialize();
  hemocell.initializeCellfield();
  hemocell.addCellType<WbcHighOrderModel>(cell_type, WBC_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(cell_type, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  // Micrometer, not LU!
  hemocell.setInitialMinimumDistanceFromSolid(cell_type, 1);
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  // Assign desired output fields for fluid and particle fields.
  hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_FORCE, OUTPUT_BOUNDARY});
  hemocell.setOutputs(cell_type,
                      {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE,
                       OUTPUT_FORCE_VOLUME, OUTPUT_FORCE_BENDING,
                       OUTPUT_FORCE_LINK, OUTPUT_FORCE_INNER_LINK,
                       OUTPUT_FORCE_AREA, OUTPUT_FORCE_VISC});

  // Start fresh or from checkpoint.
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  if (hemocell.iter == 0) {
    auto warmup = (*cfg)["parameters"]["warmup"].read<unsigned>();
    pcout << "(Capillary) fresh start: warming up cell-free fluid domain for "
          << warmup << " iterations..." << endl;
    for (unsigned i = 0; i < warmup; ++i) {
      hemocell.lattice->collideAndStream();
    }
  }

  auto t_max = (*cfg)["sim"]["tmax"].read<unsigned int>();
  auto t_measurement = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  auto t_checkpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  pcout << "(Capillary) Starting simulation..." << endl;

  while (hemocell.iter < t_max) {
    hemocell.iterate();

    // Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      driving_force);

    if (hemocell.iter % t_measurement == 0) {
      // Write output and log basic statistics.
      auto cell_count =
          CellInformationFunctionals::getTotalNumberOfCells(&hemocell);

      pcout << "(main) Stats. @ "
            << hemocell.iter
            << " (" << hemocell.iter * param::dt << " s):" << endl;
      pcout << "\t # of cells: " << cell_count;

      auto finfo = FluidInfo::calculateVelocityStatistics(&hemocell);
      auto toMpS = param::dx / param::dt;
      pcout << "\t Velocity  -  max.: " << finfo.max * toMpS
            << " m/s, mean: " << finfo.avg * toMpS
            << " m/s, rel. app. viscosity: "
            << (param::u_lbm_max * 0.5) / finfo.avg << endl;
      auto pinfo = ParticleInfo::calculateForceStatistics(&hemocell);
      auto topN = param::df * 1.0e12;
      pcout << "\t Force  -  min.: " << pinfo.min * topN
            << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max
            << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

      hemocell.writeOutput();
    }
    if (hemocell.iter % t_checkpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  pcout << "(main) Simulation finished :) " << endl;
  return 0;
}
