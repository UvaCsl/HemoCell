#include "gtest/gtest.h"
#include "helper/cellInfo.h"
#include "helper/hemoCellStretch.h"
#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"
#include "rbcHighOrderModel.h"

/// Detailed validation test of the cell stretch problem.
TEST(Validation, StretchCell) {
  char *args[] = {(char *)"test", (char *)"path", NULL};
  char *inp = (char *)"validation/stretch_cell/config_stretch_cell.xml";

  hemo::HemoCell hemocell(inp, 0, args, hemo::HemoCell::MPIHandle::External);
  hemo::Config *cfg = hemocell.cfg;
  hemo::param::lbm_base_parameters(*cfg);
  hemo::param::ef_lbm =
      (*cfg)["parameters"]["stretchForce"].read<T>() * 1e-12 / hemo::param::df;
  hemo::param::printParameters();

  plint nz = 13 * (1e-6 / (*cfg)["domain"]["dx"].read<T>());
  plint nx = 2 * nz;
  plint ny = nz;
  plint extendedEnvelopeWidth = 2;

  hemocell.lattice = new plb::MultiBlockLattice3D<T, DESCRIPTOR>(
      plb::defaultMultiBlockPolicy3D().getMultiBlockManagement(
          nx, ny, nz, extendedEnvelopeWidth),
      plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
      plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
      plb::defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new plb::GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 /
                                                          hemo::param::tau));

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);

  auto *boundaryCondition = plb::createLocalBoundaryCondition3D<T, DESCRIPTOR>();
  boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice);
  setBoundaryVelocity(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      plb::Array<T, 3>(0., 0., 0.));

  hemocell.latticeEquilibrium(1., hemo::Array<T, 3>({0., 0., 0.}));
  hemocell.lattice->initialize();

  hemocell.initializeCellfield();
  hemocell.addCellType<hemo::RbcHighOrderModel>("validation/stretch_cell/stretch_RBC", RBC_FROM_SPHERE);
  hemocell.loadParticles();
  auto cellfield = (*hemocell.cellfields)["validation/stretch_cell/stretch_RBC"];

  // Setting up the stretching
  unsigned int n_forced_lsps = 1 + 6;
  hemo::HemoCellStretch cellStretch(*cellfield, n_forced_lsps, hemo::param::ef_lbm);

  auto tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();

  // Get undeformed values
  auto to_micro_meter = 1e-6 / hemo::param::dx;
  auto initial_volume_lbm = cellfield->meshmetric->getVolume();
  auto initial_volume = initial_volume_lbm / pow(to_micro_meter, 3);

  while (hemocell.iter < tmax) {
    // The cell force has to be applied manully.
    cellStretch.applyForce();
    hemocell.iterate();

    // There can be only one stretching cell in the domain.
    auto count = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
    ASSERT_EQ(count, 1);

    if (hemocell.iter < 100)
      // Allow some setup iterations before the assertions are actually
      // asserted. This is inherited from the previous test scripts that only
      // enforce these assertions every 100 iterations, however, here we _do_
      // enforce them every iteration after the first 100 and not only every
      // 100th iteration.
      continue;

    hemo::CellInformationFunctionals::calculateCellVolume(&hemocell);
    hemo::CellInformationFunctionals::calculateCellArea(&hemocell);
    hemo::CellInformationFunctionals::calculateCellStretch(&hemocell);

    auto volume = (hemo::CellInformationFunctionals::info_per_cell[0].volume);
    auto surface = (hemo::CellInformationFunctionals::info_per_cell[0].area);
    auto stretch = (hemo::CellInformationFunctionals::info_per_cell[0].stretch);

    // convert all to µ-meter
    volume /= std::pow(to_micro_meter, 3);
    surface /= std::pow(to_micro_meter, 2);
    stretch /= to_micro_meter;

    // cell stretch upper bound 9.6µm
    ASSERT_LE(stretch, 9.6);

    // cell surface area upper and lower bounds in µm^2
    ASSERT_LE(surface, 133.04);
    ASSERT_GE(surface, 129.34);

    // cell volume upper and lower bounds in µm^3
    ASSERT_LE(volume, 81.19);
    ASSERT_GE(volume, 81.12);

    // relative volume change between 100.0% and 100.1%
    ASSERT_LE(volume / initial_volume, 1.001);
    ASSERT_GT(volume / initial_volume, 1.000);

    // clear the gathered information
    hemo::CellInformationFunctionals::clear_list();
  }
}
