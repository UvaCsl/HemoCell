#include "helper/cellInfo.h"
#include "helper/hemoCellStretch.h"
#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"
#include "rbcHighOrderModel.h"
#include "gtest/gtest.h"

// For these tests we probe the force-displacement curve after 10_000
// iterations, where the curves are approaching their steady-state values.
const unsigned max_iteration = 10000;
// The number of points subjected to the stretch force in the tests.
const unsigned n_forced_lsps = 7;
const double to_micro_meter = 1e-6;

// ForcedDiameter is an helper class to setup parameterised tests to assert the
// resulting axial and transverse deformations of an RBC subjected to the
// corresponding stretch force. For each direction we provide an approximate
// lower and upper bound as expected value.
template <typename T> class ForcedDiameter {
public:
  // The applied stretch force
  T force;
  // The approximate bound on the transverse and axial dimensions
  std::tuple<T, T> transverse;
  std::tuple<T, T> axial;

  ForcedDiameter(T force, T t_lower, T t_upper, T a_lower, T a_upper)
      : force{force}, transverse{std::make_tuple(t_lower, t_upper)},
        axial{std::make_tuple(a_lower, a_upper)} {}

  // Ensure `gtest` knows how to format the parameter class on failure.
  template <typename U>
  friend std::ostream &operator<<(std::ostream &, const ForcedDiameter<U> &);
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const ForcedDiameter<T> &fd) {
  os << "stretch force: " << fd.force << ", "
     << "transverse: (" << std::get<0>(fd.transverse) << ", "
     << std::get<1>(fd.transverse) << "), "
     << "axial: (" << std::get<0>(fd.axial) << ", " << std::get<1>(fd.axial)
     << ")\n";
  return os;
}

template<typename T>
T axial_displacement(hemo::Array<T, 6> bbox) {
  auto bbox_in_si = bbox / (to_micro_meter / hemo::param::dx);
  return bbox_in_si[1] - bbox_in_si[0];
}

template<typename T>
T transverse_displacement(hemo::Array<T, 6> bbox) {
  auto bbox_in_si = bbox / (to_micro_meter / hemo::param::dx);
  return bbox_in_si[3] - bbox_in_si[2];
}

class ValidationForceDisplacement
    : public ::testing::TestWithParam<ForcedDiameter<double>> {};

TEST_P(ValidationForceDisplacement, StretchCell) {
  auto parameters = GetParam();
  char *args[] = {(char *)"test", (char *)"path", NULL};
  char *inp = (char *)"validation/stretch_cell/config_stretch_cell.xml";

  hemo::HemoCell hemocell(inp, 0, args, hemo::HemoCell::MPIHandle::External);
  hemo::Config *cfg = hemocell.cfg;
  hemo::param::lbm_base_parameters(*cfg);
  hemo::param::ef_lbm = parameters.force * (1e-12 / hemo::param::df);
  hemo::param::printParameters();

  plint nz = 13 * (to_micro_meter / (*cfg)["domain"]["dx"].read<T>());
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

  auto *boundaryCondition =
      plb::createLocalBoundaryCondition3D<T, DESCRIPTOR>();
  boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice);
  setBoundaryVelocity(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      plb::Array<T, 3>(0., 0., 0.));

  hemocell.latticeEquilibrium(1., hemo::Array<T, 3>({0., 0., 0.}));
  hemocell.lattice->initialize();

  hemocell.initializeCellfield();
  hemocell.addCellType<hemo::RbcHighOrderModel>("validation/stretch_cell/stretch_RBC", RBC_FROM_SPHERE);
  hemocell.loadParticles();

  auto cellfield = (*hemocell.cellfields)["validation/stretch_cell/stretch_RBC"];
  hemo::HemoCellStretch cellStretch(*cellfield, n_forced_lsps, hemo::param::ef_lbm);

  auto initial_volume_lbm = cellfield->meshmetric->getVolume();
  auto initial_volume = initial_volume_lbm / pow(to_micro_meter, 3);

  while (hemocell.iter < max_iteration) {
    cellStretch.applyForce(); // cell force should be applied manually
    hemocell.iterate();
    hemo::CellInformationFunctionals::calculateCellBoundingBox(&hemocell);

    auto bbox = hemo::CellInformationFunctionals::info_per_cell[0].bbox;
    auto axial = axial_displacement(bbox);
    auto transverse = transverse_displacement(bbox);

    if (hemocell.iter == 1 || hemocell.iter % 1000 == 0)
      std::cout << "Iteration: " << hemocell.iter
                << "\ttransverse:\t" << transverse
                << "\taxial:\t" << axial
                << "\t" << parameters << std::endl;

    hemo::CellInformationFunctionals::clear_list();
  }

  // There can be only one cell in the domain.
  auto cell_count = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
  ASSERT_EQ(cell_count, 1);

  hemo::CellInformationFunctionals::calculateCellBoundingBox(&hemocell);
  auto bbox = hemo::CellInformationFunctionals::info_per_cell[0].bbox;

  // Enforce the expected deformation of the bounds
  ASSERT_GE(transverse_displacement(bbox), std::get<0>(parameters.transverse));
  ASSERT_LE(transverse_displacement(bbox), std::get<1>(parameters.transverse));
  ASSERT_GE(axial_displacement(bbox), std::get<0>(parameters.axial));
  ASSERT_LE(axial_displacement(bbox), std::get<1>(parameters.axial));

  // Approximate bounds on the volume change.
  hemo::CellInformationFunctionals::calculateCellVolume(&hemocell);
  auto volume = hemo::CellInformationFunctionals::info_per_cell[0].volume;
  volume /= std::pow(to_micro_meter, 3);
  ASSERT_LE(volume / initial_volume, 1.020);
  ASSERT_GT(volume / initial_volume, 0.980);

  // FIXME: the static nature of the forced points in the cell stretch array
  // require us to reset its content here, as the destructor will not clean up
  // the contents before the end of the complete execution.
  cellStretch.lower_lsps.clear();
  cellStretch.upper_lsps.clear();
}

// Setup the parameterised test suite to assert the axial and transverse
// diameters of a single RBC when subjected to different levels of the external
// stretching force. This asserts the force-displacement curves from Figure 4
// as published in: https://doi.org/10.3389/fphys.2017.00563.
INSTANTIATE_TEST_SUITE_P(
    Validation, ValidationForceDisplacement,
    ::testing::Values(ForcedDiameter<double>(25, 7.3, 7.9, 9.2, 9.7),
                      ForcedDiameter<double>(75, 7.0, 7.5, 11, 12),
                      ForcedDiameter<double>(125, 6.5, 7.0, 12.25, 12.75)));
