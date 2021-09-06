#include "gtest/gtest.h"
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include "palabos3D.h"
#include "palabos3D.hh"

const unsigned warmup_iterations = 100;
const unsigned max_iteration = 1000;
const auto geometry_file = "../examples/pipeflow/tube.stl";

/// Detailed validation test of the cell stretch problem.
TEST(Validation, Pipeflow) {
  char *args[] = {(char *)"test", (char *)"path", NULL};
  char *inp = (char *)"validation/pipeflow/config_pipeflow.xml";

  hemo::HemoCell hemocell(inp, 0, args, hemo::HemoCell::MPIHandle::External);
  hemo::Config * cfg = hemocell.cfg;

  std::auto_ptr<plb::MultiScalarField3D<int>> flagMatrix;
  std::auto_ptr<hemo::VoxelizedDomain3D<T>> voxelizedDomain;

  hemo::getFlagMatrixFromSTL(geometry_file,
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                       (*cfg)["domain"]["refDirN"].read<int>(),
                       (*cfg)["domain"]["refDir"].read<int>(),
                       voxelizedDomain, flagMatrix,
                       (*cfg)["domain"]["blockSize"].read<int>(),
                       (*cfg)["domain"]["particleEnvelope"].read<int>());

  hemo::param::lbm_pipe_parameters((*cfg), flagMatrix.get());
  hemo::param::printParameters();

  hemocell.lattice = new plb::MultiBlockLattice3D<T, DESCRIPTOR>(
            voxelizedDomain.get()->getMultiBlockManagement(),
            plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
            plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
            plb::defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new plb::GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/hemo::param::tau));

  defineDynamics(*hemocell.lattice, *flagMatrix.get(), (*hemocell.lattice).getBoundingBox(), new hemo::BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1., {0., 0., 0.});

  //Driving Force
  auto poiseuilleForce =  8 * hemo::param::nu_lbm * (hemo::param::u_lbm_max * 0.5) / hemo::param::pipe_radius / hemo::param::pipe_radius;
  auto driving_force = plb::Array<T, 3> {poiseuilleForce, 0., 0.};

  hemocell.lattice->initialize();
  hemocell.initializeCellfield();

  hemocell.addCellType<hemo::RbcHighOrderModel>("validation/pipeflow/RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("validation/pipeflow/RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("validation/pipeflow/RBC", 0.5);

  hemocell.addCellType<hemo::PltSimpleModel>("validation/pipeflow/PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("validation/pipeflow/PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);
  hemocell.loadParticles();

  // Enable the external force from the start.
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    driving_force);

  // Perform warmup iterations on the fluid field only.
  for (unsigned i = 0; i < warmup_iterations; ++i)
    hemocell.lattice->collideAndStream();

  while (hemocell.iter < max_iteration ) {
    hemocell.iterate();
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                driving_force);

    // After initialisation only 42 partiles (either RBC or PLT) remain within
    // the pipe for the currently specified resolution. NOTE: the initial .pos
    // files do contain significantly more cells (> 200).
    auto count = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
    ASSERT_EQ(count, 42);

    if (hemocell.iter < 100)
      // let the pipeflow run for some iterations before assertions
      continue;

    // Enforce bounds on the viscosity observed in the pipe.
    auto fluid_info = hemo::FluidInfo::calculateVelocityStatistics(&hemocell);
    auto viscosity = 0.5*hemo::param::u_lbm_max / fluid_info.avg;
    ASSERT_LT(viscosity, 3.0);
    ASSERT_GT(viscosity, 1.03);

    // Enforce average force on the particles.
    auto particle_info = hemo::ParticleInfo::calculateForceStatistics(&hemocell);
    auto average_force = particle_info.avg * hemo::param::df * 1.0e12;
    ASSERT_LT(average_force, 4.0);
  }
}
