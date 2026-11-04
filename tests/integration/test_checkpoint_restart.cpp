#include "gtest/gtest.h"
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "palabos3D.h"
#include "palabos3D.hh"

#include <memory>
#include <string>
#include <vector>

// Checkpoint -> restart round trip on the small pipeflow validation case.
// Asserts that the state which cannot be re-derived from the input files is
// written to and restored from the checkpoint.

namespace {

const auto pipeflow_config = "tests/validation/pipeflow/config_pipeflow.xml";
const auto pipeflow_geometry = "examples/pipeflow/tube.stl";
const auto rbc_name = "tests/validation/pipeflow/RBC";
const auto plt_name = "tests/validation/pipeflow/PLT";

struct PipeflowSetup {
  std::auto_ptr<plb::MultiScalarField3D<int>> flagMatrix;
  std::auto_ptr<hemo::VoxelizedDomain3D<T>> voxelizedDomain;
  plb::Array<T, 3> drivingForce;
};

// Same domain, lattice and cell types as Validation.Pipeflow, without loading
// particles (that differs between a fresh start and a restart).
void setupPipeflow(hemo::HemoCell & hemocell, PipeflowSetup & setup) {
  hemo::Config * cfg = hemocell.cfg;
  hemo::getFlagMatrixFromSTL(pipeflow_geometry,
                             (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                             (*cfg)["domain"]["refDirN"].read<int>(),
                             (*cfg)["domain"]["refDir"].read<int>(),
                             setup.voxelizedDomain, setup.flagMatrix,
                             (*cfg)["domain"]["blockSize"].read<int>(),
                             (*cfg)["domain"]["particleEnvelope"].read<int>());
  hemo::param::lbm_pipe_parameters((*cfg), setup.flagMatrix.get());

  hemocell.lattice = new plb::MultiBlockLattice3D<T, DESCRIPTOR>(
      setup.voxelizedDomain.get()->getMultiBlockManagement(),
      plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
      plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
      plb::defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new plb::GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/hemo::param::tau));
  defineDynamics(*hemocell.lattice, *setup.flagMatrix.get(), (*hemocell.lattice).getBoundingBox(),
                 new hemo::BounceBack<T, DESCRIPTOR>(1.), 0);
  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1., {0., 0., 0.});

  const T poiseuilleForce = 8 * hemo::param::nu_lbm * (hemo::param::u_lbm_max * 0.5)
                            / hemo::param::pipe_radius / hemo::param::pipe_radius;
  setup.drivingForce = plb::Array<T, 3>(poiseuilleForce, 0., 0.);

  hemocell.lattice->initialize();
  hemocell.initializeCellfield();
  hemocell.addCellType<hemo::RbcHighOrderModel>(rbc_name, RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(rbc_name, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid(rbc_name, 0.5);
  hemocell.addCellType<hemo::PltSimpleModel>(plt_name, ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(plt_name, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());
  hemocell.setSystemPeriodicity(0, true);
}

void iterate(hemo::HemoCell & hemocell, const PipeflowSetup & setup, unsigned int steps) {
  for (unsigned int i = 0; i < steps; i++) {
    hemocell.iterate();
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, setup.drivingForce);
  }
}

} // namespace

TEST(Checkpoint, RestartRestoresState) {
  char *args[] = {(char *)"test", NULL};

  unsigned int savedIteration = 0;
  unsigned int savedLastOutput = 0;
  int savedNumberOfCells = 0;
  int savedCellCount = 0;
  std::string checkpointFile;

  {
    hemo::HemoCell hemocell((char *)pipeflow_config, 0, args, hemo::HemoCell::MPIHandle::External);
    PipeflowSetup setup;
    setupPipeflow(hemocell, setup);
    hemocell.loadParticles();
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, setup.drivingForce);

    iterate(hemocell, setup, 20);
    hemocell.writeOutput();
    iterate(hemocell, setup, 10);
    hemocell.saveCheckPoint();

    savedIteration = hemocell.iter;
    savedLastOutput = hemocell.getLastOutputIteration();
    savedNumberOfCells = hemocell.cellfields->number_of_cells;
    savedCellCount = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
    checkpointFile = hemo::global.checkpointDirectory + "checkpoint.xml";
  }
  ASSERT_EQ(savedIteration, 30u);
  ASSERT_EQ(savedLastOutput, 20u);
  ASSERT_GT(savedCellCount, 0);

  // The checkpoint stores the values explicitly.
  {
    hemo::Config checkpointCfg(checkpointFile);
    ASSERT_TRUE(checkpointCfg.checkpointed);
  }

  std::vector<char> checkpointPath(checkpointFile.begin(), checkpointFile.end());
  checkpointPath.push_back('\0');
  {
    hemo::HemoCell hemocell(checkpointPath.data(), 0, args, hemo::HemoCell::MPIHandle::External);
    ASSERT_TRUE(hemocell.cfg->checkpointed);
    PipeflowSetup setup;
    setupPipeflow(hemocell, setup);
    hemocell.loadCheckPoint();

    EXPECT_EQ(hemocell.iter, savedIteration);
    // Previously reset to 0, so the next output added `iter` to every residence time.
    EXPECT_EQ(hemocell.getLastOutputIteration(), savedLastOutput);
    // Restored from the checkpoint instead of re-derived from the .pos files.
    EXPECT_TRUE(hemocell.cellfields->number_of_cells_from_checkpoint);
    EXPECT_EQ(hemocell.cellfields->number_of_cells, savedNumberOfCells);
    EXPECT_EQ(hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell), savedCellCount);

    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, setup.drivingForce);
    iterate(hemocell, setup, 10);
    EXPECT_EQ(hemocell.iter, savedIteration + 10);
    EXPECT_EQ(hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell), savedCellCount);
  }
}
