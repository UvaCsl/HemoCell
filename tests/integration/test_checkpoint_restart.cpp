#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"

#include <memory>
#include <string>
#include <vector>

// Checkpoint -> restart round trip on the small pipeflow validation case.
// Asserts that the state which cannot be re-derived from the input files is
// written to and restored from the checkpoint.

using namespace pipeflow_test;

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
