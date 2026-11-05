#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"

#include <sys/stat.h>

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Checkpoint -> restart round trip on the small pipeflow validation case.
// Asserts that the state which cannot be re-derived from the input files is
// written to and restored from the checkpoint.

using namespace pipeflow_test;

namespace {

// Modification time of a file in nanoseconds (0 if it does not exist).
long long modificationTime(const std::string & path) {
  struct stat st;
  if (stat(path.c_str(), &st) != 0) { return 0; }
  return (long long)st.st_mtim.tv_sec * 1000000000LL + st.st_mtim.tv_nsec;
}

// Run the pipeflow case for 10 iterations and save a checkpoint; returns the
// path of its checkpoint.xml.
std::string writeCheckpoint() {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.loadParticles();
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 10);
  hemocell.saveCheckPoint();
  return hemo::global.checkpointDirectory + "checkpoint.xml";
}

// Replace `from` by `to` in a file (rank 0), then synchronise.
void replaceInFile(const std::string & path, const std::string & from, const std::string & to) {
  if (plb::global::mpi().isMainProcessor()) {
    std::ifstream in(path.c_str());
    std::stringstream text;
    text << in.rdbuf();
    std::string content = text.str();
    const std::string::size_type at = content.find(from);
    ASSERT_NE(at, std::string::npos) << from;
    content.replace(at, from.size(), to);
    std::ofstream out(path.c_str());
    out << content;
  }
  plb::global::mpi().barrier();
}

// Restarting from `checkpointFile` must end in a fatal error mentioning `expected`.
void expectRestartRejected(const std::string & checkpointFile, const std::string & expected) {
  char * args[] = {(char *)"test", NULL};
  std::vector<char> path(checkpointFile.begin(), checkpointFile.end());
  path.push_back('\0');
  hemo::FatalErrorsThrow guard;
  hemo::HemoCell hemocell(path.data(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  try {
    hemocell.loadCheckPoint();
    ADD_FAILURE() << "checkpoint accepted";
  } catch (const hemo::FatalError & e) {
    EXPECT_NE(std::string(e.what()).find(expected), std::string::npos) << e.what();
  }
}

} // namespace

TEST(Checkpoint, RestartRestoresState) {
  char * args[] = {(char *)"test", NULL};

  unsigned int savedIteration = 0;
  unsigned int savedLastOutput = 0;
  int savedNumberOfCells = 0;
  int savedCellCount = 0;
  std::string checkpointFile;

  {
    hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
    PipeflowSetup setup;
    setupPipeflow(hemocell, setup);
    hemocell.loadParticles();
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, setup.drivingForce);

    iterate(hemocell, setup, 20);
    hemocell.writeOutput();
    iterate(hemocell, setup, 10);
    hemocell.saveCheckPoint();

    // checkpoint.xml is written last, under a temporary name.
    const std::string dir = hemo::global.checkpointDirectory;
    EXPECT_GE(modificationTime(dir + "checkpoint.xml"), modificationTime(dir + "lattice.dat"));
    EXPECT_GE(modificationTime(dir + "checkpoint.xml"), modificationTime(dir + "particleField.dat"));
    EXPECT_EQ(modificationTime(dir + "checkpoint.xml.tmp"), 0);

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

TEST(Checkpoint, MismatchingCellTypesAreRejected) {
  const std::string checkpointFile = writeCheckpoint();
  replaceInFile(checkpointFile, "RBC:", "WBC:");
  expectRestartRejected(checkpointFile, "cell types");
}

TEST(Checkpoint, MismatchingParticleLayoutIsRejected) {
  const std::string checkpointFile = writeCheckpoint();
  replaceInFile(checkpointFile, "record=", "record=1");
  expectRestartRejected(checkpointFile, "particle records");
}

// A checkpoint records the compile-time configuration it was written with
// (config/constant_defaults.h), for the record; a restart only warns if it
// differs (the particle layout is checked separately).
TEST(Checkpoint, RecordsBuildConfiguration) {
  const std::string checkpointFile = writeCheckpoint();
  std::ifstream in(checkpointFile.c_str());
  std::stringstream text;
  text << in.rdbuf();
  EXPECT_NE(text.str().find("<BuildConfiguration>"), std::string::npos);
  EXPECT_NE(text.str().find(HEMOCELL_BUILD_CONFIGURATION), std::string::npos);
}
