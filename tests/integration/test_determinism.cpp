#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "parallelism/mpiManager.h"

#include <fstream>
#include <iterator>
#include <string>
#include <vector>

// Two identical runs (same configuration, ranks and decomposition) write
// byte-identical checkpoint data: the particle exchanges unpack in a fixed
// order, not in MPI arrival order. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

namespace {

// Run the pipeflow case for 40 iterations with output every 10, save a
// checkpoint and return its directory.
std::string runAndCheckpoint() {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.loadParticles();
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  for (int i = 0; i < 4; i++) {
    iterate(hemocell, setup, 10);
    hemocell.writeOutput();
  }
  hemocell.saveCheckPoint();
  return hemo::global.checkpointDirectory;
}

std::vector<char> contents(const std::string & path) {
  std::ifstream in(path, std::ios::binary);
  return std::vector<char>(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

} // namespace

TEST(Determinism, IdenticalRunsWriteIdenticalCheckpoints) {
  const std::string first = runAndCheckpoint();
  const std::string second = runAndCheckpoint();
  ASSERT_NE(first, second);
  plb::global::mpi().barrier();
  if (plb::global::mpi().isMainProcessor()) {
    for (const char * file : {"particleField.dat", "lattice.dat"}) {
      const std::vector<char> a = contents(first + file);
      const std::vector<char> b = contents(second + file);
      ASSERT_FALSE(a.empty()) << first + file;
      EXPECT_TRUE(a == b) << file << " differs between " << first << " and " << second;
    }
  }
}
