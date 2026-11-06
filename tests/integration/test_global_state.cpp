#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"
#include "parallelism/mpiManager.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

// HemoCell keeps process-wide state (the Parameters, `global`, static
// pointers and statistics). Destroying the HemoCell object must put it back
// to its initial values, so a next simulation in the same process starts as
// in a new process. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

namespace {

// The pipeflow test configuration with cellsDeletedInfo enabled.
std::string writeConfig() {
  return hemo_test::testConfig(pipeflow_config, [](std::string & xml) {
    hemo_test::replaceOnce(xml, "<hemocell>", "<hemocell>\n<verbose><cellsDeletedInfo>1</cellsDeletedInfo></verbose>");
  });
}

} // namespace

TEST(GlobalState, IsResetWhenTheSimulationEnds) {
  const std::string config = writeConfig();
  std::string outputDirectory;
  {
    char * args[] = {(char *)"test", NULL};
    hemo::HemoCell hemocell(config, 0, args, hemo::HemoCell::MPIHandle::External);
    outputDirectory = plb::global::directories().getOutputDir();
    EXPECT_TRUE(hemo::global.hemoCellInitialized);
    EXPECT_TRUE(hemo::global.cellsDeletedInfo);
    PipeflowSetup setup;
    setupPipeflow(hemocell, setup);
    hemocell.loadParticles();
    iterate(hemocell, setup, 2);
    hemo::CellInformationFunctionals::calculateCellVolume(&hemocell);
    EXPECT_GT(hemo::param::dx, 0.);
    EXPECT_GT(hemo::param::pipe_radius, 0.);
    EXPECT_EQ(hemo::HemoCellParticleField::cellFields, hemocell.cellfields);
    EXPECT_NE(hemo::global.checkpointDirectory, "./checkpoint/");
  }

  EXPECT_FALSE(hemo::global.hemoCellInitialized);
  EXPECT_FALSE(hemo::global.cellsDeletedInfo);
  EXPECT_EQ(hemo::global.checkpointDirectory, "./checkpoint/");
  EXPECT_EQ(hemo::global.statistics.elapsed().count(), 0);
  EXPECT_EQ(hemo::param::dx, 0.);
  EXPECT_EQ(hemo::param::dt, 0.);
  EXPECT_EQ(hemo::param::df, 0.);
  EXPECT_EQ(hemo::param::tau, 0.);
  EXPECT_EQ(hemo::param::nu_lbm, 0.);
  EXPECT_EQ(hemo::param::u_lbm_max, 0.);
  EXPECT_EQ(hemo::param::pipe_radius, 0.);
  EXPECT_EQ(hemo::param::f_limit, 0.);
  EXPECT_EQ(hemo::HemoCellParticleField::cellFields, nullptr);
  EXPECT_TRUE(hemo::CellInformationFunctionals::info_per_cell.empty());

  plb::global::mpi().barrier();
  if (plb::global::mpi().isMainProcessor()) { std::filesystem::remove_all(outputDirectory); }
  std::remove(config.c_str());
}
