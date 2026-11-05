#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "parallelism/mpiManager.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

// HemoCell's start-up checks only require the settings the library uses:
// <sim><tmax> is read by the drivers that use it and is optional for
// sanityCheck(). Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

TEST(SanityCheck, DoesNotNeedTmax) {
  std::ifstream in(pipeflow_config);
  std::stringstream content;
  content << in.rdbuf();
  std::string xml = content.str();
  const size_t begin = xml.find("<tmax>");
  const size_t end = xml.find("</tmax>");
  ASSERT_NE(begin, std::string::npos);
  xml.erase(begin, end + 7 - begin);
  const std::string path =
      testing::TempDir() + "hemocell_no_tmax_r" + std::to_string(plb::global::mpi().getRank()) + ".xml";
  {
    std::ofstream out(path);
    out << xml;
  }

  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(path, 0, args, hemo::HemoCell::MPIHandle::External);
  EXPECT_THROW((*hemocell.cfg)["sim"]["tmax"], std::invalid_argument);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.loadParticles();
  EXPECT_NO_THROW(hemocell.iterate()); // the first iterate() runs sanityCheck()
  EXPECT_EQ(hemocell.iter, 1u);
  std::remove(path.c_str());
}
