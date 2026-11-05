#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "error.h"

#include <fstream>
#include <sstream>
#include <string>

// The compile-time configuration (config/constant_defaults.h) of a case must
// match that of the library it links; HemoCell checks it once at
// construction. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

TEST(BuildConfiguration, MismatchWithTheLibraryIsFatal) {
  const std::string other = std::string(HEMOCELL_BUILD_CONFIGURATION) + ",INTERIOR_VISCOSITY";
  char * args[] = {(char *)"test", NULL};
  hemo::FatalErrorsThrow guard;
  try {
    hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External, other.c_str());
    FAIL() << "a case with another configuration was accepted";
  } catch (const hemo::FatalError & e) {
    const std::string message = e.what();
    EXPECT_NE(message.find(other), std::string::npos) << message;
    EXPECT_NE(message.find("library variant"), std::string::npos) << message;
  }
  EXPECT_FALSE(hemo::global.hemoCellInitialized); // nothing was set up
}

TEST(BuildConfiguration, MatchingConfigurationIsAccepted) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  EXPECT_TRUE(hemo::global.hemoCellInitialized);
}
