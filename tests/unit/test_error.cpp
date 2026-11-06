#include "gtest/gtest.h"

#include "config.h"
#include "error.h"
#include "parallelism/mpiManager.h"
#include "test_support.h"

#include <fstream>
#include <string>

// fatal() ends the MPI job; under a FatalErrorsThrow guard it throws instead.

TEST(Error, FatalThrowsUnderGuardWithStreamedMessage) {
  hemo::FatalErrorsThrow guard;
  try {
    hemo::fatal("cell type ", 3, " of ", 2.5, " does not exist");
    FAIL() << "fatal() returned";
  } catch (const hemo::FatalError & e) { EXPECT_EQ(std::string(e.what()), "cell type 3 of 2.5 does not exist"); }
}

TEST(Error, GuardsNest) {
  hemo::FatalErrorsThrow outer;
  { hemo::FatalErrorsThrow inner; }
  // Still throwing after the inner guard is gone.
  EXPECT_THROW(hemo::fatal("still guarded"), hemo::FatalError);
}

TEST(Error, MissingConfigFileIsFatal) {
  hemo::FatalErrorsThrow guard;
  EXPECT_THROW(hemo::Config("/nonexistent/hemocell/config.xml"), hemo::FatalError);
}

TEST(Error, MalformedConfigFileIsFatal) {
  hemo::FatalErrorsThrow guard;
  // The rank is part of the name, so concurrent ranks do not write the same file.
  const std::string path = hemo_test::outputRoot() + "hemocell_error_test_malformed_r" +
                           std::to_string(plb::global::mpi().getRank()) + ".xml";
  {
    std::ofstream file(path);
    file << "<hemocell><domain></hemocell>";
  }
  try {
    hemo::Config cfg(path);
    FAIL() << "malformed XML was accepted";
  } catch (const hemo::FatalError & e) { EXPECT_NE(std::string(e.what()).find(path), std::string::npos); }
}
