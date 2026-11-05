#include "gtest/gtest.h"
#include "palabos3D.h"

// The test runner for all tests defined under `hemocell/tests/`. The code
// evalutaes all tests sequentially. The custom main file is provided rather
// than linking with `gtest_main`'s default main file to ensure we have explicit
// control on the MPI environment.
//
// Each initialisation of the HemoCell class invokes the `plb::plbInit` function
// which performs `MPI_Init` on its construction and `MPI_Finalize` on
// destruction. However, when running multiple tests we cannot invoke these
// functions multiples times and MPI will terminate if done so anyways.
//
// Therefore, we initialse the global MPI manager before running all the tests,
// ensuring that the MPI environment is present, afterwhich all tests are
// evaluated. In each tests, it is important to initialise HemoCell by passing
// the `ManagedMPI::External` flag, to indicate we (i.e. gtest) is managing the
// MPI instance and `plb::plbInit` should not be called.
//
// Finally, when leaving this function, the MPI environment will be teared down,
// invoking `MPI_Finalize` only once.
//
// The simulations of the tests write into one directory per run in the
// GoogleTest temporary directory (hemo_test::outputRoot(), see
// test_support.h), removed afterwards unless HEMOCELL_KEEP_TEST_OUTPUT is set.
#include "test_support.h"

#include <cstdlib>
#include <filesystem>
#include <unistd.h>

int argc;
char ** argv;

int main(int t_argc, char * t_argv[]) {
  argc = t_argc;
  argv = t_argv;
  ::testing::InitGoogleTest(&argc, argv);

  plb::plbInit(&(argc), &(argv));

  // One directory per run (the process ID of rank 0), the same on all ranks.
  std::string root = testing::TempDir() + "hemocell_tests_" + std::to_string(getpid()) + "/";
  plb::global::mpi().bCast(root);
  hemo_test::outputRoot() = root;
  if (plb::global::mpi().isMainProcessor()) { std::filesystem::create_directories(root); }
  plb::global::mpi().barrier();

  const int result = RUN_ALL_TESTS();

  plb::global::mpi().barrier();
  if (plb::global::mpi().isMainProcessor() && !std::getenv("HEMOCELL_KEEP_TEST_OUTPUT")) {
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
  }
  return result;
}
