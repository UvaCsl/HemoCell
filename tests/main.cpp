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
int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    plb::plbInit(&(argc), &(argv));
    return RUN_ALL_TESTS();
}
