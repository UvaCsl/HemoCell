#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"

// Output and statistics must not change the simulation: a run that writes
// output and takes statistics every few steps ends bit for bit in the same
// state as one that does not. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

namespace {

// Run the pipeflow case for `steps` iterations; with `readEvery` > 0, write
// output and take statistics after every `readEvery` iterations.
State run(unsigned int steps, unsigned int readEvery) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.loadParticles();
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  for (unsigned int step = 1; step <= steps; step++) {
    iterate(hemocell, setup, 1);
    if (readEvery > 0 && step % readEvery == 0) {
      EXPECT_GT(hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell), 0u);
      EXPECT_GT(hemo::CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, rbc_name), 0u);
      hemo::CellInformationFunctionals::calculateCellVolume(&hemocell);
      hemo::CellInformationFunctionals::calculateCellStretch(&hemocell);
      hemo::CellInformationFunctionals::clear_list();
      hemo::FluidInfo::calculateVelocityStatistics(&hemocell);
      hemo::ParticleInfo::calculateForceStatistics(&hemocell);
      hemocell.writeOutput();
    }
  }
  return captureState(hemocell);
}

} // namespace

TEST(Output, DoesNotChangeTrajectory) {
  // stepParticleEvery is 2 in this configuration: reading after step 7 and 21
  // falls between particle updates, after 14 and 28 on one.
  const State quiet = run(30, 0);
  const State busy = run(30, 7);
  ASSERT_FALSE(quiet.particles.empty());
  ASSERT_FALSE(quiet.fluid.empty());
  EXPECT_EQ(firstDifference(quiet.particles, busy.particles), -1)
      << quiet.particles.size() << " vs " << busy.particles.size() << " particle values";
  EXPECT_EQ(firstDifference(quiet.fluid, busy.fluid), -1);
}
