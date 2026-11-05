#include "gtest/gtest.h"
#include "pipeflow_setup.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>

using namespace pipeflow_test;

namespace {

// Largest |force_repulsion| of all particles, over all ranks.
T maxRepulsionForce(hemo::HemoCell & hemocell) {
  T localMax = 0;
  for (plint blockId : hemocell.cellfields->immersedParticles->getLocalInfo().getBlocks()) {
    for (const hemo::HemoCellParticle & particle :
         hemocell.cellfields->immersedParticles->getComponent(blockId).particles) {
      const hemo::Array<T, 3> & f = particle.sv.force_repulsion;
      localMax = std::max(localMax, std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]));
    }
  }
  T globalMax = 0;
  MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return globalMax;
}

} // namespace

// Wall repulsion without cell-cell repulsion used to add to force_repulsion at
// every update without ever resetting it, so the wall force kept growing (and
// stayed applied after a vertex moved away from the wall).
TEST(Repulsion, WallRepulsionAloneDoesNotAccumulate) {
  static_assert(sizeof(T) == sizeof(double), "MPI_DOUBLE below assumes T == double");
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup, 0.0); // allow vertices next to the wall
  hemocell.loadParticles();

  const unsigned int interval = (*hemocell.cfg)["ibm"]["stepParticleEvery"].read<unsigned int>();
  hemocell.enableBoundaryParticles(5.0, 1.0, interval); // k 5 pN, cutoff 1 um
  ASSERT_FALSE(hemocell.repulsionEnabled);

  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 2 * interval); // includes wall-repulsion updates

  // The update iterate() performs, applied twice at the same positions: the
  // result must not depend on how often it was applied.
  hemocell.cellfields->applyBoundaryRepulsionForce(!hemocell.repulsionEnabled);
  const T once = maxRepulsionForce(hemocell);
  hemocell.cellfields->applyBoundaryRepulsionForce(!hemocell.repulsionEnabled);
  const T twice = maxRepulsionForce(hemocell);

  ASSERT_GT(once, 0.0); // some vertices are within the cutoff of the wall
  EXPECT_DOUBLE_EQ(twice, once);
}

#ifdef FORCE_LIMIT
// The repulsion constant is given in pN, and the repulsion force of a vertex
// (particle and wall repulsion together) is capped at FORCE_LIMIT: without
// the cap, k * r_c / d grows without bound as a vertex approaches.
TEST(Repulsion, ConstantInPicoNewtonAndForceCapped) {
  static_assert(sizeof(T) == sizeof(double), "MPI_DOUBLE below assumes T == double");
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup, 0.0); // allow vertices next to the wall
  hemocell.loadParticles();

  const unsigned int interval = (*hemocell.cfg)["ibm"]["stepParticleEvery"].read<unsigned int>();
  hemocell.setRepulsion(5.0, 0.5); // 5 pN, 0.5 um
  EXPECT_DOUBLE_EQ(hemocell.cellfields->repulsionConstant, 5.0e-12 / hemo::param::df);
  hemocell.setRepulsionTimeScaleSeperation(interval);
  hemocell.enableBoundaryParticles(1e6, 1.0, interval); // 1e6 pN: far above the cap within the cutoff
  EXPECT_DOUBLE_EQ(hemocell.cellfields->boundaryRepulsionConstant, 1e6 * 1e-12 / hemo::param::df);

  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 1); // the first step updates (and caps) the repulsion

  const T maxForce = maxRepulsionForce(hemocell);
  EXPECT_NEAR(maxForce, hemo::param::f_limit, 1e-12 * hemo::param::f_limit);
}
#endif
