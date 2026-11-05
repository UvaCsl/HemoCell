#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"
#include "fluidInfo.h"

#include <cmath>
#include <limits>

// A vertex with a non-finite force must not be spread into the fluid: its
// cell is removed and the lattice stays finite. Run on 1 and 2 ranks.

using namespace pipeflow_test;

TEST(NonFiniteForce, CellIsRemovedAndFluidStaysFinite) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell((char *)pipeflow_config, 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.loadParticles();
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 10);
  const int before = hemo::CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, rbc_name);
  ASSERT_GT(before, 0);

  // Corrupt the force of one vertex of one RBC, on every block holding a copy
  // of it, right before its force is spread.
  hemo::HemoCellFields & fields = *hemocell.cellfields;
  const unsigned int rbcType = fields[rbc_name]->ctype;
  int target = -1;
  for (plb::plint bid : fields.immersedParticles->getLocalInfo().getBlocks()) {
    for (const hemo::HemoCellParticle & p : fields.immersedParticles->getComponent(bid).particles) {
      if (p.sv.celltype == rbcType) { target = std::max(target, fields.base_cell_id(p.sv.cellId)); }
    }
  }
  int globalTarget = target;
  MPI_Allreduce(&target, &globalTarget, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  for (plb::plint bid : fields.immersedParticles->getLocalInfo().getBlocks()) {
    for (hemo::HemoCellParticle & p : fields.immersedParticles->getComponent(bid).particles) {
      if (fields.base_cell_id(p.sv.cellId) == globalTarget && p.sv.vertexId == 0) {
        p.sv.force[0] = std::numeric_limits<T>::quiet_NaN();
      }
    }
  }
  iterate(hemocell, setup, 1);

  const hemo::FluidStatistics fluid = hemo::FluidInfo::calculateVelocityStatistics(&hemocell);
  EXPECT_TRUE(std::isfinite(fluid.max));
  EXPECT_TRUE(std::isfinite(fluid.avg));
  // Counting deletes incomplete cells, the corrupted cell among them.
  EXPECT_EQ(hemo::CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, rbc_name), before - 1);
}
