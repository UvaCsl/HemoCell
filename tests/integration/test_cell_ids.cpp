#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"

// Cell IDs grow by a multiple of the number of cells at every periodic
// crossing and are never folded back, so long runs reach IDs beyond the int
// range. Shifting all IDs past 2^32 at the start must leave the trajectory
// unchanged (IDs are only keys, and a common shift keeps their order). Run on
// 1 and 2 ranks (unit, unit_mpi): the envelope synchronisation exchanges IDs.

using namespace pipeflow_test;

namespace {

struct Result {
  State state;
  pluint cells = 0;
  plb::plint shift = 0;
};

// With `shiftIds`, add a multiple of the number of cells larger than 2^32 to
// every cell ID right after loading, as after that many periodic laps.
Result run(bool shiftIds) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell((char *)pipeflow_config, 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.loadParticles();
  Result result;
  const plb::plint cells = hemocell.cellfields->number_of_cells;
  EXPECT_GT(cells, 0);
  if (shiftIds && cells > 0) {
    result.shift = ((plb::plint(1) << 32) / cells + 1) * cells;
    plb::MultiParticleField3D<hemo::HemoCellParticleField> & field = *hemocell.cellfields->immersedParticles;
    for (plb::plint bid : field.getLocalInfo().getBlocks()) {
      hemo::HemoCellParticleField & block = field.getComponent(bid);
      for (hemo::HemoCellParticle & p : block.particles) {
        const plb::plint base = hemocell.cellfields->base_cell_id(p.sv.cellId);
        p.sv.cellId += result.shift;
        EXPECT_EQ(hemocell.cellfields->base_cell_id(p.sv.cellId), base);
      }
      block.invalidate_ppc();
      block.invalidate_lpc();
    }
  }
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 30);
  result.cells = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
  result.state = captureState(hemocell, result.shift);
  return result;
}

} // namespace

TEST(CellIds, IdsBeyondIntRangeGiveTheSameTrajectory) {
  const Result normal = run(false);
  const Result shifted = run(true);
  ASSERT_GT(normal.cells, 0u);
  ASSERT_GT(shifted.shift, plb::plint(1) << 32);
  EXPECT_EQ(shifted.cells, normal.cells);
  EXPECT_EQ(firstDifference(normal.state.particles, shifted.state.particles), -1);
  EXPECT_EQ(firstDifference(normal.state.fluid, shifted.state.fluid), -1);
}
