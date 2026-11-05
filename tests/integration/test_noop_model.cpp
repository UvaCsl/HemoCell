#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "NoOp.h"
#include "parallelism/mpiManager.h"

// The null model (NoOp) can be registered like any other model and exerts no
// membrane forces: with repulsion disabled, the vertices of its cells carry
// no force, while those of the other cell type (PLT) do. Run on 1 and 2 ranks
// (unit, unit_mpi).

using namespace pipeflow_test;

TEST(Mechanics, NoOpCellsCarryNoMembraneForce) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflow_config, 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow<hemo::NoOp>(hemocell, setup);
  hemocell.loadParticles();
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 10);

  const pluint rbc = (*hemocell.cellfields)[rbc_name]->ctype;
  long counts[3] = {0, 0, 0}; // NoOp vertices, of which with a force, PLT vertices with a force
  plb::MultiParticleField3D<hemo::HemoCellParticleField> & field = *hemocell.cellfields->immersedParticles;
  for (plb::plint bid : field.getLocalInfo().getBlocks()) {
    for (const hemo::HemoCellParticle & p : field.getComponent(bid).particles) {
      const bool hasForce = p.sv.force[0] != 0 || p.sv.force[1] != 0 || p.sv.force[2] != 0;
      if (p.sv.celltype == rbc) {
        counts[0]++;
        counts[1] += hasForce;
      } else {
        counts[2] += hasForce;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  EXPECT_GT(counts[0], 0);
  EXPECT_EQ(counts[1], 0);
  EXPECT_GT(counts[2], 0);
}
