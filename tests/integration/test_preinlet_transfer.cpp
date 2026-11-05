#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "preInlet.h"
#include "parallelism/mpiManager.h"

#include <mpi.h>

#include <fstream>
#include <string>

// The pre-inlet passes its cells into the main domain: a cell placed in the
// pre-inlet is copied to the domain ranks at the first applyPreInlet() (the
// pre-inlet is shorter than the transfer region), with its cell ID and shifted
// along the flow direction by the pre-inlet length. Set up as in
// examples/pipeflow_with_preinlet. Needs 2 ranks (unit_mpi): one for the
// pre-inlet, one for the domain.

using namespace pipeflow_test;

namespace {

void copyFile(const std::string & from, const std::string & to) {
  std::ifstream in(from);
  std::ofstream out(to);
  out << in.rdbuf();
}

} // namespace

TEST(PreInlet, CellsAreTransferredIntoTheDomain) {
  const int size = plb::global::mpi().getSize();
  if (size < 2) { GTEST_SKIP() << "a pre-inlet needs at least 2 MPI processes"; }

  // One RBC in the middle of the pre-inlet (positions are relative to the
  // pre-inlet's lower corner): 5.5 um along x, on the pipe axis.
  const std::string rbc = hemo_test::outputRoot() + "preinlet_RBC";
  const std::string plt = hemo_test::outputRoot() + "preinlet_PLT";
  if (plb::global::mpi().isMainProcessor()) {
    copyFile(std::string(rbc_name) + ".xml", rbc + ".xml");
    copyFile(std::string(plt_name) + ".xml", plt + ".xml");
    std::ofstream(rbc + ".pos") << "1\n5.5 13 13 0 0 0\n";
    std::ofstream(plt + ".pos") << "0\n";
  }
  plb::global::mpi().barrier();
  const std::string config = hemo_test::testConfig(pipeflow_config, [](std::string & xml) {
    hemo_test::replaceOnce(
        xml, "<hemocell>",
        "<hemocell>\n<preInlet><parameters><lengthN> 20 </lengthN><Re> 0.5 </Re></parameters></preInlet>");
  });

  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(config, 0, args, hemo::HemoCell::MPIHandle::External);
  hemo::Config & cfg = *hemocell.cfg;
  std::unique_ptr<plb::VoxelizedDomain3D<T>> voxelizedDomain;
  std::unique_ptr<plb::MultiScalarField3D<int>> flagMatrix;
  hemo::getFlagMatrixFromSTL(pipeflow_geometry, cfg["domain"]["fluidEnvelope"].read<int>(),
                             cfg["domain"]["refDirN"].read<int>(), cfg["domain"]["refDir"].read<int>(), voxelizedDomain,
                             flagMatrix, cfg["domain"]["blockSize"].read<int>(),
                             cfg["domain"]["particleEnvelope"].read<int>());
  hemo::param::lbm_base_parameters(cfg);
  hemocell.preInlet = new hemo::PreInlet(&hemocell, flagMatrix.get());
  plb::Box3D slice = flagMatrix->getBoundingBox();
  slice.x0 = slice.x1;
  hemocell.preInlet->preInletFromSlice(Direction::Xpos, slice);
  hemocell.initializeLattice(voxelizedDomain->getMultiBlockManagement());
  if (!hemocell.partOfpreInlet) { hemocell.lattice->periodicity().toggleAll(false); }
  hemocell.preInlet->initializePreInlet();
  hemo::boundaryFromFlagMatrix(hemocell.lattice, flagMatrix.get(), hemocell.partOfpreInlet);
  hemocell.preInlet->createBoundary();
  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.latticeEquilibrium(1., {0., 0., 0.});
  hemocell.preInlet->calculateDrivingForce();
  hemocell.lattice->initialize();
  hemocell.initializeCellfield();
  hemocell.addCellType<hemo::RbcHighOrderModel>(rbc, hemo::RBC_FROM_SPHERE);
  hemocell.addCellType<hemo::PltSimpleModel>(plt, hemo::ELLIPSOID_FROM_SPHERE);
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(1);
  if (!hemocell.partOfpreInlet) {
    // A plane, as addPressureBoundary0N() requires (the example passes three
    // layers, which only a build without PLB_DEBUG accepts).
    plb::Box3D bb = hemocell.lattice->getBoundingBox();
    plb::Box3D outlet(bb.x0, bb.x0, bb.y0, bb.y1, bb.z0, bb.z1);
    plb::OnLatticeBoundaryCondition3D<T, DESCRIPTOR> * boundary =
        new plb::BoundaryConditionInstantiator3D<T, DESCRIPTOR, plb::WrappedZouHeBoundaryManager3D<T, DESCRIPTOR>>();
    boundary->addPressureBoundary0N(outlet, *hemocell.lattice, plb::boundary::density);
    setBoundaryDensity(*hemocell.lattice, outlet, 1.0);
    delete boundary;
  }
  hemocell.loadParticles();

  // Particles on the ranks of each part before and after one step.
  auto vertices = [&hemocell](bool preInletPart) {
    long local = 0;
    if (hemocell.partOfpreInlet == preInletPart) {
      for (plb::plint bid : hemocell.cellfields->immersedParticles->getLocalInfo().getBlocks()) {
        local += hemocell.cellfields->immersedParticles->getComponent(bid).particles.size();
      }
    }
    long total = 0;
    MPI_Allreduce(&local, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return total;
  };
  const long inPreInlet = vertices(true);
  EXPECT_GT(inPreInlet, 0) << "the cell was not placed in the pre-inlet";
  EXPECT_EQ(vertices(false), 0);

  hemocell.iterate();
  if (hemocell.partOfpreInlet) { hemocell.preInlet->setDrivingForce(); }
  hemocell.preInlet->applyPreInlet();

  // The domain now holds copies (with the cell's ID 0) within its particle field.
  long local = 0, wrongId = 0;
  if (!hemocell.partOfpreInlet) {
    for (plb::plint bid : hemocell.cellfields->immersedParticles->getLocalInfo().getBlocks()) {
      for (const hemo::HemoCellParticle & p : hemocell.cellfields->immersedParticles->getComponent(bid).particles) {
        local++;
        if (hemocell.cellfields->base_cell_id(p.sv.cellId) != 0) { wrongId++; }
      }
    }
  }
  long counts[2] = {local, wrongId};
  MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  EXPECT_GT(counts[0], 0) << "no particle reached the domain";
  EXPECT_EQ(counts[1], 0);
}
