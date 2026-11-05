#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "preInlet.h"
#include "parallelism/mpiManager.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

// With a pre-inlet, the first ranks hold the pre-inlet and the others the
// domain. When the configuration splits the domain into more atomic blocks
// (mABx * mABy * mABz) than it has ranks, every block must still belong to a
// domain rank; blocks without an explicit owner used to get the rank equal to
// their block ID (a pre-inlet rank or no rank at all). Needs 2 ranks (unit_mpi).

using namespace pipeflow_test;

namespace {

void replace(std::string & text, const std::string & from, const std::string & to) {
  const size_t at = text.find(from);
  ASSERT_NE(at, std::string::npos) << from;
  text.replace(at, from.size(), to);
}

// The pipeflow test configuration with a pre-inlet and the domain split into
// three blocks along x, writing its output to the GoogleTest temporary directory.
std::string writeConfig() {
  std::ifstream in(pipeflow_config);
  std::stringstream content;
  content << in.rdbuf();
  std::string xml = content.str();
  replace(xml, "<hemocell>",
          "<hemocell>\n<preInlet><parameters><lengthN> 20 </lengthN><Re> 0.5 </Re></parameters></preInlet>");
  replace(xml, "<domain>", "<domain>\n<mABx> 3 </mABx><mABy> 1 </mABy><mABz> 1 </mABz>");
  replace(xml, "<outputDirectory>tmp</outputDirectory>",
          "<outputDirectory>" + testing::TempDir() + "hemocell_preinlet_blocks</outputDirectory>");
  const std::string path =
      testing::TempDir() + "hemocell_preinlet_blocks_r" + std::to_string(plb::global::mpi().getRank()) + ".xml";
  std::ofstream out(path);
  out << xml;
  return path;
}

// Number of blocks of `management` that are not owned by a rank in [first, end).
int blocksOutside(plb::MultiBlockManagement3D & management, int first, int end) {
  int outside = 0;
  for (const auto & bulk : management.getSparseBlockStructure().getBulks()) {
    const int rank = management.getThreadAttribution().getMpiProcess(bulk.first);
    if (rank < first || rank >= end) { outside++; }
  }
  return outside;
}

} // namespace

TEST(PreInlet, EveryBlockBelongsToARankOfItsLattice) {
  const int size = plb::global::mpi().getSize();
  if (size < 2) { GTEST_SKIP() << "a pre-inlet needs at least 2 MPI processes"; }
  const std::string config = writeConfig();
  std::string outputDirectory;
  {
    char * args[] = {(char *)"test", NULL};
    hemo::HemoCell hemocell(config, 0, args, hemo::HemoCell::MPIHandle::External);
    outputDirectory = plb::global::directories().getOutputDir();
    hemo::Config & cfg = *hemocell.cfg;
    PipeflowSetup setup;
    hemo::getFlagMatrixFromSTL(pipeflow_geometry, cfg["domain"]["fluidEnvelope"].read<int>(),
                               cfg["domain"]["refDirN"].read<int>(), cfg["domain"]["refDir"].read<int>(),
                               setup.voxelizedDomain, setup.flagMatrix, cfg["domain"]["blockSize"].read<int>(),
                               cfg["domain"]["particleEnvelope"].read<int>());
    hemo::param::lbm_base_parameters(cfg);
    hemocell.preInlet = new hemo::PreInlet(&hemocell, setup.flagMatrix.get());
    plb::Box3D slice = setup.flagMatrix->getBoundingBox();
    slice.x0 = slice.x1;
    hemocell.preInlet->preInletFromSlice(Direction::Xpos, slice);
    hemocell.initializeLattice(setup.voxelizedDomain->getMultiBlockManagement());

    const int preInletRanks = hemocell.preInlet->nProcs;
    ASSERT_GE(preInletRanks, 1);
    ASSERT_LT(preInletRanks, size);
    EXPECT_EQ(hemocell.partOfpreInlet, plb::global::mpi().getRank() < preInletRanks);
    EXPECT_EQ(hemocell.domain_lattice_management->getSparseBlockStructure().getNumBlocks(), 3);
    EXPECT_EQ(blocksOutside(*hemocell.domain_lattice_management, preInletRanks, size), 0);
    EXPECT_EQ(blocksOutside(*hemocell.preinlet_lattice_management, 0, preInletRanks), 0);
    // Each rank holds blocks of its own lattice.
    EXPECT_FALSE(hemocell.lattice->getLocalInfo().getBlocks().empty());
  }
  plb::global::mpi().barrier();
  if (plb::global::mpi().isMainProcessor()) { std::filesystem::remove_all(outputDirectory); }
  std::remove(config.c_str());
}
