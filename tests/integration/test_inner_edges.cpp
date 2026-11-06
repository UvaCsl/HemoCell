#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "error.h"
#include "parallelism/mpiManager.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

// The inner edges of a material file are vertex pairs of one particular mesh.
// A pair outside the mesh of the cell type (another shape or minNumTriangles)
// or an unreadable pair is a fatal error when the cell type is added, instead
// of an out-of-bounds access later. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

namespace {

// The PLT material of the pipeflow test (a 66-vertex mesh) with `extraEdge`
// appended to its inner edges, written to the GoogleTest temporary directory.
// Returns the cell type name (the file name without `.xml`).
std::string writeMaterial(const std::string & extraEdge) {
  std::ifstream in(std::string(plt_name) + ".xml");
  std::stringstream content;
  content << in.rdbuf();
  std::string xml = content.str();
  const size_t at = xml.find("</InnerEdges>");
  EXPECT_NE(at, std::string::npos);
  xml.insert(at, "<Edge>" + extraEdge + "</Edge>\n");
  const std::string name =
      hemo_test::outputRoot() + "hemocell_inner_edges_r" + std::to_string(plb::global::mpi().getRank());
  std::ofstream out(name + ".xml");
  out << xml;
  return name;
}

// The message of the fatal error raised when adding a PLT type with `extraEdge`.
std::string addCellTypeError(const std::string & extraEdge) {
  const std::string name = writeMaterial(extraEdge);
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  std::string message;
  {
    hemo::FatalErrorsThrow guard;
    try {
      hemocell.addCellType<hemo::PltSimpleModel>(name, hemo::ELLIPSOID_FROM_SPHERE);
    } catch (const hemo::FatalError & e) { message = e.what(); }
  }
  std::remove((name + ".xml").c_str());
  return message;
}

} // namespace

TEST(CellMaterial, InnerEdgeOutsideTheMeshIsFatal) {
  const std::string message = addCellTypeError(" 3 66 ");
  EXPECT_NE(message.find("inner edge 3 66"), std::string::npos) << message;
  EXPECT_NE(message.find("66-vertex"), std::string::npos) << message;
}

TEST(CellMaterial, UnreadableInnerEdgeIsFatal) {
  const std::string message = addCellTypeError(" 3 ");
  EXPECT_NE(message.find("is not a pair of vertex indices"), std::string::npos) << message;
}
