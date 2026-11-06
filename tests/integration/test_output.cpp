#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"
#include "genericFunctions.h"

#include <hdf5.h>
#include <hdf5_hl.h>
#include <dirent.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// HDF5 and XDMF output of the small pipeflow case: every cell is written once,
// the fluid files hold every lattice node once, and the XDMF files list the
// written blocks. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

namespace {

struct Dataset {
  std::vector<hsize_t> dims;
  std::vector<float> values;
};

Dataset readFloats(hid_t file, const char * name) {
  Dataset d;
  hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) { return d; }
  hid_t space = H5Dget_space(dataset);
  d.dims.resize(H5Sget_simple_extent_ndims(space));
  H5Sget_simple_extent_dims(space, d.dims.data(), NULL);
  hsize_t n = 1;
  for (hsize_t s : d.dims) { n *= s; }
  d.values.resize(n);
  if (n > 0) { H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d.values.data()); }
  H5Sclose(space);
  H5Dclose(dataset);
  return d;
}

// Files in `directory` starting with `prefix`.
std::vector<std::string> listFiles(const std::string & directory, const std::string & prefix) {
  std::vector<std::string> found;
  DIR * dir = opendir(directory.c_str());
  if (!dir) { return found; }
  while (struct dirent * entry = readdir(dir)) {
    const std::string name = entry->d_name;
    if (name.compare(0, prefix.size(), prefix) == 0) { found.push_back(name); }
  }
  closedir(dir);
  return found;
}

std::string readText(const std::string & path) {
  std::ifstream in(path.c_str());
  std::stringstream text;
  text << in.rdbuf();
  return text.str();
}

int countOccurrences(const std::string & text, const std::string & pattern) {
  int n = 0;
  for (std::string::size_type at = text.find(pattern); at != std::string::npos; at = text.find(pattern, at + 1)) {
    n++;
  }
  return n;
}

} // namespace

TEST(Output, CellsOnceFluidBulkAndXdmf) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  hemocell.outputInSiUnits = false; // lattice units: positions in nodes, cell density in vertices per node
  hemocell.setOutputs(rbc_name,
                      {hemo::OUTPUT_POSITION, hemo::OUTPUT_TRIANGLES, hemo::OUTPUT_CELL_ID, hemo::OUTPUT_VERTEX_ID});
  hemocell.setOutputs(plt_name, {hemo::OUTPUT_POSITION, hemo::OUTPUT_TRIANGLES, hemo::OUTPUT_CELL_ID});
  hemocell.setFluidOutputs({hemo::OUTPUT_VELOCITY, hemo::OUTPUT_DENSITY, hemo::OUTPUT_CELL_DENSITY});
  hemocell.loadParticles();
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, 20);
  hemocell.writeOutput();

  // Collective: all ranks.
  const int nRBC = hemo::CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, rbc_name);
  const int nPLT = hemo::CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, plt_name);
  plb::global::mpi().barrier();
  if (plb::global::mpi().getRank() != 0) { return; }
  ASSERT_GT(nRBC, 0);

  const std::map<plb::plint, plb::Box3D> & bulks = hemocell.lattice->getSparseBlockStructure().getBulks();
  hsize_t latticeNodes = 0;
  for (const auto & bulk : bulks) { latticeNodes += bulk.second.nCells(); }

  const std::string directory = plb::global::directories().getOutputDir() + "/";
  const std::string step = hemo::zeroPadNumber(hemocell.iter);
  const std::string stepDirectory = directory + "hdf5/" + step + "/";

  // Every cell once, with all its vertices in vertex order.
  struct CellType {
    const char * name;
    int expectedCells;
    int numVertex;
  } types[] = {{"RBC", nRBC, (*hemocell.cellfields)[rbc_name]->numVertex},
               {"PLT", nPLT, (*hemocell.cellfields)[plt_name]->numVertex}};
  hsize_t rbcVertices = 0;
  for (const CellType & type : types) {
    std::map<int, int> verticesPerCell;
    int blocksWithCells = 0;
    for (const std::string & name : listFiles(stepDirectory, std::string(type.name) + ".")) {
      hid_t file = H5Fopen((stepDirectory + name).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(file, 0) << name;
      Dataset position = readFloats(file, "Position");
      Dataset cellId = readFloats(file, "Cell Id");
      Dataset triangles = readFloats(file, "Triangles");
      H5Fclose(file);
      ASSERT_EQ(position.dims.size(), 2u);
      const hsize_t rows = position.dims[0];
      ASSERT_EQ(rows % type.numVertex, 0u) << name;
      ASSERT_EQ(cellId.values.size(), rows) << name;
      if (rows > 0) { blocksWithCells++; }
      for (hsize_t r = 0; r < rows; r++) {
        verticesPerCell[int(cellId.values[r])]++;
        // Rows are grouped by cell.
        EXPECT_EQ(cellId.values[r], cellId.values[r - r % type.numVertex]) << name;
      }
      for (float vertex : triangles.values) { EXPECT_LT(hsize_t(vertex), rows) << name; }
      { // Cell IDs are stored as 64-bit integers (exact, unlike float).
        hid_t f = H5Fopen((stepDirectory + name).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t dataset = H5Dopen2(f, "Cell Id", H5P_DEFAULT);
        hid_t type = H5Dget_type(dataset);
        EXPECT_EQ(H5Tget_class(type), H5T_INTEGER) << name;
        EXPECT_EQ(H5Tget_size(type), 8u) << name;
        H5Tclose(type);
        H5Dclose(dataset);
        H5Fclose(f);
      }
      if (std::string(type.name) == "RBC") {
        rbcVertices += rows;
        hid_t f = H5Fopen((stepDirectory + name).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        Dataset vertexId = readFloats(f, "Vertex Id");
        H5Fclose(f);
        for (hsize_t r = 0; r < rows; r++) { EXPECT_EQ(int(vertexId.values[r]), int(r % type.numVertex)) << name; }
      }
    }
    EXPECT_EQ(int(verticesPerCell.size()), type.expectedCells) << type.name;
    for (const auto & cell : verticesPerCell) { EXPECT_EQ(cell.second, type.numVertex) << type.name; }

    // The step file lists every block with cells, the series file includes the step.
    const std::string xdmf = readText(directory + type.name + "." + step + ".xmf");
    EXPECT_EQ(countOccurrences(xdmf, "<Grid Name=\"block "), blocksWithCells) << type.name;
    EXPECT_NE(xdmf.find("CollectionType=\"Spatial\""), std::string::npos);
    if (blocksWithCells > 0) {
      EXPECT_NE(xdmf.find("NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">"), std::string::npos) << type.name;
    }
    const std::string series = readText(directory + type.name + ".xmf");
    EXPECT_NE(series.find("href=\"" + std::string(type.name) + "." + step + ".xmf\""), std::string::npos);
  }

  // Fluid: the bulk nodes of every block, each node once.
  hsize_t fluidNodes = 0;
  double cellDensity = 0;
  const std::vector<std::string> fluidFiles = listFiles(stepDirectory, "Fluid.");
  for (const std::string & name : fluidFiles) {
    hid_t file = H5Fopen((stepDirectory + name).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0) << name;
    Dataset velocity = readFloats(file, "Velocity");
    Dataset density = readFloats(file, "CellDensity_RBC");
    H5Fclose(file);
    ASSERT_EQ(velocity.dims.size(), 4u);
    EXPECT_EQ(velocity.dims[3], 3u);
    fluidNodes += velocity.dims[0] * velocity.dims[1] * velocity.dims[2];
    for (float v : density.values) { cellDensity += v; }
  }
  EXPECT_EQ(fluidNodes, latticeNodes);
  EXPECT_EQ(fluidFiles.size(), bulks.size());
  // Every vertex lies in the bulk of exactly one block (a vertex beyond a
  // periodic boundary through its image).
  EXPECT_EQ(cellDensity, double(rbcVertices));

  const std::string fluidXdmf = readText(directory + "Fluid." + step + ".xmf");
  EXPECT_EQ(countOccurrences(fluidXdmf, "<Grid Name=\"block "), int(bulks.size()));
  // Velocity, Density and one CellDensity_<type> per cell type.
  EXPECT_EQ(countOccurrences(fluidXdmf, "Center=\"Cell\""), 4 * int(bulks.size()));
}
