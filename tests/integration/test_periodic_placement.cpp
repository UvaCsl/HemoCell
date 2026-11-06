#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"
#include "parallelism/mpiManager.h"

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <vector>

// A cell whose centre lies in the domain but whose surface crosses a periodic
// boundary (here x = 0 of the pipeflow domain, periodic in x) is skipped and
// counted at load by default; with setPlaceCellsAcrossPeriodicBoundaries() it
// is placed as images on both sides and stays one complete cell while the
// simulation runs. Run on 1 and 2 ranks (unit, unit_mpi).

using namespace pipeflow_test;

namespace {

void copyFile(const std::string & from, const std::string & to) {
  std::ifstream in(from);
  std::ofstream out(to);
  out << in.rdbuf();
}

struct Options {
  bool placeImages = false;
  std::string rbcPositions; // contents of the .pos file of the RBCs (micrometres, degrees)
  T forceFactor = 1;        // times the pipeflow driving force
  T initialVelocity = 0;    // of the fluid along the pipe [lattice units]
  unsigned int steps = 20;
};

struct Result {
  long cellsAfterLoad = 0;
  T relativeVolumeError = -1; // max over cells at load (images assembled into one shape)
  long cellsAfterSteps = 0;
  long distinctCellsAfterSteps = 0; // distinct base cell IDs of the counted cells
  long homeImagesWithShiftedId = 0; // counted cells whose ID is not the base ID (centre crossed)
  T volumeErrorAfterSteps = -1;     // max relative deviation from the reference volume
  long skippedCrossing = 0;         // main rank
};

// Volume errors of the cells whose centre lies on this rank (relative to the
// reference volume), and their IDs.
void measureCells(hemo::HemoCell & hemocell, const std::string & type, T & maxError, std::vector<long> & ids) {
  hemo::CellInformationFunctionals::calculateCellVolume(&hemocell);
  hemo::CellInformationFunctionals::calculateCellPosition(&hemocell);
  const T reference = (*hemocell.cellfields)[type]->meshmetric->getVolume();
  T error = -1;
  std::vector<long> local;
  for (const auto & entry : hemo::CellInformationFunctionals::info_per_cell) {
    if (!entry.second.centerLocal) { continue; }
    error = std::max(error, std::fabs(entry.second.volume - reference) / reference);
    local.push_back(long(entry.first));
  }
  hemo::CellInformationFunctionals::clear_list();
  MPI_Allreduce(&error, &maxError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  // Gather the IDs of all ranks.
  int count = local.size(), size = plb::global::mpi().getSize();
  std::vector<int> counts(size), offsets(size, 0);
  MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  for (int r = 1; r < size; r++) { offsets[r] = offsets[r - 1] + counts[r - 1]; }
  ids.resize(offsets[size - 1] + counts[size - 1]);
  MPI_Allgatherv(local.data(), count, MPI_LONG, ids.data(), counts.data(), offsets.data(), MPI_LONG, MPI_COMM_WORLD);
}

Result run(const Options & options) {
  const std::string rbc = hemo_test::outputRoot() + "periodic_RBC";
  const std::string plt = hemo_test::outputRoot() + "periodic_PLT";
  if (plb::global::mpi().isMainProcessor()) {
    copyFile(std::string(rbc_name) + ".xml", rbc + ".xml");
    copyFile(std::string(plt_name) + ".xml", plt + ".xml");
    std::ofstream(rbc + ".pos") << options.rbcPositions;
    std::ofstream(plt + ".pos") << "0\n";
  }
  plb::global::mpi().barrier();

  Result result;
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup, 0.5, rbc, plt);
  if (options.placeImages) { hemocell.setPlaceCellsAcrossPeriodicBoundaries(true); }
  if (options.initialVelocity != 0) { hemocell.latticeEquilibrium(1., {options.initialVelocity, 0., 0.}); }
  hemocell.loadParticles();
  result.skippedCrossing = (*hemocell.cellfields)[rbc]->cellsCrossingPeriodicBoundary;
  result.cellsAfterLoad = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
  std::vector<long> ids;
  measureCells(hemocell, rbc, result.relativeVolumeError, ids);

  for (int d = 0; d < 3; d++) { setup.drivingForce[d] *= options.forceFactor; }
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    setup.drivingForce);
  iterate(hemocell, setup, options.steps);
  result.cellsAfterSteps = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
  measureCells(hemocell, rbc, result.volumeErrorAfterSteps, ids);
  std::set<long> bases;
  for (const long id : ids) {
    bases.insert(long(hemocell.cellfields->base_cell_id(id)));
    if (hemocell.cellfields->base_cell_id(id) != id) { result.homeImagesWithShiftedId++; }
  }
  result.distinctCellsAfterSteps = bases.size();
  return result;
}

} // namespace

// One RBC centred 0.5 um from the x = 0 face, on the pipe axis (13 um).
const char * oneCellAcrossTheBoundary = "1\n0.5 13 13 0 0 0\n";

TEST(CellLoading, CellAcrossPeriodicBoundaryIsSkippedAndCountedByDefault) {
  Options options;
  options.rbcPositions = oneCellAcrossTheBoundary;
  const Result result = run(options);
  EXPECT_EQ(result.cellsAfterLoad, 0);
  if (plb::global::mpi().isMainProcessor()) { EXPECT_EQ(result.skippedCrossing, 1); }
}

TEST(CellLoading, CellAcrossPeriodicBoundaryIsPlacedOnRequest) {
  Options options;
  options.placeImages = true;
  options.rbcPositions = oneCellAcrossTheBoundary;
  const Result result = run(options);
  EXPECT_EQ(result.cellsAfterLoad, 1);
  EXPECT_GE(result.relativeVolumeError, 0); // measured
  EXPECT_LT(result.relativeVolumeError, 1e-9);
  EXPECT_EQ(result.cellsAfterSteps, 1);
  EXPECT_EQ(result.skippedCrossing, 0);
}

// Cells that cross the periodic boundary and (on 2 ranks) the block boundary
// in the middle of the pipe while the simulation runs stay complete: one RBC
// 0.7 um before the x = L face, whose centre crosses it (its home image then
// carries a shifted ID), and one in the middle. The fluid starts at 0.02
// lattice units per step along the pipe and the driving force is raised
// tenfold, so that the cells move a few lattice nodes in 100 steps.
TEST(CellLoading, CellsCrossingPeriodicAndBlockBoundariesAreConserved) {
  Options options;
  options.placeImages = true;
  options.rbcPositions = "2\n50.8 13 13 0 0 0\n25.75 13 13 0 0 0\n";
  options.forceFactor = 10;
  options.initialVelocity = 0.02;
  options.steps = 100;
  const Result result = run(options);
  ASSERT_EQ(result.cellsAfterLoad, 2);
  EXPECT_EQ(result.cellsAfterSteps, 2);
  EXPECT_EQ(result.distinctCellsAfterSteps, 2);
  EXPECT_EQ(result.homeImagesWithShiftedId, 1); // the first cell crossed x = L
  EXPECT_LT(result.volumeErrorAfterSteps, 0.01);
}
