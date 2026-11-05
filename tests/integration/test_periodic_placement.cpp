#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "cellInfo.h"
#include "parallelism/mpiManager.h"

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

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

struct Result {
  long cellsAfterLoad = 0;
  T relativeVolumeError = -1; // of the placed cell at load (images assembled into one shape)
  long cellsAfterSteps = 0;
  long skippedCrossing = 0; // main rank
};

Result run(bool placeImages) {
  // One RBC centred 0.5 um from the x = 0 face, on the pipe axis (13 um), and no platelets.
  const std::string rbc = testing::TempDir() + "hemocell_periodic_RBC";
  const std::string plt = testing::TempDir() + "hemocell_periodic_PLT";
  if (plb::global::mpi().isMainProcessor()) {
    copyFile(std::string(rbc_name) + ".xml", rbc + ".xml");
    copyFile(std::string(plt_name) + ".xml", plt + ".xml");
    std::ofstream(rbc + ".pos") << "1\n0.5 13 13 0 0 0\n";
    std::ofstream(plt + ".pos") << "0\n";
  }
  plb::global::mpi().barrier();

  Result result;
  {
    char * args[] = {(char *)"test", NULL};
    hemo::HemoCell hemocell(pipeflow_config, 0, args, hemo::HemoCell::MPIHandle::External);
    PipeflowSetup setup;
    setupPipeflow(hemocell, setup, 0.5, rbc, plt);
    if (placeImages) { hemocell.setPlaceCellsAcrossPeriodicBoundaries(true); }
    hemocell.loadParticles();
    result.cellsAfterLoad = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
    // The volume of the complete cell (on the rank that holds its centre)
    // equals the reference volume only if its images line up.
    hemo::CellInformationFunctionals::calculateCellVolume(&hemocell);
    hemo::CellInformationFunctionals::calculateCellPosition(&hemocell);
    T error = -1;
    for (const auto & entry : hemo::CellInformationFunctionals::info_per_cell) {
      if (!entry.second.centerLocal) { continue; }
      const T reference = (*hemocell.cellfields)[rbc]->meshmetric->getVolume();
      error = std::fabs(entry.second.volume - reference) / reference;
    }
    hemo::CellInformationFunctionals::clear_list();
    MPI_Allreduce(&error, &result.relativeVolumeError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    result.skippedCrossing = (*hemocell.cellfields)[rbc]->cellsCrossingPeriodicBoundary;
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, setup.drivingForce);
    iterate(hemocell, setup, 20);
    result.cellsAfterSteps = hemo::CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
  }
  plb::global::mpi().barrier();
  if (plb::global::mpi().isMainProcessor()) {
    for (const std::string & name : {rbc, plt}) {
      std::remove((name + ".xml").c_str());
      std::remove((name + ".pos").c_str());
    }
  }
  return result;
}

} // namespace

TEST(CellLoading, CellAcrossPeriodicBoundaryIsSkippedAndCountedByDefault) {
  const Result result = run(false);
  EXPECT_EQ(result.cellsAfterLoad, 0);
  if (plb::global::mpi().isMainProcessor()) { EXPECT_EQ(result.skippedCrossing, 1); }
}

TEST(CellLoading, CellAcrossPeriodicBoundaryIsPlacedOnRequest) {
  const Result result = run(true);
  EXPECT_EQ(result.cellsAfterLoad, 1);
  EXPECT_GE(result.relativeVolumeError, 0); // measured
  EXPECT_LT(result.relativeVolumeError, 1e-9);
  EXPECT_EQ(result.cellsAfterSteps, 1);
  EXPECT_EQ(result.skippedCrossing, 0);
}
