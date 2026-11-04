#include "gtest/gtest.h"

#include "hemoCellParticle.h"
#include "logfile.h"
#include "parallelism/mpiManager.h"

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

TEST(Array, InitializerListIsBounded) {
  hemo::Array<double,3> a = {1.0, 2.0, 3.0, 4.0}; // extra value used to overflow
  EXPECT_EQ(a[0], 1.0);
  EXPECT_EQ(a[1], 2.0);
  EXPECT_EQ(a[2], 3.0);

  a = {7.0, 8.0}; // assigns the given elements only
  EXPECT_EQ(a[0], 7.0);
  EXPECT_EQ(a[1], 8.0);
  EXPECT_EQ(a[2], 3.0);
}

TEST(Array, AddingADot3DAddsOnce) {
  hemo::Array<double,3> a = {1.0, 2.0, 3.0};
  a += plb::Dot3D(10, 20, 30); // used to add the offset three times
  EXPECT_EQ(a[0], 11.0);
  EXPECT_EQ(a[1], 22.0);
  EXPECT_EQ(a[2], 33.0);
}

namespace {
hemo::HemoCellParticle makeParticle() {
  hemo::HemoCellParticle p(hemo::Array<T,3>({1.0, 2.0, 3.0}), 7, 5, 0);
  p.sv.force = {0.5, 0.0, 0.0};
  p.force_total = {9.0, 9.0, 9.0};
  p.kernelWeights = {0.25, 0.75};
  return p;
}
}

TEST(HemoCellParticle, CopiesAndMovesKeepForcePointersOnOwnForce) {
  hemo::HemoCellParticle a = makeParticle();

  hemo::HemoCellParticle copy(a);
  EXPECT_EQ(copy.force_area, &copy.sv.force);
  EXPECT_EQ(copy.force_total[0], 9.0);
  EXPECT_EQ(copy.kernelWeights.size(), 2u);

  hemo::HemoCellParticle moved(std::move(copy));
  EXPECT_EQ(moved.force_volume, &moved.sv.force);
  EXPECT_EQ(moved.kernelWeights.size(), 2u);

  hemo::HemoCellParticle assigned(hemo::Array<T,3>({0.0, 0.0, 0.0}), 1, 1, 0);
  assigned = a;
  EXPECT_EQ(assigned.force_link, &assigned.sv.force);
  EXPECT_EQ(assigned.force_total[0], 9.0); // was not copied by operator=
  EXPECT_EQ(assigned.sv.cellId, 7);

  hemo::HemoCellParticle moveAssigned(hemo::Array<T,3>({0.0, 0.0, 0.0}), 1, 1, 0);
  moveAssigned = std::move(assigned);
  EXPECT_EQ(moveAssigned.force_bending, &moveAssigned.sv.force);
  EXPECT_EQ(moveAssigned.sv.force[0], 0.5);
}

TEST(HemoCellParticle, CopiesShareSeparatedOutputForces) {
  hemo::HemoCellParticle a = makeParticle();
  hemo::Array<T,3> separate = {0.0, 0.0, 0.0};
  a.force_area = &separate; // as during output (separateForceVectors)
  a.force_volume = &separate;

  hemo::HemoCellParticle copy(a);
  EXPECT_EQ(copy.force_area, &separate);

  std::vector<hemo::HemoCellParticle> particles;
  particles.push_back(a);
  particles.push_back(makeParticle());
  particles.reserve(100); // reallocation moves the particles
  EXPECT_EQ(particles[0].force_area, &separate);
  EXPECT_EQ(particles[1].force_area, &particles[1].sv.force);
}

TEST(Logfile, MessagesBeforeOpenAreWrittenToTheFile) {
  hemo::Logfile log;
  log.mainProcess = true;
  log << "early message " << 42 << std::endl;

  const std::string path = testing::TempDir() + "hemocell_logfile_test_r"
                         + std::to_string(plb::global::mpi().getRank());
  ASSERT_TRUE(log.open(path));
  log << "late message" << std::endl;
  log.close();

  std::ifstream file(path);
  std::stringstream content;
  content << file.rdbuf();
  EXPECT_EQ(content.str(), "early message 42\nlate message\n");
}

TEST(Logfile, NonMainProcessDropsMessages) {
  hemo::Logfile log; // mainProcess defaults to false
  log << "dropped" << std::endl;
  EXPECT_EQ(log.pending.str(), "");
}
