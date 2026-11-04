#include "gtest/gtest.h"

#include "particleGrid.h"

#include <cmath>
#include <vector>

namespace {

using hemo::HemoCellParticle;

HemoCellParticle makeParticle(T x, T y, T z, plint cellId) {
  return HemoCellParticle(hemo::Array<T, 3>({x, y, z}), cellId, 0, 0);
}

T norm3(const hemo::Array<T, 3> & v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

std::vector<HemoCellParticle> & resetRepulsion(std::vector<HemoCellParticle> & particles) {
  for (HemoCellParticle & particle : particles) { particle.sv.force_repulsion = {0., 0., 0.}; }
  return particles;
}

} // namespace

// The previous grid stored at most 10 particles per lattice node without a
// bounds check (heap overflow). The grid must hold any number per node.
TEST(ParticleGrid, HoldsMoreThanTenParticlesPerNode) {
  std::vector<HemoCellParticle> particles;
  for (int i = 0; i < 25; i++) {
    particles.push_back(makeParticle(2.0 + 0.01 * i, 3.0, 4.0, i)); // all nearest to node (2,3,4)
  }
  particles.push_back(makeParticle(0.0, 0.0, 0.0, 100));  // node (0,0,0)
  particles.push_back(makeParticle(-3.0, 0.0, 0.0, 101)); // outside the block: not binned
  particles.push_back(makeParticle(5.6, 0.0, 0.0, 102));  // rounds to x == nx: not binned

  hemo::ParticleGrid grid;
  grid.build(particles, plb::Dot3D(0, 0, 0), 6, 6, 6);

  const unsigned int node = grid.index(2, 3, 4);
  ASSERT_EQ(grid.count(node), 25u);
  unsigned int expected = 0;
  for (const unsigned int * i = grid.begin(node); i != grid.end(node); ++i) {
    EXPECT_EQ(*i, expected++); // ascending particle order within a node
  }
  EXPECT_EQ(grid.count(grid.index(0, 0, 0)), 1u);
  EXPECT_EQ(*grid.begin(grid.index(0, 0, 0)), 25u);

  unsigned int total = 0;
  for (unsigned int n = 0; n < 6 * 6 * 6; n++) { total += grid.count(n); }
  EXPECT_EQ(total, 26u);
}

TEST(ParticleGrid, BinsRelativeToBlockLocation) {
  std::vector<HemoCellParticle> particles;
  particles.push_back(makeParticle(12.2, 7.4, 3.6, 0)); // local (2.2, 2.4, 0.6) -> node (2,2,1)

  hemo::ParticleGrid grid;
  grid.build(particles, plb::Dot3D(10, 5, 3), 4, 4, 4);
  EXPECT_EQ(grid.count(grid.index(2, 2, 1)), 1u);

  // Rebuilding with an empty particle set clears the grid.
  particles.clear();
  grid.build(particles, plb::Dot3D(10, 5, 3), 4, 4, 4);
  EXPECT_EQ(grid.count(grid.index(2, 2, 1)), 0u);
}

// Pairs inside one lattice node used to be visited twice (double force),
// pairs in neighbouring nodes once. Both must now get the same force.
TEST(ParticleGrid, RepulsionIsIndependentOfNodeAssignment) {
  const T r_const = 1.0;
  const T r_cutoff = 0.9;
  const T separation = 0.4;

  std::vector<HemoCellParticle> particles;
  particles.push_back(makeParticle(2.9, 3.0, 3.0, 0));              // node (3,3,3)
  particles.push_back(makeParticle(2.9 + separation, 3.0, 3.0, 1)); // node (3,3,3)
  particles.push_back(makeParticle(7.3, 3.0, 3.0, 2));              // node (7,3,3)
  particles.push_back(makeParticle(7.3 + separation, 3.0, 3.0, 3)); // node (8,3,3)

  hemo::ParticleGrid grid;
  grid.build(particles, plb::Dot3D(0, 0, 0), 12, 8, 8);
  ASSERT_EQ(grid.count(grid.index(3, 3, 3)), 2u);
  ASSERT_EQ(grid.count(grid.index(7, 3, 3)), 1u);
  ASSERT_EQ(grid.count(grid.index(8, 3, 3)), 1u);

  hemo::computePairRepulsion(resetRepulsion(particles), grid, r_const, r_cutoff);

  const T expected = r_const * r_cutoff / separation;
  EXPECT_NEAR(norm3(particles[0].sv.force_repulsion), expected, 1e-9);
  EXPECT_NEAR(norm3(particles[2].sv.force_repulsion), expected, 1e-9);

  // Equal and opposite, directed away from the partner.
  for (int d = 0; d < 3; d++) {
    EXPECT_NEAR(particles[0].sv.force_repulsion[d] + particles[1].sv.force_repulsion[d], 0.0, 1e-12);
    EXPECT_NEAR(particles[2].sv.force_repulsion[d] + particles[3].sv.force_repulsion[d], 0.0, 1e-12);
  }
  EXPECT_LT(particles[0].sv.force_repulsion[0], 0.0);
  EXPECT_LT(particles[2].sv.force_repulsion[0], 0.0);
}

TEST(ParticleGrid, NoRepulsionWithinACellOrBeyondCutoff) {
  std::vector<HemoCellParticle> particles;
  particles.push_back(makeParticle(3.0, 3.0, 3.0, 5));
  particles.push_back(makeParticle(3.3, 3.0, 3.0, 5)); // same cell: no force
  particles.push_back(makeParticle(4.4, 3.0, 3.0, 6)); // other cell, 1.1 away: beyond cutoff

  hemo::ParticleGrid grid;
  grid.build(particles, plb::Dot3D(0, 0, 0), 8, 8, 8);
  hemo::computePairRepulsion(resetRepulsion(particles), grid, 1.0, 0.9);

  for (const HemoCellParticle & particle : particles) { EXPECT_EQ(norm3(particle.sv.force_repulsion), 0.0); }
}
