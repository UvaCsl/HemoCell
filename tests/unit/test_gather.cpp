#include "gtest/gtest.h"

#include "hemoCellFunctional.h"
#include "parallelism/mpiManager.h"

#include <map>

// HemoCellGatheringFunctional combines the maps of all ranks; cell-keyed maps
// use 64-bit keys (cell IDs grow beyond the int range).

namespace {
const plb::plint largeKey = plb::plint(1) << 40;
}

TEST(Gather, AllRanksReceiveEntriesWith64BitKeys) {
  const int rank = plb::global::mpi().getRank();
  const int size = plb::global::mpi().getSize();
  std::map<plb::plint, double> values;
  values[largeKey + rank] = 0.5 * rank;
  hemo::HemoCellGatheringFunctional<double>::gather(values);
  ASSERT_EQ(values.size(), size_t(size));
  for (int r = 0; r < size; r++) {
    ASSERT_EQ(values.count(largeKey + r), 1u);
    EXPECT_EQ(values[largeKey + r], 0.5 * r);
  }
}

TEST(Gather, RootReceivesEntriesWithIntKeys) {
  const int rank = plb::global::mpi().getRank();
  const int size = plb::global::mpi().getSize();
  std::map<int, long> values;
  values[rank] = 10L * rank;
  hemo::HemoCellGatheringFunctional<long>::gatherToRoot(values);
  if (rank == 0) {
    ASSERT_EQ(values.size(), size_t(size));
    for (int r = 0; r < size; r++) { EXPECT_EQ(values[r], 10L * r); }
  } else {
    EXPECT_EQ(values.size(), 1u);
  }
}
