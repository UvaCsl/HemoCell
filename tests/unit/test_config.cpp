#include "gtest/gtest.h"

#include "config.h"
#include "parallelism/mpiManager.h"

#include <fstream>
#include <stdexcept>
#include <string>

namespace {

// Write an XML file into the GoogleTest temporary directory. The MPI rank is
// part of the name, so concurrent ranks do not write the same file.
std::string writeTempXml(const std::string & name, const std::string & content) {
  const std::string path =
      testing::TempDir() + "hemocell_config_test_r" + std::to_string(plb::global::mpi().getRank()) + "_" + name;
  std::ofstream file(path);
  file << content;
  return path;
}

} // namespace

TEST(Config, RootWithoutPrologIsNotACheckpoint) {
  hemo::Config cfg(writeTempXml("noprolog.xml", "<hemocell><domain><dx> 5e-7 </dx></domain></hemocell>"));
  EXPECT_FALSE(cfg.checkpointed);
  EXPECT_DOUBLE_EQ(cfg["domain"]["dx"].read<double>(), 5e-7);
}

TEST(Config, LeadingCommentIsNotACheckpoint) {
  hemo::Config cfg(writeTempXml(
      "comment.xml", "<?xml version=\"1.0\"?>\n<!-- comment -->\n<hemocell><sim><tmax>10</tmax></sim></hemocell>\n"));
  EXPECT_FALSE(cfg.checkpointed);
  EXPECT_EQ(cfg["sim"]["tmax"].read<int>(), 10);
}

TEST(Config, CheckpointRootIsDetected) {
  hemo::Config cfg(writeTempXml("checkpoint.xml",
                                "<?xml version=\"1.0\"?>\n<Checkpoint><hemocell><sim><tmax>10</tmax></sim></hemocell>"
                                "<General><Iteration>5</Iteration></General></Checkpoint>\n"));
  EXPECT_TRUE(cfg.checkpointed);
  EXPECT_EQ(cfg["sim"]["tmax"].read<int>(), 10); // read below <Checkpoint><hemocell>
}

TEST(Config, ReadRejectsInvalidValues) {
  hemo::Config cfg(
      writeTempXml("values.xml", "<hemocell><a>1e6</a><b></b><c>abc</c><d> 12 </d><e>3.5x</e></hemocell>"));
  EXPECT_THROW(cfg["a"].read<int>(), std::invalid_argument); // used to silently become 1
  EXPECT_DOUBLE_EQ(cfg["a"].read<double>(), 1e6);
  EXPECT_THROW(cfg["b"].read<int>(), std::invalid_argument);    // empty element used to crash
  EXPECT_THROW(cfg["c"].read<double>(), std::invalid_argument); // used to silently become 0
  EXPECT_EQ(cfg["d"].read<int>(), 12);
  EXPECT_THROW(cfg["e"].read<double>(), std::invalid_argument);
  EXPECT_THROW(cfg["missing"].read<int>(), std::invalid_argument);
}

TEST(Config, ReadBoolAndString) {
  hemo::Config cfg(writeTempXml("boolstring.xml", "<hemocell><t>true</t><f> 0 </f><one>1</one><x>yes</x>"
                                                  "<s>  output dir/with space  </s></hemocell>"));
  EXPECT_TRUE(cfg["t"].read<bool>());
  EXPECT_FALSE(cfg["f"].read<bool>());
  EXPECT_TRUE(cfg["one"].read<bool>());
  EXPECT_THROW(cfg["x"].read<bool>(), std::invalid_argument);
  EXPECT_EQ(cfg["s"].read<std::string>(), "output dir/with space");
}

namespace {
// The message of the std::invalid_argument thrown by `access`.
template <typename F>
std::string errorOf(F access) {
  try {
    access();
  } catch (const std::invalid_argument & e) { return e.what(); }
  return "";
}
} // namespace

// Errors name the file, the element path and its line, so a missing or
// invalid setting can be found without knowing which code read it.
TEST(Config, ErrorsNameFileElementAndLine) {
  const std::string path = writeTempXml("errors.xml", "<hemocell>\n<sim>\n<tmax></tmax>\n<tmeas>x</tmeas>\n</sim>\n"
                                                      "</hemocell>\n");
  hemo::Config cfg(path);

  const std::string missing = errorOf([&] { cfg["sim"]["tcheckpoint"]; });
  EXPECT_NE(missing.find(path), std::string::npos) << missing;
  EXPECT_NE(missing.find("line 2: hemocell/sim has no element <tcheckpoint>"), std::string::npos) << missing;

  const std::string empty = errorOf([&] { cfg["sim"]["tmax"].read<int>(); });
  EXPECT_NE(empty.find("line 3: hemocell/sim/tmax is empty"), std::string::npos) << empty;

  const std::string invalid = errorOf([&] { cfg["sim"]["tmeas"].read<int>(); });
  EXPECT_NE(invalid.find("line 4: hemocell/sim/tmeas: cannot convert \"x\""), std::string::npos) << invalid;

  const std::string root = errorOf([&] { cfg["domain"]; });
  EXPECT_NE(root.find("line 1: hemocell has no element <domain>"), std::string::npos) << root;
}
