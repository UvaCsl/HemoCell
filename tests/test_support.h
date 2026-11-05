// Helpers shared by the tests.
#ifndef HEMOCELL_TESTS_TEST_SUPPORT_H
#define HEMOCELL_TESTS_TEST_SUPPORT_H

#include "gtest/gtest.h"
#include "parallelism/mpiManager.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>

namespace hemo_test {

/// Directory for everything the tests' simulations write (output, logs,
/// checkpoints) and for the configuration copies of testConfig(); set by
/// tests/main.cpp (one per run, in the GoogleTest temporary directory) and
/// removed at the end of the run unless HEMOCELL_KEEP_TEST_OUTPUT is set.
inline std::string & outputRoot() {
  static std::string root;
  return root;
}

/**
 * @brief A copy of the configuration file @p path whose output directory lies
 * in outputRoot(), so that the tests do not write into the repository.
 *
 * The copy is private to this rank (its name contains the rank and a running
 * number); the output directory is the same on all ranks, as HemoCell needs.
 * Input paths in the file stay relative to the working directory.
 * @param[in] path configuration file, relative to the working directory
 * @param[in] edit applied to the text of the copy (e.g. to change a setting)
 * @return path of the copy
 */
inline std::string testConfig(const std::string & path, const std::function<void(std::string &)> & edit = nullptr) {
  static int copies = 0;
  std::ifstream in(path);
  EXPECT_TRUE(in.good()) << path;
  std::stringstream content;
  content << in.rdbuf();
  std::string xml = content.str();
  const std::string directory = "<outputDirectory>" + outputRoot() + "run</outputDirectory>";
  const size_t begin = xml.find("<outputDirectory>");
  const size_t end = xml.find("</outputDirectory>", begin);
  if (begin != std::string::npos && end != std::string::npos) {
    xml.replace(begin, end + 18 - begin, directory);
  } else if (xml.find("<parameters>") != std::string::npos) {
    xml.replace(xml.find("<parameters>"), 12, "<parameters>" + directory);
  } else {
    xml.replace(xml.find("<hemocell>"), 10, "<hemocell><parameters>" + directory + "</parameters>");
  }
  if (edit) { edit(xml); }
  const std::string copy = outputRoot() + "config_" + std::to_string(copies++) + "_r" +
                           std::to_string(plb::global::mpi().getRank()) + "_" +
                           std::filesystem::path(path).filename().string();
  std::ofstream out(copy);
  out << xml;
  return copy;
}

/// Replace the first `from` in `text` by `to` (a failed expectation if absent).
inline void replaceOnce(std::string & text, const std::string & from, const std::string & to) {
  const size_t at = text.find(from);
  EXPECT_NE(at, std::string::npos) << from;
  if (at != std::string::npos) { text.replace(at, from.size(), to); }
}

} // namespace hemo_test

#endif
