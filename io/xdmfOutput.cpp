/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "xdmfOutput.h"
#include "genericFunctions.h"
#include "logfile.h"

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"

#include <dirent.h>
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace hemo {
namespace {

// Block records travel to rank 0 as text: one block per line, fields
// separated by tabs (identifiers and dataset names contain neither tabs nor
// newlines). Doubles are written with 17 digits, so they arrive unchanged.
std::string serialize(const std::vector<XdmfBlockFile> & blocks) {
  std::ostringstream out;
  out << std::setprecision(17);
  for (const XdmfBlockFile & b : blocks) {
    out << b.identifier << '\t' << b.file << '\t' << b.blockId << '\t' << (b.structured ? 1 : 0);
    for (int i = 0; i < 3; i++) { out << '\t' << b.nodes[i]; }
    for (int i = 0; i < 3; i++) { out << '\t' << b.origin[i]; }
    for (int i = 0; i < 3; i++) { out << '\t' << b.spacing[i]; }
    out << '\t' << b.vertices << '\t' << b.triangles << '\t' << b.innerLinks << '\t' << b.attributes.size();
    for (const XdmfAttribute & a : b.attributes) {
      out << '\t' << a.name << '\t' << a.components << '\t' << a.numberType << '\t' << a.precision;
    }
    out << '\n';
  }
  return out.str();
}

std::vector<XdmfBlockFile> parse(const std::string & text) {
  std::vector<XdmfBlockFile> blocks;
  std::istringstream lines(text);
  std::string line;
  while (std::getline(lines, line)) {
    if (line.empty()) { continue; }
    std::vector<std::string> f;
    std::istringstream fields(line);
    std::string field;
    while (std::getline(fields, field, '\t')) { f.push_back(field); }
    std::size_t k = 0;
    XdmfBlockFile b;
    b.identifier = f[k++];
    b.file = f[k++];
    b.blockId = std::stol(f[k++]);
    b.structured = f[k++] == "1";
    for (int i = 0; i < 3; i++) { b.nodes[i] = std::stol(f[k++]); }
    for (int i = 0; i < 3; i++) { b.origin[i] = std::stod(f[k++]); }
    for (int i = 0; i < 3; i++) { b.spacing[i] = std::stod(f[k++]); }
    b.vertices = std::stol(f[k++]);
    b.triangles = std::stol(f[k++]);
    b.innerLinks = std::stol(f[k++]);
    const std::size_t nAttributes = std::stoul(f[k++]);
    for (std::size_t a = 0; a < nAttributes; a++) {
      XdmfAttribute attribute;
      attribute.name = f[k++];
      attribute.components = std::stoi(f[k++]);
      attribute.numberType = f[k++];
      attribute.precision = std::stoi(f[k++]);
      b.attributes.push_back(attribute);
    }
    blocks.push_back(b);
  }
  return blocks;
}

// Concatenation of the strings of all ranks, on rank 0 (empty elsewhere).
std::string gatherToRoot(const std::string & local) {
  const int rank = plb::global::mpi().getRank();
  const int size = plb::global::mpi().getSize();
  int length = local.size();
  std::vector<int> lengths(rank == 0 ? size : 0);
  MPI_Gather(&length, 1, MPI_INT, lengths.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<int> offsets(lengths.size(), 0);
  int total = 0;
  for (std::size_t r = 0; r < lengths.size(); r++) {
    offsets[r] = total;
    total += lengths[r];
  }
  std::vector<char> all(total + 1);
  MPI_Gatherv(local.data(), length, MPI_CHAR, all.data(), lengths.data(), offsets.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
  return std::string(all.data(), total);
}

std::string escape(const std::string & text) {
  std::string out;
  for (char c : text) {
    switch (c) {
      case '&':
        out += "&amp;";
        break;
      case '<':
        out += "&lt;";
        break;
      case '>':
        out += "&gt;";
        break;
      case '"':
        out += "&quot;";
        break;
      default:
        out += c;
    }
  }
  return out;
}

const char * attributeType(int components) {
  switch (components) {
    case 1:
      return "Scalar";
    case 3:
      return "Vector";
    case 6:
      return "Tensor6";
    case 9:
      return "Tensor";
    default:
      return "Matrix";
  }
}

void writeHdfItem(std::ostream & out, const std::string & dimensions, const std::string & numberType,
                  const std::string & file, const std::string & dataset, int precision = 4) {
  out << "<DataItem Dimensions=\"" << dimensions << "\" NumberType=\"" << numberType << "\" Precision=\"" << precision
      << "\" Format=\"HDF\">" << escape(file) << ":/" << escape(dataset) << "</DataItem>";
}

void writeStructuredBlock(std::ostream & out, const XdmfBlockFile & b) {
  out << "      <Grid Name=\"block " << b.blockId << "\" GridType=\"Uniform\">\n";
  // A node value fills the cell around the node: one point more than nodes.
  out << "        <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << b.nodes[0] + 1 << ' ' << b.nodes[1] + 1
      << ' ' << b.nodes[2] + 1 << "\"/>\n";
  out << "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
  out << "          <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">" << b.origin[0]
      << ' ' << b.origin[1] << ' ' << b.origin[2] << "</DataItem>\n";
  out << "          <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">" << b.spacing[0]
      << ' ' << b.spacing[1] << ' ' << b.spacing[2] << "</DataItem>\n";
  out << "        </Geometry>\n";
  for (const XdmfAttribute & a : b.attributes) {
    std::ostringstream dims;
    dims << b.nodes[0] << ' ' << b.nodes[1] << ' ' << b.nodes[2] << ' ' << a.components;
    out << "        <Attribute Name=\"" << escape(a.name) << "\" AttributeType=\"" << attributeType(a.components)
        << "\" Center=\"Cell\">\n          ";
    writeHdfItem(out, dims.str(), a.numberType, b.file, a.name, a.precision);
    out << "\n        </Attribute>\n";
  }
  out << "      </Grid>\n";
}

void writeCellBlock(std::ostream & out, const XdmfBlockFile & b) {
  std::ostringstream positionDims;
  positionDims << b.vertices << " 3";
  out << "      <Grid Name=\"block " << b.blockId << "\" GridType=\"Uniform\">\n";
  if (b.triangles > 0) {
    std::ostringstream dims;
    dims << b.triangles << " 3";
    out << "        <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << b.triangles << "\">\n          ";
    writeHdfItem(out, dims.str(), "Int", b.file, "Triangles");
    out << "\n        </Topology>\n";
  } else {
    out << "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << b.vertices
        << "\" NodesPerElement=\"1\"/>\n";
  }
  out << "        <Geometry GeometryType=\"XYZ\">\n          ";
  writeHdfItem(out, positionDims.str(), "Float", b.file, "Position");
  out << "\n        </Geometry>\n";
  for (const XdmfAttribute & a : b.attributes) {
    std::ostringstream dims;
    dims << b.vertices << ' ' << a.components;
    out << "        <Attribute Name=\"" << escape(a.name) << "\" AttributeType=\"" << attributeType(a.components)
        << "\" Center=\"Node\">\n          ";
    writeHdfItem(out, dims.str(), a.numberType, b.file, a.name, a.precision);
    out << "\n        </Attribute>\n";
  }
  out << "      </Grid>\n";
  if (b.innerLinks > 0) {
    std::ostringstream dims;
    dims << b.innerLinks << " 2";
    out << "      <Grid Name=\"block " << b.blockId << " inner links\" GridType=\"Uniform\">\n";
    out << "        <Topology TopologyType=\"Polyline\" NumberOfElements=\"" << b.innerLinks
        << "\" NodesPerElement=\"2\">\n          ";
    writeHdfItem(out, dims.str(), "Int", b.file, "InnerLinks");
    out << "\n        </Topology>\n";
    out << "        <Geometry GeometryType=\"XYZ\">\n          ";
    writeHdfItem(out, positionDims.str(), "Float", b.file, "Position");
    out << "\n        </Geometry>\n";
    out << "      </Grid>\n";
  }
}

// Write to a temporary file and rename it, so a reader (e.g. ParaView during
// the run) never sees a partially written file.
void writeFile(const std::string & path, const std::string & content) {
  const std::string temporary = path + ".tmp";
  bool ok;
  {
    std::ofstream out(temporary.c_str(), std::ios::trunc);
    out << content;
    ok = bool(out);
  }
  if (!ok || std::rename(temporary.c_str(), path.c_str()) != 0) {
    hlog << "(HemoCell) (Output) WARNING: could not write " << path << std::endl;
  }
}

// Step files of `identifier` in `directory` (<identifier>.<digits>.xmf), sorted.
std::vector<std::string> findStepFiles(const std::string & directory, const std::string & identifier) {
  std::vector<std::string> found;
  DIR * dir = opendir(directory.c_str());
  if (!dir) { return found; }
  const std::string prefix = identifier + ".";
  const std::string suffix = ".xmf";
  while (struct dirent * entry = readdir(dir)) {
    const std::string name = entry->d_name;
    if (name.size() <= prefix.size() + suffix.size() || name.compare(0, prefix.size(), prefix) != 0 ||
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) {
      continue;
    }
    const std::string number = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
    if (number.find_first_not_of("0123456789") == std::string::npos) { found.push_back(name); }
  }
  closedir(dir);
  // Zero-padded iterations: lexicographic order is the order in time.
  std::sort(found.begin(), found.end());
  return found;
}

} // namespace

void writeXdmf(const std::vector<XdmfBlockFile> & blocks, plb::plint iter, double time) {
  const std::string all = gatherToRoot(serialize(blocks));
  if (plb::global::mpi().getRank() != 0) { return; }

  std::map<std::string, std::vector<XdmfBlockFile>> byIdentifier;
  for (const XdmfBlockFile & b : parse(all)) { byIdentifier[b.identifier].push_back(b); }

  std::string directory = plb::global::directories().getOutputDir();
  if (directory.empty() || directory[directory.size() - 1] != '/') { directory += "/"; }

  for (auto & entry : byIdentifier) {
    const std::string & identifier = entry.first;
    std::vector<XdmfBlockFile> & files = entry.second;
    std::sort(files.begin(), files.end(),
              [](const XdmfBlockFile & a, const XdmfBlockFile & b) { return a.blockId < b.blockId; });

    std::ostringstream step;
    // 15 digits: origins print as e.g. 1.325e-05 rather than 1.3249999999999999e-05.
    step << std::setprecision(15);
    step << "<?xml version=\"1.0\" ?>\n<Xdmf Version=\"3.0\">\n  <Domain>\n";
    step << "    <Grid Name=\"" << escape(identifier) << "\" GridType=\"Collection\" CollectionType=\"Spatial\">\n";
    step << "      <Time Value=\"" << std::setprecision(12) << time << std::setprecision(15) << "\"/>\n";
    for (const XdmfBlockFile & b : files) {
      if (b.structured) {
        writeStructuredBlock(step, b);
      } else if (b.vertices > 0) {
        writeCellBlock(step, b);
      }
    }
    step << "    </Grid>\n  </Domain>\n</Xdmf>\n";
    writeFile(directory + identifier + "." + zeroPadNumber(iter) + ".xmf", step.str());

    std::ostringstream series;
    series << "<?xml version=\"1.0\" ?>\n<Xdmf Version=\"3.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n"
           << "  <Domain>\n    <Grid Name=\"" << escape(identifier)
           << "\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
    // element(/1/1/1): the first Grid of the Domain of the included file.
    for (const std::string & name : findStepFiles(directory, identifier)) {
      series << "      <xi:include href=\"" << escape(name) << "\" xpointer=\"element(/1/1/1)\"/>\n";
    }
    series << "    </Grid>\n  </Domain>\n</Xdmf>\n";
    writeFile(directory + identifier + ".xmf", series.str());
  }
}

} // namespace hemo
