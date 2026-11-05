#include "gtest/gtest.h"
#include "pipeflow_setup.h"
#include "wbcHighOrderModel.h"
#include "rbcMalariaModel.h"
#include "parallelism/mpiManager.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Properties of the membrane forces of every cell model (CellMechanics::
// ParticleMechanics()) on a single cell built from the model's reference mesh:
// no force at the rest shape, zero net force of the internal forces, forces
// that rotate and translate with the cell, and zero net torque of all terms
// but bending. The bending force of a vertex is balanced by equal shares on
// its neighbours, which gives zero net force but not zero torque (see
// developer/known_limitations.rst); its torque is bounded at the measured
// level. Runs on 1 and 2 ranks (unit, unit_mpi); the cell lives on each rank
// separately (no communication involved).

using namespace pipeflow_test;

namespace {

const std::string tempDir() {
  return testing::TempDir();
}

// A material file of another case, copied to the temporary directory (with
// `replace` -> `with` applied); returns the cell type name.
std::string copyMaterial(const std::string & from, const std::string & name, const std::string & replace = "",
                         const std::string & with = "") {
  std::ifstream in(from);
  std::stringstream content;
  content << in.rdbuf();
  std::string xml = content.str();
  if (!replace.empty()) {
    const size_t at = xml.find(replace);
    EXPECT_NE(at, std::string::npos) << replace;
    if (at != std::string::npos) { xml.replace(at, replace.size(), with); }
  }
  const std::string type =
      tempDir() + "hemocell_mechanics_" + name + "_r" + std::to_string(plb::global::mpi().getRank());
  std::ofstream out(type + ".xml");
  out << xml;
  return type;
}

// One cell of a cell type, with its vertices at given positions.
class Cell {
public:
  Cell(hemo::HemoCellField & field_) : field(field_) {
    const plb::TriangularSurfaceMesh<T> & mesh = *field.meshElement;
    particles.reserve(mesh.getNumVertices()); // no reallocation: force pointers stay valid
    for (plb::plint i = 0; i < mesh.getNumVertices(); i++) {
      const plb::Array<T, 3> & v = mesh.getVertex(i);
      particles.emplace_back(hemo::Array<T, 3>{v[0], v[1], v[2]}, 0, i, field.ctype);
    }
  }

  std::vector<hemo::Array<T, 3>> positions() const {
    std::vector<hemo::Array<T, 3>> result;
    for (const hemo::HemoCellParticle & p : particles) { result.push_back(p.sv.position); }
    return result;
  }

  // Forces [lattice units] of the cell model at the given vertex positions;
  // without the bending term if `withoutBending`.
  std::vector<hemo::Array<T, 3>> forces(const std::vector<hemo::Array<T, 3>> & positions, bool withoutBending = false) {
    std::map<plb::plint, std::vector<hemo::HemoCellParticle *>> ppc;
    for (pluint i = 0; i < particles.size(); i++) {
      particles[i].sv.position = positions[i];
      particles[i].sv.v = {0., 0., 0.};
      particles[i].sv.force = {0., 0., 0.}; // the force components alias it
      ppc[0].push_back(&particles[i]);
    }
    const std::map<plb::plint, bool> lpc = {{0, true}};
    std::vector<hemo::Array<T, 3>> discarded(particles.size(), hemo::Array<T, 3>{0., 0., 0.});
    if (withoutBending) {
      for (pluint i = 0; i < particles.size(); i++) { particles[i].force_bending = &discarded[i]; }
    }
    field.mechanics->ParticleMechanics(ppc, lpc, field.ctype);
    for (hemo::HemoCellParticle & p : particles) { p.force_bending = &p.sv.force; }
    std::vector<hemo::Array<T, 3>> result;
    for (const hemo::HemoCellParticle & p : particles) { result.push_back(p.sv.force); }
    return result;
  }

private:
  hemo::HemoCellField & field;
  std::vector<hemo::HemoCellParticle> particles;
};

T norm3(const hemo::Array<T, 3> & a) {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

T maxNorm(const std::vector<hemo::Array<T, 3>> & v) {
  T m = 0;
  for (const auto & a : v) { m = std::max(m, norm3(a)); }
  return m;
}

hemo::Array<T, 3> centroid(const std::vector<hemo::Array<T, 3>> & x) {
  hemo::Array<T, 3> c = {0., 0., 0.};
  for (const auto & a : x) { c += a; }
  return c / T(x.size());
}

// A deformed shape: stretched along x, compressed along y, plus a fixed
// pseudo-random displacement of each vertex (2 % of the mean edge length).
std::vector<hemo::Array<T, 3>> deform(const std::vector<hemo::Array<T, 3>> & rest, T edge) {
  const hemo::Array<T, 3> c = centroid(rest);
  std::vector<hemo::Array<T, 3>> x;
  unsigned long state = 12345;
  auto random = [&state]() {
    state = state * 6364136223846793005UL + 1442695040888963407UL;
    return T((state >> 11) & 0xFFFFF) / T(0xFFFFF) - 0.5;
  };
  for (const auto & r : rest) {
    hemo::Array<T, 3> d = r - c;
    d[0] *= 1.05;
    d[1] *= 0.97;
    x.push_back(c + d + hemo::Array<T, 3>{random(), random(), random()} * (0.04 * edge));
  }
  return x;
}

// Rotation about an axis through the centroid, then a translation.
hemo::Array<T, 3> rotate(const hemo::Array<T, 3> & a) {
  const T ax = 0.3, ay = -1.1, az = 0.7; // Euler angles [rad], X then Y then Z
  hemo::Array<T, 3> b = {a[0], std::cos(ax) * a[1] - std::sin(ax) * a[2], std::sin(ax) * a[1] + std::cos(ax) * a[2]};
  b = {std::cos(ay) * b[0] + std::sin(ay) * b[2], b[1], -std::sin(ay) * b[0] + std::cos(ay) * b[2]};
  return {std::cos(az) * b[0] - std::sin(az) * b[1], std::sin(az) * b[0] + std::cos(az) * b[1], b[2]};
}

struct Measurement {
  T restForce;               // max |F| at rest / max |F| deformed
  T netForce;                // |sum F| / sum |F|, deformed
  T netTorque;               // |sum (x - c) x F| / sum |x - c| |F|, deformed
  T netTorqueWithoutBending; // the same without the bending term
  T rotationError;           // max |F(R x + t) - R F(x)| / max |F(x)|
};

// |sum (x - c) x F| / sum |x - c| |F|
T relativeTorque(const std::vector<hemo::Array<T, 3>> & x, const std::vector<hemo::Array<T, 3>> & f) {
  const hemo::Array<T, 3> c = centroid(x);
  hemo::Array<T, 3> sum = {0., 0., 0.};
  T sumAbs = 0;
  for (pluint i = 0; i < x.size(); i++) {
    const hemo::Array<T, 3> r = x[i] - c;
    sum += hemo::Array<T, 3>{r[1] * f[i][2] - r[2] * f[i][1], r[2] * f[i][0] - r[0] * f[i][2],
                             r[0] * f[i][1] - r[1] * f[i][0]};
    sumAbs += norm3(r) * norm3(f[i]);
  }
  return norm3(sum) / sumAbs;
}

Measurement measure(hemo::HemoCellField & field) {
  Cell cell(field);
  const std::vector<hemo::Array<T, 3>> rest = cell.positions();
  const T edge = field.mechanics->cellConstants.edge_mean_eq;
  const std::vector<hemo::Array<T, 3>> restForces = cell.forces(rest);
  const std::vector<hemo::Array<T, 3>> x = deform(rest, edge);
  const std::vector<hemo::Array<T, 3>> f = cell.forces(x);
  const T scale = maxNorm(f);

  Measurement m;
  m.restForce = maxNorm(restForces) / scale;

  hemo::Array<T, 3> sumF = {0., 0., 0.};
  T sumAbsF = 0;
  for (pluint i = 0; i < x.size(); i++) {
    sumF += f[i];
    sumAbsF += norm3(f[i]);
  }
  m.netForce = norm3(sumF) / sumAbsF;
  m.netTorque = relativeTorque(x, f);
  m.netTorqueWithoutBending = relativeTorque(x, cell.forces(x, true));

  std::vector<hemo::Array<T, 3>> y;
  const hemo::Array<T, 3> shift = {3.25, -1.5, 7.0};
  for (const auto & a : x) { y.push_back(rotate(a) + shift); }
  const std::vector<hemo::Array<T, 3>> g = cell.forces(y);
  T error = 0;
  for (pluint i = 0; i < x.size(); i++) { error = std::max(error, norm3(g[i] - rotate(f[i]))); }
  m.rotationError = error / scale;
  return m;
}

struct Models {
  hemo::HemoCell * hemocell = nullptr;
  std::map<std::string, std::string> names; // model -> cell type name
};

// The pipeflow setup (RBC: RbcHighOrderModel, PLT: PltSimpleModel) plus a WBC
// (WbcHighOrderModel, cases/pipeflow_wbc) and a malaria-infected RBC
// (RbcMalariaModel, cases/pipeflowMalaria, its STL mesh).
void addModels(hemo::HemoCell & hemocell, Models & models) {
  models.hemocell = &hemocell;
  models.names["RbcHighOrderModel"] = rbc_name;
  models.names["PltSimpleModel"] = plt_name;
  const std::string wbc = copyMaterial("cases/pipeflow_wbc/WBC.xml", "WBC");
  hemocell.addCellType<hemo::WbcHighOrderModel>(wbc, hemo::WBC_SPHERE);
  models.names["WbcHighOrderModel"] = wbc;
  char cwd[4096];
  const std::string root = getcwd(cwd, sizeof(cwd)) ? std::string(cwd) : std::string(".");
  const std::string malaria =
      copyMaterial("cases/pipeflowMalaria/RBC_MALARIA.xml", "MALARIA", "<StlFile>vRBC_uniform.stl</StlFile>",
                   "<StlFile>" + root + "/cases/pipeflowMalaria/vRBC_uniform.stl</StlFile>");
  hemocell.addCellType<hemo::RbcMalariaModel>(malaria, hemo::MESH_FROM_STL);
  models.names["RbcMalariaModel"] = malaria;
}

} // namespace

TEST(Mechanics, MembraneForcesOfEveryModel) {
  char * args[] = {(char *)"test", NULL};
  hemo::HemoCell hemocell(pipeflowConfig(), 0, args, hemo::HemoCell::MPIHandle::External);
  PipeflowSetup setup;
  setupPipeflow(hemocell, setup);
  Models models;
  addModels(hemocell, models);
  ASSERT_EQ(models.names.size(), 4u);
  for (const auto & model : models.names) {
    SCOPED_TRACE(model.first);
    const Measurement m = measure(*(*hemocell.cellfields)[model.second]);
    // Round-off level (measured: 1e-15 to 2e-13).
    EXPECT_LT(m.restForce, 1e-10);
    // Exact up to round-off (measured: about 1e-17).
    EXPECT_LT(m.netForce, 1e-12);
    // Exact up to round-off (measured: about 1e-13).
    EXPECT_LT(m.rotationError, 1e-10);
    EXPECT_LT(m.netTorqueWithoutBending, 1e-12);
    // Bending is not torque-free (measured: 1.5e-5 to 3e-5 for this
    // deformation); a regression guard, not a property.
    EXPECT_LT(m.netTorque, 1e-4);
  }
}
