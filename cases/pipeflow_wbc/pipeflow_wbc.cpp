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
// Pipe flow with RBCs, PLTs and one WBC, driven by the pressure gradient of a
// given Newtonian wall shear rate, with a running check of the numerical
// stability of the cells (see README.md).
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "wbcHighOrderModel.h"
#include "fluidInfo.h"

#include "palabos3D.h"
#include "palabos3D.hh"

#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

typedef double T;

using namespace hemo;
using namespace plb;
using namespace std;

namespace {

// Shape and force metrics of the cells of one type.
struct TypeMetrics {
  double cells = 0, vertices = 0, cappedVertices = 0, nonFinite = 0;
  double maxForce = 0;          // |F| per vertex [lattice force]
  double maxEdgeStrain = 0;     // L / L0 - 1 of membrane links
  double minEdgeStrain = 0;     //
  double maxAreaChange = 0;     // |A / A0 - 1| per triangle
  double maxVolumeChange = 0;   // |V / V0 - 1| per cell
  double maxElongation = 0;     // largest vertex distance from the centroid / that of the reference mesh
  double minCompression = 1e30; // smallest vertex distance from the centroid / that of the reference mesh
  // Summed over the cells (used for the single WBC): centroid and mean x velocity [lattice].
  double centroid[3] = {0, 0, 0};
  double velocityX = 0;
};

// Largest and smallest vertex distance from the centroid of the reference mesh.
void referenceExtent(HemoCellField & field, double & rmax, double & rmin) {
  plb::TriangularSurfaceMesh<T> & mesh = *field.meshElement;
  const plint n = mesh.getNumVertices();
  plb::Array<T, 3> c(0., 0., 0.);
  for (plint i = 0; i < n; i++) { c += mesh.getVertex(i); }
  c /= T(n);
  rmax = 0;
  rmin = 1e30;
  for (plint i = 0; i < n; i++) {
    const plb::Array<T, 3> d = mesh.getVertex(i) - c;
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    rmax = std::max(rmax, r);
    rmin = std::min(rmin, r);
  }
}

// Collective: metrics of every cell type. Between particle updates a block
// keeps only a 3-node envelope, so a cell spanning blocks is complete on no
// rank; every rank therefore contributes the vertices in its bulk (each vertex
// is in exactly one bulk, a vertex beyond the periodic boundary through its
// image) and the cells are reassembled on all ranks. `buffer` is reused.
std::vector<TypeMetrics> measureCells(HemoCell & hemocell, const std::vector<double> & rmax0,
                                      const std::vector<double> & rmin0, T periodX, std::vector<double> & buffer,
                                      std::vector<int> & types) {
  HemoCellFields & fields = *hemocell.cellfields;
  const unsigned int nTypes = fields.size();
  const unsigned int nCells = fields.number_of_cells;
  unsigned int maxVertices = 0;
  for (unsigned int k = 0; k < nTypes; k++) { maxVertices = std::max<unsigned int>(maxVertices, fields[k]->numVertex); }
  // Per vertex: x, y, z, |F|, v_x, present; per cell: cell type + 1 (0: no vertex).
  const unsigned int stride = 6;
  buffer.assign(size_t(nCells) * maxVertices * stride, 0.0);
  types.assign(nCells, 0);
  for (plint bid : fields.immersedParticles->getLocalInfo().getBlocks()) {
    HemoCellParticleField & pf = fields.immersedParticles->getComponent(bid);
    for (const HemoCellParticle & particle : pf.particles) {
      if (!pf.isContainedABS(particle.sv.position, pf.localDomain)) { continue; }
      const plint cell = fields.base_cell_id(particle.sv.cellId);
      double * v = &buffer[(size_t(cell) * maxVertices + particle.sv.vertexId) * stride];
      v[0] = particle.sv.position[0];
      v[1] = particle.sv.position[1];
      v[2] = particle.sv.position[2];
      v[3] = norm(particle.sv.force);
      v[4] = particle.sv.v[0];
      v[5] = 1.0;
      types[cell] = particle.sv.celltype + 1;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, buffer.data(), int(buffer.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, types.data(), int(types.size()), MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  std::vector<TypeMetrics> m(nTypes);
  std::vector<hemo::Array<T, 3>> p;
  for (unsigned int cell = 0; cell < nCells; cell++) {
    if (types[cell] == 0) { continue; }
    const unsigned int ctype = types[cell] - 1;
    const double * v = &buffer[size_t(cell) * maxVertices * stride];
    unsigned int present = 0;
    for (unsigned int i = 0; i < maxVertices; i++) { present += v[i * stride + 5] > 0 ? 1 : 0; }
    if (int(present) != fields[ctype]->numVertex) { continue; } // incomplete: counted as lost
    const CommonCellConstants & cc = fields[ctype]->mechanics->cellConstants;
    TypeMetrics & t = m[ctype];
    t.cells += 1;
    p.resize(present);
    hemo::Array<T, 3> c = {0., 0., 0.};
    double vx = 0;
    for (unsigned int i = 0; i < present; i++) {
      p[i] = {v[i * stride], v[i * stride + 1], v[i * stride + 2]};
      // Unwrap periodic images relative to vertex 0.
      while (p[i][0] - p[0][0] > 0.5 * periodX) { p[i][0] -= periodX; }
      while (p[0][0] - p[i][0] > 0.5 * periodX) { p[i][0] += periodX; }
      c += p[i];
      vx += v[i * stride + 4];
      const double f = v[i * stride + 3];
      if (!(std::isfinite(f) && std::isfinite(p[i][0]) && std::isfinite(p[i][1]) && std::isfinite(p[i][2]))) {
        t.nonFinite += 1;
      }
      t.maxForce = std::max(t.maxForce, f);
      if (f > param::f_limit) { t.cappedVertices += 1; }
    }
    t.vertices += present;
    c /= T(present);

    for (unsigned int e = 0; e < cc.edge_list.size(); e++) {
      const double s = norm(p[cc.edge_list[e][1]] - p[cc.edge_list[e][0]]) / cc.edge_length_eq_list[e] - 1.0;
      t.maxEdgeStrain = std::max(t.maxEdgeStrain, s);
      t.minEdgeStrain = std::min(t.minEdgeStrain, s);
    }
    double volume = 0;
    for (unsigned int k = 0; k < cc.triangle_list.size(); k++) {
      const hemo::Array<T, 3> a = p[cc.triangle_list[k][0]] - c;
      const hemo::Array<T, 3> b = p[cc.triangle_list[k][1]] - c;
      const hemo::Array<T, 3> d = p[cc.triangle_list[k][2]] - c;
      const double area = 0.5 * norm(crossProduct(b - a, d - a));
      t.maxAreaChange = std::max(t.maxAreaChange, std::fabs(area / cc.triangle_area_eq_list[k] - 1.0));
      volume += dot(a, crossProduct(b, d)) / 6.0;
    }
    t.maxVolumeChange = std::max(t.maxVolumeChange, std::fabs(volume / cc.volume_eq - 1.0));

    double rmax = 0, rmin = 1e30;
    for (const hemo::Array<T, 3> & q : p) {
      const double r = norm(q - c);
      rmax = std::max(rmax, r);
      rmin = std::min(rmin, r);
    }
    t.maxElongation = std::max(t.maxElongation, rmax / rmax0[ctype]);
    t.minCompression = std::min(t.minCompression, rmin / rmin0[ctype]);
    for (int d = 0; d < 3; d++) { t.centroid[d] += c[d]; }
    t.velocityX += vx / present;
  }
  return m;
}

template <typename V>
V readOptional(Config & cfg, const char * section, const char * key, V fallback) {
  try {
    return cfg[section][key].read<V>();
  } catch (std::invalid_argument &) { return fallback; }
}

} // namespace

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  std::unique_ptr<MultiScalarField3D<int>> flagMatrix;
  std::unique_ptr<VoxelizedDomain3D<T>> voxelizedDomain;
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(), (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                       (*cfg)["domain"]["refDirN"].read<int>(), (*cfg)["domain"]["refDir"].read<int>(), voxelizedDomain,
                       flagMatrix, (*cfg)["domain"]["blockSize"].read<int>(),
                       (*cfg)["domain"]["particleEnvelope"].read<int>());

  // Driving force for a given wall shear rate. For a Newtonian fluid the wall
  // shear stress of Poiseuille flow is tau_w = G R / 2 (G: force density); in
  // steady state it follows from the force balance on the whole pipe, with or
  // without cells, so with a cell-free layer at the wall the wall shear rate is
  // close to that of the cell-free flow.
  param::lbm_base_parameters(*cfg);
  Box3D section = flagMatrix->getBoundingBox();
  section.x1 = section.x0;
  param::pipe_radius = std::sqrt(T(computeSum<int>(*flagMatrix, section)) / PI);
  const T wallShearRate = (*cfg)["domain"]["wallShearRate"].read<T>(); // [1/s]
  const T shear_lbm = wallShearRate * param::dt;
  const T drivingForce = 2.0 * param::nu_lbm * shear_lbm / param::pipe_radius;
  param::u_lbm_max = shear_lbm * param::pipe_radius / 2.0;
  param::re = param::u_lbm_max * 2.0 * param::pipe_radius / param::nu_lbm;
  hlog << "(PipeflowWBC) pipe radius " << param::pipe_radius << " LU (" << param::pipe_radius * param::dx * 1e6
       << " um), wall shear rate " << wallShearRate << " 1/s, cell-free centreline velocity "
       << param::u_lbm_max * param::dx / param::dt << " m/s, Re " << param::re << endl;
  param::printParameters();

  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
      voxelizedDomain.get()->getMultiBlockManagement(), defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));
  defineDynamics(*hemocell.lattice, *flagMatrix.get(), (*hemocell.lattice).getBoundingBox(),
                 new BounceBack<T, DESCRIPTOR>(1.), 0);
  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1., plb::Array<T, 3>(0., 0., 0.));
  hemocell.lattice->initialize();

  hemocell.initializeCellfield();
  const int stepMaterial = (*cfg)["ibm"]["stepMaterialEvery"].read<int>();
  const int stepParticle = (*cfg)["ibm"]["stepParticleEvery"].read<int>();
  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.addCellType<WbcHighOrderModel>("WBC", WBC_SPHERE);
  for (const char * name : {"RBC", "PLT", "WBC"}) {
    hemocell.setMaterialTimeScaleSeparation(name, stepMaterial);
    hemocell.setInitialMinimumDistanceFromSolid(name, 0.5); // um
  }
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(stepParticle);
  if ((*cfg)["domain"]["repulsion"].read<int>()) {
    hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
    hemocell.setRepulsionTimeScaleSeperation(stepMaterial);
  }
  hemocell.setOutputs("RBC", {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE, OUTPUT_VELOCITY});
  hemocell.setOutputs("PLT", {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE, OUTPUT_VELOCITY});
  hemocell.setOutputs("WBC", {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE, OUTPUT_VELOCITY, OUTPUT_INNER_LINKS});
  hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_SHEAR_RATE});
  hemocell.setSystemPeriodicity(0, true);

  // Run length and intervals.
  const T wallStrain = (*cfg)["sim"]["wallStrain"].read<T>();
  const plint tcheck = (*cfg)["sim"]["tcheck"].read<plint>();
  if (tcheck <= 0 || tcheck % stepMaterial != 0) {
    hlog << "(PipeflowWBC) Error: tcheck must be a positive multiple of stepMaterialEvery" << endl;
    return 1;
  }
  const plint minIterations = (*cfg)["sim"]["minIterations"].read<plint>();
  plint tmax = std::max<plint>(minIterations, plint(std::ceil(wallStrain / shear_lbm)));
  tmax = std::min<plint>((tmax + tcheck - 1) / tcheck * tcheck, (*cfg)["sim"]["tmax"].read<plint>());
  const T outputStrain = (*cfg)["sim"]["outputStrain"].read<T>();
  plint tout =
      outputStrain > 0 ? std::max<plint>(tcheck, plint(std::round(outputStrain / shear_lbm / tcheck)) * tcheck) : 0;
  // Stability limits: beyond the singular points of the area (30 %) and volume
  // (10 %) terms the restoring force decreases, so the mesh is broken.
  const T limitEdgeStrain = readOptional<T>(*cfg, "stability", "maxEdgeStrain", 2.0);
  const T limitAreaChange = readOptional<T>(*cfg, "stability", "maxAreaChange", std::sqrt(MaxCellSurfaceAreaChange));
  const T limitVolumeChange = readOptional<T>(*cfg, "stability", "maxVolumeChange", std::sqrt(MaxCellVolumetricChange));
  const T limitVelocity = readOptional<T>(*cfg, "stability", "maxLatticeVelocity", 0.3);
  hlog << "(PipeflowWBC) " << tmax << " iterations (wall strain " << tmax * shear_lbm << "), check every " << tcheck
       << ", output every " << tout << endl;

  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }
  if (hemocell.iter == 0) {
    setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      plb::Array<T, DESCRIPTOR<T>::d>(drivingForce, 0.0, 0.0));
    for (plint i = 0; i < (*cfg)["parameters"]["warmup"].read<plint>(); ++i) { hemocell.lattice->collideAndStream(); }
  }

  const unsigned int nTypes = hemocell.cellfields->size();
  std::vector<double> rmax0(nTypes), rmin0(nTypes);
  for (unsigned int k = 0; k < nTypes; k++) {
    referenceExtent(*(*hemocell.cellfields)[k], rmax0[k], rmin0[k]);
    const CommonCellConstants & cc = (*hemocell.cellfields)[k]->mechanics->cellConstants;
    if (!cc.inner_edge_list.empty()) {
      T lmin = 1e30, lmax = 0;
      for (T l : cc.inner_edge_length_eq_list) {
        lmin = std::min(lmin, l);
        lmax = std::max(lmax, l);
      }
      hlog << "(PipeflowWBC) " << (*hemocell.cellfields)[k]->name << ": " << cc.inner_edge_list.size()
           << " inner links, rest length / reference diameter " << lmin / (2 * rmax0[k]) << " - "
           << lmax / (2 * rmax0[k]) << endl;
    }
  }
  const plb::Box3D box = hemocell.lattice->getBoundingBox();
  const T axisY = 0.5 * (box.y0 + box.y1), axisZ = 0.5 * (box.z0 + box.z1);
  const T periodX = box.getNx();
  std::vector<double> gatherBuffer;
  std::vector<int> gatherTypes;
  const int wbc = int((*hemocell.cellfields)["WBC"]->ctype);

  const bool root = global::mpi().getRank() == 0;
  std::ofstream csv;
  if (root) {
    csv.open((global::directories().getOutputDir() + "/stability.csv").c_str());
    csv << "iteration,time_s,wall_strain,fluid_umax_lbm,fluid_umean_mps";
    for (unsigned int k = 0; k < nTypes; k++) {
      const std::string n = (*hemocell.cellfields)[k]->name;
      for (const char * c : {"cells", "max_force_pN", "capped_fraction", "max_edge_strain", "min_edge_strain",
                             "max_area_change", "max_volume_change", "max_elongation", "min_compression"}) {
        csv << "," << n << "_" << c;
      }
    }
    csv << ",WBC_radial_offset_um,WBC_velocity_mps\n";
  }
  std::vector<double> initialCells;
  std::vector<TypeMetrics> worst(nTypes);
  const T topN = param::df * 1e12;
  std::string failure;

  while (hemocell.iter < tmax) {
    hemocell.iterate();
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      plb::Array<T, DESCRIPTOR<T>::d>(drivingForce, 0.0, 0.0));

    // Forces are fresh (not yet capped by spreadParticleForce) one iteration
    // after a material update.
    if ((hemocell.iter - 1) % tcheck == 0) {
      std::vector<TypeMetrics> m = measureCells(hemocell, rmax0, rmin0, periodX, gatherBuffer, gatherTypes);
      const FluidStatistics fluid = FluidInfo::calculateVelocityStatistics(&hemocell);
      if (initialCells.empty()) {
        for (const TypeMetrics & t : m) { initialCells.push_back(t.cells); }
      }
      std::ostringstream reason;
      if (!std::isfinite(fluid.max) || fluid.max > limitVelocity) {
        reason << "fluid velocity " << fluid.max << " (lattice)";
      }
      for (unsigned int k = 0; k < nTypes; k++) {
        const TypeMetrics & t = m[k];
        const std::string n = (*hemocell.cellfields)[k]->name;
        if (t.nonFinite > 0) { reason << n << ": non-finite position or force; "; }
        if (t.cells < initialCells[k]) { reason << n << ": " << initialCells[k] - t.cells << " cell(s) lost; "; }
        if (t.maxEdgeStrain > limitEdgeStrain) { reason << n << ": edge strain " << t.maxEdgeStrain << "; "; }
        if (t.maxAreaChange > limitAreaChange) { reason << n << ": area change " << t.maxAreaChange << "; "; }
        if (t.maxVolumeChange > limitVolumeChange) { reason << n << ": volume change " << t.maxVolumeChange << "; "; }
        TypeMetrics & w = worst[k];
        w.maxForce = std::max(w.maxForce, t.maxForce);
        w.cappedVertices = std::max(w.cappedVertices, t.vertices > 0 ? t.cappedVertices / t.vertices : 0.);
        w.maxEdgeStrain = std::max(w.maxEdgeStrain, t.maxEdgeStrain);
        w.minEdgeStrain = std::min(w.minEdgeStrain, t.minEdgeStrain);
        w.maxAreaChange = std::max(w.maxAreaChange, t.maxAreaChange);
        w.maxVolumeChange = std::max(w.maxVolumeChange, t.maxVolumeChange);
        w.maxElongation = std::max(w.maxElongation, t.maxElongation);
        w.minCompression = std::min(w.minCompression, t.minCompression);
      }
      const TypeMetrics & cell = m[wbc];
      const T wy = cell.cells > 0 ? cell.centroid[1] - axisY : 0, wz = cell.cells > 0 ? cell.centroid[2] - axisZ : 0;
      const T offset = std::sqrt(wy * wy + wz * wz) * param::dx * 1e6;
      const T wbcVelocity = cell.velocityX * param::dx / param::dt;
      if (root) {
        csv << hemocell.iter << "," << hemocell.iter * param::dt << "," << hemocell.iter * shear_lbm << "," << fluid.max
            << "," << fluid.avg * param::dx / param::dt;
        for (const TypeMetrics & t : m) {
          csv << "," << t.cells << "," << t.maxForce * topN << ","
              << (t.vertices > 0 ? t.cappedVertices / t.vertices : 0.) << "," << t.maxEdgeStrain << ","
              << t.minEdgeStrain << "," << t.maxAreaChange << "," << t.maxVolumeChange << "," << t.maxElongation << ","
              << t.minCompression;
        }
        csv << "," << offset << "," << wbcVelocity << std::endl;
      }
      if (!reason.str().empty()) {
        failure = reason.str();
        break;
      }
    }
    if (tout > 0 && hemocell.iter % tout == 0) { hemocell.writeOutput(); }
  }

  std::ostringstream summary;
  summary << std::setprecision(4);
  for (unsigned int k = 0; k < nTypes; k++) {
    const TypeMetrics & w = worst[k];
    summary << " | " << (*hemocell.cellfields)[k]->name << ": max force " << w.maxForce * topN
            << " pN, capped fraction " << w.cappedVertices << ", edge strain " << w.minEdgeStrain << " .. "
            << w.maxEdgeStrain << ", area change " << w.maxAreaChange << ", volume change " << w.maxVolumeChange
            << ", elongation " << w.maxElongation << ", compression " << w.minCompression;
  }
  if (!failure.empty()) {
    hlog << "(PipeflowWBC) RESULT UNSTABLE at iteration " << hemocell.iter << " (wall strain "
         << hemocell.iter * shear_lbm << "): " << failure << summary.str() << endl;
    hemocell.writeOutput();
    return 2;
  }
  hlog << "(PipeflowWBC) RESULT STABLE for " << hemocell.iter << " iterations (wall strain "
       << hemocell.iter * shear_lbm << ")" << summary.str() << endl;
  return 0;
}
