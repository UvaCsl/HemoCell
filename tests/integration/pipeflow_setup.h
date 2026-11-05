// Shared setup of the small pipeflow case (same domain, lattice and cell
// types as Validation.Pipeflow) for integration tests.
#ifndef HEMOCELL_TESTS_PIPEFLOW_SETUP_H
#define HEMOCELL_TESTS_PIPEFLOW_SETUP_H

#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "palabos3D.h"
#include "palabos3D.hh"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace pipeflow_test {

constexpr const char * pipeflow_config = "tests/validation/pipeflow/config_pipeflow.xml";
constexpr const char * pipeflow_geometry = "examples/pipeflow/tube.stl";
constexpr const char * rbc_name = "tests/validation/pipeflow/RBC";
constexpr const char * plt_name = "tests/validation/pipeflow/PLT";

struct PipeflowSetup {
  std::unique_ptr<plb::MultiScalarField3D<int>> flagMatrix;
  std::unique_ptr<plb::VoxelizedDomain3D<T>> voxelizedDomain;
  plb::Array<T, 3> drivingForce;
};

// Same domain, lattice and cell types as Validation.Pipeflow, without loading
// particles (that differs between a fresh start and a restart). `RbcModel` is
// the constitutive model of the RBC type; `rbc` and `plt` name the cell types
// (their .xml and .pos files).
template <class RbcModel = hemo::RbcHighOrderModel>
void setupPipeflow(hemo::HemoCell & hemocell, PipeflowSetup & setup, T rbcMinimumDistanceFromSolid = 0.5,
                   const std::string & rbc = rbc_name, const std::string & plt = plt_name) {
  hemo::Config * cfg = hemocell.cfg;
  hemo::getFlagMatrixFromSTL(pipeflow_geometry, (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                             (*cfg)["domain"]["refDirN"].read<int>(), (*cfg)["domain"]["refDir"].read<int>(),
                             setup.voxelizedDomain, setup.flagMatrix, (*cfg)["domain"]["blockSize"].read<int>(),
                             (*cfg)["domain"]["particleEnvelope"].read<int>());
  hemo::param::lbm_pipe_parameters((*cfg), setup.flagMatrix.get());

  hemocell.initializeLatticeFromFlagMatrix(setup.voxelizedDomain->getMultiBlockManagement(), *setup.flagMatrix);

  const T poiseuilleForce =
      8 * hemo::param::nu_lbm * (hemo::param::u_lbm_max * 0.5) / hemo::param::pipe_radius / hemo::param::pipe_radius;
  setup.drivingForce = plb::Array<T, 3>(poiseuilleForce, 0., 0.);

  hemocell.lattice->initialize();
  hemocell.initializeCellfield();
  hemocell.addCellType<RbcModel>(rbc, hemo::RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(rbc, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid(rbc, rbcMinimumDistanceFromSolid);
  hemocell.addCellType<hemo::PltSimpleModel>(plt, hemo::ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(plt, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());
  hemocell.setSystemPeriodicity(0, true);
}

inline void iterate(hemo::HemoCell & hemocell, const PipeflowSetup & setup, unsigned int steps) {
  for (unsigned int i = 0; i < steps; i++) {
    hemocell.iterate();
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, setup.drivingForce);
  }
}

// State of a run on this rank, for exact comparisons between runs with the
// same decomposition.
struct State {
  // Particles of every local block in storage order (the order sets the
  // summation order of force spreading), without the residence time, which
  // output brings up to date; cell IDs minus the `idShift` of captureState().
  std::vector<T> particles;
  // Populations and external force of the bulk nodes of every local block.
  std::vector<T> fluid;
};

inline State captureState(hemo::HemoCell & hemocell, plb::plint idShift = 0) {
  State state;
  plb::MultiParticleField3D<hemo::HemoCellParticleField> & field = *hemocell.cellfields->immersedParticles;
  for (plb::plint bid : field.getLocalInfo().getBlocks()) {
    for (const hemo::HemoCellParticle & p : field.getComponent(bid).particles) {
      for (const hemo::Array<T, 3> * a : {&p.sv.v, &p.sv.position, &p.sv.force, &p.sv.force_repulsion}) {
        state.particles.insert(state.particles.end(), a->begin(), a->end());
      }
      state.particles.push_back(T(p.sv.cellId - idShift));
      state.particles.push_back(T(p.sv.vertexId));
      state.particles.push_back(T(p.sv.celltype));
    }
  }
  plb::MultiBlockLattice3D<T, DESCRIPTOR> & lattice = *hemocell.lattice;
  for (plb::plint bid : lattice.getLocalInfo().getBlocks()) {
    plb::SmartBulk3D bulk(lattice.getMultiBlockManagement(), bid);
    const plb::Box3D box = bulk.toLocal(bulk.getBulk());
    plb::BlockLattice3D<T, DESCRIPTOR> & block = lattice.getComponent(bid);
    for (plb::plint x = box.x0; x <= box.x1; x++) {
      for (plb::plint y = box.y0; y <= box.y1; y++) {
        for (plb::plint z = box.z0; z <= box.z1; z++) {
          plb::Cell<T, DESCRIPTOR> & cell = block.get(x, y, z);
          for (plb::plint i = 0; i < DESCRIPTOR<T>::q; i++) { state.fluid.push_back(cell[i]); }
          const T * force = cell.getExternal(DESCRIPTOR<T>::ExternalField::forceBeginsAt);
          state.fluid.insert(state.fluid.end(), force, force + 3);
        }
      }
    }
  }
  return state;
}

// Index of the first difference (exact comparison), -1 if identical.
inline long firstDifference(const std::vector<T> & a, const std::vector<T> & b) {
  if (a.size() != b.size()) { return long(std::min(a.size(), b.size())); }
  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] != b[i]) { return long(i); }
  }
  return -1;
}

} // namespace pipeflow_test

#endif
