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

#include <memory>

namespace pipeflow_test {

constexpr const char * pipeflow_config = "tests/validation/pipeflow/config_pipeflow.xml";
constexpr const char * pipeflow_geometry = "examples/pipeflow/tube.stl";
constexpr const char * rbc_name = "tests/validation/pipeflow/RBC";
constexpr const char * plt_name = "tests/validation/pipeflow/PLT";

struct PipeflowSetup {
  std::auto_ptr<plb::MultiScalarField3D<int>> flagMatrix;
  std::auto_ptr<hemo::VoxelizedDomain3D<T>> voxelizedDomain;
  plb::Array<T, 3> drivingForce;
};

// Same domain, lattice and cell types as Validation.Pipeflow, without loading
// particles (that differs between a fresh start and a restart).
inline void setupPipeflow(hemo::HemoCell & hemocell, PipeflowSetup & setup,
                          T rbcMinimumDistanceFromSolid = 0.5) {
  hemo::Config * cfg = hemocell.cfg;
  hemo::getFlagMatrixFromSTL(pipeflow_geometry,
                             (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                             (*cfg)["domain"]["refDirN"].read<int>(),
                             (*cfg)["domain"]["refDir"].read<int>(),
                             setup.voxelizedDomain, setup.flagMatrix,
                             (*cfg)["domain"]["blockSize"].read<int>(),
                             (*cfg)["domain"]["particleEnvelope"].read<int>());
  hemo::param::lbm_pipe_parameters((*cfg), setup.flagMatrix.get());

  hemocell.lattice = new plb::MultiBlockLattice3D<T, DESCRIPTOR>(
      setup.voxelizedDomain.get()->getMultiBlockManagement(),
      plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
      plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
      plb::defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
      new plb::GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/hemo::param::tau));
  defineDynamics(*hemocell.lattice, *setup.flagMatrix.get(), (*hemocell.lattice).getBoundingBox(),
                 new hemo::BounceBack<T, DESCRIPTOR>(1.), 0);
  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1., {0., 0., 0.});

  const T poiseuilleForce = 8 * hemo::param::nu_lbm * (hemo::param::u_lbm_max * 0.5)
                            / hemo::param::pipe_radius / hemo::param::pipe_radius;
  setup.drivingForce = plb::Array<T, 3>(poiseuilleForce, 0., 0.);

  hemocell.lattice->initialize();
  hemocell.initializeCellfield();
  hemocell.addCellType<hemo::RbcHighOrderModel>(rbc_name, RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(rbc_name, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid(rbc_name, rbcMinimumDistanceFromSolid);
  hemocell.addCellType<hemo::PltSimpleModel>(plt_name, ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(plt_name, (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
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

} // namespace pipeflow_test

#endif
