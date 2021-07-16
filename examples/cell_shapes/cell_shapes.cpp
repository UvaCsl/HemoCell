#include "cellInfo.h"
#include "fluidInfo.h"
#include "hemocell.h"
#include "particleInfo.h"
#include "pltSimpleModel.h"
#include "rbcHighOrderModel.h"
#include "wbcHighOrderModel.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config *cfg = hemocell.cfg;

  param::lbm_pipe_parameters((*cfg), 50);
  param::printParameters();

  T poiseuilleForce = 8 * param::nu_lbm * (param::u_lbm_max * 0.5) /
                      param::pipe_radius / param::pipe_radius;

  auto management =
      defaultMultiBlockPolicy3D().getMultiBlockManagement(50, 50, 50, 2);

  hemocell.initializeLattice(management);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0, true);

  Box3D topChannel(0, 49, 0, 49, 49, 49);
  Box3D bottomChannel(0, 49, 0, 49, 0, 0);
  Box3D backChannel(0, 49, 49, 49, 0, 49);
  Box3D frontChannel(0, 49, 0, 0, 0, 49);

  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR>);
  defineDynamics(*hemocell.lattice, bottomChannel,
                 new BounceBack<T, DESCRIPTOR>);
  defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR>);
  defineDynamics(*hemocell.lattice, frontChannel,
                 new BounceBack<T, DESCRIPTOR>);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.latticeEquilibrium(1., plb::Array<double, 3>(0., 0., 0.));
  hemocell.lattice->initialize();

  auto tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  auto tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();

  hemocell.initializeCellfield();

  // Red blood cell from sphere using RbcHighOrderModel
  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", 20);

  // Red blood cell using a custom shape read from STL directly
  hemocell.addCellType<RbcHighOrderModel>("RBC_FROM_STL", MESH_FROM_STL);
  hemocell.setMaterialTimeScaleSeparation("RBC_FROM_STL", 20);

  // Platelet from ellipsoid
  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", 20);

  // Platelet from ellipsoid using higher resolution discretisation
  hemocell.addCellType<PltSimpleModel>("PLT_HO", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT_HO", 20);

  // White blood cell from a Spherical model
  hemocell.addCellType<WbcHighOrderModel>("WBC_HO", WBC_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("WBC_HO", 20);

  // Only update the integrated velocity every X timesteps.
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(5);

  // Apply outputs to fluid field and all cell types
  hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_BOUNDARY, OUTPUT_FORCE});
  hemocell.setOutputs("RBC_HO", {OUTPUT_POSITION, OUTPUT_TRIANGLES});
  hemocell.setOutputs("RBC_FROM_STL", {OUTPUT_POSITION, OUTPUT_TRIANGLES});
  hemocell.setOutputs("PLT", {OUTPUT_POSITION, OUTPUT_TRIANGLES});
  hemocell.setOutputs("PLT_HO", {OUTPUT_POSITION, OUTPUT_TRIANGLES});
  hemocell.setOutputs("WBC_HO", {OUTPUT_POSITION, OUTPUT_TRIANGLES});

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);

  // Load the particles from all the *.pos files
  hemocell.loadParticles();

  while (hemocell.iter < tmax) {
    // Advance the fluid field and cellfields one tick.
    hemocell.iterate();

    // Set driving force as required after each iteration
    setExternalVector(
        *hemocell.lattice, hemocell.lattice->getBoundingBox(),
        DESCRIPTOR<T>::ExternalField::forceBeginsAt,
        plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    if (hemocell.iter % tmeas == 0) {
      hlog << "Iteration : " << hemocell.iter << " # of cells: "
           << CellInformationFunctionals::getTotalNumberOfCells(&hemocell)
           << std::endl;
      hemocell.writeOutput();
    }
  }
  return 0;
}
