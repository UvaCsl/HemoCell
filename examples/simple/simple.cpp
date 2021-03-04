#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

int main (int argc, char * argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1],argc,argv);
  Config * cfg = hemocell.cfg;

  param::lbm_pipe_parameters((*cfg), 50);
  param::printParameters();

  T poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / param::pipe_radius / param::pipe_radius;

  MultiBlockManagement3D management = defaultMultiBlockPolicy3D().getMultiBlockManagement(50, 50, 50, 2);

  hemocell.initializeLattice(management);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);

  Box3D topChannel(0, 49, 0, 49, 49, 49);
  Box3D bottomChannel( 0, 49, 0, 49, 0, 0);
  Box3D backChannel( 0, 49, 49, 49, 0, 49);
  Box3D frontChannel( 0, 49, 0, 0, 0, 49);

  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, frontChannel, new BounceBack<T, DESCRIPTOR> );

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));
  hemocell.lattice->initialize();

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();

  hemocell.initializeCellfield();

  // Add a particleType to the simulation, the template argument refers to the
  // corresponding mechanics in the mechanics/ folder
  // The first argument must correspont with the CELL.xml and CELL.pos present in
  // the directory (where CELL is the string input).
  // The second argument defines how a cell is build up. see
  // config/constant_defaults.h for options.
  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);

  // Only update the forces resulting from the mechanical deformation every X
  // timesteps, recalculating this is the most costly step and since our
  // timestep is so small it can be done intermittently
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", 20);

  // Only update the integrated velocity (from the fluid field to the particles)
  // every X timesteps.
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(5);

  hemocell.setFluidOutputs( { OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_FORCE,
                              OUTPUT_SHEAR_RATE, OUTPUT_STRAIN_RATE,
                              OUTPUT_SHEAR_STRESS, OUTPUT_BOUNDARY, OUTPUT_OMEGA,
                              OUTPUT_CELL_DENSITY } );

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);

  //Load the particles from all the *.pos files
  // hemocell.loadParticles();

  while (hemocell.iter < tmax ) {
    //Advance the fluid field and cellfields one tick.
    hemocell.iterate();

    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    // When we want to save
    if (hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();
    }
  }
  return 0;
}
