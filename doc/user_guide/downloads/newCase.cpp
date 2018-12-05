// Include most of the interface offered by the hemocell library
#include <hemocell.h>
// This is the mechanical model for the cells that we want to use later on,
// alternatives can be found in the mechanics folder
#include <rbcHighOrderModel.h>
// These are functions found in the helpers folder, they are not in the core of
// hemoCell but can be handy nonetheless
#include <cellInfo.h>
#include <fluidInfo.h>

int main (int argc, char * argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  // The first argument is the config.xml location, the second and third argument
  // are necessary as a passthrough for the palabos initialization
  HemoCell hemocell(argv[1],argc,argv);

  // Calculate and load in the lattice boltzmann parameters from the config file
  // that will be used later on. Pretend that we are calculating the parameters
  // for a pipe, to get an acceptable maximum velocity.
  param::lbm_pipe_parameters((*hemocell.cfg),50);
  // Also print the parameters so we have visual confirmation.
  param::printParameters();

  // Although we are not creating a pipe, we still must define a driving force,
  // We pretend that this is a pipe, therefore the resulting velocity will be higher,
  // but acceptable. It is possible to analytically solve this correctly if you
  // want.
  T poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / param::pipe_radius / param::pipe_radius;

  // First we create a Palabos management object 
  // The first three arguments are the number of fluid cells in x,y and z
  // direction, so this is a 50x50x50 block, the fourth argument is the fluid
  // envelope size and must be two
  MultiBlockManagement3D management = defaultMultiBlockPolicy3D().getMultiBlockManagement(50, 50, 50, 2);

  // Initialize the fluid lattice within hemocell
  hemocell.initializeLattice(management);

  // Just to be sure disable all periodicity. Afterwards enable it in the
  // x-direction
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);

  //After we set up the fluid, it is time to set up the particles in the
  //simulation
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

  // Request outputs from the simulation, here we have requested all of the
  // possible outputs!
  hemocell.setOutputs("RBC_HO", { OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE, 
                                  OUTPUT_FORCE_VOLUME, OUTPUT_FORCE_BENDING, OUTPUT_FORCE_REPULSION, 
                                  OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC,
                                  OUTPUT_INNER_LINKS, OUTPUT_CELL_ID, OUTPUT_VERTEX_ID } );
  hemocell.setFluidOutputs( { OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_FORCE,
                              OUTPUT_SHEAR_RATE, OUTPUT_STRAIN_RATE,
                              OUTPUT_SHEAR_STRESS, OUTPUT_BOUNDARY, OUTPUT_OMEGA, 
															OUTPUT_CELL_DENSITY } );

  //Load the particles from all the *.pos files
  hemocell.loadParticles();


  // Load some basic values from the config.xml file that define how long the
  // simulation must run and when we want to save output
  unsigned int tmax = (*hemocell.cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*hemocell.cfg)["sim"]["tmeas"].read<unsigned int>();


  //This is the main running loop, run for tmax iterations.
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
