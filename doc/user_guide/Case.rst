Creating your own HemoCell Case
===============================

This tutorial assumes that you have already compiled the HemoCell library,
either by following :ref:`singularity` or by following :ref:`from_source`.

To create a new hemoCell case it is the easiest to create a new folder within
the ``cases`` directory with the name of your case. When this folder is created
you can run ``make cmakefiles`` in the cases directory to create the cMakeFile
in the new directory.

.. code-block:: console

    vikko@the9:~/HemoCell/cases$ ls
    CMakeLists_template.txt  oneCellShear  PLT_template.xml     stretchCell
    Makefile                 pipeflow      RBC_HO_template.xml
    vikko@the9:~/HemoCell/cases$ mkdir newCase
    vikko@the9:~/HemoCell/cases$ make cmakefiles 
    cp  CMakeLists_template.txt newCase/CMakeLists.txt
    sed -i 's/FOLDER_NAME__/newCase/g' newCase/CMakeLists.txt
    vikko@the9:~/HemoCell/cases$ ls newCase/
    CMakeLists.txt

Within this directory there must be a ``.cpp`` file with the same name, so:

.. code-block:: console
    
    vikko@the9:~/HemoCell/cases/newCase$ touch newCase.cpp

Editing your newCase.cpp
------------------------

Firstly we have to include the headers we need in ``newCase.cpp``.

.. code-block:: c++

  // Include most of the interface offered by the HemoCell library
  #include <hemocell.h> 
  // This is the mechanical model for the cells that we want to use later on,
  // alternatives can be found in the mechanics folder
  #include <rbcHighOrderModel.h>
  // These are functions found in the helpers folder, they are not in the core of
  // hemoCell but can be handy nonetheless
  #include <cellInfo.h>
  #include <fluidInfo.h>

After that we have to define a main function as this will become the main
program. it is useful to include the ``config.xml`` file as first argument,
and also check for it:

.. code-block:: c++

  int main (int argc, char * argv[]) {
    if(argc < 2) {
      cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
      return -1;
    }
  }

The first thing that should be done in a hemoCell case is initializing the
HemoCell object:

.. code-block:: c++

  // The first argument is the config.xml location, the second and third argument
  // are necessary as a passthrough for the palabos initialization
  HemoCell hemocell(argv[1],argc,argv);


Afterwards we must define the parameters for the lattice boltzmann simulation.
These are read in from the ``config.xml`` file:

.. code-block:: c++
  
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


Since we want to create the simplest possible case we do not load in any stl
file but just create a cube with one periodic direction. An example of how to load in a stl file
can be found in ``pipeflow.cpp`` within the :ref:`pipeflow` case.

.. code-block:: c++

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
  // Set up bounceback boundaries in the other directions
  Box3D topChannel(0, 49, 0, 49, 49, 49);
  Box3D bottomChannel( 0, 49, 0, 49, 0, 0);
  Box3D backChannel( 0, 49, 49, 49, 0, 49);
  Box3D frontChannel( 0, 49, 0, 0, 0, 49);

  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, frontChannel, new BounceBack<T, DESCRIPTOR> );
  //Disable statistics to run faster
  hemocell.lattice->toggleInternalStatistics(false);
  //Equilibrate everything
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));
  //Finalize everything
  hemocell.lattice->initialize();

Then we set up the rest of the simulation, the comments should explain
everything:

.. code-block:: c++

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

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);

  //Load the particles from all the *.pos files
  hemocell.loadParticles();


  // Load some basic values from the config.xml file that define how long the
  // simulation must run and when we want to save output
  unsigned int tmax = (*hemocell.cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*hemocell.cfg)["sim"]["tmeas"].read<unsigned int>();


Finally we come to the main running loop, this case is very simple and has no
checkpointing etc. built in, these features can be found in the other example
cases:

.. code-block:: c++

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


You can download this file from :download:`here <downloads/newCase.cpp>`

Creating a bare config.xml
--------------------------

For this case we have minimalized the values read from the config.xml file. This
means that the following config file is enough to run our newCase.

.. code-block:: xml

  <?xml version="1.0" ?>
  <hemocell>

  <domain>
      <rhoP> 1025 </rhoP>   <!--Density of the surrounding fluid, Physical units [kg/m^3]-->
      <nuP> 1.1e-6 </nuP>   <!-- Kinematic viscosity of blood plasma, physical units [m^2/s]-->
      <dx> 5e-7 </dx> <!--Physical length of 1 Lattice Unit -->
      <dt> 1e-7 </dt> <!-- Time step for the LBM system. A negative value will set Tau=1 and calc. the corresponding time-step. -->
      <kBT> 4.100531391e-21 </kBT> <!-- in SI, m2 kg s-2 (or J) for T=300 -->
      <Re> 1.5 </Re>   <!--Reynolds number-->
      <particleEnvelope> 25 </particleEnvelope>
  </domain>

  <sim>
      <tmax> 50000 </tmax> <!-- total number of iterations -->
      <tmeas>  500 </tmeas> <!-- interval after which data is written -->
  </sim>

  </hemocell>

Now there is only one more xml file missing, namely the RBC_HO.xml file.
Fortunately this file is included in the cases folder, you can copy it to the
newCase as following:

.. code-block:: console

    vikko@the9:~/HemoCell/cases$ cp RBC_HO_template.xml newCase/RBC_HO.xml
    vikko@the9:~/HemoCell/cases$ ls newCase/
    CMakeLists.txt  config.xml  newCase.cpp  RBC_HO.xml

Creating the initial positions for the Cells
--------------------------------------------

As a final touch we must create an RBC_HO.pos file which contains the positions
of the RBC's that we want in our simulation. For this we use the tool that is
described in :ref:`packcells`. Run packCells with the following command to
create only RBC in a 25x25x25 domain:

.. code-block:: console

    vikko@the9:~/HemoCell/tools/packCells$ ./packCells
    Insufficient arguments.

    USAGE: packCells sX sY sZ [OPTIONAL ARGUMENTS ...]

    OPTIONAL ARGUMENTS:
      --hematocrit <0-1.0>                 -h The hematocrit of the solution
      --plt_ratio <ratio>                     The ratio of PLT per RBC, default=0.07
      --rbc <n>                               Number of Red Blood Cells
      --plt <n>                               Number of Platelets
      --wbc <n>                               Number of White Blood Cells
      --vrbc <n>                              Number of Stage V gametocytes
      --cell <name> <n> <e1, e2, diameter>    Custom Celltype described by ellipsoid
      --allowRotate                        -r Allow for rotation of ellipsoids
      --scale <ratio>                         Scales the neighbourhood grid (only change this if you know what you are doing!)
      --maxiter <n>                           Maximum number of iterations
      --help                                  Print this
    OUTPUT:
      <Cell>.pos for every celltype. First line is the number of cells.
      The rest of the lines is the cells in "Location<X Y Z> Rotation<X Y Z>" format.
      Cells.pov for visualization in, for example, povray

    NOTE:
      sX, sY and sZ are the domain size
      sX, sY, sZ and output are in micrometers[µm]
      --hematocrit and --RBC are mutually exclusive
      --hematocrit and --PLT are mutually exclusive
      --PLT-ratio is an No-Op without --hematocrit
    vikko@the9:~/HemoCell/tools/packCells$ ./packCells  25 25 25 --plt_ratio 0 --hematocrit 0.3 -r
    Loaded parameters, we found:
      Domain Size (µm): ( 25.000000 , 25.000000 , 25.000000 )
      Maximum Iterations : 2147483547
      Scale              : 0.250000
      Rotation           : 1
      Hematocrit    : 0.300000
      PLT/RBC Ratio : 0.000000
    We have found the following Cells:
      RBC
        No   : 48
        Sizes: (8.400000 , 4.400000 , 8.400000 )

    Nominal requested volume fraction: 0.499380

         Steps     Actual       Nominal        Inner         Outer             Force
                  density       density       diameter      diameter       per particle

        68764  0.1604380013  0.1604380013  1.2985355219  1.2985355219  0.000000000000000 PACKING DONE 
    vikko@the9:~/HemoCell/tools/packCells$ cp RBC.pos ../../cases/newCase/RBC_HO.pos

With the RBC_HO.pos file present in the newCase directory all the pieces should
be there to run our first newly created case!

Running our newly created case
------------------------------

Finally everything should be in place! confirm this by executing the following
command and checking if you get similar output:

.. code-block:: console

    vikko@the9:~/HemoCell/cases$ ls newCase/
    CMakeLists.txt  config.xml  newCase.cpp  RBC_HO.pos RBC_HO.xml

Compile our case by executing the folling commands, replace X by the number of
cores you want to run on:

.. code-block:: console 

    vikko@the9:~/HemoCell/cases/newCase$ mkdir build
    vikko@the9:~/HemoCell/cases/newCase$ cd build
    vikko@the9:~/HemoCell/cases/newCase/build$ cmake ../
    vikko@the9:~/HemoCell/cases/newCase/build$ make -j4
    vikko@the9:~/HemoCell/cases/newCase/build$ cd ../
    vikko@the9:~/HemoCell/cases/newCase/$ mpirun -n X ./newCase config.xml

Finally the output should be stored in ``tmp/``. see :ref:`read_output` on how
to parse this output.
