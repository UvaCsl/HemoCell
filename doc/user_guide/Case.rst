Creating your own HemoCell Case
===============================

This tutorial assumes that you have already compiled the HemoCell library
following :ref:`from_source`.

To add a new HemoCell case, you can follow the following steps:

- Duplicate the template directory ``hemocell/examples/template`` with the name
  of your case::

    user@local: ~/hemocell/$ cp -r ./examples/template ./examples/newCase

  This provides a directory ``case`` with an already defined ``CMakeLists.txt``.
- Register your case for compilation in ``hemocell/examples/CMakeLists.txt`` by
  appending a line containing ``add_subdirectory("newCase")``. This ensures that
  ``CMake`` picks up the directory of your example during the compilation
  process and links it to the required libraries::

    user@local: ~/hemocell/examples$ cp -r template newCase
    user@local: ~/hemocell/examples$ echo 'add_subdirectory("newCase")' >> CMakeLists.txt

- Within this directory the main file should be defined. Here, we assume the
  file is given the same name as the directory::

    user@local:~/hemocell/examples/case$ touch newCase.cpp

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
  // HemoCell but can be handy nonetheless
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

The first thing that should be done in a HemoCell case is initializing the
HemoCell object:

.. code-block:: c++

  // The first argument is the config.xml location, the second and third argument
  // are necessary as a passthrough for the Palabos initialization
  HemoCell hemocell(argv[1],argc,argv);


Afterwards we must define the parameters for the lattice Boltzmann simulation.
These are read in from the ``config.xml`` file:

.. code-block:: c++

  // Calculate and load in the lattice Boltzmann parameters from the config file
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
can be found in ``pipeflow.cpp`` within the :ref:`cases/pipeflow:Pipe flow` case.

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
  // The first argument must correspond with the CELL.xml and CELL.pos present in
  // the directory (where CELL is the string input).
  // The second argument defines how a cell is build up. see
  // config/constant_defaults.h for options.
  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);

  // Only update the forces resulting from the mechanical deformation every X
  // timesteps, recalculating this is the most costly step and since our
  // timestep is so small it can be done intermittently
  hemocell.setMaterialTimeScaleSeparation("RBC", 20);

  // Only update the integrated velocity (from the fluid field to the particles)
  // every X timesteps.
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(5);

  // Request outputs from the simulation, here we have requested all of the
  // possible outputs!
  hemocell.setOutputs("RBC", { OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_FORCE,
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

Now there is only one more xml file missing, namely the RBC.xml file.
Fortunately this file is included in the examples folder, you can copy it to the
newCase as following:

.. code-block:: console

    user@local:~/hemocell/examples$ cp RBC_template.xml newCase/RBC.xml
    user@local:~/hemocell/examples$ ls newCase/
    CMakeLists.txt  config.xml  newCase.cpp  RBC.xml

Creating the initial positions for the Cells
--------------------------------------------

As a final touch we must create an RBC.pos file which contains the positions
of the RBC's that we want in our simulation. For this we use the tool that is
described in :ref:`packcells`. Run packCells with the following command to
create only RBC in a 25x25x25 domain:

.. code-block:: console

    user@local:~/hemocell/tools/packCells$ ./packCells
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
    user@local:~/hemocell/tools/packCells$ ./packCells  25 25 25 --plt_ratio 0 --hematocrit 0.3 -r
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
    user@local:~/hemocell/tools/packCells$ cp RBC.pos ../../examples/newCase/RBC.pos

With the RBC.pos file present in the newCase directory all the pieces should
be there to run our first newly created case!

Running our newly created case
------------------------------

Finally everything should be in place! confirm this by executing the following
command and checking if you get similar output::

    user@local:~/hemocell/examples$ ls newCase/
    CMakeLists.txt  config.xml  newCase.cpp  RBC.pos RBC.xml

To compile the case, creating the ``hemocell/examples/newCase/newCase``
executable, run the following commands from ``hemocell/``::

    user@local:~/hemocell/ mkdir build
    user@local:~/hemocell/ cd build
    user@local:~/hemocell/ cmake ..
    user@local:~/hemocell/ cmake --build . --parallel $(nproc) --target newCase

Or alternatively evaluate the ``./compile.sh`` script from within the
``hemocell/examples/newCase`` directory which automates those steps.

Afterwards, the executable ``newCase`` should be located at
``hemocell/examples/newCase/``. The simulation can then be started as::

    user@local:~/hemocell/examples/newCase/$ mpirun -n $(nproc) ./newCase config.xml

where the number of used cores can be set by the ``-n`` flag. Once the
simulation completes, its output will be stored in ``newCase/tmp/``. See
:ref:`read_output` on how to parse and interpret the generated output.

Compiling examples with special features
----------------------------------------

By default the examples are linked to the default HemoCell library
(``libhemocell.a``). However, some examples require additional features to be
enabled in the library. Currently, these consist of interior viscosity,
solidification mechanics, and load-balancing through ``Parmetis``. To enable
these features, you are required to change the compilation target of your
example. This can be done by modifying the ``target_link_libraries`` command in
the ``CMakeLists.txt`` located in your example directory. To enable either of
those options, please change the line accordingly:

- Default: ``target_link_libraries(${EXEC} ${PROJECT_NAME})``
- Interior viscosity: ``target_link_libraries(${EXEC} ${PROJECT_NAME}_interior_viscosity)``
- Solidify mechanics: ``target_link_libraries(${EXEC} ${PROJECT_NAME}_solidify_mechanics)``
- Load-balancing: ``target_link_libraries(${EXEC} ${PROJECT_NAME}_parmetis)``

This changes the target library of your example to the corresponding library with
the required features enable.

Note to enable load-balancing features, the optional dependency ``Parmetis``
should be present on the system (see :ref:`from_source`).
