Configuration Files
===================

.. _config.xml:

Config.xml
----------

The configuration that is used at runtime by the case. The configuration file 
should be valid xml. The first tag should always be ``<hemocell>``. Within the 
hemocell you can specify any tag that can be used in your simulation after you 
initialized the config file like this:: 
   
  (*hemocell.cfg)["domain"]["<your attribute>"].read<type>() 

Most tags are used within the hemocell framework to configure options. However
some are used from the **case.cpp** file and thus are not present in all cases.
These options are denoted with a bold **case.cpp**. 


* ``<hemocell>`` outer tag, used for identification

  * ``<verbose>`` options related to verbosity

    * ``<cellsDeletedInfo>`` option to print the location of deleted cells,
      because this option needs to aggregate location data with an if statement it
      has a small performance impact and is disabled by default

  * ``<parameters>`` parameters regarding the simulation, some options might be
    moved but are still here because of legacy reasons.

    * ``<warmup>`` **case.cpp** files to let the fluid field warm
      up (usually develop a parabolic profile) before placing the cells into the
      flow. the argument is an integer and that many iterations of the fluid field
      are usually performed
    * ``<prctForced>`` **case.cpp** Percentage (between 0-1) of points of the
      cell over which a stretchingforce is applied. Only used in stretchCell.cpp
    * ``<stretchForce>`` **case.cpp** Total stretching force (in pN) applied to
      the cell. Only used in stretchCell.cpp
    * ``<outputDirectory>`` The base directory where output is saved, _x is
      appended if the directory already exists.
    * ``<checkpointDirectory>`` A relative directory (to the output directory)
      where the checkpoints (if any are requested) are saved.
    * ``<logDirectory>`` A directory relative to the output directory where the
      logfiles are saved
    * ``<logFile>`` The name of a logfile, if such a name exists then .x is
      appended (useful for restarting from a checkpoint)

  * ``<ibm>``

    * ``<radius>`` **case.cpp** Used only in oneCellShear as the original radius
      of the cell currently being sheared.
    * ``<stepMaterialEvery>`` **case.cpp** Update the particle material model after this many fluid time steps
    * ``<stepParticleEvery>`` **case.cpp** Update particle velocity after this many fluid time steps

  * ``<domain>``

    * ``<fluidEnvelope>`` **case.cpp** Legacy option, must be 2 if used
    * ``<geometry>`` **case.cpp** Used within the pipeflow case to denote the
      location of the stl file which is used to create the boundaries
    * ``<rhoP>`` Density of the fluid in SI units (kg/m³)
    * ``<nuP>``  Viscosity of the fluid (specifically the blood plasma in the case of HemoCell) in SI units (m²/s)
    * ``<dx>`` The length of a lattice unit in SI (m)
    * ``<dt>`` The duration of one timestep in LBM in SI units (s)
    * ``<refDir>`` **case.cpp** Used for determining reference direction of system when created from stl-file
    * ``<refDirN>`` **case.cpp** The number of lattice nodes in the refDir direction. used in
      conjunction with refDir. And for bodyforce calculations from ``<Re>`` as well
    * ``<blockSize>`` **case.cpp** Used to set a desired edge-size of an atomic block. Usefull in combination with load balancing
    * ``<kBT>`` the boltzmann constant times the temperature. in SI (m² kg s¯² (or J) for T=300)
    * ``<Re>`` Used for calculation of a bodyforce if used. **Note:** calculated
      bodyforce must still be applied within **case.cpp**, otherwise this has no
      effect
    * ``<particleEnvelope>`` This option is denoted in ``<dx>``. Should be a bit larger than the longest stretch of
      a particle in the current simulation. otherwise particles will be deleted.
      Usually a value of 25 is used, otherwise a warning is displayed.
    * ``<kRep>`` **case.cpp** Repulsion constant used for repulsion force. 
      Uncomment line in pipeflow.cpp if you want to use this.
    * ``<RepCutoff>`` **case.cpp** Cutoff distance in **micrometer!** for the
      repulsion force. Uncomment line in pipeflow.cpp if you want to use this.

  * ``<sim>``

    * ``<tmax>`` **case.cpp** Total number of iterations to run simulation
    * ``<tmeas>`` **case.cpp** Interval after wich data is written
    * ``<tcheckpoint>`` **case.cpp** Interval after which data is checkpointed
    * ``<tbalance>`` **case.cpp** Interval after which atomic blocks are balanced over processors, only in combination with load-balancing library


CELL.xml and CELL.pos
---------------------

The configuration that defines a celltype (.xml) and where those cells should be 
positioned (.pos). The ``<Cell>.xml`` file should be a valid xml file. The outer tag 
should be ``<hemocell>`` and within it should be a tag called ``<MaterialModel>``.
Within ``<Materialmodel>`` The following tags are used, also depending on the
material model used (defined as template in the addCellType function within a **case.cpp** file). For
example: The red blood cell material model will ignore inner edges.

  * **kBend** Bending force modulus for membrane + cytoskeleton (in **kBT** units)
  * **kVolume** Volume conservation coefficient (dimensionless)
  * **kArea** Local area conservation coefficient (dimensionless)
  * **kLink** Link force coefficient (dimensionless)
  * **minNumTriangles** Minimum number of triangles to create when not loading
    stl file, final number can be larger
  * **radius** Radius of the cell in SI units (m). used for scaling in hemocell 
  * **Volume** Volume of a cell in µm, only used for density output.
  * **enableInteriorViscosity** [0,1] use enable viscosity, should be used in
    combination with **viscosityRatio**
  * **viscosityRatio** ratio between interior and exterior viscosity
  * **eta_m** membrane viscosity, currently not used
  * **InnerEdges** contains **Edge** which contains two integers denoting which
    vertices in the model should have an inner edge between them.


The ``<Cell>.pos`` file should contain the number of cells (and thus number of following
lines) on the first line. Then each following line should contain 6 floats.
the first three deterimine the place in µm in X,Y,Z respectively and the last three determine
the rotation in degrees in X,Y,Z respectively

.. _constants:

config/constant_defaults.h
--------------------------

This file is used to define compile time constants for the hemocell library. The
hemocell library (in build/hemocell) is used to link all the **case.cpp** files
against. For example, whenever you want to use hemocell with interior viscosity
you must uncomment ``#define INTERIOR_VISCOSITY`` such that it is enabled. Below
we listed the options present in this file and when you can use them

* ``SOLIDIFY_MECHANICS`` Used for thrombus formation, relevant examples and
  source code is not yet available in V2.0
* ``INTERIOR_VISCOSITY`` Enable if you want to run cases with interior
  viscosity, adds two vectors to the hemocellparticle class, and thus has a
  measurable performance impact (don't enable when not needed)
* ``HEMOCELL_MATERIAL_INTEGRATION`` Defines how the velocity of the fluid is
  integrated to the particles. Euler [1] or Adams-Bashforth [2]. See
  ``src/hemoCellParticle.h`` for implementation details
* ``DESCRIPTOR`` The collision operator and dimensionality of the underlying
  lattice boltzmann fluid. This collision operator is only used in the palabos
  part of hemocell, find more information about it on `palabos.org`_.
* ``FORCE_LIMIT`` Limits the force the particles can exert on the fluid field.
  This means that a particle can deform more, but in return the fluid field
  stays stable. The force is in picoNewton.
* ``HEMOCELL_PARTICLE_FILED`` This has been built in and can't change anymore,
  it used to be interchangable with palabos particle fields very early on.
* ``OUTPUT_XXX`` Defines the outputs that can be requested in a **case.cpp** for
  either the fluid or the cells.
* ``T`` It is possible to define ``T`` as float instead of double, this
  decreases accuracy but increases speed, also decreases memory footprint.
* ``constructMeshElement``

  * ``RBC_FROM_SPHERE`` create the RBC model from mathematical equations (see
    hemocell paper), accepts a minimum number of to be created vertices.
  * ``ELLIPSOID_FROM_SPHERE`` mathimatically create a discretization which is
    mainly used for the PLT model.
  * ``STRING_FROM_VERTEXES`` legacy, not used anymore
  * ``WBC_SPHERE`` mathematically create a sphere in the form of a white blood
    cell
  * ``MESH_FROM_STL`` load the vertices from a stl file defined in the CELL.xml
    file 

.. _palabos.org: http://palabos.org
