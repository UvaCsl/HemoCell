Structure of a HemoCell Case
===================

Inside a case folder you will find various files that make up a specific case.


Case config.xml
---------------

The configuration that is used at runtime by the case. The configuration file
should be valid xml. The first tag should always be ``<hemocell>``. Within the
hemocell you can specify any tag that can be used in your simulation after you
initialized the config file like this::
  
  (*hemocell.cfg)["domain"]["<your attribute>"].read<type>()

A few tags have a predefined meaning and can be used by HemoCell to determine
other parameters.
  
  * **ibm**

    * **stepMaterialEvery** Update the particle material model after this many fluid time steps
    * **stepParticleEvery** Update particle velocity after this many fluid time steps
  
  * **domain**

    * **fluidEnvelope** Palabos specific, envelope for the fluid, should be 1
    * **rhoP** Density of the fluid in SI units (kg/m³)
    * **nuP**  Viscosity of the fluid (specifically the blood plasma in the case of HemoCell) in SI units (m²/s)
    * **dx** The length of a lattice unit in SI (m)
    * **dt** The duration of one timestep in LBM in SI units (s)
    * **refDir** Used for determining reference direction of system when created from stl-file
    * **refDirN** The number of lattice nodes in the refDir direction. used in conjunction with refDir. And for bodyforce calculations from **Re** as well
    * **blockSize** Used to set a desired edge-size of an atomic block. Usefull in combination with load balancing
    * **kBT** the boltzmann constant times the temperature. in SI (m² kg s¯² (or J) for T=300) 
    * **Re** Used for calculation of a bodyforce if used
    * **particleEnvelope** Should be a bit larger than the longest stretch of a particle in the current simulation. otherwise particles will be deleted
    * **kRep** Repulsion constant used for repulsion force
    * **RepCutoff** Cutoff distance in **micrometer!** for the repulsion force

  * **sim** 

    * **tmax** Total number of iterations to run simulation
    * **tmeas** Interval after wich data is written
    * **tcheckpoint** Interval after which data is checkpointed
    * **tbalance** Interval after which atomic blocks are balanced over processors, only in combination with load-balancing library

Case <Cell>.xml and <Cell>.pos
---------------

The configuration that defines a celltype (.xml) and where those cells should be 
positioned (.pos). The ``<Cell>.xml`` file should be a valid xml file. The outer tag 
should be ``<hemocell>`` and within it should be a tag called ``<MaterialModel>``.
Within ``<Materialmodel>`` The following tags are usually used:

  * **kBend** Bending force modulus for membrane + cytoskeleton (in **kBT** units)
  * **kVolume** Volume conservation coefficient (dimensionless)
  * **kArea** Local area conservation coefficient (dimensionless)
  * **kLink** Link force coefficient (dimensionless)
  * **minNumTriangles** Minimum number of triangles to divide stl or mesh in, final number can be larger
  * **radius** Radius of the cell in SI units (m). used for scaling in hemocell 
  * **Volume** Volume of a cell in µm, only used for density output.

The ``<Cell>.pos`` file should contain the number of cells (and thus number of following
lines) on the first line. Then each following line should contain 6 floats.
the first three deterimine the place in µm in X,Y,Z respectively and the last three determine
the rotation in degrees in X,Y,Z respectively

Case output folder
---------------

The output of a case is usually written to the ``<case>/tmp`` folder. The
checkpoints are the ``.xml`` and ``.dat`` files. When a new checkpoint is
created they are moved to ``.xml.old and ``.dat.old``. The hdf5 output is stored
per timestep in ``tmp/hdf5`` and the csv output in ``tmp/csv``. See
:any:`read_output` and :any:`bpp` for more info.

Resuming from a checkpoint
--------------------------

To resume from a checkpoint you should run the executable from the directory you
ran it originally from (so the directory with the ``.xml`` and ``.pos`` files
visible. The first argument should be ``tmp/checkpoint.xml`` instead of
``config.xml``. HemoCell should then automatically resume from the last saved
checkpoint.

.. note::
  
  The number of processors on which you run the case doesn't need to be the
  same!
