.. _packcells:

Initial cell positions (packCells)
==================================

``tools/packCells`` places cells at random positions and orientations in a
box, so that no two cells overlap, and writes one ``<type>.pos`` file per cell
type for :ref:`loadParticles() <reference/hemocell_api:Setting up the cells>` to read. Each
cell is represented by an enclosing ellipsoid, somewhat larger than the cell
itself; the ellipsoids are packed with a force-biased algorithm
:cite:`bargiel2014geometrical`, which pushes overlapping ellipsoids apart and
turns them while it slowly shrinks them, until none overlap.

Building
--------

``packCells`` has its own CMake project (it does not need Palabos). From the
repository root::

  cmake -S tools/packCells -B tools/packCells/build
  cmake --build tools/packCells/build
  tools/packCells/packCells --help

The binary is written next to the sources, to ``tools/packCells/packCells``
(``-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=<dir>`` writes it elsewhere; give each of
several build trees its own). OpenMP, if found, parallelises the contact
search; the result does not depend on it or on the number of threads
(``OMP_NUM_THREADS``). ``-DPACKCELLS_NATIVE_ARCH=OFF`` builds for other CPUs.

A first packing
---------------

::

  ./packCells 40 40 40 --hematocrit 0.3 --seed 1

packs RBCs at a hematocrit of 0.3, and platelets (0.07 per RBC), into a
periodic box of 40 x 40 x 40 µm, and writes ``RBC.pos`` and ``PLT.pos`` to the
current directory. ``--rbc <n>``, ``--plt <n>``, ``--wbc <n>`` and the other
type options give numbers of cells instead; ``--help`` lists them. The same
``--seed`` gives the same packing; without it a random seed is used and
printed.

The output ends with the final scale of the ellipsoids and what it means::

  Final scale of the ellipsoids: 0.9998 (start 1.0000, minimum 0.9310)
  The ellipsoids do not overlap at their nominal size: the cells start out of each other's interaction range.

The progress lines before it show, for every tenth iteration, the volume
fraction of the ellipsoids and their scale, both at the *outer* scale (the
current size of the ellipsoids, which shrinks) and at the *inner* one (the
scale at which the closest overlapping pair would touch). The packing is done
when the two are equal, i.e. nothing overlaps. The exit status is 0 when the
packing succeeded, 1 for invalid arguments, 2 if it was not done within
``--maxiter`` iterations and 3 if it stopped at the minimum scale (see below);
with 2 and 3 the files are written, but the cells may overlap.

The ellipsoids and their scale
------------------------------

.. list-table::
   :header-rows: 1

   * - Option
     - Cells (file)
     - Ellipsoid [µm]
     - Cell mesh [µm]
     - Mesh fit
   * - ``--rbc``
     - human RBCs (``RBC.pos``)
     - 8.4 x 4.4 x 8.4
     - 7.82 x 2.29 x 7.82
     - 0.931
   * - ``--plt``
     - human platelets (``PLT.pos``)
     - 2.7 x 1.17 x 2.7
     - 2.50 x 1.09 x 2.50
     - 0.929
   * - ``--wbc``
     - white blood cells (``WBC.pos``)
     - 8.4 x 8.4 x 8.4
     - 8.0
     - 0.952
   * - ``--vrbc``
     - stage V gametocytes (``vRBC.pos``)
     - 3.5 x 6.0 x 11.0
     - 2.94 x 4.83 x 10.8
     - 1.089
   * - ``--rbc_m``
     - mouse RBCs (``RBC_m.pos``)
     - 6.23 x 3.26 x 6.23
     - 5.80 x 1.70 x 5.80
     - 0.931
   * - ``--plt_m``
     - mouse wild-type platelets (``PLT_m.pos``)
     - 2.05 x 1.17 x 2.05
     - 1.90 x 1.08 x 1.90
     - 0.927
   * - ``--plt_mko``
     - mouse knockout platelets (``PLT_mko.pos``)
     - 3.68
     - 3.42
     - 0.929

The sizes are the full diameters along the cell's own x, y and z axes; the
meshes are those HemoCell builds from ``tools/cell_templates`` (RBC and PLT
lie with their thin axis along y). The ellipsoids are larger than the meshes:
they stand for the cells' effective size, which includes the interaction range
of the immersed boundary method (cells interact from about one kernel width
apart). ``--cell <name> <n> <dx> <dy> <dz>`` adds a type of your own with an
ellipsoid of these diameters.

The ellipsoids start at ``--start_scale`` (default 1) times this nominal size
and shrink until they no longer overlap. A final scale below 1 means that
some cells start within each other's interaction range: a relaxation that
makes high volume fractions reachable. Below the *mesh fit*, the scale at
which an ellipsoid still contains its cell's mesh, the meshes themselves may
intersect, so the packing stops there, with exit status 3: the default
``--min_scale`` is the largest mesh fit of the types packed (1 for ``--cell``
types). The vRBC ellipsoid is smaller than the gametocyte mesh, so vRBCs start
at 1.14.

The ellipsoids shrink by about 0.5 % of the start scale every 1000 iterations,
so the run time is bounded: from 1 to the RBC's 0.93 takes about 14,000
iterations. On 24 threads a 40 µm cube at Ht 0.3 packs at scale 1.00 in
about 100 iterations (well under a second), at Ht 0.45 at 0.96 in 0.7 s,
and at Ht 0.5 it reaches 0.93 and stops; a 100 x 50 x 50 µm box at Ht 0.45
takes about 3 s. A larger ``--start_scale`` leaves more room between the cells
but takes correspondingly longer.

Hematocrit
----------

``--hematocrit <Ht>`` adds as many human RBCs as fill that fraction of the box
(or of the cylinder, see below), counting each RBC with its effective volume
of 97 µm³, the 81 µm³ of its mesh plus the interaction range; this effective
size, not the mesh, matters for the collisions between cells and so for the
viscosity. Platelets are added with ``--plt_ratio`` per RBC (default 0.07). The
output lists, per type, the number of cells and their volume fraction by the
mesh volume, and for RBCs the hematocrit by the effective volume. Other types
(WBCs, mouse cells) are added with their numbers and do not enter the
hematocrit.

Using the packing in a case
---------------------------

* **File names.** Each file must be named after the cell type it is for, the
  name given to ``addCellType()``; rename them where they differ (e.g.
  ``RBC.pos`` to ``RBC_HO.pos``).
* **Positions** are in µm, measured from lattice node 0 of the case (in a
  pre-inlet, from the pre-inlet's lower corner). Each line holds the centre and
  the rotation in degrees about x, then y, then z, about the cell's centre, as
  HemoCell applies it to the cell's mesh.
* **Periodic axes.** Along a periodic axis the cell centres lie in the box,
  but cells may cross its faces. Make the box the periodic length of the case
  (number of nodes times ``dx``) and call
  ``hemocell.setPlaceCellsAcrossPeriodicBoundaries(true)`` before
  ``loadParticles()``; otherwise HemoCell skips the cells that cross the
  boundary (and says how many in the log). With a box of another size those
  cells would overlap cells at the other end.
* **Walls.** HemoCell does not place a cell that reaches a wall: any of its
  vertices within the ``setInitialMinimumDistanceFromSolid()`` distance of a
  wall node. Pack for the walls, as below, or many cells are lost and with them
  part of the hematocrit. The log says how many cells were placed (``placed N
  of M cells``).

Walls and pipes
~~~~~~~~~~~~~~~

``--walls <axes>`` closes the box along these axes (they are then not
periodic), ``--cylinder <axis> <R>`` confines the cells to a pipe of radius R
along that axis, through the centre of the box, and ``--offset <x> <y> <z>``
is added to the positions written, to move the box onto the fluid region of
the case. A cell touches a wall as it touches another cell. With ``dx`` in µm
and a minimum distance from solid of ``k`` nodes (the distance divided by
``dx``, rounded down):

* Flat walls at the nodes ``w0`` and ``w1`` of an axis: ``--walls`` with that
  axis, a box of ``(w1 - w0 - 1 - 2k) dx`` along it and the offset
  ``(w0 + 0.5 + k) dx``. For ``examples/parallelplanes`` (walls at the z
  nodes 0 and 199, ``dx`` 0.5 µm, ``k`` 1)::

    ./packCells 50 50 98 --hematocrit 0.3 --walls z --offset 0 0 0.75

* A pipe whose axis lies at the nodes ``(yc, zc)``, with the nearest wall node
  ``r`` nodes from it: ``--cylinder`` with the radius
  ``R = (r - (k + 0.5) sqrt(2)) dx``, a box of ``2R`` across it and the offset
  ``yc dx - R``, ``zc dx - R``. For ``examples/pipeflow`` (axis at the nodes
  26, 26, ``r`` 25, 103 nodes long, ``k`` 1)::

    ./packCells 51.5 22.9 22.9 --hematocrit 0.3 --cylinder x 11.45 --offset 0 1.55 1.55

With these, HemoCell placed every cell of both examples (except those across
the periodic boundary, which these drivers skip); periodic packings of the
whole domain lost 8 % of the cells at the walls of ``parallelplanes`` and
70 % of the RBCs in the pipe. The hematocrit is that of the packing region,
so the thin layer kept free at the walls is not counted.

Checking a packing
------------------

``tools/packCells/check_packing.py`` (Python 3 with numpy) checks the files
independently: for every pair of cells, also through the periodic images, the
scale at which their ellipsoids touch (Perram-Wertheim contact function), and
for the walls and the cylinder the scale at which an ellipsoid touches them.
Give it the box and walls as to ``packCells``::

  python3 tools/packCells/check_packing.py 51.5 22.9 22.9 --cylinder x 11.45 --offset 0 1.55 1.55 RBC.pos PLT.pos

It reports the smallest scales and fails (exit status 1) if one lies below the
minimum scale (``--min_scale``, by default the largest mesh fit of the
types). ``--cell <name> <dx> <dy> <dz>`` gives the ellipsoid of a type of your
own or of a renamed file. The CI packs and checks a periodic box, a slab and a
pipe with every change.
