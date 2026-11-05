Output fields
=============

You select what HemoCell writes to disk with two calls from your case:

* ``hemocell.setOutputs("RBC", { ... })`` — per-cell fields, **once per cell
  type**, and
* ``hemocell.setFluidOutputs({ ... })`` — fluid (lattice) fields.

Each takes a list of ``OUTPUT_*`` constants defined in
``config/constant_defaults.h``. Requesting only what you need keeps the output
smaller and the simulation faster. By default the values are written in SI units
(set ``hemocell.outputInSiUnits = false`` for lattice units). The fields end up
in compressed HDF5 files under ``tmp/hdf5``, described by XDMF files for
ParaView (see `Files`_), and cell scalar data additionally in ``tmp/csv`` (see
:ref:`read_output`).

Files
-----

Each ``writeOutput()`` writes, in the output directory (``tmp/`` by default):

``hdf5/<iteration>/<name>.<iteration>.p.<block>.h5``
   One file per atomic block for every cell type (``<name>`` is the cell type,
   without a directory part), the fluid (``Fluid``) and the CEPAC field
   (``CEPAC``); ``_PRE`` is appended in a pre-inlet. Datasets are single
   precision (the triangle and inner-link lists are integers), byte-shuffled
   and deflate-compressed. The attributes hold ``dx``, ``dt`` and the
   ``iteration``.
``<name>.<iteration>.xmf``
   XDMF (version 3) description of that step: all blocks with the time of the
   step.
``<name>.xmf``
   All steps of ``<name>`` in the directory, with their time (seconds, or
   iterations for lattice-unit output). Open this file in ParaView to get a
   time series.

The XDMF files are rewritten at every output, so a running simulation can be
inspected. ParaView reads them with its XDMF3 readers and with the legacy
XDMF reader. For output of HemoCell versions that did not write XDMF files,
generate them with :ref:`batchPostProcess.sh <bpp>`.

Cell files hold the vertices of the cells a block writes: every cell is
written once, by the block whose bulk contains the cell's vertex 0, with all
its vertices in vertex order (``numVertex`` rows per cell, cells in ascending
cell ID). A cell crossing a periodic boundary is written once, at the position
of its vertex 0, and may therefore extend past the domain boundary.
``Triangles`` and ``InnerLinks`` index the rows of the block's ``Position``.

Fluid files hold the values of the block's bulk lattice nodes as
``(z, y, x, component)`` arrays (attribute ``subdomainSize`` is ``(nz, ny,
nx)``). Every node is the centre of one cell of size ``dx`` of a regular grid
whose lower corner is ``relativePosition`` (``z, y, x``, the first node minus
``dx/2``), so the blocks tile the domain without overlap and every node is
written once.

.. note::

   Older versions wrote every cell from each block that held part of it (and
   both images of a cell crossing a periodic boundary), and wrote one extra
   node layer around each fluid block. Analysis scripts that removed
   duplicate cells or the extra layer are not needed for current output.

Per-cell outputs
----------------

Pass these to ``setOutputs(name, {...})``.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Constant
     - Field written
   * - ``OUTPUT_POSITION``
     - Position of each membrane vertex.
   * - ``OUTPUT_TRIANGLES``
     - The triangulated surface mesh (connectivity), needed to render the cell
       surface.
   * - ``OUTPUT_VELOCITY``
     - Velocity of each membrane vertex.
   * - ``OUTPUT_FORCE``
     - Total force on each vertex: the membrane force plus the repulsion force
       that will be applied to the fluid in the next iteration. For cell types
       with a material time-scale separation > 1 this is the force from the
       last material update, while the individual components below are
       evaluated for the current shape, so they need not add up exactly.
   * - ``OUTPUT_FORCE_VOLUME``
     - Volume-conservation force contribution.
   * - ``OUTPUT_FORCE_AREA``
     - Local area-conservation force contribution.
   * - ``OUTPUT_FORCE_BENDING``
     - Bending force contribution.
   * - ``OUTPUT_FORCE_LINK``
     - Link (stretching) force contribution.
   * - ``OUTPUT_FORCE_VISC``
     - Viscous force contribution.
   * - ``OUTPUT_FORCE_INNER_LINK``
     - Inner-link force contribution (models with internal structure).
   * - ``OUTPUT_FORCE_REPULSION``
     - Repulsion force contribution (requires repulsion enabled), as computed
       at the last repulsion update.
   * - ``OUTPUT_INNER_LINKS``
     - The inner-edge connectivity of the cell.
   * - ``OUTPUT_CELL_ID``
     - The unique id of the cell each vertex belongs to.
   * - ``OUTPUT_VERTEX_ID``
     - The id of each vertex.
   * - ``OUTPUT_RES_TIME``
     - Residence time of the cell in the domain.

.. note::

   ``OUTPUT_BINDING_SITES`` and ``OUTPUT_INTERIOR_POINTS`` are also defined and
   relate to the binding/solidification and interior-viscosity features
   respectively. They are only meaningful when those features are enabled; verify
   against the model you use before relying on them.

Fluid outputs
-------------

Pass these to ``setFluidOutputs({...})``.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Constant
     - Field written
   * - ``OUTPUT_VELOCITY``
     - Fluid velocity field.
   * - ``OUTPUT_DENSITY``
     - Fluid density field (kg/m\ :sup:`3` in SI output; for the CEPAC field
       the transported concentration, unscaled).
   * - ``OUTPUT_FORCE``
     - Total force field acting on the fluid: the external (body) force plus
       the force spread from the immersed cells. Writing it leaves the
       lattice force unchanged, so a body force set before ``writeOutput()``
       does not need to be re-applied.
   * - ``OUTPUT_SHEAR_STRESS``
     - Shear stress field.
   * - ``OUTPUT_SHEAR_RATE``
     - Shear rate field.
   * - ``OUTPUT_STRAIN_RATE``
     - Strain-rate field.
   * - ``OUTPUT_OMEGA``
     - Relaxation parameter (``omega = 1/tau``, dimensionless) field.
   * - ``OUTPUT_BOUNDARY``
     - Boundary/solid mask of the domain.
   * - ``OUTPUT_CELL_DENSITY``
     - Number density of cells mapped onto the fluid grid, one dataset
       ``CellDensity_<name>`` per cell type: the number of vertices nearest to
       each node, times the volume fraction of a vertex in SI output.

.. tip::

   The full list of ``OUTPUT_*`` constants lives in
   :ref:`config/constant_defaults.h <constants>`. If a constant is not listed
   above, it is either experimental or specific to an optional feature.
