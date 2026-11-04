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
in compressed HDF5 files under ``tmp/hdf5`` (convert them for ParaView with
:ref:`batchPostProcess.sh <bpp>`), and cell scalar data additionally in
``tmp/csv`` (see :ref:`read_output`).

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
     - Fluid density field.
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
     - Relaxation parameter (``omega = 1/tau``) field.
   * - ``OUTPUT_BOUNDARY``
     - Boundary/solid mask of the domain.
   * - ``OUTPUT_CELL_DENSITY``
     - Number density of cells mapped onto the fluid grid.

.. tip::

   The full list of ``OUTPUT_*`` constants lives in
   :ref:`config/constant_defaults.h <constants>`. If a constant is not listed
   above, it is either experimental or specific to an optional feature.
