Mechanical (cell) models
========================

A *mechanical model* defines the constitutive (force) behaviour of a cell type.
You select one as the template argument to
``hemocell.addCellType<Model>(name, constructType)``. All models live in the
``mechanics/`` directory, derive from ``CellMechanics``, and implement a
``ParticleMechanics`` method that computes the membrane forces from the cell's
current shape (see :ref:`concepts/cell_model:The cell model`).

Available models
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Class (header)
     - Typical use
     - Description
   * - ``RbcHighOrderModel`` (``rbcHighOrderModel.h``)
     - Red blood cells
     - The validated, optimised RBC model. Uses the link, area, bending, and
       volume force terms (``kLink``, ``kArea``, ``kBend``, ``kVolume``). The
       recommended starting point for RBCs :cite:`Zavodszky:2017,tarksalooyeh2019optimizing`.
   * - ``PltSimpleModel`` (``pltSimpleModel.h``)
     - Platelets (PLT)
     - A simpler model for the small, stiff platelet, usually constructed with
       ``ELLIPSOID_FROM_SPHERE``. Shares the same four force terms but is tuned
       (and allowed a much larger bending angle) for a near-rigid ellipsoid.
   * - ``WbcHighOrderModel`` (``wbcHighOrderModel.h``)
     - White blood cells / nucleated cells
     - Extends the RBC-style membrane with an inner rigid core
       (``k_inner_rigid``, ``k_cytoskeleton``, ``core_radius``), modelling a
       generic single-nucleus cell with a stiff nucleus.
   * - ``RbcMalariaModel`` (``rbcMalariaModel.h``)
     - Infected / stiffened RBCs
     - An RBC variant with an additional inner-link term (``k_inner_link``) used
       to study malaria-infected red blood cells :cite:`czaja2020influence`.
   * - ``NoOp`` (``NoOp.h``)
     - Passive tracers
     - A null model: ``ParticleMechanics`` does nothing, so the cell experiences
       no internal forces. Its vertices follow the fluid like tracers and the
       cell does not keep its shape (it is not rigid). The cell type still
       needs its ``<name>.xml``, for the mesh.

Force coefficients
------------------

The coefficients shared by the deformable models are set per cell type in the
``CELL.xml`` file (see :ref:`xml_files:CELL.xml and CELL.pos`):

``kLink``, ``kArea``, ``kBend``, ``kVolume``, and (where used) ``k_inner_link`` /
``k_inner_rigid``. The hard limits on each deformation mode are compiled in as
the ``cellgroup`` constants in ``constant_defaults.h``:

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Constant
     - Meaning
   * - ``MaxCellVolumetricChange``
     - Maximum volumetric change allowed (quasi-incompressibility).
   * - ``MaxCellSurfaceAreaChange``
     - Maximum surface-area change allowed.
   * - ``MaxCellBendingAngle``
     - Maximum bending angle for red/white blood cells.
   * - ``MaxPLTBendingAngle``
     - Maximum bending angle for platelets (much larger, ~near rigid).
   * - ``MaxCellPersistenceLength``
     - Maximum persistence length, allowing up to ~300 % stretching.

.. note::

   These constants were tuned for validated, stable behaviour and are stored as
   squared values for efficiency. Changing them affects stability; see the
   commentary in ``constant_defaults.h`` before adjusting.

Writing your own model
----------------------

To create a new constitutive model, duplicate an existing model's ``.h`` and
``.cpp`` (``rbcHighOrderModel`` is the usual template), rename the class, and
re-implement ``ParticleMechanics`` with your force law. Any response that can be
expressed as a force as a function of deformation will work. Then use it via
``addCellType<YourModel>(...)``. See the :ref:`FAQ <FAQ:Frequently Asked
Questions>` entries on adding constitutive models and internal structural
mechanics, and :ref:`advanced_cases/custom_cells:Custom cells and reading cells
from STL`.
