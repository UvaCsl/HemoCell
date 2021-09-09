One stretching cell
-------------------

The example in ``examples/stretchCell`` presents a validation example for the
material models used in ``HemoCell``. A single cell is initialised in a
periodic domain, i.e. all external surfaces are subjected to periodic boundary
conditions.

The setup mimics the optical-tweezer stretching measurement.
A pair of external forces are applied on the outer points of the cell pointing
in opposite directions. These forces *stretch* the cell, where its deformation
can be validated with respect to experimental data.

After :ref:`compilation<compilation>`, the example can be run using single core as:

.. code::

   # run the simulation from the `examples/stretchCell/`
   mpirun -n 1 ./stretchCell config.xml

   # generate Paraview compatible output files
   ../../scripts/batchPostProcess.sh

.. note::
   The stretching cell examples should be run with just a single processor, so
   either ``mpirun -n 1 ./stretchCell config.xml`` as above, or directly
   ``./stretchCell config.xml``, as the helper function used to specify the stretching forces on the cell only
   supports sequential operation.

The outcome files are generated in ``tmp/``, where the flow field and particle
fields can be visualised separately by viewing the ``tmp/Fluid.*.xmf`` and
``tmp/RBC.*.xmf`` files in `Paraview`_.

.. figure:: ../_static/cases/one-stretched-cell.png
   :alt: A cell subjected to stretching forces.
   :align: center
   :figwidth: 90%

   A deformed red blood cell after being subjected to opposing forces at either
   side of the RBC, where the arrows visualise the applied forcing.

.. note::
   The `stretchCell` executable also collects the axial and transverse
   dimensions of the cell for each iteration for the used stretch force ``$f``.
   These are collection in `./stretch-$f.log` and can be used for validation
   purposes.

Configuration
=============

There is one parameter that can be modified in this example, i.e. the
applied stretching force in pico Newton:

* ``<parameters><stretchForce>``: the applied stretching force in pico Newton.


Validation
==========

The script ``validation.sh`` helps to perform basic validation of the cell
stretch model as originally presented in :cite:`Zavodszky:2017`. The script runs
through various evaluations of ``stretchCell`` with different values for the
considered stretch force ``stretchFoce`` (:ref:`cases/stretchcell:Configuration`). The resulting
log files ``stretch-*log`` contain the axial and transverse dimensions of the
RBC for the different stretch forces. From here, Figure 4
(:cite:`Zavodszky:2017`) can be recreated. If `gnuplot`_ is available on the
system, the figure is automatically generated under
``stretchCell/validation/validation.png``.

After :ref:`compilation<compilation>` ``stretchCell`` the validation can be
evaluated through:

.. code::

   ./validation.sh
   open validation/validation.png

.. note::
   In ``validation/`` a number of reference values are provided that correspond
   to the curves presented in :cite:`Zavodszky:2017` which provide a clear
   reference point when investigating cell stretch models.

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/
