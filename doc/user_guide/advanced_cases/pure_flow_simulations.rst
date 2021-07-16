Running a pure fluid flow
=========================

With HemoCell it is possible to run a pure fluid flow, without any cells present
within the domain. This is technically the same as just running a fluid
simulation with Palabos. The simplest way to remove all the cells is to set the
number of cells in the ``CELL.pos`` files to zero (first line). Or to remove the
``addCellType`` and all subsequent functions that request that specific
cell-type (such as ``setOutputs``) from the ``case.cpp`` file.

However, this still leaves a small overhead as the function for the material
model are still called (although little work is done inside them as there are no
cells). This can be circumvented by replacing

.. code-block:: c++

   hemocell.iterate()

with

.. code-block:: c++

   hemocell.lattice->collideAndStream();
   hemocell.iter++;

within the ``case.cpp`` file

.. _other_topics:csv-output:
