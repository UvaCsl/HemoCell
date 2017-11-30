.. role:: cpplus(code)
     :language: cpp

Creating your own HemoCell case
===============================

To get started with HemoCell you can use one of the many cases that are readily
available in the ``hemocell/cases`` folder

.. toctree::
   :maxdepth: 1

   hemocell_cases

HemoCell Interface
------------------

HemoCell offers a high level C++ interface for defining and running simulations.

To get started with a HemoCell case the :cpplus:`#include "hemocell.h"`
directive
should always be present on the top of the file.
this includes the ``HemoCell`` class which you want one of per simulation.

.. doxygenclass:: HemoCell
  :project: hemocell

