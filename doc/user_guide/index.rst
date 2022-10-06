.. HemoCell documentation master file, created by
   sphinx-quickstart on Mon Nov 27 16:26:36 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HemoCell
========

HemoCell is a parallel computing framework for simulation of dense deformable
capsule suspensions, with special emphasis on blood flows and blood related
vesicles (cells). The library implements validated mechanical models for
Red Blood Cells (RBCs) and is capable of reproducing emergent transport
characteristics of such complex cellular systems :cite:`Zavodszky:2017`.
HemoCell is capable of handling large simulation domain sizes and high shear-rate
flows providing a virtual environment to evaluate a wide palette of microfluidic
scenarios :cite:`zavodszky2019red,czaja2018cell,van2020haemodynamic`.

For the simulation of dense flows, HemoCell employs the Immersed Boundary
Method (IBM) to couple the immersed vesicles, e.g. RBCs, platelets (PLTs),
leukocytes, or other custom cell types to the fluid (e.g., plasma). The particles
are tracked using a Lagrangian discrete element approach, while the flow field is
implemented using the lattice Boltzmann method (LBM). The current implementation
is based on the `Palabos`_ library. HemoCell manages all required data structures,
such as materials and cell models (particles), their interactions within the flow field,
load-balancing, and communication among the processors.

The library provides validated and highly optimised mechanical models for the
simulation of red blood cells :cite:`tarksalooyeh2019optimizing,Alowayyed:2018`.
Furthermore, the library is extensible and allows to implement different mechanical
(cell) models :cite:`czaja2020influence,Haan:2018` and cell-binding techniques
:cite:`van2019identifying` to study numerous application-specific behaviour.

.. figure:: _static/hemocell-structure.png
   :alt: Code structure of HemoCell
   :align: center
   :figwidth: 90%

The code is implemented in C/C++ with parallelism achieved through `MPI
<https://en.wikipedia.org/wiki/Message_Passing_Interface>`_, although only for
advanced use-cases the users are required to interact with the parallelism,
which is otherwise hidden from the user. The system is build with `CMake
<https://cmake.org/>`_ and runs on a variety of systems and HPC clusters (see
:ref:`getting started<QuickStart:HemoCell Getting Started>`).

Multiple examples are provided to illustrate typical use-cases of HemoCell:

* Developing mechanical models for different cell types and their interaction.
  These are typically quick running simulations on single processors (order of
  seconds to minutes) that aim to investigate/validate different formulations of
  mechanical models for the immersed cells, e.g.
  :ref:`shearing<cases/stretchcell:One stretching cell>`,
  :ref:`stretching<cases/onecellshear:One shearing cell>`, or
  :ref:`"parachuting"<cases/parachuting:A parachuting cell>` of a single cell.
  Additionally, one might want to study the interaction between colliding
  particles, e.g. :ref:`cases/cellCollision_interior_viscosity:Colliding cells
  with interior viscosity`.

  .. image:: _static/cases/rbc-plt-trajectory.png
     :width: 49%
  .. image:: _static/cases/parachuting-sideview.png
     :width: 49%

* Studying large simulation domains with large number of immersed particles.
  These simulations are typically derived from straight channel flow conditions,
  where the domain size, number of immersed particles, and flow conditions are varied.
  These simulations can vary from quick running simulations on small hardware
  (desktop/workstation) to long lasting simulations on large HPC compute
  clusters with thousands of cores. Examples of smaller pipe flow cases are
  presented in :ref:`cases/pipeflow:Pipe flow` and
  :ref:`cases/pipeflow_with_preinlet:Pipe flow with periodic inflow`.

  .. image:: _static/cases/pipeflow-initial.png
     :width: 49%
  .. image:: _static/cases/pipeflow-large.jpg
     :width: 49%

*When using HemoCell please cite the corresponding HemoCell paper(s)*
:cite:`Zavodszky:2017`.

User Guide
==========

.. toctree::
   :maxdepth: 2

   QuickStart
   hemocell_cases
   Case
   xml_files
   helper_scripts
   advanced_cases
   visualization
   utilities
   common_mistakes

.. toctree::
   :maxdepth: 1

   Downloads

.. toctree::
   :maxdepth: 1

   FAQ

.. toctree::
   :maxdepth: 2

   hemocell_doxygen

Acknowledgments
================

HemoCell is developed and maintained by the following persons, where any
questions can be directed towards: info@hemocell.eu.

.. list-table::
  :header-rows: 0

  * - Gábor Závodszky
    - Developer and Co-PI
    - G.Zavodszky at uva.nl
  * - Alfons Hoekstra
    - Co-PI
    -
  * - Ben Czaja
    - Developer
    - B.E.Czaja at uva.nl
  * - Christian Spieker
    - Developer
    - C.J.Spieker at uva.nl
  * - Mark Wijzenbroek
    - Package maintainer
    -
  * - Max van der Kolk
    - Former developer
    -
  * - Lampros Mountrakis
    - Former developer
    -
  * - Victor Azizi
    - Former developer
    -
  * - Britt van Rooij
    - Former developer
    -
  * - Saad Allowayyed
    - Former developer
    -
  * - Daan van Ingen
    - Former contributor
    -
  * - Hendrik Cornelisse
    - Former contributor
    -
  * - Mike de Haan
    - Former contributor
    -
  * - Kevin de Vries
    - Former contributor
    -
  * - Jonathan de Bouter
    - Former contributor
    -
  * - Roland Joo-Kovacs
    - Former contributor
    -

Citing HemoCell
---------------

*When using HemoCell please cite the HemoCell paper:*

.. code-block:: text

  @article{Zavodszky:2017,
    author={Závodszky, Gábor and van Rooij, Britt and Azizi, Victor and Hoekstra, Alfons},
    title={Cellular Level In-silico Modeling of Blood Rheology with An Improved Material Model for Red Blood Cells},
    journal={Frontiers in Physiology},
    volume={8},
    pages={563},
    year={2017},
    url={https://www.frontiersin.org/article/10.3389/fphys.2017.00563},
    doi={10.3389/fphys.2017.00563},
    issn={1664-042X},
  }

HemoCell related publications
=============================

.. bibliography:: refs.bib
   :style: unsrt
   :all:

.. _Palabos: https://palabos.unige.ch/

Contributions
=============

Before you contribute
---------------------

* Please make sure that your contribution falls under HemoCell license.
* If you want to resolve a bug : make sure that it still exists. 
  You can build the latest master branch and verify that the error is reproducable.
* Make sure that the bug you want to report is not already reported in our Github issues
  and that no one is working on it.
* If you have any questions about the software or if you are facing any issues using it
  feel free to open a Github issue or email at info@hemocell.eu.

Code contribution
-----------------

* Create a fork of HemoCell repository.
* Create a new branch from the develop branch for the issue you want to work on. 
  Please give a name to your branch that is relevant to the issue.
* Modify/add code to your branch.
* Before pushing make sure you have not included unrelated changes and that the project builds properly.
* After you are done, push your changes and create a pull request from your branch to the develop branch.
