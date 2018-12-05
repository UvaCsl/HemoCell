.. HemoCell documentation master file, created by
   sphinx-quickstart on Mon Nov 27 16:26:36 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HemoCell's User guide
====================================

*When using HemoCell please* **cite** *the HemoCell paper* [#HC]_.

HemoCell implements an immersed boundary method with a lattice Boltzmann method
and langrangian structure points for the cells. The lattice Boltzmann method is
implemented with the help of the palabos library (`palabos.org`_). On top of that
HemoCell manages the structures (material, cells) and the interactions between
these and the fluid. The update of the forces of the mechanical model of a cell
are highly optimized as is the interaction between the fluid and the mechanical
model. 


.. toctree::
   :maxdepth: 2

   QuickStart
   hemocell_cases
   Case
   xml_files
   Scripts
   other_topics

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
----------------

HemoCell is developed and maintained by the following persons. Any questions
about hemocell can be directed toward: info@hemocell.eu.

.. list-table::
  :header-rows: 0

  * - Victor azizi  
    - Lead Programmer
    - V.W.AziziTarksalooyeh at uva.nl
  * - Gábor Zavodszky 
    - Developer and user 
    - G.Zavodszky at uva.nl
  * - Britt van Rooij 
    - Developer and user 
    - B.J.M.vanRooij at uva.nl

*When using HemoCell please* **cite** *the HemoCell paper* [#HC]_.


..
  Indices and tables
  ==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`

.. rubric:: References
.. [#HC] `Zavodszky, G., van Rooij, B., Azizi, V., Alowayyed, S., & Hoekstra, A. (2017).  Cellular Level In-silico Modeling of Blood Rheology with An Improved Material Model for Red Blood Cells, Fronties in physiology, 8, 563. <https://doi.org/10.3389/fphys.2017.00563>`_

.. _palabos.org: http://palabos.org
