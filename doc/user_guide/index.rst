.. HemoCell documentation master file, created by
   sphinx-quickstart on Mon Nov 27 16:26:36 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HemoCell's User guide
====================================

*When using HemoCell please cite the corresponding HemoCell paper(s)*
:cite:`Zavodszky:2017`.

HemoCell implements an immersed boundary method with a lattice Boltzmann method
and langrangian structure points for the cells. The lattice Boltzmann method is
implemented with the help of the `Palabos`_ library.
On top of that HemoCell manages the structures (material, cells) and the
interactions between these and the fluid. The update of the forces of the
mechanical model of a cell are highly optimized as is the interaction between
the fluid and the mechanical model. 


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
  * - Max van der Kolk
    - Developer and package maintainer
    - M.vanderkolk at uva.nl  
  * - Mark Wijzenbroek
    - Package maintainer
    - 
  * - Victor Azizi  
    - Former developer
    - 
  * - Britt van Rooij 
    - Former developer 
    - 
  * - Daan van Ingen
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

HemoCell related papers
-----------------------

.. bibliography:: refs.bib
   :style: unsrt
   :all:

.. _Palabos: https://palabos.unige.ch/
