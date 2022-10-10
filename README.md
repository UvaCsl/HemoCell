     _   _   ____   __  __   _____    ___   ____   __     __   
    ( )_( ) ( ___) (  \/  ) (  _  )  / __) ( ___) (  )   (  )  
     ) _ (   )__)   )    (   )(_)(  ( (__   )__)   )(__   )(__ 
    (_) (_) (____) (_/\/\_) (_____)  \___) (____) (____) (____)   
    
          HighpErformance MicrOscopic CELlular Library


About HemoCell
==============

`HemoCell` is a parallel computing framework for simulation of dense deformable capsule suspensions, with special emphasis on blood flows and blood related vesicles (cells). The library implements validated mechanical models for Red Blood Cells (RBCs) and is capable of reproducing emergent transport characteristics of such complex cellular systems. HemoCell is capable of handling large simulation domain sizes and high shear-rate flows providing a virtual environment to evaluate a wide palette of microfluidic scenarios.

For the simulation of dense flows, HemoCell employs the Immersed Boundary Method (IBM) to couple the immersed vesicles, e.g. RBCs, platelets (PLTs), leukocytes, or other custom cell types to the fluid (e.g., plasma). The particles are tracked using a Lagrangian discrete element approach, while the flow field is implemented using the lattice Boltzmann method (LBM). The current implementation is based on the Palabos library (http://www.palabos.org). HemoCell manages all required data structures, such as materials and cell models (particles), their interactions within the flow field, load-balancing, and communication among the processors.

The library provides validated and highly optimised mechanical models for the simulation of red blood cells. Furthermore, the library is extensible and allows to implement different mechanical (cell) models and cell-binding techniques to study numerous application-specific behaviour.

Documentation
=============

Everything you need to know about how to get started with HemoCell can be found
on [hemocell.eu](https://hemocell.eu/user_guide).


License
=======

Please see the file `LICENSE` in the root directory. 


Contributions
=====================

For contributions please read our contribution guidelines on HemoCell documentation.




