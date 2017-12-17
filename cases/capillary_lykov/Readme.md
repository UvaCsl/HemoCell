#Microcontraction Lykov

This simulation based on this article:

>[Lykov et al. 2017 Probing eukaryotic cell mechanics via mesoscopic simulations.](http://doi.org/10.1371/journal.pcbi.1005726)
>The geometry parameters are the same as in the [Fig3](http://doi.org/10.1371/journal.pcbi.1005726.g003)

 The gaps between the obstacles are 10, 12 and 15 um, scaled with the size of the wbc diameter (now: wbc diameter 8um, gap is 6um).
 
 The gap should be larger than the wbc rigid core diameter!!

 The length is 25 um, the larger height is 10 um and the smaller one is 3.4 um. The spacing between the rows of obstacles is 60 um
 Pressure drop: 0.67Pa/μm -> in cpp dpdz = 6.7e5 

       ^    |¯   -   _                                         |¯   -   _
       |    |             ¯    _                               |
       |    |                        |   ^                     |
     10 um  |                        |   | 3.4 um              |
       |    |                        |   |                     |
       |    |             _    ¯         v                     |
       v    |_   -  ¯                                          |_   -  ¯
            ^
            |
          10~15 um                    <---------60um---------->
            |
            v
            |¯   -   _                                         |¯   -   _
            |             ¯    _                               |
            |                        |                         |
            |                        |                         |
            |                        |                         |
            |             _    ¯                               |
            |_   -  ¯                                          |_   -  ¯
            <------------25um------->