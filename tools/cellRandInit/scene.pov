#version  3.7;
#include "colors.inc"
global_settings{assumed_gamma 1.0}
#default{ finish{ ambient 0.1 diffuse 0.65 phong 0.05 phong_size 50 roughness 0.1 }}
#default{ pigment { color rgb< 187, 96, 96>} }// { color Red }}

 //background { color MediumBlue }
 background { color White }

 camera {
    location <1600, 800, 850>
    look_at  <512.000000, 226.000000, 256.000000>
  }

light_source { <1500, 0, 0> color White}
light_source { <-100, 0, 0> color White}
light_source { <0, 1500, 0> color White}
light_source { <0, -100, 0> color White}
light_source { <0, 0, 1500> color White}
light_source { <0, 0, -100> color White}

// Use RBC.pov for biconcave shape
#include "RBC.pov"
//#declare RBC = sphere {<0, 0, 0>, 1}
#declare PLT = sphere {<0, 0, 0>, 1}
#declare WBC = sphere {<0, 0, 0>, 1}

#include "cells.pov"