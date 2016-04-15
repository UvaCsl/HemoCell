#version  3.7;
#include "colors.inc"
global_settings{assumed_gamma 1.0}
#default{ finish{ ambient 0.1 diffuse 0.9 }}
#default{ pigment { color Red }}
 
 background { color MediumBlue }
 camera {
    location <110, 46, 58>
    look_at  <27, 27,  27>
  }

light_source { <100, 100, -100> color White}
light_source { <-100, 100, 100> color White}
light_source { <100, -100, 100> color White}
light_source { <150, 60, 60> color White}

#include "ellipsoids.pov"