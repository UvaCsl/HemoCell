/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CONSTANT_DEFAULTS_H
#define CONSTANT_DEFAULTS_H
// ============== Compile time options - Set these

/*
 * Which version are we at
 */
#define VERSION_MAJOR 1
#define VERSION_MINOR 4

/*
Choose kernel. 
Phi1 [1], Phi2 [2] - Default, Phi3 [3], Phi4 [4],  Phi4c [5]
*/
#ifndef HEMOCELL_KERNEL
#define HEMOCELL_KERNEL 2
#endif

/*
Choose material model.
1 - Dao/Suresh model (+Fedosov 2010 improvements)
2 - HO model (Default - this is the validated model)
*/
#ifndef HEMOCELL_MATERIAL_MODEL
#define HEMOCELL_MATERIAL_MODEL 2
#endif

/*
Choose to use precalculated curvature angles for equilibrium cell shapes.
*/
#ifndef PRECALCULATED_ANGLES
#define PRECALCULATED_ANGLES
#endif

/*
Choose material integration method.
Euler [1], Adams-Bashforth [2]
*/
#ifndef HEMOCELL_MATERIAL_INTEGRATION
#define HEMOCELL_MATERIAL_INTEGRATION 1
#endif

/*
Select BGK dynamics with Guo's forcing method.
*/
#ifndef DESCRIPTOR
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#define DESCRIPTOR plb::descriptors::ForcedD3Q19Descriptor
#endif

#ifndef CEPAC_DESCRIPTOR
#include "latticeBoltzmann/advectionDiffusionLattices.h"
#define CEPAC_DESCRIPTOR plb::descriptors::AdvectionDiffusionD3Q19Descriptor
#endif

/*
Force unconditional stability on the material model. Note: it will not make the model magically correct, only stable!
FORCE_LIMIT sets the allowed maximal force coming from the constitutive model (in LBM units).
*/
// in [pN / surface particle].
#ifndef FORCE_LIMIT
#define FORCE_LIMIT 50.0 
#endif
/*
Choose bending force implementation. This is a numerical modelling choice.
[1] - Local only bending implementation acting on the two opposing vertices of the surfaces with the bending angle.
[2] - Distributed bending using all four vertices with wieghting.
[3] - Distributed bending using all four vertices.
Note: 	[1] is advised for cases where structural rigidity is needed. 
		[2] is an intermediate modell between [1] and [3] using four vertices, however, weighting them. So far seem to be the best option.
		[3] is useful for having increased numerical stability (req.: dx<= 0.5 um, otherwise oscillates). 
*/
#ifndef HEMOCELL_MEMBRANE_BENDING
#define HEMOCELL_MEMBRANE_BENDING 1
#endif

/*
 * Particle Field Type:
 * Correct options are:
 * DenseParticleField3D [Palabos]
 * LightParticleField3D [Palabos]
 * HemoParticleField3D  [HemoCell]
 */
#ifndef HEMOCELL_PARTICLE_FIELD
#define HEMOCELL_PARTICLE_FIELD class HemoCellParticleField
#endif

/*
 * Used for defining the shape of the cells through constructMeshElement function
 */
#define RBC_FROM_SPHERE 1
#define ELLIPSOID_FROM_SPHERE 6
#define STRING_FROM_VERTEXES 7
#define WBC_SPHERE 0
#define MESH_FROM_STL 2

/*
 * Defines for desired output per celltype, can save space/time etc. etc. also works for specific fluid ones (vel, force etc)
 */
#define OUTPUT_POSITION 1
#define OUTPUT_FORCE 2
#define OUTPUT_FORCE_VOLUME 21
#define OUTPUT_FORCE_BENDING 22
#define OUTPUT_FORCE_AREA 23
#define OUTPUT_FORCE_LINK 24
#define OUTPUT_FORCE_VISC 25
#define OUTPUT_FORCE_INNER_LINK 26
#define OUTPUT_FORCE_REPULSION 27
#define OUTPUT_TRIANGLES 3
#define OUTPUT_VELOCITY 4
#define OUTPUT_DENSITY 5
#define OUTPUT_VERTEX_ID 7
#define OUTPUT_CELL_ID 8
#define OUTPUT_CELL_DENSITY 9
#define OUTPUT_SHEAR_STRESS 10
#define OUTPUT_INNER_LINKS 11


//==================== Not really an option but a nice shortcut
#define param Parameters

#ifndef T
typedef double T;
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

//Define signed and unsigned ints
#ifndef plint
typedef long int plint;
#endif
#ifndef pluint
typedef long unsigned int pluint;
#endif

#endif
