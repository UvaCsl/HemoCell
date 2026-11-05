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

#include <iomanip>

#ifndef CONSTANT_DEFAULTS_H
#define CONSTANT_DEFAULTS_H
// ============== Compile time options - Set these

/*
 * Which version are we at
 */
#define VERSION_MAJOR 2
#define VERSION_MINOR "8"

// Use for solidifying, performance impact so only enable when using
#ifndef SOLIDIFY_MECHANICS
#define SOLIDIFY_MECHANICS
#endif

// Use for Interior Viscosity mechanics, performance impact so only enable when using
#ifndef INTERIOR_VISCOSITY
//#define INTERIOR_VISCOSITY
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
FORCE_LIMIT sets the allowed maximal force per surface particle, in [pN]. It is
converted to lattice units in Parameters (param::f_limit) and applied to the
total membrane force when spreading it to the fluid; a quarter of it caps the
membrane-viscosity force per edge. Define FORCE_LIMIT differently (or remove
this default) at compile time to change or disable the cap.
*/
#ifndef FORCE_LIMIT
#define FORCE_LIMIT 50.0
#endif

namespace hemo {
/*
 * Shape of the reference mesh of a cell type (HemoCell::addCellType(),
 * constructMeshElement()). Constants in namespace hemo (formerly macros).
 */
constexpr int WBC_SPHERE = 0;
constexpr int RBC_FROM_SPHERE = 1;
constexpr int MESH_FROM_STL = 2;
constexpr int ELLIPSOID_FROM_SPHERE = 6;
constexpr int STRING_FROM_VERTEXES = 7;

/*
 * Output fields per cell type (HemoCell::setOutputs()) and of the fluid
 * (HemoCell::setFluidOutputs(), setCEPACOutputs()); see
 * doc/user_guide/reference/outputs.rst. Constants in namespace hemo
 * (formerly macros).
 */
constexpr int OUTPUT_POSITION = 1;
constexpr int OUTPUT_FORCE = 2;
constexpr int OUTPUT_FORCE_VOLUME = 21;
constexpr int OUTPUT_FORCE_BENDING = 22;
constexpr int OUTPUT_FORCE_AREA = 23;
constexpr int OUTPUT_FORCE_LINK = 24;
constexpr int OUTPUT_FORCE_VISC = 25;
constexpr int OUTPUT_FORCE_INNER_LINK = 26;
constexpr int OUTPUT_FORCE_REPULSION = 27;
constexpr int OUTPUT_TRIANGLES = 3;
constexpr int OUTPUT_VELOCITY = 4;
constexpr int OUTPUT_DENSITY = 5;
constexpr int OUTPUT_VERTEX_ID = 7;
constexpr int OUTPUT_CELL_ID = 8;
constexpr int OUTPUT_CELL_DENSITY = 9;
constexpr int OUTPUT_SHEAR_STRESS = 10;
constexpr int OUTPUT_INNER_LINKS = 11;
constexpr int OUTPUT_OMEGA = 12;
constexpr int OUTPUT_BOUNDARY = 13;
constexpr int OUTPUT_BINDING_SITES = 14;
constexpr int OUTPUT_INTERIOR_POINTS = 15;
constexpr int OUTPUT_SHEAR_RATE = 16;
constexpr int OUTPUT_STRAIN_RATE = 17;
constexpr int OUTPUT_RES_TIME = 18;
} // namespace hemo

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

/** @defgroup cellgroup Cell mechanics constants
 * A set of constants defining the cell stiffness for various deformation modes.
 * The presentation and derivation of these constants are presented in the
 * validation paper: https://doi.org/10.3389/fphys.2017.00563.
 *
 * Note that these parameters are not unique, i.e. multiple combinations of the
 * constants can result in validated cell behaviour.
 *
 * The current set of parameters is slightly different as presented in the
 * original publication, specifically \p MaxCellBendingAngle has been reduced to
 * achieve improved stability when RBCs are subjected to high shear flows.
 * Although this results in slightly stiffer cell mechanics, the obtained
 * force-displacement curves are still within the range of the validation data.
 * For reference, please refer to `examples/stretchCell`.
 *
 * All parameters are expressed as their squared value to prevent computation of
 * the square at each evaluation.
 *
 * @{
 */

/** @brief The maximum volumetric change of an cell to maintain
 * quasi-incompressibility.
 */
constexpr T MaxCellVolumetricChange = 0.01;
/** @brief The maximum change of the cells surface area is limited to 30%. */
constexpr T MaxCellSurfaceAreaChange = 0.09;
/** @brief The maximum bending angle for red/white blood cells.
 * The bending angle of the RBC/WBC are changed to allow for increased stability
 * when encountering high shear flows during the simulations.
 */
constexpr T MaxCellBendingAngle = 0.0555;
/** @brief The maximum bending angle for platelets.
 * Approximately (pi/2)^2 as derived from the smallest representable curvature
 * of the triangulated surface mesh.
 */
constexpr T MaxPLTBendingAngle = 2.467;
/** @brief The maximum link force derived from persistence lengths.
 * This maximum cell persistence length allows for at most 300% stretching.
 */
constexpr T MaxCellPersistenceLength = 9.0;

/** @} */ // end of the cell group constants

namespace hemo {
///Used to circumvent buffer initialization of characters
struct NoInitChar {
  char value;
  NoInitChar() {
    // do nothing
    static_assert(sizeof *this == sizeof value, "invalid size");
    static_assert(__alignof * this == __alignof value, "invalid alignment");
  }
  ~NoInitChar(){};
};
}

#endif
