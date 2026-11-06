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

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "packing.h"

namespace {

// The built-in cell types. A cell is represented by an enclosing ellipsoid,
// given by its three full diameters [um] along the cell's own x, y and z axes
// (HemoCell's RBC and PLT meshes have their thin axis along y). The ellipsoids
// are larger than the meshes HemoCell builds from the templates in
// tools/cell_templates: the margin stands for the IBM interaction range (cells
// interact from about one kernel width apart), the effective size of the
// cells. meshFit is the scale of the ellipsoid at which it just contains the
// mesh (measured on the generated meshes): below it, packed cells may
// intersect.
struct BuiltinType {
  const char * option;      // --<option> <n>
  const char * name;        // output file <name>.pos
  double dx, dy, dz;        // diameters of the enclosing ellipsoid [um]
  double meshFit;           // see above
  const char * description; // for --help
};
const BuiltinType builtinTypes[] = {
    // RBC_FROM_SPHERE, radius 3.91 um: 7.82 x 2.29 x 7.82 um.
    {"rbc", "RBC", 8.4, 4.4, 8.4, 0.931, "human RBCs"},
    // ELLIPSOID_FROM_SPHERE, radius 1.25 um, aspect ratio 0.435: 2.50 x 1.09 x
    // 2.50 um (the RBC's margin; it was 2.4 x 1.05 x 2.4, smaller than the
    // mesh).
    {"plt", "PLT", 2.7, 1.17, 2.7, 0.929, "human platelets"},
    // WBC_SPHERE, radius 4.0 um: 8.0 um.
    {"wbc", "WBC", 8.4, 8.4, 8.4, 0.952, "white blood cells"},
    // MESH_FROM_STL of cases/pipeflowMalaria/vRBC_uniform.stl at radius 5.4 um:
    // 2.94 x 4.83 x 10.8 um. The ellipsoid does not contain the mesh (fit
    // 1.089), so these pack from a larger start; with the RBC's margin the
    // ellipsoid would be 3.85 x 5.3 x 14.5 um.
    {"vrbc", "vRBC", 3.5, 6.0, 11.0, 1.089, "stage V gametocytes"},
    // RBC_FROM_SPHERE, radius 2.9 um: 5.80 x 1.70 x 5.80 um (the human RBC's
    // ellipsoid scaled by 2.9 / 3.91; it was 5.8 x 3.4 x 5.8).
    {"rbc_m", "RBC_m", 6.23, 3.26, 6.23, 0.931, "mouse RBCs"},
    // ELLIPSOID_FROM_SPHERE, radius 0.95 um, aspect ratio 0.57: 1.90 x 1.08 x
    // 1.90 um (it was 1.84 x 1.05 x 1.84).
    {"plt_m", "PLT_m", 2.05, 1.17, 2.05, 0.927, "mouse wild-type PLTs"},
    // ELLIPSOID_FROM_SPHERE, radius 1.71 um, aspect ratio 1: 3.42 um (it was
    // 1.71 um, the radius taken as the diameter).
    {"plt_mko", "PLT_mko", 3.68, 3.68, 3.68, 0.929, "mouse knockout PLTs"},
};
const int numBuiltinTypes = sizeof(builtinTypes) / sizeof(builtinTypes[0]);

// Effective volume of a human RBC [um^3] for --hematocrit: the 81 um^3 of the
// mesh plus the IBM interaction range. The effective size, not the mesh, is
// what matters for the collisions, and so for the viscosity.
const double rbcEffectiveVolume = 97.0;

enum OptionId {
  OPT_HEMATOCRIT = 1000,
  OPT_PLT_RATIO,
  OPT_CELL,
  OPT_NOROTATE,
  OPT_START_SCALE,
  OPT_MIN_SCALE,
  OPT_SCALE_REMOVED,
  OPT_MAXITER,
  OPT_SEED,
  OPT_HELP,
  OPT_BUILTIN = 2000 // + index into builtinTypes
};

void printHelp(std::ostream & out) {
  out << "USAGE: packCells sX sY sZ [OPTIONS]\n"
         "\n"
         "Places cells at random positions and orientations in a box of sX x sY x sZ um,\n"
         "periodic in all three directions, so that their enclosing ellipsoids do not\n"
         "overlap, and writes <type>.pos for every cell type.\n"
         "\n"
         "The ellipsoids are larger than the cells: they stand for the cells' effective\n"
         "size, which includes the IBM interaction range. They start at --start_scale\n"
         "times this nominal size and shrink until they no longer overlap. Below scale 1,\n"
         "cells may start within each other's interaction range: a relaxation that makes\n"
         "high volume fractions reachable. Below the scale at which an ellipsoid still\n"
         "contains its cell's mesh (\"mesh fit\"), the meshes may intersect, so the\n"
         "packing stops there (--min_scale).\n"
         "\n"
         "CELL TYPES (the number of cells of a type; 0 writes a file without cells):\n";
  char line[200];
  snprintf(line, sizeof(line), "  %-14s  %-21s %-12s %-21s  %s\n", "", "cells", "file", "ellipsoid [um]", "mesh fit");
  out << line;
  for (int i = 0; i < numBuiltinTypes; i++) {
    const BuiltinType & t = builtinTypes[i];
    snprintf(line, sizeof(line), "  --%-8s <n>  %-21s %-12s %5.2f x %5.2f x %5.2f  %.2f\n", t.option, t.description,
             (std::string(t.name) + ".pos").c_str(), t.dx, t.dy, t.dz, t.meshFit);
    out << line;
  }
  out << "  --cell <name> <n> <dx> <dy> <dz>\n"
         "                 a type of your own (<name>.pos): the full diameters [um] of its\n"
         "                 enclosing ellipsoid along the cell's x, y and z axes; its mesh\n"
         "                 fit is taken as 1\n"
         "  --hematocrit <Ht>, -h <Ht>\n"
         "                 RBCs for this hematocrit of the box (0 < Ht < 1), counted with the\n"
         "                 effective RBC volume of 97 um^3 (the mesh plus the IBM interaction\n"
         "                 range), and PLTs by --plt_ratio; not with --rbc or --plt\n"
         "  --plt_ratio <r>  platelets per RBC with --hematocrit (default 0.07)\n"
         "\n"
         "PACKING:\n"
         "  --start_scale <s>\n"
         "                 starting scale of the ellipsoids (default 1, or 5 % above\n"
         "                 --min_scale if that is larger); a larger start leaves more room\n"
         "                 between the cells but takes longer: the ellipsoids shrink by about\n"
         "                 0.5 % of the start every 1000 iterations (1 to 0.93: 14,000)\n"
         "  --min_scale <s>\n"
         "                 smallest scale accepted (default: the largest mesh fit of the\n"
         "                 types packed); lower it to accept intersecting meshes\n"
         "  --noRotate       keep the orientation of the cells fixed (all angles 0)\n"
         "  --maxiter <n>    stop after n iterations\n"
         "  --seed <n>       seed of the random numbers (default: a random one, printed): the\n"
         "                 same seed gives the same packing, with any number of threads\n"
         "                 (OMP_NUM_THREADS)\n"
         "  --help           print this\n"
         "\n"
         "OUTPUT (in the current directory; existing files are overwritten):\n"
         "  <type>.pos: the number of cells, then one line per cell, \"x y z rx ry rz\": the\n"
         "  centre [um] in [0, sX) x [0, sY) x [0, sZ), and the rotation [degrees] about x,\n"
         "  then y, then z, as HemoCell's loader applies it. The file names must match the\n"
         "  cell type names of the case (addCellType); rename them where they differ.\n"
         "  The progress shows the volume fraction of the ellipsoids and their scale, both\n"
         "  at the outer scale (their current size) and the inner one (at which the\n"
         "  closest overlapping pair touches; equal to the outer once none overlap).\n"
         "\n"
         "NOTES:\n"
         "  Sizes are in um, not lattice units.\n"
         "  The packing is periodic in all three directions: cell centres lie in the box,\n"
         "  but cells may cross its faces. Placing those cells in HemoCell needs\n"
         "  setPlaceCellsAcrossPeriodicBoundaries(true), and is correct only if the box\n"
         "  equals the periodic domain of the case.\n"
         "\n"
         "EXIT STATUS:\n"
         "  0 packed; 1 invalid arguments; 2 not finished within --maxiter; 3 stopped below\n"
         "  --min_scale (lower the hematocrit or enlarge the box). With 2 and 3 the files\n"
         "  are written, but the cells may overlap.\n";
}

[[noreturn]] void usageError(const std::string & message) {
  std::cerr << "packCells: " << message << "\nRun packCells --help for the options." << std::endl;
  std::exit(1);
}

// The whole argument as a finite number.
double parseReal(const char * text, const std::string & what) {
  if (text) {
    errno = 0;
    char * end = nullptr;
    const double value = std::strtod(text, &end);
    if (end != text && *end == '\0' && errno == 0 && std::isfinite(value)) { return value; }
  }
  usageError(what + ": expected a number, got '" + (text ? text : "") + "'");
}

// The whole argument as an integer in [lo, hi].
long parseInteger(const char * text, const std::string & what, long lo, long hi) {
  if (text) {
    errno = 0;
    char * end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end != text && *end == '\0' && errno == 0 && value >= lo && value <= hi) { return value; }
  }
  usageError(what + ": expected an integer from " + std::to_string(lo) + " to " + std::to_string(hi) + ", got '" +
             (text ? text : "") + "'");
}

double parsePositive(const char * text, const std::string & what) {
  const double value = parseReal(text, what);
  if (!(value > 0)) { usageError(what + ": must be positive, got '" + std::string(text) + "'"); }
  return value;
}

} // namespace

int main(int argc, char * argv[]) {
  const int maxCount = INT_MAX / 2;
  double startScale = -1.0; // default: from minScale
  double minScale = -1.0;   // default: the largest mesh fit
  double plt_ratio = 0.07;
  bool plt_ratio_set = false;
  int maxIter = INT_MAX - 100; // Nrot_step = maxIter + 1 must not overflow
  double hematocrit = 0.0;
  bool hematocrit_set = false;
  bool rbc_plt_set = false; // --rbc or --plt
  bool doRotate = true;
  long seed = long(std::random_device()() & 0x7fffffff);
  vector<CellType> cellTypes;

  if (argc == 1) {
    printHelp(std::cerr);
    return 1;
  }

  std::vector<option> longOpts = {{"hematocrit", required_argument, nullptr, OPT_HEMATOCRIT},
                                  {"plt_ratio", required_argument, nullptr, OPT_PLT_RATIO},
                                  {"cell", required_argument, nullptr, OPT_CELL},
                                  {"noRotate", no_argument, nullptr, OPT_NOROTATE},
                                  {"start_scale", required_argument, nullptr, OPT_START_SCALE},
                                  {"min_scale", required_argument, nullptr, OPT_MIN_SCALE},
                                  {"scale", required_argument, nullptr, OPT_SCALE_REMOVED},
                                  {"maxiter", required_argument, nullptr, OPT_MAXITER},
                                  {"seed", required_argument, nullptr, OPT_SEED},
                                  {"help", no_argument, nullptr, OPT_HELP}};
  for (int i = 0; i < numBuiltinTypes; i++) {
    longOpts.push_back({builtinTypes[i].option, required_argument, nullptr, OPT_BUILTIN + i});
  }
  longOpts.push_back({nullptr, 0, nullptr, 0});

  int val;
  while ((val = getopt_long_only(argc, argv, "h:s:", longOpts.data(), nullptr)) != -1) {
    switch (val) {
      case 'h':
      case OPT_HEMATOCRIT:
        hematocrit = parseReal(optarg, "--hematocrit");
        hematocrit_set = true;
        break;
      case OPT_PLT_RATIO:
        plt_ratio = parseReal(optarg, "--plt_ratio");
        plt_ratio_set = true;
        break;
      case OPT_CELL: {
        // --cell <name> <n> <dx> <dy> <dz>: optarg is the name, the other four
        // follow it.
        if (optind + 4 > argc) { usageError("--cell needs <name> <n> <dx> <dy> <dz>"); }
        const std::string name = optarg;
        const std::string what = "--cell " + name;
        const int n = int(parseInteger(argv[optind], what + ": the number of cells", 0, maxCount));
        double d[3];
        for (int k = 0; k < 3; k++) { d[k] = parsePositive(argv[optind + 1 + k], what + ": diameter"); }
        cellTypes.push_back({name, d[0], d[1], d[2], n, 0.0});
        optind += 4;
        break;
      }
      case OPT_NOROTATE:
        doRotate = false;
        break;
      case OPT_START_SCALE:
        startScale = parsePositive(optarg, "--start_scale");
        break;
      case OPT_MIN_SCALE:
        minScale = parsePositive(optarg, "--min_scale");
        break;
      case 's':
      case OPT_SCALE_REMOVED:
        usageError("--scale (-s) no longer exists: the ellipsoids start at --start_scale (default 1) times "
                   "their nominal size, and the collision-check grid is set automatically");
      case OPT_MAXITER:
        maxIter = int(parseInteger(optarg, "--maxiter", 1, INT_MAX - 100));
        break;
      case OPT_SEED:
        seed = parseInteger(optarg, "--seed", 0, LONG_MAX);
        break;
      case OPT_HELP:
        printHelp(std::cout);
        return 0;
      default:
        if (val >= OPT_BUILTIN && val < OPT_BUILTIN + numBuiltinTypes) {
          const BuiltinType & t = builtinTypes[val - OPT_BUILTIN];
          const int n = int(parseInteger(optarg, std::string("--") + t.option, 0, maxCount));
          cellTypes.push_back({t.name, t.dx, t.dy, t.dz, n, t.meshFit});
          if (std::string(t.option) == "rbc" || std::string(t.option) == "plt") { rbc_plt_set = true; }
          break;
        }
        // getopt has printed the reason
        std::cerr << "Run packCells --help for the options." << std::endl;
        return 1;
    }
  }

  // Check all arguments before packing.
  if (argc - optind != 3) {
    usageError("expected the three domain sizes sX sY sZ [um], got " + std::to_string(argc - optind) +
               " further arguments");
  }
  const double sX = parsePositive(argv[optind], "sX");
  const double sY = parsePositive(argv[optind + 1], "sY");
  const double sZ = parsePositive(argv[optind + 2], "sZ");

  if (hematocrit_set) {
    if (rbc_plt_set) { usageError("--hematocrit sets the RBCs and PLTs; it cannot be combined with --rbc or --plt"); }
    if (!(hematocrit > 0 && hematocrit < 1)) { usageError("--hematocrit: must lie between 0 and 1"); }
    if (!(plt_ratio >= 0)) { usageError("--plt_ratio: must not be negative"); }
    const BuiltinType & rbc = builtinTypes[0];
    const BuiltinType & plt = builtinTypes[1];
    const double nRBC = std::round(hematocrit * sX * sY * sZ / rbcEffectiveVolume);
    if (nRBC < 1) { usageError("the box holds no RBC at this hematocrit"); }
    if (nRBC > maxCount) { usageError("too many RBCs for this hematocrit and box"); }
    cellTypes.push_back({rbc.name, rbc.dx, rbc.dy, rbc.dz, int(nRBC), rbc.meshFit});
    const int nPLT = int(std::round(nRBC * plt_ratio));
    if (nPLT > 0) { cellTypes.push_back({plt.name, plt.dx, plt.dy, plt.dz, nPLT, plt.meshFit}); }
  } else if (plt_ratio_set) {
    usageError("--plt_ratio needs --hematocrit (give the platelets with --plt otherwise)");
  }

  std::set<std::string> names;
  long totalCells = 0;
  double largestMeshFit = 0; // of the types with cells; 1 for --cell types
  for (const CellType & cell : cellTypes) {
    if (cell.name.empty() || cell.name.find('/') != std::string::npos) {
      usageError("--cell: invalid name '" + cell.name + "'");
    }
    if (!names.insert(cell.name).second) {
      usageError("the cell type " + cell.name + " is given twice (both would be written to " + cell.name + ".pos)");
    }
    totalCells += cell.number;
    if (cell.number > 0) { largestMeshFit = std::max(largestMeshFit, cell.meshFit > 0 ? cell.meshFit : 1.0); }
  }
  if (totalCells == 0) {
    usageError("no cells to pack: give at least one cell (--rbc, --plt, --wbc, --cell, ...) or --hematocrit");
  }
  if (totalCells > maxCount) { usageError("too many cells"); }

  if (minScale < 0) { minScale = largestMeshFit; }
  if (startScale < 0) { startScale = std::max(1.0, 1.05 * minScale); }
  if (startScale < minScale) {
    usageError("--start_scale " + std::to_string(startScale) + " is below the minimum scale " +
               std::to_string(minScale) + " (--min_scale)");
  }

  cout.setf(ios::fixed, ios::floatfield);
  cout << "Parameters:" << endl;
  cout << "  Domain Size [µm]: ( " << sX << " , " << sY << " , " << sZ << " )" << endl;
  cout << "  Maximum Iterations : " << maxIter << endl;
  cout << "  Start scale        : " << startScale << endl;
  cout << "  Minimum scale      : " << minScale << endl;
  cout << "  Rotation           : " << doRotate << endl;
  cout << "  Seed               : " << seed << endl;
  if (hematocrit_set) {
    cout << "  Hematocrit    : " << hematocrit << endl;
    cout << "  PLT/RBC Ratio : " << plt_ratio << endl;
  }
  cout << "The following cell types are defined:" << endl;
  for (CellType & cell : cellTypes) {
    cout << "  " << cell.name << endl;
    cout << "    No   : " << cell.number << endl;
    cout << "    Sizes: (" << cell.dx << " , " << cell.dy << " , " << cell.dz << " )" << endl;
    cout << "    Ellipsoid volume: " << 4.0 / 3.0 * PI * cell.dx * cell.dy * cell.dz / 8.0 << " [µm3]." << endl;
    if (cell.meshFit > 0) { cout << "    Mesh fit: " << cell.meshFit << endl; }
  }

  Random::init(seed);

  Packing pack;

  pack.setRndRotation(doRotate); // This needs to be set first!

  pack.initBlood(sX, sY, sZ, maxIter, startScale, minScale, cellTypes);

  const Packing::Result result = pack.execute();

  pack.saveBloodCellPositions();

  const double finalScale = pack.getFinalScale();
  cout.precision(4);
  cout << "Final scale of the ellipsoids: " << finalScale << " (start " << startScale << ", minimum " << minScale << ")"
       << endl;
  switch (result) {
    case Packing::Packed:
      // Within 0.1 % of the nominal size counts as nominal: the scale shrinks
      // a little in every iteration, so a start at 1 ends just below it.
      if (finalScale >= 0.999) {
        cout << "The ellipsoids do not overlap at their nominal size: the cells start out of each other's "
                "interaction range."
             << endl;
      } else if (finalScale >= largestMeshFit) {
        cout << "The ellipsoids do not overlap at this scale: cells may start within each other's interaction "
                "range, but their meshes do not intersect."
             << endl;
      } else {
        cout << "WARNING: the ellipsoids do not overlap at this scale, which is below the mesh fit of the cells ("
             << largestMeshFit << "): some meshes may intersect." << endl;
      }
      return 0;
    case Packing::NotFinished:
      cout << "WARNING: not finished within --maxiter " << maxIter << ": the cells may overlap. Rerun with --seed "
           << seed << " to reproduce it." << endl;
      return 2;
    case Packing::BelowMinimumScale:
      cout << "WARNING: the ellipsoids shrank below the minimum scale " << minScale
           << " before they stopped overlapping: the cells may overlap. Lower the hematocrit or enlarge the box"
           << (minScale > 1.0 ? "." : ", or accept a smaller scale with --min_scale.") << endl;
      return 3;
  }
  return 0;
}
