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

#ifndef PACKCELLS_H
#define PACKCELLS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>

#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>

#include "geometry.h"
#include "ellipsoid.h"
#include "rnd_utils.h"


using namespace std;

class Packing {
  bool isFinished;

  int NumParts, NumSpecies;
  double Epsilon, Epsilon_scl, Epsilon_rot;
  Ellipsoid ** particles;
  Species ** species;

  double Force_step;
  // Scale of the ellipsoids relative to their nominal size (1): Douter is the
  // current scale, which starts at Douter0 and only shrinks; Dinner the scale
  // at which the closest overlapping pair touches (Douter once no pair
  // overlaps). The packing stops below MinScale.
  double Dinner, Douter0, Douter, Douter2, MinScale;
  // Volume fractions of the ellipsoids: nomPackDens at scale 1,
  // DensityNominal at Douter and Pactual at Dinner.
  double Pactual, nomPackDens, DensityNominal;

  int No_cells_x, No_cells_y, No_cells_z;
  int No_cells, Ncell_min;
  int Max_steps, NumSteps = 0, Ntau;
  double Rc_max, Relax;
  vector3 Box, Half;
  // Edge lengths of the collision-check cells: Box / No_cells per axis, at
  // least 1 (the interaction range in scaled units), so that the box itself is
  // exactly the requested size and periodic with that size.
  vector3 CellSize;
  vector3 RequestedBox; // requested domain size (scaled)

  int *Link_head, *Link_list;
  int *Ncell_bound_x, *Ncell_bound_y, *Ncell_bound_z;
  double *Pbc_x, *Pbc_y, *Pbc_z;

  bool doExit;
  int Nfig, Nsf;
  int printEveryNSteps;
  int Nrot_step;


  double Sizing = 1.0;

  bool rndRotation = true;

  // An overlapping pair (i, j), j < i, found for particle i, and the forces it
  // gives. The contacts are searched in parallel (each thread writes only the
  // list of its own particle i) and then applied in the order of i, so the
  // result does not depend on the number of threads: it equals the result of
  // one thread with the earlier code, which applied each pair as it was found.
  struct Contact {
    int j;
    double f_AB_scl;    // squared scale at which the pair touches
    vector3 f;          // translation force on j (-f on i)
    vector3 fu_i, fu_j; // rotation forces on i (subtracted) and j (zero if not rotated)
  };
  std::vector<std::vector<Contact>> contacts; // per particle i
  std::vector<double> forceMagnitude;         // |force| per particle, for Force_step

  void init();
  void iterate();
  void calc_forces();
  void calc_motion();
  // Finds the contacts of particle ipart with the particles j < ipart (the
  // collision-check cell lists must be complete).
  void collectContacts(int ipart);
  // The contact of ei and ej (ej at ei + rij); false if they do not overlap.
  bool evaluatePair(const Ellipsoid & ei, const Ellipsoid & ej, const vector3 & rij, Contact & contact) const;
  void outputHeader();
  void output();
  double zeroin(BinaryEllipsoidSystem & ell, double ax, double bx) const;
  // Collision-check cell of a position in the box (a position on the upper
  // face belongs to the last cell).
  int cellIndex(const vector3 & pos) const {
    const int x = std::min(No_cells_x - 1, std::max(0, int(pos[0] / CellSize[0])));
    const int y = std::min(No_cells_y - 1, std::max(0, int(pos[1] / CellSize[1])));
    const int z = std::min(No_cells_z - 1, std::max(0, int(pos[2] / CellSize[2])));
    return No_cells_z * (No_cells_y * x + y) + z;
  }

public:
  enum Result {
    Packed,           // no overlap at the final scale, which is at least the minimum
    NotFinished,      // the maximum number of iterations was reached first
    BelowMinimumScale // the scale fell below the minimum before the overlaps were resolved
  };

  Packing();
  ~Packing();
  // Packs the cells.
  Result execute();
  // The scale of the ellipsoids at the end (relative to their nominal size).
  double getFinalScale() const { return Douter; }

  // Sets up the packing of ctypes in a periodic box of sizeX x sizeY x sizeZ
  // [um] (the cell diameters in ctypes are in um as well). The ellipsoids start
  // at startScale times their nominal size and shrink until they no longer
  // overlap; the packing stops at maxSteps iterations or below minScale.
  void initBlood(double sizeX, double sizeY, double sizeZ, int maxSteps, double startScale, double minScale,
                 vector<CellType> & ctypes);
  // void initSuspension(vector<int> nPartsPerComponent, vector<vector3> diametersPerComponent, vector<int> domainSize, double nominalPackingDensity, int maxSteps, double sizing);
  void savePov(const char * fileName, int wbcNumber);
  void saveBloodCellPositions();
  void getOutput(vector<vector<vector3>> & positions, vector<vector<vector3>> & angles);
  void setRndRotation(bool rndRotation_) { rndRotation = rndRotation_; }
};


// The Perram-Wertheim contact function of the pair: f_AB_scl is the square of
// the scale at which the two nominal ellipsoids touch; they overlap at the
// current scale Douter if f_AB = f_AB_scl / Douter^2 < 1. Reads the particles
// only, so it may run on several threads at once.
bool Packing::evaluatePair(const Ellipsoid & ei, const Ellipsoid & ej, const vector3 & rij, Contact & contact) const {
  BinaryEllipsoidSystem e(ei, ej, rij);
  double lambda = zeroin(e, 0., 1.);
  e.get_lambda() = lambda;

  if (fabs(lambda) < 1e-10) {
#pragma omp critical(packing_messages)
    {
      cout << "e1" << endl << ei << endl;
      cout << "e2" << endl << ej << endl;
      cout << "lambda too small, returning!" << endl;
    }
    return false;
  }
  double f_AB_scl = 4 * e.f_AB(lambda);
  double f_AB = f_AB_scl / Douter2;
  vector3 n = e.count_n();
  if (f_AB >= 1) return false;
  if (f_AB <= 0) {
#pragma omp critical(packing_messages)
    {
      cout << "fab " << f_AB << endl;
      cout << ei << endl;
      cout << ej << endl;
      cout << e << endl;
      cout << "-----------------------------------------------" << endl;
    }
  }
  contact.f_AB_scl = f_AB_scl;
  //	shift
  contact.f = (1 - f_AB) * norm(n);

  //	rotation
  contact.fu_i = vector3();
  contact.fu_j = vector3();
  if (Epsilon_rot < 1e-10 || NumSteps % Nrot_step) return true;

  if (!ei.getSpecies()->getIsSphere()) contact.fu_i = (1 - f_AB) * norm(prod(e.count_rac(), n));
  if (!ej.getSpecies()->getIsSphere()) contact.fu_j = (1 - f_AB) * norm(prod(e.count_rbc(), n));
  return true;
}


Packing::Packing() {
  Dinner = 0;
  isFinished = 0;
  doExit = 0;
  Nfig = 16; //8;
  Nsf = 1;
  cout.setf(ios::fixed, ios::floatfield);
}

Packing::~Packing() {
  delete[] Link_head;
  delete[] Link_list;
  delete[] Ncell_bound_x;
  delete[] Ncell_bound_y;
  delete[] Ncell_bound_z;
  delete[] Pbc_x;
  delete[] Pbc_y;
  delete[] Pbc_z;
  for (int i = 0; i < NumParts; i++) delete particles[i];
  for (int i = 0; i < NumSpecies; i++) delete species[i];
  delete[] particles;
  delete[] species;
}

Packing::Result Packing::execute() {
  init();
  calc_forces();

  Pactual = nomPackDens * Dinner * Dinner * Dinner;
  outputHeader();

  // The scale only shrinks: once below the minimum, the packing cannot succeed.
  for (NumSteps = 0; NumSteps < Max_steps && !isFinished && Douter >= MinScale; NumSteps++) {
    iterate();
    if (NumSteps % printEveryNSteps == 0) output();
  }

  if (Douter < MinScale) {
    output();
    cout << " PACKING STOPPED(!) " << endl;
    return BelowMinimumScale;
  }
  if (!isFinished) {
    cout << " PACKING TERMINATED(!) " << endl;
    return NotFinished;
  }
  cout << " PACKING DONE " << endl; // A configuration close to equilibrium found
  return Packed;
}


// Takes dimensions in um!
void Packing::initBlood(double sizeX, double sizeY, double sizeZ, int maxSteps, double startScale, double minScale,
                        vector<CellType> & ctypes) {
  NumParts = 0;
  NumSpecies = 0;
  Epsilon = 0.1;
  Epsilon_rot = 3.0;
  Douter0 = startScale;
  MinScale = minScale;

  // Internal length unit: 1 / 1.2 of the largest diameter. The collision-check
  // cells are about that size (at least one unit), and the translation step
  // (Epsilon) and the rotation step (Epsilon_rot) are in this unit. This was
  // the --scale option (default 1.2 / the largest diameter of the most numerous
  // type), which also set the starting scale of the ellipsoids.
  double maxDiameter = 0;
  for (const CellType & ctype : ctypes) {
    maxDiameter = std::max(maxDiameter, std::max(ctype.dx, std::max(ctype.dy, ctype.dz)));
  }
  Sizing = 1.2 / maxDiameter;
  for (CellType & ctype : ctypes) {
    NumParts += ctype.number;
    NumSpecies++;
    ctype.dx *= Sizing;
    ctype.dy *= Sizing;
    ctype.dz *= Sizing;
  }

  sizeX *= Sizing;
  sizeY *= Sizing;
  sizeZ *= Sizing;

  // Collision-check cells of at least unit size (in the same quantity as the
  // cell diameters); the box keeps its exact size. (It used to be rounded up to
  // whole cells, so the packing was periodic with a larger box than requested
  // and placed cells beyond the requested size.)
  No_cells_x = std::max(1, int(floor(sizeX)));
  No_cells_y = std::max(1, int(floor(sizeY)));
  No_cells_z = std::max(1, int(floor(sizeZ)));
  RequestedBox = vector3(sizeX, sizeY, sizeZ);

  cout << endl
       << "Collision-check grid: " << No_cells_x << " x " << No_cells_y << " x " << No_cells_z << " cells of about "
       << 1. / Sizing << " um" << endl;

  Ntau = 102400;

  Max_steps = maxSteps; // if force-free configuration is not possible, still stop calculation at some point

  // Set up species
  species = new Species *[NumSpecies];
  int counter = 0;
  for (CellType & ctype : ctypes) {
    species[counter] = new Species(ctype.name, ctype.number, vector3(ctype.dx, ctype.dy, ctype.dz));
    counter++;
  }

  // Calc nominal packing density
  double domainVol = sizeX * sizeY * sizeZ;
  double cellVol = 0;
  for (CellType & ctype : ctypes) {
    cellVol += ((4. / 3. * PI * ctype.dx / 2. * ctype.dy / 2. * ctype.dz / 2.) * ctype.number);
  }

  // Volume fraction of the ellipsoids at their nominal size
  nomPackDens = cellVol / domainVol;

  cout << endl << "Volume fraction of the nominal ellipsoids: " << nomPackDens << endl;

  if (nomPackDens > 0.7) {
    cout << "*** WARNING ***" << endl;
    cout << "The enclosing ellipsoids fill " << nomPackDens << " of the box at their nominal size." << endl;
    cout << "Highest possible dense packing is around 0.77!!! You might want to reconsider..." << endl;
  }

  if (rndRotation)
    Nrot_step = 1; // Execute rotation every n-th step
  else
    Nrot_step = maxSteps + 1; // a.k.a. never

  // Output properties
  printEveryNSteps = 10;
}


void Packing::init() {
  int imult;
  int ndim_x, ndim_y, ndim_z;
  No_cells = No_cells_x * No_cells_y * No_cells_z;
  Ncell_min = (No_cells_x < No_cells_y) ? No_cells_x : No_cells_y;
  Ncell_min = (No_cells_z < Ncell_min) ? No_cells_z : Ncell_min;
  Box = RequestedBox;
  Half = 0.5 * Box;
  CellSize = vector3(Box[0] / No_cells_x, Box[1] / No_cells_y, Box[2] / No_cells_z);
  double corr = 0;
  // Rc_max: the largest diameter of all types, along any axis (it used to be
  // the x diameter only, which missed contacts of cells elongated along y or z,
  // e.g. vRBC).
  Rc_max = 0;
  for (int is = 0; is < NumSpecies; is++) { Rc_max = std::max(Rc_max, species[is]->getMaxDiameter()); }
  // The ellipsoids start at the given scale (Douter0, from initBlood()). It
  // used to be (scale * the mean of dx * dy * dz)^(1/3), so the --scale option
  // set it: 0.97 for RBCs by default, less for mixtures with small cells.
  Douter = Douter0;
  DensityNominal = nomPackDens * Douter * Douter * Douter;
  Relax = 0.5 * Douter / Ntau;
  Douter2 = Douter * Douter;
  // Rotation step from the largest diameter of all types (it was the x
  // diameter of the first type given, so it depended on the order of the
  // options; unchanged for the human cell types).
  double alpha = Rc_max;
  if (Epsilon_rot > .1) Epsilon_rot = 50. / (alpha * alpha * alpha);

  // Random positions and, unless rotation is off, uniformly random orientations
  // (spheres keep theirs). All cells used to start with the same orientation
  // and turned only when they collided, so at low volume fractions most kept
  // it (at Ht 0.1, 19 of 66 RBCs came out unrotated, lying mostly flat in the
  // xz-plane).
  particles = new Ellipsoid *[NumParts];
  int ipart = 0;
  for (int is = 0; is < NumSpecies; is++) {
    for (int ip = 0; ip < species[is]->getNum(); ip++) {
      Ellipsoid * e = new Ellipsoid(species[is], Box);
      if (rndRotation && !species[is]->getIsSphere()) e->get_q().random();
      particles[ipart++] = e;
    }
  }

  Link_head = new int[No_cells];
  Link_list = new int[NumParts];
  contacts.resize(NumParts);
  forceMagnitude.resize(NumParts);

  ndim_x = 3 * No_cells_x - 1;
  ndim_y = 3 * No_cells_y - 1;
  ndim_z = 3 * No_cells_z - 1;

  Ncell_bound_x = new int[ndim_x];
  Ncell_bound_y = new int[ndim_y];
  Ncell_bound_z = new int[ndim_z];

  Pbc_x = new double[ndim_x];
  Pbc_y = new double[ndim_y];
  Pbc_z = new double[ndim_z];

  // Cell index offsets and periodic image shifts (Pbc_*, in box lengths) for
  // a search index running over three periods: [0, N) is the image below,
  // [N, 2N) the box itself, [2N, 3N) the image above.
  imult = 1;
  corr = -Box[0];
  for (int i = 1; i < ndim_x; i++) {
    Ncell_bound_x[i] = imult * No_cells_y * No_cells_z;
    Pbc_x[i] = corr;
    if (++imult >= No_cells_x) {
      imult = 0;
      corr += Box[0];
    }
  }
  imult = 1;
  corr = -Box[1];
  for (int i = 1; i < ndim_y; i++) {
    Ncell_bound_y[i] = imult * No_cells_z;
    Pbc_y[i] = corr;
    if (++imult >= No_cells_y) {
      imult = 0;
      corr += Box[1];
    }
  }
  imult = 1;
  corr = -Box[2];
  for (int i = 1; i < ndim_z; i++) {
    Ncell_bound_z[i] = imult;
    Pbc_z[i] = corr;
    if (++imult >= No_cells_z) {
      imult = 0;
      corr += Box[2];
    }
  }
}

void Packing::iterate() {
  int iprec;
  Douter -= Relax;
  Douter2 = Douter * Douter;
  calc_motion();
  calc_forces();
  DensityNominal = nomPackDens * Douter * Douter * Douter;
  Pactual = nomPackDens * Dinner * Dinner * Dinner;

  // Stop criterion
  isFinished = (Force_step < 1e-15);

  if (isFinished) {
    output();
    return;
  }

  if (doExit) return;
  iprec = (int)(-log10(DensityNominal - Pactual));
  if (iprec >= Nsf) {
    Relax *= 0.5;
    if (++Nsf >= Nfig) doExit = true;
  }
}

// The overlap forces of all pairs. The contact search used to run in parallel
// while the cell lists were being built (each particle was inserted after its
// own search), with the contact range in shared members and the forces of both
// particles of a pair updated without synchronisation: which pairs a thread
// saw depended on timing. Now the cell lists are complete before the search,
// every thread writes only the contacts of its own particles, and the forces
// are applied afterwards in the order of the particles (as one thread did).
void Packing::calc_forces() {
  Dinner = Douter * Douter;

  for (int i = 0; i < No_cells; i++) Link_head[i] = -1;
  for (int i = 0; i < NumParts; i++) {
    particles[i]->set_force();
    particles[i]->set_forceu();
  }
  // Each cell list runs from the highest particle index down.
  for (int i = 0; i < NumParts; i++) {
    const int icell = cellIndex(particles[i]->get_pos());
    Link_list[i] = Link_head[icell];
    Link_head[icell] = i;
  }

  // Small chunks: the work per particle varies with its contacts, and a few
  // hundred particles in chunks of 16 left most of 24 threads idle (Ht 0.45 in
  // a 40 um cube: 74 instead of 263 us per iteration).
#pragma omp parallel for schedule(dynamic, 4)
  for (int ipart = 0; ipart < NumParts; ipart++) collectContacts(ipart);

  for (int i = 0; i < NumParts; i++) {
    Ellipsoid & pi = *particles[i];
    for (const Contact & c : contacts[i]) {
      Ellipsoid & pj = *particles[c.j];
      if (c.f_AB_scl < Dinner) Dinner = c.f_AB_scl;
      pi.get_f() -= c.f;
      pj.get_f() += c.f;
      pi.get_fu() -= c.fu_i;
      pj.get_fu() += c.fu_j;
    }
  }
  Dinner = sqrt(Dinner);
}

void Packing::collectContacts(int ipart) {
  std::vector<Contact> & found = contacts[ipart];
  found.clear();
  const Ellipsoid & pi = *particles[ipart];
  const vector3 pos = pi.get_pos();
  const double diameter_i = pi.getSpecies()->getMaxDiameter();
  // Contact range: half the sum of the largest diameters, plus 10 %.
  const double rcut = 0.55 * Douter * (diameter_i + Rc_max);

  Contact contact;
  // Pairs whose bounding spheres (radius: half the largest diameter at the
  // current scale) are more than 1 % apart cannot overlap: skip the contact
  // function (the expensive part) for them. The result does not change.
  auto consider = [&](int jpart, const vector3 & rij) {
    const Ellipsoid & pj = *particles[jpart];
    const double reach = 0.505 * Douter * (diameter_i + pj.getSpecies()->getMaxDiameter());
    if (rij * rij > reach * reach) return;
    if (evaluatePair(pi, pj, rij, contact)) {
      contact.j = jpart;
      found.push_back(contact);
    }
  };

  const double minCellSize = std::min(CellSize[0], std::min(CellSize[1], CellSize[2]));
  if ((int)(2.0 * rcut / minCellSize) < Ncell_min - 1) {
    // The cells within the contact range, in a search index over three periods
    // (see init()); the partners j < ipart, as each pair is found once.
    const int low_cell_x = (int)((pos[0] - rcut) / CellSize[0] + No_cells_x);
    const int low_cell_y = (int)((pos[1] - rcut) / CellSize[1] + No_cells_y);
    const int low_cell_z = (int)((pos[2] - rcut) / CellSize[2] + No_cells_z);
    const int lim_cell_x = (int)((pos[0] + rcut) / CellSize[0] + No_cells_x);
    const int lim_cell_y = (int)((pos[1] + rcut) / CellSize[1] + No_cells_y);
    const int lim_cell_z = (int)((pos[2] + rcut) / CellSize[2] + No_cells_z);
    for (int leap_x = low_cell_x; leap_x <= lim_cell_x; leap_x++) {
      const int icell_x = Ncell_bound_x[leap_x];
      for (int leap_y = low_cell_y; leap_y <= lim_cell_y; leap_y++) {
        const int icell_xy = Ncell_bound_y[leap_y] + icell_x;
        for (int leap_z = low_cell_z; leap_z <= lim_cell_z; leap_z++) {
          const int icell = Ncell_bound_z[leap_z] + icell_xy;
          int jpart = Link_head[icell];
          while (jpart >= ipart) jpart = Link_list[jpart]; // the lists run downwards; -1 ends them
          // The neighbour cell may lie in a periodic image of the box: shift
          // the neighbour into it.
          const vector3 shift(Pbc_x[leap_x], Pbc_y[leap_y], Pbc_z[leap_z]);
          for (; jpart != -1; jpart = Link_list[jpart]) { consider(jpart, particles[jpart]->get_pos() - pos + shift); }
        }
      }
    }
  } else {
    // The contact range reaches across the box: all partners, nearest image.
    for (int jpart = 0; jpart < ipart; jpart++) {
      vector3 rij = particles[jpart]->get_pos() - pos;
      rij.pbc_diff(Half, Box);
      consider(jpart, rij);
    }
  }
}

void Packing::calc_motion() {
  Force_step = 0;
  Epsilon_scl = Epsilon * Douter0;

#pragma omp parallel for schedule(static)
  for (int i = 0; i < NumParts; i++) {
    Ellipsoid * p = particles[i];
    forceMagnitude[i] = sqrt(p->get_f() * p->get_f());
    //	translation
    vector3 buff = p->get_pos() + p->get_f() * Epsilon_scl;
    buff.pbc(Box);
    p->get_pos() = buff;

    //	rotation
    if ((NumSteps + 1) % Nrot_step) continue;
    if (fabs(Epsilon_rot) < 1e-5) continue;
    if (p->getSpecies()->getIsSphere()) continue; //	sphere
    double len = Epsilon_rot * sqrt(p->get_fu() * p->get_fu());
    double cp = 1 / sqrt(1 + len * len);
    double sp = len * cp;
    double cp2 = sqrt(0.5 * (1 + cp));
    double sp2 = 0.5 * sp / cp2;
    Quaternion q(cp2, sp2 * p->get_fu());
    p->rotate(q);
  }

  // Summed in a fixed order, so that the stop criterion does not depend on the
  // number of threads (it was summed by all threads without synchronisation).
  for (int i = 0; i < NumParts; i++) Force_step += forceMagnitude[i];
  Force_step /= NumParts;
}

// Columns: the volume fraction of the ellipsoids at the inner and the outer
// scale, the inner scale (at which the closest overlapping pair touches; the
// outer scale once nothing overlaps), the outer scale (the current size of the
// ellipsoids, 1 = nominal) and the mean overlap force per cell.
void Packing::outputHeader() {
  cout << endl;
  cout << "    Steps  Vol.fraction  Vol.fraction         Scale         Scale              Force" << endl;
  cout << "                (inner)       (outer)       (inner)       (outer)       per particle" << endl;
  cout << endl;
}
void Packing::output() {
  cout.precision(10);
  cout << "\r";
  cout << setw(9) << NumSteps;
  cout << setw(14) << Pactual;
  cout << setw(14) << DensityNominal;
  cout << setw(14) << Dinner;
  cout << setw(14) << Douter;
  cout.precision(15);
  cout << setw(19) << Force_step;
  cout << std::flush;
}

void Packing::getOutput(vector<vector<vector3>> & positions, vector<vector<vector3>> & angles) {
  positions.resize(NumSpecies);
  angles.resize(NumSpecies);

  int counter = 0;
  for (int ns = 0; ns < NumSpecies; ns++) {
    int numCells = species[ns]->getNum();

    positions[ns].resize(numCells);
    angles[ns].resize(numCells);

    for (int nc = 0; nc < numCells; nc++) {
      Ellipsoid * pi = particles[counter];

      vector3 pos = pi->get_pos();
      matrix33 Q = pi->get_q().countQ();
      vector3 euler(atan2(Q(1, 2), Q(2, 2)), -asin(Q(0, 2)), atan2(Q(0, 1), Q(0, 0)));
      //euler *= 180 / PI; //Rad to Deg

      positions[ns][nc] = vector3(pos[0], pos[1], pos[2]) * (1. / Sizing);
      angles[ns][nc] = vector3(euler[0], euler[1], euler[2]); // Euler angles in RAD

      counter++;
    }
  }
}

//Legacy, needs rewriting
void Packing::savePov(const char * fileName, int wbcNumber) {
  ofstream povf(fileName);
  povf.setf(ios::fixed, ios::floatfield);
  povf.precision(6);


  // Write out ellipsoids
  for (int i = 0; i < NumParts; i++) {
    Ellipsoid * pi = particles[i];
    vector3 pos = pi->get_pos() * (1. / Sizing);
    Species * k = pi->getSpecies();
    vector3 rad = 0.5 * Dinner * k->getSize() * (1. / Sizing);
    bool sph = k->getIsSphere();

    if (i < k->getNum())
      povf << "object { RBC" << endl;
    else if (k->getNum() == wbcNumber)
      povf << "object { WBC" << endl;
    else
      povf << "object { PLT" << endl;

    povf << "scale <" << setw(12) << rad[0] << ", " << setw(12) << rad[1] << ", " << setw(12) << rad[2] << "> " << endl;
    if (!sph) { //	sphere
      matrix33 Q = pi->get_q().countQ();
      vector3 euler(atan2(Q(1, 2), Q(2, 2)), -asin(Q(0, 2)), atan2(Q(0, 1), Q(0, 0)));
      euler *= 180 / PI;
      povf << "rotate " << setw(12) << euler[0] << "*x" << endl;
      povf << "rotate " << setw(12) << euler[1] << "*y" << endl;
      povf << "rotate " << setw(12) << euler[2] << "*z" << endl;
    }
    povf << "translate <" << setw(12) << pos[0] << ", " << setw(12) << pos[1] << ", " << setw(12) << pos[2] << "> "
         << endl;

    if (NumSpecies > 1) {
      povf << "pigment { color ";
      if (i < k->getNum())
        povf << "Red";
      else if (k->getNum() == wbcNumber)
        povf << "White";
      else
        povf << "Yellow";
      povf << " }" << endl;
    }

    povf << "}" << endl;
  }
  povf.close();
}

void Packing::saveBloodCellPositions() {
  int speciesCounter = 0;

  for (int j = 0; j < NumSpecies; j++) {
    ofstream cellsFile(species[j]->name + ".pos");

    //cellsFile << No_cells_x << " " << No_cells_y << " " << No_cells_z << endl; // Dimensions

    cellsFile << species[j]->getNum() << endl; // Num. of cells of this type

    for (int i = speciesCounter; i < (speciesCounter + species[j]->getNum()); i++) {
      Ellipsoid * pi = particles[i];

      vector3 pos = pi->get_pos() * (1. / Sizing);
      matrix33 Q = pi->get_q().countQ();
      vector3 euler(atan2(Q(1, 2), Q(2, 2)), -asin(Q(0, 2)), atan2(Q(0, 1), Q(0, 0)));
      euler *= 180 / PI; //Rad to Deg

      //if(i < species[0]->getn())
      cellsFile << pos[0] << " " << pos[1] << " " << pos[2] << " " << euler[0] << " " << euler[1] << " " << euler[2]
                << endl;
      //else
      //    pltFile << pos[0] << " " << pos[1] << " " << pos[2] << " " << euler[0] << " " << euler[1] << " " << euler[2] << endl;
    }

    cellsFile.close();
    speciesCounter += species[j]->getNum();
  }
}

/*
 ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

double Packing::zeroin(BinaryEllipsoidSystem & ell, double ax, double bx) const {
  double a, b, c;
  double fa, fb, fc;
  double EPSILON = 1e-15;
  double tol = 0.0;

  a = ax;
  b = bx;
  c = a;
  fa = ell.f_AB_d(a);
  fb = ell.f_AB_d(b);
  fc = fa;

  for (;;) {
    double prev_step = b - a;
    double tol_act;
    double p;
    double q;
    double new_step;
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol_act = 2 * EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;

    if (fabs(new_step) <= tol_act || fb == (double)0) return b;

    if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
      register double t1, cb, t2;
      cb = c - b;
      if (a == c) {
        t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        q = fa / fc;
        t1 = fb / fc;
        t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > 0.)
        q = -q;
      else
        p = -p;

      if (p < (0.75 * cb * q - fabs(tol_act * q) / 2) && p < fabs(prev_step * q / 2)) new_step = p / q;
    }

    if (fabs(new_step) < tol_act) {
      if (new_step > 0.)
        new_step = tol_act;
      else
        new_step = -tol_act;
    }
    a = b;
    fa = fb;
    b += new_step;
    fb = ell.f_AB_d(b);
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
    }
  }
}


#endif /* PACKCELLS_H */
