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
  double Dinner, Douter0, Douter, Douter2;
  double Pactual, nomPackDens, DensityNominal;

  int No_cells_x, No_cells_y, No_cells_z;
  int No_cells, Ncell_min;
  int Max_steps, NumSteps = 0, Ntau;
  double Sphere_vol, Diam_dens;
  double Rc_max, Rcut, Rcut2, Relax;
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

  void init();
  void iterate();
  void calc_forces();
  void calc_motion();
  void force_part(int ipart_p);
  void force_all(int ipart_p);
  void outputHeader();
  void output();
  double zeroin(BinaryEllipsoidSystem & ell, double ax, double bx);
  // Collision-check cell of a position in the box (a position on the upper
  // face belongs to the last cell).
  int cellIndex(const vector3 & pos) const {
    const int x = std::min(No_cells_x - 1, std::max(0, int(pos[0] / CellSize[0])));
    const int y = std::min(No_cells_y - 1, std::max(0, int(pos[1] / CellSize[1])));
    const int z = std::min(No_cells_z - 1, std::max(0, int(pos[2] / CellSize[2])));
    return No_cells_z * (No_cells_y * x + y) + z;
  }

public:
  Packing();
  ~Packing();
  // Packs the cells; returns whether an overlap-free packing was reached
  // within the maximum number of iterations.
  bool execute();

  double calc_forces(BinaryEllipsoidSystem & e, Ellipsoid & ei, Ellipsoid & ej);

  // Initialise blood suspension with RBCs and platelets
  void initBlood(double sizeX, double sizeY, double sizeZ, int maxSteps, double sizing, vector<CellType> & ctypes);
  // void initSuspension(vector<int> nPartsPerComponent, vector<vector3> diametersPerComponent, vector<int> domainSize, double nominalPackingDensity, int maxSteps, double sizing);
  void savePov(const char * fileName, int wbcNumber);
  void saveBloodCellPositions();
  void getOutput(vector<vector<vector3>> & positions, vector<vector<vector3>> & angles);
  void setRndRotation(bool rndRotation_) { rndRotation = rndRotation_; }
};


double Packing::calc_forces(BinaryEllipsoidSystem & e, Ellipsoid & ei, Ellipsoid & ej) {
  double lambda = zeroin(e, 0., 1.);
  e.get_lambda() = lambda;

  if (fabs(lambda) < 1e-10) {
    cout << "e1" << endl << ei << endl;
    cout << "e2" << endl << ej << endl;
    cout << "lambda too small, returning!" << endl;
    return 0;
  }
  double f_AB_scl = 4 * e.f_AB(lambda);
  double f_AB = f_AB_scl / Douter2;
  vector3 n = e.count_n();
  if (f_AB >= 1) return f_AB;
  if (f_AB <= 0) {
    cout << "fab " << f_AB << endl;
    cout << ei << endl;
    cout << ej << endl;
    cout << e << endl;
    cout << "-----------------------------------------------" << endl;
  }
  //	shift
  if (f_AB_scl < Dinner) Dinner = f_AB_scl;
  vector3 f = (1 - f_AB) * norm(n);
  ei.get_f() -= f;
  ej.get_f() += f;

  //	rotation
  if (Epsilon_rot < 1e-10 || NumSteps % Nrot_step) return f_AB;

  if (!ei.getSpecies()->getIsSphere()) ei.get_fu() -= (1 - f_AB) * norm(prod(e.count_rac(), n));
  if (!ej.getSpecies()->getIsSphere()) ej.get_fu() += (1 - f_AB) * norm(prod(e.count_rbc(), n));
  return f_AB;
}


Packing::Packing() {
  Dinner = 0;
  isFinished = 0;
  doExit = 0;
  Nfig = 16; //8;
  Nsf = 1;
  Sphere_vol = PI / 6.;
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

bool Packing::execute() {
  init();
  calc_forces();

  Pactual = Dinner * Dinner * Dinner / Diam_dens;
  outputHeader();

  for (NumSteps = 0; NumSteps < Max_steps && !isFinished; NumSteps++) {
    iterate();
    if (NumSteps % printEveryNSteps == 0) output();
  }

  if (isFinished) {
    cout << " PACKING DONE " << endl; // A configuration close to equilibrium found
  } else {
    // No equilibrium within the maximum number of iterations
    cout << " PACKING TERMINATED(!) " << endl;
    cout << "WARNING: no overlap-free packing within " << Max_steps
         << " iterations: the cells may overlap (inner diameter < outer diameter)." << endl;
  }
  return isFinished;
}


// Takes dimensions in um!
void Packing::initBlood(double sizeX, double sizeY, double sizeZ, int maxSteps, double sizing,
                        vector<CellType> & ctypes) {
  NumParts = 0;
  NumSpecies = 0;
  Epsilon = 0.1;
  Epsilon_rot = 3.0;

  Sizing = sizing;
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

  cout << endl << "Sizing parameter: " << Sizing << endl;
  cout << "Resolution of collision check: " << No_cells_x << " x " << No_cells_y << " x " << No_cells_z << endl;

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

  // Get nominal volume ratio
  nomPackDens = cellVol / domainVol;

  cout << endl << "Nominal requested volume fraction: " << nomPackDens << endl;

  if (nomPackDens > 0.7) {
    cout << "*** WARNING ***" << endl;
    cout << "For the requested cells the enclosing ellipsoids yield a nominal volume ratio of: " << nomPackDens << endl;
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
  Diam_dens = Box[0] * Box[1] * Box[2] / Sphere_vol;
  double corr = 0;
  // Rc_max: the largest diameter of all types, along any axis (it used to be
  // the x diameter only, which missed contacts of cells elongated along y or z,
  // e.g. vRBC).
  Rc_max = 0;
  for (int is = 0; is < NumSpecies; is++) {
    corr += species[is]->getNum();
    Rc_max = std::max(Rc_max, species[is]->getMaxDiameter());
  }
  Diam_dens /= corr;
  DensityNominal = nomPackDens;
  Douter0 = pow(Diam_dens * nomPackDens, 1. / 3.);
  Douter = Douter0;
  Relax = 0.5 * Douter / Ntau;
  Douter2 = Douter * Douter;
  // Rotation step from the largest diameter of all types (it was the x
  // diameter of the first type given, so it depended on the order of the
  // options; unchanged for the human cell types).
  double alpha = Rc_max;
  if (Epsilon_rot > .1) Epsilon_rot = 50. / (alpha * alpha * alpha);

  particles = new Ellipsoid *[NumParts];
  int ipart = 0;
  for (int is = 0; is < NumSpecies; is++)
    for (int ip = 0; ip < species[is]->getNum(); ip++) particles[ipart++] = new Ellipsoid(species[is], Box);

  Link_head = new int[No_cells];
  Link_list = new int[NumParts];

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
  DensityNominal = Douter * Douter * Douter / Diam_dens;
  Pactual = Dinner * Dinner * Dinner / Diam_dens;

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

void Packing::calc_forces() {
  Dinner = Douter * Douter;
  int i;

#pragma omp parallel for private(i)
  for (i = 0; i < No_cells; i++) Link_head[i] = -1;

#pragma omp parallel for private(i)
  for (i = 0; i < NumParts; i++) {
    particles[i]->set_force();
    particles[i]->set_forceu();
  }

  int ipart;

#pragma omp parallel for private(ipart) schedule(static)
  for (ipart = 0; ipart < NumParts; ipart++) {
    Species * k = particles[ipart]->getSpecies();
    // Contact range: half the sum of the largest diameters, plus 10 %.
    Rcut = 0.55 * Douter * (k->getMaxDiameter() + Rc_max);
    Rcut2 = Rcut * Rcut;

    const double minCellSize = std::min(CellSize[0], std::min(CellSize[1], CellSize[2]));
    if ((int)(2.0 * Rcut / minCellSize) < Ncell_min - 1)
      force_part(ipart);
    else
      force_all(ipart);
  }
  Dinner = sqrt(Dinner);
}

void Packing::calc_motion() {
  vector3 buff;
  Force_step = 0;
  Epsilon_scl = Epsilon * Douter0;
  int i;

#pragma omp parallel for private(i, buff)
  for (i = 0; i < NumParts; i++) {
    Ellipsoid * p = particles[i];
    Force_step += sqrt(p->get_f() * p->get_f());
    //	translation
    buff = p->get_pos() + p->get_f() * Epsilon_scl;
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

  Force_step /= NumParts;
}

void Packing::force_part(int ipart_p) {
  int jpart, ifirst;
  int icell_x, icell_xy, icell;
  int leap_x, leap_y, leap_z;
  Ellipsoid *pi = particles[ipart_p], *pj;
  vector3 rij;
  vector3 pos = pi->get_pos();

  int low_cell_x = (int)((pos[0] - Rcut) / CellSize[0] + No_cells_x);
  int low_cell_y = (int)((pos[1] - Rcut) / CellSize[1] + No_cells_y);
  int low_cell_z = (int)((pos[2] - Rcut) / CellSize[2] + No_cells_z);
  int lim_cell_x = (int)((pos[0] + Rcut) / CellSize[0] + No_cells_x);
  int lim_cell_y = (int)((pos[1] + Rcut) / CellSize[1] + No_cells_y);
  int lim_cell_z = (int)((pos[2] + Rcut) / CellSize[2] + No_cells_z);

#pragma omp parallel
  for (leap_x = low_cell_x; leap_x <= lim_cell_x; leap_x++) {
    icell_x = Ncell_bound_x[leap_x];
    for (leap_y = low_cell_y; leap_y <= lim_cell_y; leap_y++) {
      icell_xy = Ncell_bound_y[leap_y] + icell_x;
#pragma omp for private(leap_z)
      for (leap_z = low_cell_z; leap_z <= lim_cell_z; leap_z++) {
        icell = Ncell_bound_z[leap_z] + icell_xy;
        jpart = Link_head[icell];
        while (jpart != -1) {
          pj = particles[jpart];
          // The neighbour cell may lie in a periodic image of the box: shift
          // the neighbour into it (the shift was computed but not applied, so
          // cells were not repelled across the box faces and overlapped there).
          rij = pj->get_pos() - pi->get_pos() + vector3(Pbc_x[leap_x], Pbc_y[leap_y], Pbc_z[leap_z]);
          BinaryEllipsoidSystem eij(*pi, *pj, rij);
          calc_forces(eij, *pi, *pj);
          jpart = Link_list[jpart];
        }
      }
    }
  }

  icell = cellIndex(pos);
  ifirst = Link_head[icell];
  Link_list[ipart_p] = ifirst;
  Link_head[icell] = ipart_p;
}

void Packing::force_all(int ipart_p) {
  Ellipsoid *pi = particles[ipart_p], *pj;
  vector3 pos = pi->get_pos();
  vector3 rij;
  int jpart;

#pragma omp parallel for private(jpart)
  for (jpart = 0; jpart < ipart_p; jpart++) {
    pj = particles[jpart];
    rij = pj->get_pos() - pi->get_pos();
    rij.pbc_diff(Half, Box);
    BinaryEllipsoidSystem eij(*pi, *pj, rij);
    calc_forces(eij, *pi, *pj);
  }
  int icell = cellIndex(pos);
  int ifirst = Link_head[icell];
  Link_list[ipart_p] = ifirst;
  Link_head[icell] = ipart_p;
}

void Packing::outputHeader() {
  cout << endl;
  cout << "     Steps";
  cout << "     Actual";
  cout << "       Nominal";
  cout << "        Inner";
  cout << "         Outer";
  cout << "             Force";
  cout << endl;
  cout << "         ";
  cout << "     density";
  cout << "       density";
  cout << "       diameter";
  cout << "      diameter";
  cout << "       per particle";
  cout << endl;
  cout << endl;
  return;
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

double Packing::zeroin(BinaryEllipsoidSystem & ell, double ax, double bx) {
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
