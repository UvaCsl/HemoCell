/* 
 * File:   packCells.h
 * Author: vikko
 *
 * Created on 11 December 2017, 10:06
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

// Oversize ellipsoids by 10% (only for init suspension)
const double OVERSIZE = 1.1;	// This helps to avoid too close membranes -> problematic overlaps for IBM

using namespace std;

class Packing {
  bool isFinished;

  int NumParts, NumSpecies;
  double Epsilon, Epsilon_scl, Epsilon_rot;
  Ellipsoid **particles;
  Species **species;

  double Force_step;
  double Dinner, Douter0, Douter, Douter2;
  double Pactual, nomPackDens, DensityNominal;

  int No_cells_x, No_cells_y, No_cells_z;
  int No_cells, Ncell_min;
  int Max_steps, NumSteps=0, Ntau;
  double Sphere_vol, Diam_dens;
  double Rc_max, Rcut, Rcut2, Relax;
  vector3 Box, Half;

  int *Link_head, *Link_list;
  int *Ncell_bound_x, *Ncell_bound_y, *Ncell_bound_z;
  double *Pbc_x, *Pbc_y, *Pbc_z;

  bool doExit;
  int Nfig, Nsf;
  int printEveryNSteps;
  int Nrot_step;

  // Standard diameters for bounDinnerg ellipsoids of RBCs and platelets in um [= 1e-6 m]
  //double rbcA = 8.0, rbcB = 2.5, rbcC = 8.0;
  //double plateletA = 2.5 , plateletB = 0.75, plateletC = 2.5;

  //double rbcA = 9.0, rbcB = 4.4, rbcC = 9.0; // Inreased to cover biconcave shape of RBCs

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
  double zeroin(BinaryEllipsoidSystem& ell, double ax, double bx);
public:
  Packing ();
  ~Packing();
  void execute();

  double calc_forces(BinaryEllipsoidSystem& e, Ellipsoid& ei, Ellipsoid& ej);

  // Initialise blood suspension with RBCs and platelets
  void initBlood(float sizeX, float sizeY, float sizeZ, int maxSteps, double sizing, vector<CellType> & ctypes);
  void initSuspension(vector<int> nPartsPerComponent, vector<vector3> diametersPerComponent, vector<int> domainSize, double nominalPackingDensity, int maxSteps, double sizing);
  void savePov(const char * fileName, int wbcNumber);
  void saveBloodCellPositions();
  void getOutput(vector<vector<vector3> > &positions, vector<vector<vector3> > &angles);
  void setRndRotation(bool rndRotation_) {rndRotation = rndRotation_;}
};



double Packing::calc_forces(BinaryEllipsoidSystem& e, Ellipsoid& ei, Ellipsoid& ej) {

  double lambda = zeroin(e, 0., 1.);
  e.get_lambda() = lambda;

  if (fabs(lambda) < 1e-10) {
    cout << "e1" << endl << ei << endl;
    cout << "e2" << endl << ej << endl;
    cout << "lambda too small, returning!"<<endl;
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

  if (!ei.getSpecies()->getIsSphere())
          ei.get_fu() -= (1 - f_AB) * norm(prod(e.count_rac(), n));
  if (!ej.getSpecies()->getIsSphere())
          ej.get_fu() += (1 - f_AB) * norm(prod(e.count_rbc(), n));
  return f_AB;
}


Packing::Packing () {
	Dinner = 0;
	isFinished = 0;
	doExit = 0;
	Nfig = 16;//8;
	Nsf = 1;
	Sphere_vol = PI / 6.;
	Random::init();
	cout.setf (ios::fixed, ios::floatfield);

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

void Packing::execute() {
	//cout << "execute begin" << endl;
	init();
	//cout << "init end" << endl;
	calc_forces();
	//cout << "calc_forces end" << endl;

	Pactual = Dinner*Dinner*Dinner / Diam_dens;
	outputHeader();
	
	for (NumSteps=0; NumSteps < Max_steps && !isFinished; NumSteps++) {
		iterate();
		if (NumSteps % printEveryNSteps == 0) output();
	}
	if(NumSteps == Max_steps-1)
        cout << "WARNING: Force-free packing was NOT achieved! Potential overlaps might still exist." << endl;
	
	if (isFinished) cout << " PACKING DONE " << endl;		// A configuration close to equilibrium found
	else cout << " PACKING TERMINATED(!) " << endl;	// No equilibrium yet, but used up maximal iterations, or there is no better configuration
}


// Takes dimensions in um!
void Packing::initBlood(float sizeX, float sizeY, float sizeZ, int maxSteps, double sizing, vector<CellType> & ctypes)
{
  NumParts = 0;
  NumSpecies = 0;
  Epsilon = 0.1;
  Epsilon_rot = 3.0;

  Sizing = sizing;
  for (CellType & ctype:ctypes) {
    NumParts += ctype.number;
    NumSpecies++;
    ctype.dx *= Sizing;
    ctype.dy *= Sizing;
    ctype.dz *= Sizing;
  }

  sizeX *= Sizing;
  sizeY *= Sizing;
  sizeZ *= Sizing;

  No_cells_x = ceil(sizeX); // in the same quantity as cell diameters (umeter)
  No_cells_y = ceil(sizeY);
  No_cells_z = ceil(sizeZ);

  Ntau = 102400;

  Max_steps = maxSteps; // if force-free configuration is not possible, still stop calculation at some point

  // Set up species
  species = new Species*[NumSpecies];
  int counter = 0;
  for (CellType & ctype : ctypes) {
    species[counter] = new Species(ctype.name, ctype.number,vector3(ctype.dx,ctype.dy,ctype.dz));
    counter++;
  }

  // Calc nominal packing density
  double domainVol = sizeX * sizeY * sizeZ;
  double cellVol = 0;
  for (CellType & ctype : ctypes) {
    cellVol += ((4./3. * PI * ctype.dx/2. * ctype.dy/2. * ctype.dz/2.) * ctype.number);
  }

  // Get nominal volume ratio
  nomPackDens = cellVol / domainVol;

  if (nomPackDens > 0.7){
      cout << "*** WARNING ***" << endl;
      cout << "For the requested hematocrit the enclosing ellipsoids yield a nominal volume ratio of: " << nomPackDens << endl;
      cout << "Highest achieved packing (lots of iterations) is around 0.77!!! You might want to reconsider..." << endl;
  }

  if(rndRotation)
      Nrot_step = 1;	// Execute rotation every n-th step
  else
      Nrot_step = maxSteps+1; // a.k.a. never

  // Output properties
  printEveryNSteps = 10;
}

void Packing::initSuspension(vector<int> nPartsPerComponent, vector<vector3> diametersPerComponent, vector<int> domainSize, double nominalPackingDensity, int maxSteps = 25000, double sizing = 1.0)
{
    NumParts = 0;
    for(unsigned int i = 0; i < nPartsPerComponent.size(); i++)
        NumParts += nPartsPerComponent[i];

    NumSpecies = nPartsPerComponent.size();
    Epsilon = 0.1;//0.1;
    Epsilon_rot = 3.0;

    Sizing = sizing;

    No_cells_x = ceil(domainSize[0] * Sizing); // in the same quantity as cell diameters
    No_cells_y = ceil(domainSize[1] * Sizing);
    No_cells_z = ceil(domainSize[2] * Sizing);
    
    Ntau = 102400;

    if(rndRotation)
    	Nrot_step = 1;	// Execute rotation every n-th step
    else
    	Nrot_step = maxSteps+1; // a.k.a. never

    Max_steps = maxSteps; // if force-free configuration is not possible, still stop calculation at some point

    // Set up species
    species = new Species*[NumSpecies];
    for(int i = 0; i < NumSpecies; i++)
        species[i] = new Species(nPartsPerComponent[i],
                                 diametersPerComponent[i] * OVERSIZE * Sizing); // Inflate by 10%

    // Get nominal volume ratio
    nomPackDens = nominalPackingDensity * OVERSIZE;

    // Output properties
    printEveryNSteps = 250;

    cout << "Number of maximal iterations: " << Max_steps << endl;
    cout << "Rotation happens every nth step: " << Nrot_step << endl;
    cout << "Number of bins for neighbour cutoff: " << No_cells_x << " x " << No_cells_y << " x " << No_cells_z << endl;
    cout << "Number of cells to pack:  " << endl;
    for(int i = 0; i < NumSpecies; i++)
        cout << "    Type: " << i << " - number: " << nPartsPerComponent[i] << " - diameters: " << 	diametersPerComponent[i][0] << "x" << diametersPerComponent[i][1] << "x" << diametersPerComponent[i][2] << endl;    

}

void Packing::init() {
	int imult;
	int ndim_x, ndim_y, ndim_z;
	No_cells = No_cells_x * No_cells_y * No_cells_z;
	Ncell_min = (No_cells_x < No_cells_y) ?  No_cells_x : No_cells_y;
	Ncell_min = (No_cells_z < Ncell_min) ? No_cells_z : Ncell_min;
	Box = vector3(No_cells_x, No_cells_y, No_cells_z);
	Half = 0.5 * Box;
	Diam_dens = No_cells / Sphere_vol;
	double corr = 0;
	vector3 r;
	Rc_max = 0;
	for (int is = 0; is < NumSpecies; is++) {
		r = species[is]->getRot();
		corr += species[is]->getNum();// * r[0]*r[1]*r[2];
		if (r[0] > Rc_max) Rc_max = r[0];
	}
	Diam_dens /= corr;
	DensityNominal = nomPackDens;
	Douter0 = pow (Diam_dens * nomPackDens , 1./3.);
	Douter = Douter0;
	Relax = 0.5 * Douter / Ntau;
	Douter2 = Douter * Douter;
	double alpha = species[0]->getRot()[0];
	if(Epsilon_rot > .1) Epsilon_rot = 50. / (alpha*alpha*alpha);

	particles = new Ellipsoid*[NumParts];
	int ipart = 0;
	for (int is = 0; is < NumSpecies; is++)
		for (int ip = 0; ip < species[is]->getNum(); ip++)
			particles[ipart++] = new Ellipsoid (species[is], Box);

	Link_head = new int [No_cells];
	Link_list = new int [NumParts];

	ndim_x = 3 * No_cells_x - 1;
	ndim_y = 3 * No_cells_y - 1;
	ndim_z = 3 * No_cells_z - 1;

	Ncell_bound_x = new int [ndim_x];
	Ncell_bound_y = new int [ndim_y];
	Ncell_bound_z = new int [ndim_z];

	Pbc_x = new double [ndim_x];
	Pbc_y = new double [ndim_y];
	Pbc_z = new double [ndim_z];

	imult = 1;
	corr = -No_cells_x;
	for (int i = 1; i < ndim_x; i++) {
		Ncell_bound_x[i] = imult * No_cells_y * No_cells_z;
		Pbc_x[i] = corr;
		if (++imult >= No_cells_x) {
			imult = 0;
			corr += No_cells_x;
		}
	}
	imult = 1;
	corr = -No_cells_y;
	for (int i = 1; i < ndim_y; i++) {
		Ncell_bound_y[i] = imult * No_cells_z;
		Pbc_y[i] = corr;
		if (++imult >= No_cells_y) {
			imult = 0;
			corr += No_cells_y;
		}
	}
	imult = 1;
	corr = -No_cells_z;
	for (int i = 1; i < ndim_z; i++) {
		Ncell_bound_z[i] = imult;
		Pbc_z[i] = corr;
		if (++imult >= No_cells_z) {
			imult = 0;
			corr += No_cells_z;
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
	iprec = (int) (-log10 (DensityNominal - Pactual));
	if (iprec >= Nsf) {
		Relax *= 0.5;
		if (++Nsf >= Nfig) doExit = true;
	}
}

void Packing::calc_forces() {

	Dinner = Douter*Douter;
	int i;
	
	#pragma omp parallel  for private(i)
	for (i = 0; i < No_cells; i++)
		Link_head[i] = -1;
	
	#pragma omp parallel  for private(i)
	for (i = 0; i < NumParts; i++) {
		particles[i]->set_force();
		particles[i]->set_forceu();
	}
	
	int ipart;

	#pragma omp parallel  for private(ipart) schedule(static)
	for (ipart = 0; ipart < NumParts; ipart++) {
		Species* k = particles[ipart]->getSpecies();
		Rcut = 0.55 * Douter * ((k->getRot())[0] + Rc_max);
		Rcut2 = Rcut * Rcut;

		if ((int) (2.0 * Rcut) < Ncell_min - 1) force_part(ipart);
		else force_all(ipart);
	}
	Dinner = sqrt(Dinner);
}

void Packing::calc_motion() {
	vector3 buff;
	Force_step = 0;
	Epsilon_scl = Epsilon * Douter0;
	int i;

	#pragma omp parallel for private(i,buff)
	for (i = 0; i < NumParts; i++) {
		Ellipsoid *p = particles[i];
		Force_step += sqrt(p->get_f()*p->get_f());
//	translation
		buff = p->get_pos() + p->get_f() * Epsilon_scl;
		buff.pbc(Box);
		p->get_pos() = buff;

//	rotation
		if ((NumSteps+1) % Nrot_step) continue;
		if (fabs(Epsilon_rot) < 1e-5) continue;
		if (p->getSpecies()->getIsSphere()) continue;	//	sphere
		double len = Epsilon_rot * sqrt(p->get_fu()*p->get_fu());
		double cp = 1 / sqrt(1 + len*len);
		double sp = len * cp;
		double cp2 = sqrt(0.5 * (1 + cp));
		double sp2 = 0.5 * sp / cp2;
		Quaternion q(cp2, sp2*p->get_fu());
		p->rotate(q);
	}

	Force_step /= NumParts;
}

void Packing::force_part(int ipart_p) {
	int	jpart, ifirst;
	int	icell_x, icell_xy, icell_y, icell_z, icell;
	int leap_x, leap_y, leap_z;
	Ellipsoid *pi = particles[ipart_p], *pj;
	vector3 rij;
	vector3 pos = pi->get_pos();

	int low_cell_x = (int) (pos[0] + No_cells_x - Rcut);
	int low_cell_y = (int) (pos[1] + No_cells_y - Rcut);
	int low_cell_z = (int) (pos[2] + No_cells_z - Rcut);
	int lim_cell_x = (int) (pos[0] + No_cells_x + Rcut);
	int lim_cell_y = (int) (pos[1] + No_cells_y + Rcut);
	int lim_cell_z = (int) (pos[2] + No_cells_z + Rcut);

	#pragma omp  parallel
    for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
        icell_x = Ncell_bound_x[leap_x];
        for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
            icell_xy = Ncell_bound_y[leap_y] + icell_x;
            #pragma omp for private(leap_z)
            for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
                icell = Ncell_bound_z[leap_z] + icell_xy;
                jpart = Link_head[icell];
                while (jpart != -1) {
                    pj = particles[jpart];
                    rij = pj->get_pos() - pi->get_pos();
                    BinaryEllipsoidSystem eij(*pi, *pj, rij);
                    calc_forces(eij, *pi, *pj);
                    jpart = Link_list[jpart];
                }
            }
        }
    }

	icell_x = (int) pos[0];
	icell_y = (int) pos[1];
	icell_z = (int) pos[2];
	icell = No_cells_z * (No_cells_y * icell_x + icell_y) + icell_z;
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
	int icell_x = (int) pos[0];
	int icell_y = (int) pos[1];
	int icell_z = (int) pos[2];
	int icell = No_cells_z * (No_cells_y * icell_x + icell_y) + icell_z;
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

void Packing::getOutput(vector<vector<vector3> > &positions, vector<vector<vector3> > &angles)
{
    positions.resize(NumSpecies);
    angles.resize(NumSpecies);

    int counter = 0;
    for(int ns = 0; ns < NumSpecies; ns++){
        int numCells = species[ns]->getNum();

        positions[ns].resize(numCells);
        angles[ns].resize(numCells);

        for(int nc = 0; nc < numCells; nc++)
        {
            Ellipsoid *pi = particles[counter];

            vector3 pos = pi->get_pos();
            matrix33 Q = pi->get_q().countQ();
            vector3 euler(atan2(Q(1,2),Q(2,2)), -asin(Q(0,2)), atan2(Q(0,1),Q(0,0)));
            //euler *= 180 / PI; //Rad to Deg

            positions[ns][nc] = vector3(pos[0], pos[1], pos[2]) * (1./Sizing);
            angles[ns][nc] = vector3(euler[0], euler[1], euler[2]);  // Euler angles in RAD

            counter++;
        }

    }
}

//Legacy, needs rewriting
void Packing::savePov(const char *fileName, int wbcNumber)
{
    ofstream povf (fileName);
    povf.setf (ios::fixed, ios::floatfield);
    povf.precision (6);
    

    // Write out ellipsoids
    for (int i = 0; i < NumParts; i++) {
        Ellipsoid *pi = particles[i];
        vector3 pos = pi->get_pos() * (1./Sizing);
        Species* k = pi->getSpecies();
        vector3 rad = 0.5 * Dinner * k->getRot()* (1./Sizing);
        bool sph = k->getIsSphere();

        if(i < k->getNum())
        	povf << "object { RBC" << endl;
        else if(k->getNum() == wbcNumber)
            povf << "object { WBC" << endl;
        else
            povf << "object { PLT" << endl;

        povf << "scale <" <<
        setw(12) << rad[0] << ", " <<
        setw(12) << rad[1] << ", " <<
        setw(12) << rad[2] << "> " << endl;
        if (!sph) {	//	sphere
            matrix33 Q = pi->get_q().countQ();
            vector3 euler(atan2(Q(1,2),Q(2,2)), -asin(Q(0,2)), atan2(Q(0,1),Q(0,0)));
            euler *= 180 / PI;
            povf << "rotate " << setw(12) << euler[0] << "*x" << endl;
            povf << "rotate " << setw(12) << euler[1] << "*y" << endl;
            povf << "rotate " << setw(12) << euler[2] << "*z" << endl;
        }
        povf << "translate <" <<
        setw(12) << pos[0] << ", " <<
        setw(12) << pos[1] << ", " <<
        setw(12) << pos[2] << "> " << endl;

        if(NumSpecies > 1)
        {
            povf << "pigment { color ";
            if(i < k->getNum())
                povf << "Red";                
            else if(k->getNum() == wbcNumber)
                povf << "White";
            else
                povf << "Yellow";
            povf << " }" << endl;
        }

        povf << "}" << endl;

    }
    povf.close();
}

void Packing::saveBloodCellPositions()
{
	int speciesCounter = 0;

	for (int j = 0; j < NumSpecies; j++){

		ofstream cellsFile (species[j]->name + ".pos");

		//cellsFile << No_cells_x << " " << No_cells_y << " " << No_cells_z << endl; // Dimensions

		cellsFile << species[j]->getNum() << endl; // Num. of cells of this type

		for(int i = speciesCounter; i < (speciesCounter+species[j]->getNum() ); i++)
		{
			Ellipsoid *pi = particles[i];

			vector3 pos = pi->get_pos() * (1./Sizing);
			matrix33 Q = pi->get_q().countQ();
			vector3 euler(atan2(Q(1,2),Q(2,2)), -asin(Q(0,2)), atan2(Q(0,1),Q(0,0)));
			euler *= 180 / PI; //Rad to Deg

			//if(i < species[0]->getn())
			cellsFile << pos[0] << " " << pos[1] << " " << pos[2] << " " << euler[0] << " " << euler[1] << " " << euler[2] << endl;
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

double Packing::zeroin(BinaryEllipsoidSystem& ell, double ax, double bx) {
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

	for(;;) {
		double prev_step = b-a;
		double tol_act;
		double p;
		double q;
		double new_step;
		if( fabs(fc) < fabs(fb) ) {
			a = b;
			b = c;
			c = a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol_act = 2*EPSILON*fabs(b) + tol/2;
		new_step = (c-b)/2;

		if( fabs(new_step) <= tol_act || fb == (double)0 )
			return b;

		if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
			register double t1, cb, t2;
			cb = c-b;
			if (a==c) {
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			} else {
				q = fa/fc;
				t1 = fb/fc;
				t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
			if (p > 0.) q = -q;
			else p = -p;

			if (p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2))
				new_step = p/q;
		}

		if (fabs(new_step) < tol_act) {
			if (new_step > 0.) new_step = tol_act;
			else new_step = -tol_act;
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

