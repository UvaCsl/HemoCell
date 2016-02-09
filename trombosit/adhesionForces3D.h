#ifndef ADHESION_FORCES_3D_H
#define ADHESION_FORCES_3D_H

//#include "palabos3D.h"
//#include "palabos3D.hh"
#include <math.h>
using namespace std;
using namespace plb;




//template<typename T>
//class CellCellForce3D;

//template<typename T, template<typename U> class Descriptor>
//class CellField3D;

// Aleksey Belyaev: This file contains classes for adhesive interactions
//                  that inherit from CellCellForce3D.
//                  - AdhesiveLennardJonesPotential - 
//                  - AdhesiveMorsePotential - 
//                  - AdhesiveFENEForce -     

template<typename T>
class AdhesiveLennardJonesPotential : public CellCellForce3D<T> {
// Lennard-Jones Potential is:
//    phi = 4*epsilonLJ*(pow(sigmaLJ*1.0/r, 12.0) - pow(sigmaLJ*1.0/r, 6.0));
// with:
//    epsilonLJ: is the characteristic energy
//    sigmaLJ: is the characteristic length
//    r_cut:   cut-off distance
public:
    AdhesiveLennardJonesPotential(T epsilonLJ_, T sigmaLJ_, T r_cut_) : epsilonLJ(epsilonLJ_), sigmaLJ(sigmaLJ_), r_cut(r_cut_) {};
    virtual ~AdhesiveLennardJonesPotential() { }
    virtual AdhesiveLennardJonesPotential<T>* clone() const { return new AdhesiveLennardJonesPotential<T>(epsilonLJ, sigmaLJ, r_cut); } ;
public:
    virtual T calculatePotential (T r) { 
        T x = pow(sigmaLJ*1.0/r, 6.0); return 4*epsilonLJ*(x*x - x); 
    };
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
        T x = pow(sigmaLJ*1.0/r, 6.0);
        T force=0.0;
        if (r<r_cut) {
           force +=  (24.0/r) * epsilonLJ * (2*x*x - x);
        }
        return   force*eij;
    };
private:
    T epsilonLJ, sigmaLJ, r_cut;
};


template<typename T>
class AdhesiveMorsePotential : public CellCellForce3D<T> {
// Morse Potential is:
//    phi = De*( exp(2*beta*(r0 - r)) - 2*exp(beta*(r0 - r)));
// with:
//    De: well depth of the potential    (0.3 kbT ) Values for 500 vertices
//    beta: well depth of the potential  (1.5e6 m)
//    r0: zero force distance            (0.3e-6 m)
//    r_cut:  cut-off distance           (0.9e-6 m)
public:
    AdhesiveMorsePotential(T De_, T beta_, T r0_, T r_cut_)
        : De(De_), beta(beta_), r0(r0_), r_cut(r_cut_) {};
    AdhesiveMorsePotential(T dx, plint numVerticesPerCell, T kBT, bool useOtherParameters)
    {
        De = 0.3*kBT * (500.0 / numVerticesPerCell); beta=1.5e6*dx; r0=0.3e-6/dx; r_cut=0.9e-6/dx;
    };
    virtual ~AdhesiveMorsePotential() { }
    virtual AdhesiveMorsePotential<T>* clone() const { return new AdhesiveMorsePotential<T>(De, beta, r0, r_cut); } ;
public:
    virtual T calculatePotential (T r) { 
         return De*( exp(2*beta*(r0 - r)) - 2*exp(beta*(r0 - r))); 
    };
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
            T x = exp( beta * (-r + r0) );
            T force=0.0;
            if (r<r_cut) {
                force=2.0*beta*De*(x*x  - x); 
            }
            return force* eij;
    };
private:
    T De, beta, r0, r_cut;
};


template<typename T>
class AdhesiveFENEForce: public CellCellForce3D<T> {
// FENE potential is:
//    phi = -0.5 * H* el_max * log(1-((r-r0)*1.0/el_max)*((r-r0)*1.0/el_max) )  +
//           + 4*epsilonLJ*(pow(sigmaLJ*1.0/r, 12.0) - pow(sigmaLJ*1.0/r, 6.0)) + epsilonLJ;
// FENE force is:
//    f  =  H*el_max*1.0/(1.0 - ((r-r0)*1.0/el_max)*((r-r0)*1.0/el_max) ) + f_LJ    
// with:
//    H:     stiffness of a bond
//    el_max:  maximal possible elongation of a bond, that means length_max = el_max + r0 
//    r0:    zero force distance
//    r_cut:  cut-off distance - we need it here, because force diverges at r=r0+el_max
//    epsilonLJ: is the characteristic energy
//    sigmaLJ: is the characteristic length
// so that potential has a repulsive contribution from Lennard Jones potential
// We might want to make this repulsion switchable (on/off):
//    flag_rep:  bool flag for turn on/off (true/false) the repulsion from LJ        
public:
    AdhesiveFENEForce(T H_, T el_max_, T r0_, T r_cut_, bool flag_rep_, T epsilonLJ_, T sigmaLJ_)
        : H(H_), el_max(el_max_), r0(r0_), r_cut(r_cut_), flag_rep(flag_rep_),
          epsilonLJ(epsilonLJ_), sigmaLJ(sigmaLJ_) 
    {};
    AdhesiveFENEForce(T dx, plint numVerticesPerCell, T kBT, bool useOtherParameters)
    {
        H = 0.3*kBT * (500.0 / numVerticesPerCell); 
        el_max= 1.0e-5/dx; 
        r0 = 0.3e-6/dx;
        r_cut = 0.9e-5/dx;
        flag_rep= false;
        epsilonLJ = 0.0;
        sigmaLJ = 0.0;
    };
    virtual ~AdhesiveFENEForce() { }
    virtual AdhesiveFENEForce<T>* clone() const { return new AdhesiveFENEForce<T>(H, el_max, r0, r_cut, flag_rep, epsilonLJ, sigmaLJ); } ;
public:
    virtual T calculatePotential (T r) { // Only the attractive contribution is accounted for the potential.
        return -0.5* H* el_max * log(1-((r-r0)*1.0/el_max)*((r-r0)*1.0/el_max) ) ;
    };
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
            T force = 0.0;
            if ( (r>r0) && (r<r_cut) && ((r-r0)<el_max) ) {
                force+= -H*(r-r0)*1.0/(1.0 - ( (r-r0)*1.0/el_max )*( (r-r0)*1.0/el_max ) );
            }
            if ( (flag_rep) && (r<r_cut) ){
                            T x = pow(sigmaLJ*1.0/r, 6.0);
                            force+= (24.0/r) * epsilonLJ * (2*x*x - x);
            } 
            return force * eij;
    };
private:
    T H, el_max, r0, r_cut;
    T epsilonLJ, sigmaLJ;
    bool flag_rep;
};



#endif  // ADHESION_FORCES_3D_H

