#ifndef CELL_CELL_FORCES_3D_H
#define CELL_CELL_FORCES_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <math.h>

template<typename T>
class CellCellForce3D;



template<typename T, template<typename U> class Descriptor>
class ComputeCellCellForces3D : public BoxProcessingFunctional3D
{
public:
    ComputeCellCellForces3D (CellCellForce3D<T> & calcForce_, T cutoffRadius_);
    virtual ~ComputeCellCellForces3D();
    ComputeCellCellForces3D(ComputeCellCellForces3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeCellCellForces3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
public:
    /*
     All of the following conditions have to be met, in order for the force to be applied:
           * Particles have to belong to different cells
           * Particle 1 has to have larger tagId that particle 2
           * Their distance has to be less that the cutoffRadius.
     */
    bool conditionsAreMet(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T& r, Array<T,3>& eij);
    void applyForce(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T r, Array<T,3> & eij);
private:
    CellCellForce3D<T> & calcForce;
    T cutoffRadius;
};



template<typename T, template<typename U> class Descriptor>
class ComputeDifferentCellForces3D : public BoxProcessingFunctional3D
{
public:
    ComputeDifferentCellForces3D (CellCellForce3D<T> & calcForce_, T cutoffRadius_);
    virtual ~ComputeDifferentCellForces3D();
    ComputeDifferentCellForces3D(ComputeDifferentCellForces3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeDifferentCellForces3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
public:
    /*
     All of the following conditions have to be met, in order for the force to be applied:
           * Particles have to belong to different cells
           * Particle 1 has to have larger tagId that particle 2
           * Their distance has to be less that the cutoffRadius.
     */
    bool conditionsAreMet(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T& r, Array<T,3>& eij);
    void applyForce(Particle3D<T,Descriptor> * p1, Particle3D<T,Descriptor> * p2, T & r, Array<T,3> & eij);
private:
    CellCellForce3D<T> & calcForce;
    T cutoffRadius;
};





template<typename T>
class CellCellForce3D {
public:
    CellCellForce3D() { } ;
    virtual ~CellCellForce3D() { }
    virtual CellCellForce3D<T>* clone() const=0;
    // CellCellForce3D X; X(r, eij) returns the force.
    virtual Array<T,3> operator() (T r, Array<T,3> & eij);
    virtual Array<T,3> operator() (Array<T,3> & x1, Array<T,3> & x2) ;
public:
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) = 0;
    virtual T calculatePotential (T r) = 0;
    virtual T calculatePotential (T r, Array<T,3> & eij);
private:
};


template<typename T>
class PowerLawForce : public CellCellForce3D<T> {
// PowerLaw Force is:
//    F = -k_int * (pow((DeltaX*1.0/r),k) - DeltaXoverRink)*eij;
// with:
//    k_int: is the characteristic energy
//    DeltaX: is the characteristic length
//    k: powerlaw exponent (>1.0, Krueger == 2.0)
//    R: is the characteristic length (= DeltaX)
public:
    PowerLawForce(T k_int_, T DeltaX_, T R_, T k_) : k_int(k_int_), DeltaX(DeltaX_), R(R_), k(k_) {
        DeltaXoverRink = pow((DeltaX*1.0/R),k);
    };
    virtual ~PowerLawForce() { }
    virtual PowerLawForce<T>* clone() const { return new PowerLawForce<T>(k_int, DeltaX, R, k); } ;
public:
    virtual T calculatePotential (T r) { return 0; }
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
        return k_int * (pow((DeltaX*1.0/r),k) - DeltaXoverRink)*eij;
    }
private:
    T k_int, DeltaX, R, k;
    T DeltaXoverRink;
};


template<typename T>
class LennardJonesPotential : public CellCellForce3D<T> {
// Lennard-Jones Potential is:
//    phi = 4*epsilonLJ*(pow(sigmaLJ*1.0/r, 12.0) - pow(sigmaLJ*1.0/r, 6.0));
// with:
//    epsilonLJ: is the characteristic energy
//    sigmaLJ: is the characteristic length
public:
    LennardJonesPotential(T epsilonLJ_, T sigmaLJ_) : epsilonLJ(epsilonLJ_), sigmaLJ(sigmaLJ_) {};
    virtual ~LennardJonesPotential() { }
    virtual LennardJonesPotential<T>* clone() const { return new LennardJonesPotential<T>(epsilonLJ, sigmaLJ); } ;
public:
    virtual T calculatePotential (T r) { T x = pow(sigmaLJ*1.0/r, 6.0); return 4*epsilonLJ*(x*x - x); }
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
        T x = pow(sigmaLJ*1.0/r, 6.0);
        return   (24.0/r) * epsilonLJ * (2*x*x - x) * eij;
    }
private:
    T epsilonLJ, sigmaLJ;
};


template<typename T>
class MorsePotential : public CellCellForce3D<T> {
// Morse Potential is:
//    phi = De*( exp(2*beta*(ro - r)) - 2*exp(beta*(r0 - r)));
// with:
//    De: well depth of the potential    (0.3 kbT ) Values for 500 vertices
//    beta: well depth of the potential  (1.5e6 m)
//    r0: zero force distance            (0.3e-6 m)
public:
    MorsePotential(T De_, T beta_, T r0_)
        : De(De_), beta(beta_), r0(r0_) {};
    MorsePotential(T dx, plint numVerticesPerCell, T kBT, bool useOtherParameters)
    {
        De = 0.3*kBT * (500.0 / numVerticesPerCell); beta=1.5e6*dx; r0=0.5e-6/dx;
    };
    virtual ~MorsePotential() { }
    virtual MorsePotential<T>* clone() const { return new MorsePotential<T>(De, beta, r0); } ;
public:
    virtual T calculatePotential (T r) { return De*( exp(2*beta*(r0 - r)) - 2*exp(beta*(r0 - r))); }
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
            T x = exp( beta * (-r + r0) );
            return 2*beta*De*(x*x  - x) * eij;
    }
private:
    T De, beta, r0;
};

template<typename T>
class JohnsonKendallRobertsForce : public CellCellForce3D<T> {
// Johnson Kendall Roberts force F is:
//    F = K*a*a*a/R - sqrt( 6 * PI * sigmaJKR * K * a*a*a );
// and the penetration depth h:
//    h = a*a/R - 2.0/3.0 * sqrt( 6 * PI * sigmaJKR * a / K);
// with:
//    R: the radius of the tip
//    a: the contact area radius
//    sigmaJKR: the adhesion work
//    K: the effective Young's modulus

public:
    JohnsonKendallRobertsForce() {};
    virtual ~JohnsonKendallRobertsForce() { }
    virtual JohnsonKendallRobertsForce<T>* clone() const { return new JohnsonKendallRobertsForce<T>(); } ;
public:
    virtual T calculatePotential (T r) {return 0.0; }
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
        return Array<T,3>(0,0,0);
    }
private:
};

#ifndef JKR
#define JKR JohnsonKendallRoberts
#endif


template<typename T>
class MorseAndPowerLawForce: public CellCellForce3D<T> {
// Only the Morse contribution is accounted for the potential.
public:
    MorseAndPowerLawForce(T De_, T beta_, T r0_, T k_int_, T DeltaX_, T R_, T k_)
        : De(De_), beta(beta_), r0(r0_),
          k_int(k_int_), DeltaX(DeltaX_), R(R_), k(k_)
    {};
    MorseAndPowerLawForce(T dx, plint numVerticesPerCell, T kBT, bool useOtherParameters)
    {
        De = 2*kBT * (500.0 / numVerticesPerCell); beta=1.5e6*dx; r0=0.3e-6/dx;
        k_int = 5e5*kBT * (500.0 / numVerticesPerCell); DeltaX=2.0; R=dx, k = 0.5*dx;
        DeltaXoverRink = pow((DeltaX*1.0/R),k);
    };
    virtual ~MorseAndPowerLawForce() { }
    virtual MorseAndPowerLawForce<T>* clone() const { return new MorseAndPowerLawForce<T>(De, beta, r0, k_int, DeltaX, R, k); } ;
public:
    virtual T calculatePotential (T r) { // Only the Morse contribution is accounted for the potential.
        return De*( exp(2*beta*(r0 - r)) - 2*exp(beta*(r0 - r)));
    }
    virtual Array<T,3> calculateForce (T r, Array<T,3> & eij) {
            T x = exp( beta * (-r + r0) );
            Array<T,3> force = 2*beta*De*(x*x  - x) * eij;
            force = force + k_int * (pow((DeltaX*1.0/r),k) - DeltaXoverRink)*eij;
            return force;
    }
private:
    T De, beta, r0;
    T k_int, DeltaX, R, k;
    T DeltaXoverRink;
};


#include "cellCellForces3D.hh"

#endif  // CELL_CELL_FORCES_3D_H

