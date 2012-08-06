#ifndef FICSIONINIT_HH
#define FICSIONINIT_HH

#include <limits>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedWallParticle3D.h"
#include "immersedWallParticle3D.hh"
#include "immersedWallParticleFunctional3D.h"
#include "immersedWallParticleFunctional3D.hh"
#include "immersedWallParticleVtk3D.h"
#include "immersedWallParticleVtk3D.hh"
#include "shellModel3D.h"
#include "shellModel3D.hh"

using namespace plb;
using namespace std;


typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

#define NMAX 50

const T pi = 4.*atan(1.);


static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
	
    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();
	
    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a))));
    }
	
    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu
	
    T deltaP = - (alpha * nu);
	
    return deltaP;
}

T poiseuilleVelocity(plint iY, plint iZ, IncomprFlowParam<T> const& parameters, plint maxN)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
	
    const T y = (T)iY - a / (T)2;
    const T z = (T)iZ - b / (T)2;
	
    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();
	
    T sum = T();
	
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
		
        sum += (cos(twoNplusOne*pi*y/a)*cosh(twoNplusOne*pi*z/a)
				/ ( pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
		
        sum -= (cos(twoNplusOne*pi*y/a)*cosh(twoNplusOne*pi*z/a)
				/ ( pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
	
    sum *= ((T)4 * alpha * a *a /pow(pi,3));
    sum += (alpha / (T)2 * (y * y - a*a / (T)4));
	
    return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
	: parameters(parameters_),
	maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = poiseuilleVelocity(iY, iZ, parameters, maxN);
        u[1] = T();
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
	: parameters(parameters_),
	maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = poiseuilleVelocity(iY, iZ, parameters, maxN);
        u[1] = T();
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};


#endif  // FICSIONINIT_HH
