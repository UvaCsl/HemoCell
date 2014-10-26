#ifndef FICSIONINIT_H
#define FICSIONINIT_H

#include <limits>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedCellParticleFunctional3D.h"
#include "shellModel3D.h"
#include "cellModel3D.h"

using namespace plb;
using namespace std;


typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

#define NMAX 50

#ifndef PI__
#define PI__
const T pi = 4.*atan(1.);
#endif  // PI__

// Position Cells inside the domain
template<typename T>
void positionCells(plint shape, T radius, plint & npar, IncomprFlowParam<T> const& parameters,
        std::vector<Array<T,3> > & centers, std::vector<T> & radii, T dx, plint flowType) ;

/* ************* Functions poiseuillePressure and poiseuilleVelocity ******************* */

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN, T Re);
T poiseuilleVelocity(plint iY, plint iZ, IncomprFlowParam<T> const& parameters, plint maxN, T Re);

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_, T Re_);
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const ;
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
    T Re;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_, T Re_);
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  ;
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
    T Re;
};

/* ************************************************************************************* */

/* ******************* iniLattice ***************************************** */
void iniLatticeOutlets( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition);

void iniLatticeSquarePoiseuille( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                 T Re);

void iniLatticePoiseuilleWithBodyForce(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition, T poiseuilleForce);

template<typename T, template<class U> class Descriptor>
void iniLatticeFullyPeriodic(MultiBlockLattice3D<T,Descriptor>& lattice, IncomprFlowParam<T> const& parameters, Array<T,3> uInit);

void iniLatticeSquareCouette( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition, T shearRate);

/* ******************* changeCouetteShearRate ***************************************** */
void changeCouetteShearRate( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
        IncomprFlowParam<T> const& parameters,
        OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition, T shearRate);

/* ******************* writeFicsionLogFile ***************************************** */
template<typename T>
void writeFicsionLogFile(IncomprFlowParam<T> const& parameters,
                  std::string const& title, T Re, T shearRate, plint flowType,
                  plint npar, T volumePerCell_LU);


/* ******************* writeVTK ***************************************** */
template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter);



/* ************* Class GetTensorFieldFromExternalVectorFunctional3D ******************* */
template<typename T, template<typename U> class Descriptor, int nDim>
class GetTensorFieldFromExternalVectorFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor, T, nDim> {
public:
    GetTensorFieldFromExternalVectorFunctional3D (int vectorStartsAt_ );
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,nDim>& tensor);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const ;
    virtual GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>* clone() const ;
private:
    int vectorStartsAt;
};


#endif  // FICSIONINIT_H
