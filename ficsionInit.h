#ifndef FICSIONINIT_H
#define FICSIONINIT_H

#include <limits>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedCellParticle3D.hh"
#include "immersedCellParticleFunctional3D.h"
#include "immersedCellParticleFunctional3D.hh"
#include "immersedCellParticleVtk3D.h"
#include "immersedCellParticleVtk3D.hh"
#include "shellModel3D.h"
#include "shellModel3D.hh"
#include "cellModel3D.h"
#include "cellModel3D.hh"

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
void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                 T Re);

/* ******************* writeFicsionLogFile ***************************************** */
template<typename T>
void writeFicsionLogFile(IncomprFlowParam<T> const& parameters,
                  std::string const& title, T Re);


/* ******************* writeVTK ***************************************** */
template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter);

/// Export the surface mesh as an ASCII STL file considering the dx.
template<typename T> // Copied from void TriangularSurfaceMesh<T>::writeAsciiSTL(std::string fname).
void writeMeshAsciiSTL(TriangleBoundary3D<T> & Cells, std::string fname, T dx=0.0);

/// Export the surface mesh as an binary STL file considering the dx.
template<typename T> // Copied from void TriangularSurfaceMesh<T>::writeBinarySTL(std::string fname).
void writeMeshBinarySTL(TriangleBoundary3D<T> & Cells, std::string fname, T dx=0.0);


#endif  // FICSIONINIT_H
