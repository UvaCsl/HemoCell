#ifndef FICSIONINIT_H
#define FICSIONINIT_H

#include <limits>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedCellParticleFunctional3D.h"
#include "immersedCellParticleVtk3D.h"
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
        std::vector<Array<T,3> > & centers, std::vector<T> & radii, plint flowType) ;

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
void iniLatticeSquarePoiseuille( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                 T Re);

void iniLatticePoiseuilleWithBodyForce(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition, T poiseuilleForce);

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
                  std::string const& title, T Re, T shearRate, plint flowType);


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

/* ******************* copyXMLreader2XMLwriter ***************************************** */
void copyXMLreader2XMLwriter(XMLreader const& reader, XMLwriter & writer);
void copyXMLreader2XMLwriter(XMLreaderProxy readerProxy, XMLwriter & writer);


#include "ficsionInit.hh"

#endif  // FICSIONINIT_H
