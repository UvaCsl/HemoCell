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
#ifndef HEMOCELLINIT_HH
#define HEMOCELLINIT_HH

#include "palabos3D.h"
#include "palabos3D.hh"

template <typename T>
class CouetteDensityAndVelocity {
public:
    CouetteDensityAndVelocity(T shearRateLU_, T zeroCoordinate_);
    void operator()(plint iX, plint iY, plint iZ, T &rho, hemo::Array<T,3>& u) const ;
private:
    T shearRateLU;
    T zeroCoordinate;
};

T couetteVelocity(plint iY, T shearRateLU, T zeroCoordinate)
{
      return shearRateLU*(iY - zeroCoordinate);
}
template <typename T>
CouetteDensityAndVelocity<T>::CouetteDensityAndVelocity(T shearRateLU_, T zeroCoordinate_)
: shearRateLU(shearRateLU_), zeroCoordinate(zeroCoordinate_)
{ }

template <typename T>
void CouetteDensityAndVelocity<T>::operator()(plint iX, plint iY, plint iZ, T &rho, hemo::Array<T,3>& u) const {
          rho = (T)1;
          u[0] = couetteVelocity(iY, shearRateLU, zeroCoordinate);
          u[1] = T();
          u[2] = T();
}

/* ************* iniLatticeSquareCouette ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeSquareCouette(plb::MultiBlockLattice3D<T,Descriptor>& lattice,
                 plint nx, plint ny, plint nz,
                 plb::OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T shearRate)
{

    plb::Box3D top = plb::Box3D(0, nx-1, 0, ny-1, nz-1, nz-1);
    plb::Box3D bottom = plb::Box3D(0, nx-1, 0, ny-1, 0, 0);

    //lattice.periodicity().toggleAll(true);
    lattice.periodicity().toggle(0, true);
    lattice.periodicity().toggle(1, true);
    lattice.periodicity().toggle(2, false);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    T vHalf = (nz-1)*shearRate*0.5;
    setBoundaryVelocity(lattice, top, plb::Array<T,3>(-vHalf,0.0,0.0));
    setBoundaryVelocity(lattice, bottom, plb::Array<T,3>(vHalf,0.0,0.0));

    // plb::Cell<double, plb::descriptors::ForcedD3Q19Descriptor> cell = lattice.get(0,0,0);
    // // plb::Dynamics<double, plb::descriptors::ForcedD3Q19Descriptor> cellDynamics = cell.getDynamics();
    // plb::Array<double, 19ul> rawPop = cell.getRawPopulations();
    // cout << rawPop[0] << std::endl;


    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, plb::Array<T,Descriptor<T>::d>(0.0,0.0,0.0));

    // For shearing test this results in a dirac delta in force -> viscous force will not like it
    //initializeAtEquilibrium(lattice, lattice.getBoundingBox(), CouetteDensityAndVelocity<T>(shearRate, (ny-1)*0.5 ));

    lattice.initialize();
}

#endif  // HEMOCELLINIT_HH
