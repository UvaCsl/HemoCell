#ifndef HEMOCELLINIT_HH
#define HEMOCELLINIT_HH

template <typename T>
class CouetteDensityAndVelocity {
public:
    CouetteDensityAndVelocity(T shearRateLU_, T zeroCoordinate_);
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const ;
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
void CouetteDensityAndVelocity<T>::operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
          rho = (T)1;
          u[0] = couetteVelocity(iY, shearRateLU, zeroCoordinate);
          u[1] = T();
          u[2] = T();
}

/* ************* iniLatticeSquareCouette ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeSquareCouette( MultiBlockLattice3D<T,Descriptor>& lattice,
                 plint nx, plint ny, plint nz,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T shearRate)
{

    Box3D top   = Box3D(0, nx-1, 0,  0,    0, nz-1);
    Box3D bottom  = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);

    //lattice.periodicity().toggleAll(true);
    lattice.periodicity().toggle(0, true);
    lattice.periodicity().toggle(1, false);
    lattice.periodicity().toggle(2, true);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    T vHalf = (ny-1)*shearRate*0.5;
    setBoundaryVelocity(lattice, top, Array<T,3>(-vHalf,0.0,0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>(vHalf,0.0,0.0));

    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(0.0,0.0,0.0));

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), CouetteDensityAndVelocity<T>(shearRate, (ny-1)*0.5 ));

    lattice.initialize();
}

#endif  // HEMOCELLINIT_HH
