#ifndef FICSIONINIT_HH
#define FICSIONINIT_HH

#include "ficsionInit.h"

// Position Cells inside the domain
template<typename T>
void positionCells(plint shape, T radius, plint & npar, IncomprFlowParam<T> const& parameters,
        std::vector<Array<T,3> > & centers, std::vector<T> & radii, T dx, plint flowType) {
    /*  Fills slices first and then the whole grid.     */

    std::vector<T> posX, posY, posZ;
    const plint nx = parameters.getNx() - 1 ;
    const plint ny = parameters.getNy()  - 1;
    const plint nz = parameters.getNz()  - 1;
    T dX = 2.05 * radius ;
    T dY = 2.05 * radius * ( (shape==1) ? 0.4 : 1 );
    T dZ = 2.05 * radius;
    plint NdX, NdY, NdZ;
    if (flowType == 2) { // If RBC disaggregation
        dY = 2.7e-6/dx;
        npar = min(npar, plint((ny-2)*0.5/dY - 1) );

        for (plint iA = 0; iA < npar; ++iA) {
            centers.push_back(Array<T,3>(0.5*nx, 0.5*ny + iA*dY, 0.5*nz));
            radii.push_back(radius);
        }


    } else {
        NdX = (nx-2)*1.0/dX;
        NdY = (ny-2)*1.0/dY;
        NdZ = (nz-2)*1.0/dZ;
        // dX = (nx-2.0)*1.0/NdX; // Needs to be re-adjusted.
        dY = (ny-2.0)*1.0/NdY;
        dZ = (nz-2.0)*1.0/NdZ;

        npar = npar<(NdX*NdY*NdZ)?npar:(NdX*NdY*NdZ);
        plint slices = npar/(NdY*NdZ);
        if (slices > 0) { dX = (nx-2.0)*1.0/slices; }
        else { dX = nx; }

        for (plint i = 0; i < slices; ++i) {
            for (plint iy = 0; iy < NdY; ++iy) {
                for (plint iz = 0; iz < NdZ; ++iz) {
                    posX.push_back((i+0.5)*dX);
                    posY.push_back((iy+0.5)*dY);
                    posZ.push_back((iz+0.5)*dZ);
                }
            }
        }
        plint lastSlice = npar%(NdY*NdZ);
        plint rows = lastSlice/NdZ;
        for (plint iy = 1; iy <= rows; ++iy) {
            for (plint iz = 0; iz < NdZ; ++iz) {
                posX.push_back((slices+0.5)*dX);
                posY.push_back(ny * iy*1.0/(rows+1 + 1.0));
                posZ.push_back((iz+0.5)*dZ);
            }
        }

        plint mods = lastSlice%NdZ;
        for (plint iz = 1; iz <= mods; ++iz) {
            posX.push_back((slices+0.5)*dX);
            posY.push_back(ny * (rows+1)*1.0/(rows+1 + 1.0));
            posZ.push_back(nz * iz*1.0/(mods + 1.0));
        }

        T addToX = 0.5;
        T addToY = 0.0;
    //    addToX = (NdX - slices) * dX * 0.5 ;

        for (pluint iA = 0; iA < posX.size(); ++iA) {
            centers.push_back(Array<T,3>(posX[iA]+addToX,posY[iA]+addToY,posZ[iA]));
            radii.push_back(radius);
        }
    }
}

/* ************* Functions poiseuillePressure and poiseuilleVelocity ******************* */
static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN, T Re)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
    const T nu = parameters.getLatticeNu();
    const T uMax = Re * nu/min(a,b);
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

T poiseuilleVelocity(plint iY, plint iZ, IncomprFlowParam<T> const& parameters, plint maxN, T Re)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
    const T y = (T)iY - a / (T)2;
    const T z = (T)iZ - b / (T)2;
    const T alpha = - poiseuillePressure(parameters,maxN, Re) / parameters.getLatticeNu();
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
SquarePoiseuilleDensityAndVelocity<T>::SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_, T Re_)
: parameters(parameters_), maxN(maxN_), Re(Re_)
{ }

template <typename T>
void SquarePoiseuilleDensityAndVelocity<T>::operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = poiseuilleVelocity(iY, iZ, parameters, maxN, Re);
        u[1] = T();
        u[2] = T();
}


template <typename T>
SquarePoiseuilleVelocity<T>::SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_, T Re_)
: parameters(parameters_), maxN(maxN_), Re(Re_)
{ }

template <typename T>
void SquarePoiseuilleVelocity<T>::operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
    u[0] = poiseuilleVelocity(iY, iZ, parameters, maxN, Re);
    u[1] = T();
    u[2] = T();
}

/* ************* iniLatticeSquarePoiseuille ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeOutlets( MultiBlockLattice3D<T,Descriptor>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D top    = Box3D(0,    nx-1, 0, ny-1, 0,    0);
    Box3D bottom = Box3D(0,    nx-1, 0, ny-1, nz-1, nz-1);

    Box3D left   = Box3D(0, nx-1, 0,    0,    0, nz-1);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);

    Box3D inlet  = Box3D(0,    0,    0,    ny-1, 0, nz-1);
    Box3D outlet = Box3D(nx-1, nx-1, 0,    ny-1, 0, nz-1);


    boundaryCondition.setVelocityConditionOnBlockBoundaries (
                                   lattice, top, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
                                   lattice, bottom, boundary::outflow );

    boundaryCondition.setVelocityConditionOnBlockBoundaries (
                                   lattice, left, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
                                   lattice, right, boundary::outflow );

    boundaryCondition.setVelocityConditionOnBlockBoundaries (
                                   lattice, inlet, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
                                   lattice, outlet, boundary::outflow );

    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(0.0,0.0,0.0));

    lattice.initialize();
}


/* ************* iniLatticeSquarePoiseuille ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeSquarePoiseuille( MultiBlockLattice3D<T,Descriptor>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T Re)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D top    = Box3D(0,    nx-1, 0, ny-1, 0,    0);
    Box3D bottom = Box3D(0,    nx-1, 0, ny-1, nz-1, nz-1);

    Box3D left   = Box3D(0, nx-1, 0,    0,    1, nz-2);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 1, nz-2);

    Box3D inlet  = Box3D(0,    0,    1,    ny-2, 1, nz-2);
    Box3D outlet = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );

    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX, Re));
    setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX, Re));

    setBoundaryVelocity(lattice, top, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>(0.0,0.0,0.0));

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX, Re));

    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(0.0,0.0,0.0));

    lattice.initialize();
}

/* ************* iniLatticeSquarePoiseuille ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticePoiseuilleWithBodyForce(MultiBlockLattice3D<T,Descriptor>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T poiseuilleForce)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D left   = Box3D(0, nx-1, 0,    0,    1, nz-2);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 1, nz-2);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );
    setBoundaryVelocity(lattice, left, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>(0.0,0.0,0.0));

    T rhoInit=1.0; Array<T,3> uInit(0.0,0.0,0.0);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rhoInit, uInit);

    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(poiseuilleForce, 0.0, 0.0));

    lattice.initialize();
}

/* ************* iniLatticeFullyPeriodic ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeFullyPeriodic(MultiBlockLattice3D<T,Descriptor>& lattice, IncomprFlowParam<T> const& parameters, Array<T,3> uInit)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    T rhoInit=1.0;
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rhoInit, uInit);

    lattice.periodicity().toggleAll(true);
    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(0.0, 0.0, 0.0));
    lattice.initialize();
}

/* ************* iniLatticeSquareCouette ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeSquareCouette( MultiBlockLattice3D<T,Descriptor>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T shearRate)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D left   = Box3D(0, nx-1, 0,    0,    0, nz-1);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);

    lattice.periodicity().toggleAll(true);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );

    T vHalf = (ny-1)*shearRate*0.5;
    setBoundaryVelocity(lattice, left, Array<T,3>(vHalf,0.0,0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>(-vHalf,0.0,0.0));

    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(0.0,0.0,0.0));

    lattice.initialize();
}

/* ************* iniLatticeSquareCouetteMeasureStress ******************* */
template<typename T, template<class U> class Descriptor>
void iniLatticeSquareCouetteMeasureStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition, T shearRate,
                 Array<plint, 3> & forceIds, plint & nMomentumExchangeCells)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D top    = Box3D(0, nx-1,    0,    0,   0, nz-1);
    Box3D bottom = Box3D(0, nx-1, ny-2, ny-2,   0, nz-1);
    Box3D lid    = Box3D(0, nx-1, ny-1, ny-1,   0, nz-1);

    lattice.periodicity().toggleAll(true);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
//    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, lid );

    T v = (ny-2.5)*shearRate;
    setBoundaryVelocity(lattice, top, Array<T,3>(v, 0.0, 0.0));
//    setBoundaryVelocity(lattice, lid, Array<T,3>(0.0, 0.0,0.0));

    defineDynamics(lattice, bottom, new MomentumExchangeBounceBack<T,Descriptor>(forceIds));
    initializeMomentumExchange(lattice, lattice.getBoundingBox());
    defineDynamics(lattice, lid, new BounceBack<T,Descriptor>);

//    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX, Re));
    nMomentumExchangeCells = bottom.nCells();
    setExternalVector( lattice, lattice.getBoundingBox(),
            Descriptor<T>::ExternalField::forceBeginsAt, Array<T,Descriptor<T>::d>(0.0,0.0,0.0));

    lattice.initialize();
}


/* ************* changeCouetteShearRate ******************* */
void changeCouetteShearRate( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
        IncomprFlowParam<T> const& parameters,
        OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition, T shearRate)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D left   = Box3D(0, nx-1, 0,    0,    1, nz-2);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 1, nz-2);
//    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
//    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );
    T vHalf = (nz-1)*shearRate*0.5;
    setBoundaryVelocity(lattice, left, Array<T,3>(vHalf,0.0,0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>(-vHalf,0.0,0.0));
//    lattice.initialize();
}


/* ************* writeFicsionLogFile ******************* */
template<typename T>
void writeFicsionLogFile(IncomprFlowParam<T> const& parameters,
                  std::string const& title, T Re, T shearRate, plint flowType,
                  plint npar, T volumePerCell_LU)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
    const T nu = parameters.getLatticeNu();
    T uMax ;
    if (flowType != 1) {
        uMax = Re * nu/min(a,b);
    } else {
        uMax = b * shearRate / 2.0;
    }
    const T dx = parameters.getDeltaX();
    const T dt = parameters.getDeltaT();
    std::string fullName = global::directories().getLogOutDir() + "plbLog.log";
    plb_ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Velocity [m/s]:              u=" << parameters.getLatticeU()*dx*1.0/dt << "\n";
    ofile << "Max Velocity [m/s]:          uMax=" << uMax*dx*1.0/dt << "\n";
    ofile << "Mach number [x100%]:         Ma=" << uMax*sqrt(3.0)*100 << "\n";
    if (flowType != 1) {
        ofile << "Reynolds number:             Re=" << Re << "\n";
        ofile << "Wall Shear rate [1/s]:       ydot=" << 4*uMax/a/dt << "\n";
    } else {
        ofile << "Reynolds number (shear):     Re=" << shearRate*pow(min(a,b),2)/nu << "\n";
        ofile << "Shear rate [1/s]:            ydot=" << shearRate*1.0/dt << "\n";
    }
    ofile << "Lattice resolution:          N=" << parameters.getResolution() << "\n";
    ofile << "Relaxation frequency:        omega=" << parameters.getOmega()*1.0 << "\n";
    ofile << "Relaxation time:             tau=" << parameters.getTau()*1.0 << "\n";
    ofile << "Grid spacing deltaX [m]:     dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:            dt=" << parameters.getDeltaT() << "\n";
    ofile << "Lattice  Nu:                 nu=" << nu << "\n";
    ofile << "Physical Nu [m2/s]:          nu_p=" << nu*dx*dx/dt << "\n";
    ofile << "Extent of the system [m]:    lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system [m]:    ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system [m]:    lz=" << parameters.getLz() << "\n";
    ofile << "Number of Cells:             Np=" << npar << "\n";
    T H = (npar*volumePerCell_LU) * (dx*dx*dx) / (parameters.getLx() * parameters.getLy() * parameters.getLz());
    ofile << "Volume fraction [x100%]:     H=" << H*100.0 << "\n";
}


/* ************* Class GetTensorFieldFromExternalVectorFunctional3D ******************* */
template<typename T, template<typename U> class Descriptor, int nDim>
GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>::GetTensorFieldFromExternalVectorFunctional3D (
    int vectorStartsAt_ ) : vectorStartsAt(vectorStartsAt_)
{
    PLB_ASSERT( vectorStartsAt+nDim <=
    Descriptor<T>::ExternalField::numScalars );
};

template<typename T, template<typename U> class Descriptor, int nDim>
void GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,nDim>& tensor) {
    Dot3D offset = computeRelativeDisplacement(lattice, tensor);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                Array<T,nDim> externalVector;

                for (plint iD=0; iD<nDim; ++iD) {
                    externalVector[iD] = *cell.getExternal(vectorStartsAt+iD);
                }
                tensor.get(oX,oY,oZ) = externalVector;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int nDim>
void GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor, int nDim>
GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>* GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>::clone() const {
    return new GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>(*this);
}



/* ************* writeVTK ******************* */
template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    MultiTensorField3D<T,3> force(lattice);
    applyProcessingFunctional(new GetTensorFieldFromExternalVectorFunctional3D<T,DESCRIPTOR,3>(
        DESCRIPTOR<T>::ExternalField::forceBeginsAt), lattice.getBoundingBox(), lattice, force);

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 8), dx);
    vtkOut.writeData<T>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<T>(*computeDensity(lattice), "density", 1.0);
    vtkOut.writeData<3,T>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,T>(force, "force",  (T) 1.0);
//    vtkOut.writeData<T>(*computeNorm(force, force.getBoundingBox()), "forceNorm",  dx/dt/dt);
    vtkOut.writeData<T>(*computeNorm(force, force.getBoundingBox()), "forceNorm LU",  1.0);

//    ImageWriter<T> imageWriter("leeloo");
//    add(force, forceScalar, force.getBoundingBox());
//    imageWriter.writeScaledPpm(scalarField, createFileName("PPM", iter, 8));

//     vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
//     vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}


template<typename T>
void writeMeshAsciiSTL(TriangleBoundary3D<T> & Cells, std::string fname, T dx=0.0)
{
    TriangularSurfaceMesh<T> mesh = Cells.getMesh();
    if (0.0 == dx) {
        dx = Cells.getDx();
    }
    // Output only from one MPI process.
    if (!global::mpi().isMainProcessor()) {
        return;
    }
    FILE *fp = fopen(fname.c_str(), "w");
    PLB_ASSERT(fp != NULL);

    char fmt1[64] = "  facet normal ";
    char fmt2[64] = "      vertex ";
    if (sizeof(T) == sizeof(long double)) {
        strcat(fmt1, "% Le % Le % Le\n");
        strcat(fmt2, "% Le % Le % Le\n");
    }
    else if (sizeof(T) == sizeof(float) ||
             sizeof(T) == sizeof(double)) {
        strcat(fmt1, "% e % e % e\n");
        strcat(fmt2, "% e % e % e\n");
    }
    else {
        PLB_ASSERT(false);
    }

    fprintf(fp, "solid surface\n");
    for (plint i = 0; i < mesh.getNumTriangles(); i++) {
        Array<T,3> n = mesh.computeTriangleNormal(i);
        Array<T,3> v;
        fprintf(fp, fmt1, n[0], n[1], n[2]);
        fprintf(fp, "    outer loop\n");
        v = dx * mesh.getVertex(i, 0);
        fprintf(fp, fmt2, v[0], v[1], v[2]);
        v = dx * mesh.getVertex(i, 1);
        fprintf(fp, fmt2, v[0], v[1], v[2]);
        v = dx * mesh.getVertex(i, 2);
        fprintf(fp, fmt2, v[0], v[1], v[2]);
        fprintf(fp, "    endloop\n");
        fprintf(fp, "  endfacet\n");
    }
    fprintf(fp, "endsolid surface\n");

    fclose(fp);
}

template<typename T>
void writeMeshBinarySTL(TriangleBoundary3D<T> & Cells, std::string fname, T dx=0.0)
{
    TriangularSurfaceMesh<T> mesh = Cells.getMesh();
    if (0.0 == dx) {
        dx = Cells.getDx();
    }
    // Output only from one MPI process.
    if (!global::mpi().isMainProcessor()) {
        return;
    }
    FILE *fp = fopen(fname.c_str(), "wb");
    PLB_ASSERT(fp != NULL);

    unsigned int nt = (unsigned int) mesh.getNumTriangles();
    unsigned short abc = 0;
    char buf[80];

    for (int i = 0; i < 80; i++)
        buf[i] = '\0';

    fwrite(buf, sizeof(char), 80, fp);
    fwrite(&nt, sizeof(unsigned int), 1, fp);
    for (plint i = 0; i < mesh.getNumTriangles(); i++) {
        Array<T,3> vertex;
        Array<T,3> normal = mesh.computeTriangleNormal(i);
        float n[3];
        n[0] = normal[0];
        n[1] = normal[1];
        n[2] = normal[2];
        fwrite((void *) n, sizeof(float), 3, fp);
        vertex = dx * mesh.getVertex(i, 0);
        float v[3];
        v[0] = vertex[0];
        v[1] = vertex[1];
        v[2] = vertex[2];
        fwrite((void *) v, sizeof(float), 3, fp);
        vertex = dx * mesh.getVertex(i, 1);
        v[0] = vertex[0];
        v[1] = vertex[1];
        v[2] = vertex[2];
        fwrite((void *) v, sizeof(float), 3, fp);
        vertex = dx * mesh.getVertex(i, 2);
        v[0] = vertex[0];
        v[1] = vertex[1];
        v[2] = vertex[2];
        fwrite((void *) v, sizeof(float), 3, fp);
        fwrite(&abc, sizeof(unsigned short), 1, fp);
    }

    fclose(fp);
}




#endif  // FICSIONINIT_HH
