#ifndef FICSION_HDF5IO_H
#define FICSION_HDF5IO_H

#include "hdfIO.h"




template<class BlockLatticeT>
void writeHDF5(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    //  Get the number of processes.
    int  p = MPI::COMM_WORLD.Get_size();
    //  Get the individual process ID.
     int id = MPI::COMM_WORLD.Get_rank();


     hsize_t   dimScalar[3];
     hsize_t   dimVector[4];
     dims[0] = DIM0;
     dims[1] = DIM1;
     dims[2] = DIM2;
     dims[3] = DIM3;


    MultiTensorField3D<T,3> force(lattice);
    applyProcessingFunctional(new GetTensorFieldFromExternalVectorFunctional3D<T,DESCRIPTOR,3>(
        DESCRIPTOR<T>::ExternalField::forceBeginsAt), lattice.getBoundingBox(), lattice, force);


    auto_ptr<MultiTensorField3D<T,3> > vel(computeVelocity(lattice));
    auto_ptr<MultiTensorField3D<T,3> > vorticity(computeVorticity(*vel));
    auto_ptr<MultiScalarField3D<T> > density(computeDensity(lattice));

    auto_ptr<MultiScalarField3D<T> > velNorm(computeNorm(vel, vel.getBoundingBox()));
    auto_ptr<MultiScalarField3D<T> > forceNorm(computeNorm(force, force.getBoundingBox()));


//    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 8), dx);
//    vtkOut.writeData<T>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
//    vtkOut.writeData<T>(*computeDensity(lattice), "density", 1.0);
//    vtkOut.writeData<3,T>(*computeVelocity(lattice), "velocity", dx/dt);
//    vtkOut.writeData<3,T>(force, "force",  (T) 1.0);
//    vtkOut.writeData<T>(*computeNorm(force, force.getBoundingBox()), "forceNorm LU",  1.0);

    hid_t     file_id;
    file_id = H5Fcreate("xdmf2d_co.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5LTmake_dataset_double(file_id,"velocityNorm", 3, dimScalar, velNorm);
    H5LTmake_dataset_double(file_id,"density", 3, dimScalar, density);
    H5LTmake_dataset_double(file_id,"forceNormLU", 3, dimScalar, forceNorm);

    H5LTmake_dataset_double(file_id,"velocity", 4, dimVector, vel);
    H5LTmake_dataset_double(file_id,"force", 4, dimVector, force);
    H5LTmake_dataset_double(file_id,"vorticity", 4, dimVector, vorticity);


    H5Fclose(file_id);
//     vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

