#ifndef FICSION_HDF5IO_HH
#define FICSION_HDF5IO_HH

#include "hdf5IO.h"
#include "ficsionInit.h"


/* ******** WriteInMultipleHDF5Files *********************************** */
template<typename T, template<typename U> class Descriptor>
WriteInMultipleHDF5Files<T,Descriptor>::WriteInMultipleHDF5Files (
        std::vector<std::string> & hdf5ContainerNames_,
        std::vector<plint> & hdf5ContainerDimensions_,
        plint iter_) :
            hdf5ContainerNames(hdf5ContainerNames_),
            hdf5ContainerDimensions(hdf5ContainerDimensions_),
            iter(iter_) {};

template<typename T, template<typename U> class Descriptor>
void WriteInMultipleHDF5Files<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size() == hdf5ContainerNames.size() );

#ifdef PLB_MPI_PARALLEL
    //  Get the number of processes.
    int  p = MPI::COMM_WORLD.Get_size();
    //  Get the individual process ID.
     int id = MPI::COMM_WORLD.Get_rank();
#else
     int p=1;
     int id = 0;
#endif

     std::string fileName = global::directories().getOutputDir() + createFileName("Fluid.p",id,3) + createFileName(".",iter,8) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t   dim[4];
//    dims[0] = DIM0;     dims[1] = DIM1;    dims[2] = DIM2;    dims[3] = 3;
    for (int i = 0; i < hdf5ContainerNames.size(); ++i) {
        if (hdf5ContainerDimensions[i] == 1) {
            ScalarField3D<T>& scalarF =
                    *dynamic_cast<ScalarField3D<T>*>(blocks[i]);
    //        H5LTmake_dataset_double(file_id, hdf5ContainerNames[i].c_str(), 3, dim, matrixScalar);
        } else {
            TensorField3D<T,3>& tensorF =
                    *dynamic_cast<TensorField3D<T,3>*>(blocks[i]);
    //        H5LTmake_dataset_double(file_id, hdf5ContainerNames[i].c_str(), 4, dim, matrixTensor);
        }
    }
    H5Fclose(file_id);
}

template<typename T, template<typename U> class Descriptor>
WriteInMultipleHDF5Files<T,Descriptor>* WriteInMultipleHDF5Files<T,Descriptor>::clone() const {
    return new WriteInMultipleHDF5Files<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void WriteInMultipleHDF5Files<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < hdf5ContainerNames.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT WriteInMultipleHDF5Files<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}




template<typename T, template<typename U> class Descriptor>
void writeHDF5(MultiBlockLattice3D<T, Descriptor>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();


    MultiTensorField3D<T,3> vel = *computeVelocity(lattice);
    MultiTensorField3D<T,3> vorticity = *computeVorticity(vel);
    MultiTensorField3D<T,3> force(lattice);
    applyProcessingFunctional(new GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,3>(
            Descriptor<T>::ExternalField::forceBeginsAt), lattice.getBoundingBox(), lattice, force);

    MultiScalarField3D<T> density = *computeDensity(lattice);
    MultiScalarField3D<T> velNorm = *computeNorm(vel, vel.getBoundingBox());
    MultiScalarField3D<T> forceNorm = *computeNorm(force, force.getBoundingBox());

    std::vector<MultiBlock3D*> hdf5ContainerArg;
    std::vector<std::string> hdf5ContainerNames;
    std::vector<plint> hdf5ContainerDimensions;
    hdf5ContainerArg.push_back(&density);   hdf5ContainerNames.push_back("density");    hdf5ContainerDimensions.push_back(1);
    hdf5ContainerArg.push_back(&velNorm);   hdf5ContainerNames.push_back("velNorm");    hdf5ContainerDimensions.push_back(1);
    hdf5ContainerArg.push_back(&forceNorm); hdf5ContainerNames.push_back("forceNorm");  hdf5ContainerDimensions.push_back(1);
    hdf5ContainerArg.push_back(&vel);       hdf5ContainerNames.push_back("velocity");   hdf5ContainerDimensions.push_back(3);
    hdf5ContainerArg.push_back(&force);     hdf5ContainerNames.push_back("force");      hdf5ContainerDimensions.push_back(3);
    hdf5ContainerArg.push_back(&vorticity); hdf5ContainerNames.push_back("vorticity");  hdf5ContainerDimensions.push_back(3);


    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new WriteInMultipleHDF5Files<T,Descriptor> (hdf5ContainerNames, hdf5ContainerDimensions, iter),
            lattice.getBoundingBox(), hdf5ContainerArg );

}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

