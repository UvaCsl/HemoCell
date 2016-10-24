#ifndef FICSION_HDF5IO_HH
#define FICSION_HDF5IO_HH

#include "hdf5IO.h"
#include "ficsionInit.h"
#include <algorithm>    // std::swap

/* ******** WriteInMultipleHDF5Files *********************************** */
template<typename T>
WriteInMultipleHDF5Files<T>::WriteInMultipleHDF5Files (
        std::vector<std::string> & hdf5ContainerNames_,
        std::vector<plint> & hdf5ContainerDimensions_,
        plint iter_, T dx_, T dt_,
        plint envelopeWidth_, bool invertXZ_for_XDMF_) :
            hdf5ContainerNames(hdf5ContainerNames_),
            hdf5ContainerDimensions(hdf5ContainerDimensions_),
            iter(iter_), dx(dx_), dt(dt_),
            envelopeWidth(envelopeWidth_), invertXZ_for_XDMF(invertXZ_for_XDMF_) {};

template<typename T>
void WriteInMultipleHDF5Files<T>::processGenericBlocks (
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
     hsize_t dim[4];
     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName("Fluid.",iter,8) + createFileName(".p.",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

     int Nx, Ny, Nz;
     Nx = domain.x1 - domain.x0 + 1;
     Ny = domain.y1 - domain.y0 + 1;
     Nz = domain.z1 - domain.z0 + 1;
     Dot3D relativePositionDot3D = blocks[0]->getLocation();
     if (invertXZ_for_XDMF)  {
    	 std::swap(Nx, Nz);
    	 std::swap(relativePositionDot3D.x, relativePositionDot3D.z);
     }
     long int relativePosition[] = {relativePositionDot3D.x + envelopeWidth, relativePositionDot3D.y + envelopeWidth, relativePositionDot3D.z + envelopeWidth};
     long int subdomainSize[] = {Nx+1, Ny+1, Nz+1};
//     long int domainSize[] = {blocks[0]->getNx(), blocks[0]->getNy(), blocks[0]->getNz()};

     H5LTset_attribute_long(file_id, "/", "relativePosition", relativePosition, 3);
     H5LTset_attribute_long(file_id, "/", "subdomainSize", subdomainSize, 3);
//     H5LTset_attribute_long(file_id, "/", "domainSize", domainSize, 3);

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     long int itrtHDF5=iter;
     long int XZ_inverted_for_XDMF=invertXZ_for_XDMF;
     H5LTset_attribute_long (file_id, "/", "iteration", &itrtHDF5, 1);
     H5LTset_attribute_long (file_id, "/", "XZ_inverted_for_XDMF", &XZ_inverted_for_XDMF, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

     dim[0] = Nx; dim[1] = Ny; dim[2] = Nz;
     dim[3] = 3;

    for (pluint i = 0; i < hdf5ContainerNames.size(); ++i) {
        int nDim = hdf5ContainerDimensions[i];
        if (nDim == 1) {
            ScalarField3D<T>& scalarF =
                    *dynamic_cast<ScalarField3D<T>*>(blocks[i]);
            int Np = Nx*Ny*Nz;
            float * matrixScalar = new float [Np];
            int iter = 0;
            if (not invertXZ_for_XDMF)  {
				for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
					for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
						for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
							matrixScalar[iter++] = scalarF.get(iX,iY,iZ);
						}
					}
				}
            } else {
				for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
					for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
						for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
							matrixScalar[iter++] = scalarF.get(iX,iY,iZ);
						}
					}
				}
            }
#ifdef NO_COMPRESSION
            H5LTmake_dataset_float(file_id, hdf5ContainerNames[i].c_str(), 3, dim, matrixScalar);
#else            
            int sid = H5Screate_simple(3,dim,NULL);
            int plist_id = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 3, dim); 
            H5Pset_deflate(plist_id, 6);
            int did = H5Dcreate2(file_id,hdf5ContainerNames[i].c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
            H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,matrixScalar);
            H5Dclose(did);
            H5Sclose(sid);
#endif
            delete [] matrixScalar;

        } else {
            TensorField3D<T,3>& tensorF =
                    *dynamic_cast<TensorField3D<T,3>*>(blocks[i]);
            int Np = Nx*Ny*Nz*3;
            float * matrixTensor = new float [Np];
            int iter = 0;
            if (not invertXZ_for_XDMF)  {
				for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
					for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
						for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
							Array<T,3> vector = tensorF.get(iX,iY,iZ);
							matrixTensor[iter++] = vector[0];
							matrixTensor[iter++] = vector[1];
							matrixTensor[iter++] = vector[2];
						}
					}
				}
            } else {
				for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
					for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
						for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
							Array<T,3> vector = tensorF.get(iX,iY,iZ);
							matrixTensor[iter++] = vector[0];
							matrixTensor[iter++] = vector[1];
							matrixTensor[iter++] = vector[2];
						}
					}
				}
			}
#ifdef NO_COMPRESSION
            H5LTmake_dataset_float(file_id, hdf5ContainerNames[i].c_str(), 4, dim, matrixTensor);
#else
            int sid = H5Screate_simple(4,dim,NULL);
            int plist_id = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 4, dim); 
            H5Pset_deflate(plist_id, 6);
            int did = H5Dcreate2(file_id,hdf5ContainerNames[i].c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
            H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,matrixTensor);
            H5Dclose(did);
            H5Sclose(sid);
#endif
            delete [] matrixTensor;
        }
    }
    H5Fclose(file_id);
}

template<typename T>
WriteInMultipleHDF5Files<T>* WriteInMultipleHDF5Files<T>::clone() const {
    return new WriteInMultipleHDF5Files<T>(*this);
}

template<typename T>
void WriteInMultipleHDF5Files<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < hdf5ContainerNames.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T>
BlockDomain::DomainT WriteInMultipleHDF5Files<T>::appliesTo () const {
    return BlockDomain::bulk;
}




template<typename T, template<typename U> class Descriptor>
void writeHDF5(MultiBlockLattice3D<T, Descriptor>& lattice,
              T dx, T dt, plint iter, bool invertXZ_for_XDMF)
{
    plint envelopeWidth = lattice.getMultiBlockManagement().getEnvelopeWidth();
//    cout << "env " << envelopeWidth << std::endl;
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
            new WriteInMultipleHDF5Files<T> (hdf5ContainerNames, hdf5ContainerDimensions, iter, dx, dt, envelopeWidth, invertXZ_for_XDMF),
            lattice.getBoundingBox(), hdf5ContainerArg );

}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

