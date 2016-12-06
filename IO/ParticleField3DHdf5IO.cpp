#ifndef FICSION_PARTICLEFIELD3D_HDF5IO_HH
#define FICSION_PARTICLEFIELD3D_HDF5IO_HH

#include "ParticleField3DHdf5IO.h"
#include "ficsionInit.h"



/* ******** WriteParticleField3DInMultipleHDF5Files *********************************** */
template<typename T, template<typename U> class Descriptor>
WriteParticleField3DInMultipleHDF5Files<T,Descriptor>::WriteParticleField3DInMultipleHDF5Files (
        plint iter_, std::string identifier_,
        T dx_, T dt_) :  iter(iter_), identifier(identifier_), dx(dx_), dt(dt_) {};

template<typename T, template<typename U> class Descriptor>
void WriteParticleField3DInMultipleHDF5Files<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size() > 0 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    int p = global::mpi().getSize();
    int id = global::mpi().getRank();
    plint Nx = particleField.getNx();
    plint Ny = particleField.getNy();
    plint Nz = particleField.getNz();

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
//    particleField.findParticles(particleField.getBoundingBox(), particles);
    plint Np=particles.size();

     /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
     
     hsize_t dimVertices[2]; 
     dimVertices[0] = Np; dimVertices[1] = 3;

     hsize_t chunk[2];  
     chunk[0] = (Np/p) > 1 ? (Np/p) : 1; chunk[1] = 3;

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     long int iterHDF5=iter;
     H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

     H5LTset_attribute_long (file_id, "/", "numberOfParticles", &Np, 1);

     /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
     /* Positions */
     float * positions = new float [3 * Np];
     Array<T,3> tmp;
     plint itr=0;
     for (plint iP = 0; iP < Np; ++iP) {
    	 tmp = particles[iP]->getPosition();
    	 positions[itr++] = tmp[0];
    	 positions[itr++] = tmp[1];
    	 positions[itr++] = tmp[2];
     }
#ifdef NO_COMPRESSION
     H5LTmake_dataset_float(file_id, "position", 2, dimVertices, positions);
#else            
     {
    int sid = H5Screate_simple(2,dimVertices,NULL);
    int plist_id = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, chunk); 
    H5Pset_deflate(plist_id, 7);
    int did = H5Dcreate2(file_id,"position",H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
    H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,positions);
    H5Dclose(did);
    H5Sclose(sid);
     }
#endif
     delete [] positions;

     /* Tags */
     int * tags= new int [Np];
     for (plint iP = 0; iP < Np; ++iP) { tags[iP] = particles[iP]->getTag(); }
#ifdef NO_COMPRESSION
     H5LTmake_dataset_int(file_id, "tag", 1, dimVertices, tags);
#else            
     {
    int sid = H5Screate_simple(1,dimVertices,NULL);
    int plist_id = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 1, chunk); 
    H5Pset_deflate(plist_id, 7);
    int did = H5Dcreate2(file_id,"tag",H5T_NATIVE_INT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
    H5Dwrite(did,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,tags);
    H5Dclose(did);
    H5Sclose(sid);
     }
#endif
     delete [] tags;

     /* Tags */
     int * ids= new int [Np];
     for (plint iP = 0; iP < Np; ++iP) { ids[iP] = particles[iP]->getId(); }
#ifdef NO_COMPRESSION
     H5LTmake_dataset_int(file_id, "id", 1, dimVertices, ids);
#else            
     {
    int sid = H5Screate_simple(1,dimVertices,NULL);
    int plist_id = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 1, chunk); 
    H5Pset_deflate(plist_id, 7);
    int did = H5Dcreate2(file_id,"id",H5T_NATIVE_INT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
    H5Dwrite(did,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ids);
    H5Dclose(did);
    H5Sclose(sid);
     }
#endif
     delete [] ids;

     /* Processor */
     int * pIds= new int [Np];
     for (plint iP = 0; iP < Np; ++iP) { pIds[iP] = id; }
#ifdef NO_COMPRESSION
     H5LTmake_dataset_int(file_id, "processor", 1, dimVertices, pIds);
#else            
     {

    int sid = H5Screate_simple(1,dimVertices,NULL);
    int plist_id = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 1, chunk); 
    H5Pset_deflate(plist_id, 7);
    int did = H5Dcreate2(file_id,"processor",H5T_NATIVE_INT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
    H5Dwrite(did,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pIds);
    H5Dclose(did);
    H5Sclose(sid);
     }
#endif
    delete [] pIds;


     H5Fclose(file_id);

}

template<typename T, template<typename U> class Descriptor>
WriteParticleField3DInMultipleHDF5Files<T,Descriptor>* WriteParticleField3DInMultipleHDF5Files<T,Descriptor>::clone() const {
    return new WriteParticleField3DInMultipleHDF5Files<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void WriteParticleField3DInMultipleHDF5Files<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT WriteParticleField3DInMultipleHDF5Files<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor>
void writeParticleField3D_HDF5(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField, T dx, T dt, plint iter, std::string identifier)
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);

    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new WriteParticleField3DInMultipleHDF5Files<T,Descriptor> (iter, identifier, dx, dt),
            particleField.getBoundingBox(), particleArg );

}





#endif  // FICSION_PARTICLEFIELD3D_HDF5IO_HH

//  Get the number of processes.

