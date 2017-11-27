#ifndef FICSION_BONDPARTICLEFIELD3D_HDF5IO_HH
#define FICSION_BONDPARTICLEFIELD3D_HDF5IO_HH

#include "BondParticleField3DHdf5IO.h"
#include "ficsionInit.h"


namespace trombocit {

/* ******** WriteBondParticleField3D *********************************** */
template<typename T, template<typename U> class Descriptor>
WriteBondParticleField3D<T,Descriptor>::WriteBondParticleField3D (
        plint iter_, std::string identifier_,
        T dx_, T dt_) :  iter(iter_), identifier(identifier_), dx(dx_), dt(dt_) {};

template<typename T, template<typename U> class Descriptor>
void WriteBondParticleField3D<T,Descriptor>::processGenericBlocks (
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
    plint Nv = Np*2;
     /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
     hsize_t dimVertices[2]; dimVertices[0] = Np; dimVertices[1] = 3;

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     long int iterHDF5=iter;
     H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

     H5LTset_attribute_long (file_id, "/", "numberOfParticles", &Np, 1);
     H5LTset_attribute_long (file_id, "/", "numberOfVertices", &Nv, 1);

//     serializer.addValue<plint>(processor);
//     serializer.addValue<T>(r);
//     serializer.addValue<T>(bondTime);
//     serializer.addValues<T,3>(eij);
// 	serializeString(serializer, uid);
//
//     for (int var = 0; var < 2; ++var) {
//     	if (particles[var] != NULL) { // Update positions etc
//     		ImmersedCellParticle3D<T,Descriptor>* pv = castParticleToICP3D(particles[var]);
// 			serializer.addValues<T,3>(pv->getPosition());
// 			serializer.addValues<T,3>(pv->get_v());
// 			serializer.addValue<plint>(pv->getMpiProcessor());
// 			serializer.addValue<plint>(pv->get_cellId());
// 			serializer.addValue<plint>(pv->getVertexId());
//     	} else {
// 			serializer.addValues<T,3>(positions[var]);
// 			serializer.addValues<T,3>(velocities[var]);
// 			serializer.addValue<plint>(processors[var]);
// 			serializer.addValue<plint>(cellId[var]);
// 			serializer.addValue<plint>(vertexId[var]);
//     	}
//     }


     /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
     /* Particles */
     float * geometry = new float [3 * 2 *Np];
     int * connections= new int [2*Np];
     plint itr2=0, itr=0;
     hemo::Array<T,3> tmp;
     for (plint iP = 0; iP < Np; ++iP) {
    	 BondParticle3D<T,Descriptor>* bp = castParticle3DToBondParticle3D(particles[iP]);
    	 tmp = bp->getPositions(0);
    	 geometry[itr++] = tmp[0];
    	 geometry[itr++] = tmp[1];
    	 geometry[itr++] = tmp[2];
    	 tmp = bp->getPositions(1);
    	 geometry[itr++] = tmp[0];
    	 geometry[itr++] = tmp[1];
    	 geometry[itr++] = tmp[2];
    	 connections[itr2] = itr2; itr2++;
    	 connections[itr2] = itr2; itr2++;
     }
     hsize_t edgeDimVertices[2]; edgeDimVertices[0] = 2*Np; edgeDimVertices[1] = 3;
     H5LTmake_dataset_float(file_id, "position", 2, edgeDimVertices, geometry);
     H5LTmake_dataset_int(file_id, "connections", 1, edgeDimVertices, connections);
     delete [] geometry;
     delete [] connections;

     /* BondTime */
     int * bondTimes= new int [Nv];
     itr=0;
     for (plint iP = 0; iP < Np; ++iP) {
    	 bondTimes[itr++] = castParticle3DToBondParticle3D(particles[iP])->getBondTime();
    	 bondTimes[itr++] = castParticle3DToBondParticle3D(particles[iP])->getBondTime();
     }
     H5LTmake_dataset_int(file_id, "bondTime", 1, edgeDimVertices, bondTimes);
     delete [] bondTimes;

     /* Processor */
     itr=0;
     int * pIds= new int [Nv];
     for (plint iP = 0; iP < Np; ++iP) {
         BondParticle3D<T,Descriptor>* bp = castParticle3DToBondParticle3D(particles[iP]);
    	 pIds[itr++] = bp->getProcessors(0);
    	 pIds[itr++] = bp->getProcessors(1);
     }
     H5LTmake_dataset_int(file_id, "processor", 1, edgeDimVertices, pIds);
     delete [] pIds;

     /* CellIds */
     itr=0;
     pIds= new int [Nv];
     for (plint iP = 0; iP < Np; ++iP) {
         BondParticle3D<T,Descriptor>* bp = castParticle3DToBondParticle3D(particles[iP]);
    	 pIds[itr++] = bp->getCellIds(0);
    	 pIds[itr++] = bp->getCellIds(1);
     }
     H5LTmake_dataset_int(file_id, "cellId", 1, edgeDimVertices, pIds);
     delete [] pIds;

     H5Fclose(file_id);

}

template<typename T, template<typename U> class Descriptor>
WriteBondParticleField3D<T,Descriptor>* WriteBondParticleField3D<T,Descriptor>::clone() const {
    return new WriteBondParticleField3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void WriteBondParticleField3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT WriteBondParticleField3D<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor>
void writeBondParticleField3D_HDF5(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField, T dx, T dt, plint iter, std::string identifier)
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particleField);

    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new WriteBondParticleField3D<T,Descriptor> (iter, identifier, dx, dt),
            particleField.getBoundingBox(), particleArg );

}


}


#endif  // FICSION_BONDPARTICLEFIELD3D_HDF5IO_H

//  Get the number of processes.

