#ifndef FICSION_HDF5IO_HH
#define FICSION_HDF5IO_HH

#include "hdf5IO.h"
#include "ficsionInit.h"


void writeCell3D_data(hid_t & file_id) {

}


/* ******** WriteCellField3DInMultipleHDF5Files *********************************** */
template<typename T, template<typename U> class Descriptor>
WriteCellField3DInMultipleHDF5Files<T,Descriptor>::WriteCellField3DInMultipleHDF5Files (
        CellField3D<T, Descriptor>& cellField3D_,
        plint iter_, std::string identifier_,
        T dx_, T dt_) :
        cellField3D(cellField3D_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_) {};

template<typename T, template<typename U> class Descriptor>
void WriteCellField3DInMultipleHDF5Files<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size() > 0 );

#ifdef PLB_MPI_PARALLEL
    //  Get the number of processes.
    int  p = MPI::COMM_WORLD.Get_size();
    //  Get the individual process ID.
     int id = MPI::COMM_WORLD.Get_rank();
#else
     int p=1;
     int id = 0;
#endif
    /************************************************************/
   /**            Fill triangle and particle lists            **/
  /************************************************************/

     std::map<plint, Cell3D<T,Descriptor> > cellIdToCell3D = cellField3D.getCellIdToCell3D();
     std::vector< plint > triangles;
     std::vector<Particle3D<T,Descriptor>* > particles;
     plint sumLocalVertices=0;

     typename std::map<plint, Cell3D<T,Descriptor> >::iterator iter;
     for (iter  = cellIdToCell3D.begin(); iter != cellIdToCell3D.end(); ++iter) {
         Cell3D<T,Descriptor> & cell3d = (iter->second);
         plint cellId = cell3d.get_cellId();
         std::vector<plint> const& cellVertices = cell3d.getVertices();
         std::vector<plint> const& cellTriangles = cell3d.getTriangles();
         for (std::vector<plint>::const_iterator iVP = cellVertices.begin(); iVP != cellVertices.end(); ++iVP)
         {
             particles.push_back( cell3D.getParticle3D(*iVP) );
         }
         std::map<plint, plint> iv = cell3d.getInvertVertices();
         for (pluint iT=0; iT < cellTriangles.size(); iT++) {
             plint iTriangle=cellTriangles[iT];
             plint t0 = iv[cell3d.getVertexId(iTriangle, 0)]+sumLocalVertices;
             plint t1 = iv[cell3d.getVertexId(iTriangle, 1)]+sumLocalVertices;
             plint t2 = iv[cell3d.getVertexId(iTriangle, 2)]+sumLocalVertices;
             triangles.push_back(t0);
             triangles.push_back(t1);
             triangles.push_back(t2);
         }
         sumLocalVertices += cellVertices.size();
     }

     pluint Np = particles.size();
     pluint Nt = triangles.size();

     /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     hsize_t dim[4];
     std::string fileName = global::directories().getOutputDir() + createFileName(identifier.c_str(),iter,8) + createFileName("_p",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

     hsize_t dimTriangles[2]; dimTriangles[0] = Nt/3; dimTriangles[1] = 3;
     hsize_t dimVertices[2]; dimVertices[0] = Np; dimVertices[1] = 3;
     H5LTmake_dataset_double(file_id, "triangles", 2, dim, &triangles[0]);

     /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/


     /*  Take care of Vectors    */
     ImmersedCellParticle3D<T,Descriptor> icParticle = castParticleToICP3D(particles[0]);
     plint vN = icParticle->getVectorsNumber();
     double * matrixTensor = new double [3 * Np];
     for (int ivN = 0; ivN < vN; ++ivN) {
         plint iter=0;
         for (int iP = 0; iP < Np; ++iP) {
            Array<T,3> vector =  castParticleToICP3D(particles[iP])->getVector(ivN);
            matrixTensor[iter++] = vector[0];
            matrixTensor[iter++] = vector[1];
            matrixTensor[iter++] = vector[2];
         }
         H5LTmake_dataset_double(file_id, icParticle->getVectorName(ivN).c_str(), 2, dimVertices, matrixTensor);
     }
     delete [] matrixTensor;

     /*  Take care of Scalars    */
     pluint Np = particles.size();
     plint sN = icParticle->getScalarsNumber();
     double * scalarTensor = new double [Np];
     for (int isN = 0; isN < sN; ++isN) {
         plint iter=0;
         for (int iP = 0; iP < Np; ++iP) {
            scalarTensor[iter++] = castParticleToICP3D(particles[iP])->getScalar(isN);
         }
         H5LTmake_dataset_double(file_id, icParticle->getScalarName(isN).c_str(), 1, dimVertices, scalarTensor);
     }
     delete [] scalarTensor;


     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     H5LTset_attribute_long (file_id, "/", "iteration", &iter, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfParticles", &Np, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfTriangles", &Nt, 1);
     H5Fclose(file_id);


}

template<typename T, template<typename U> class Descriptor>
WriteCellField3DInMultipleHDF5Files<T,Descriptor>* WriteCellField3DInMultipleHDF5Files<T,Descriptor>::clone() const {
    return new WriteCellField3DInMultipleHDF5Files<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void WriteCellField3DInMultipleHDF5Files<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT WriteCellField3DInMultipleHDF5Files<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor>
void writeCellField3D_HDF5(CellField3D<T, Descriptor>& cellField3D,
              IncomprFlowParam<T> const& parameters, plint iter, std::string identifier)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();


    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new WriteCellField3DInMultipleHDF5Files<T,Descriptor> (cellField3D, iter, identifier, dx, dt),
            cellField3D.getBoundingBox(), cellField3D.getParticleArg() );

}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

