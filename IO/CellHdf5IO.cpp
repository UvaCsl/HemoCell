#ifndef FICSION_CELL_HDF5IO_CPP
#define FICSION_CELL_HDF5IO_CPP

#include "CellHdf5IO.h"
#include "ficsionInit.h"




/* ******** WriteCellField3DInMultipleHDF5Files *********************************** */
template<typename T, template<typename U> class Descriptor>
WriteCell3DInMultipleHDF5Files<T,Descriptor>::WriteCell3DInMultipleHDF5Files (
        CellField3D<T, Descriptor>& cellField3D_,
        plint iter_, std::string identifier_,
        T dx_, T dt_) :
        cellField3D(cellField3D_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_) {};

template<typename T, template<typename U> class Descriptor>
void WriteCell3DInMultipleHDF5Files<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size() > 0 );
    int p = global::mpi().getSize();
    int id = global::mpi().getRank();
    // plint Nx = cellField3D.getParticleArg()[0]->getNx();
    // plint Ny = cellField3D.getParticleArg()[0]->getNy();
    // plint Nz = cellField3D.getParticleArg()[0]->getNz();

    /************************************************************/
   /**            Fill triangle and particle lists            **/
  /************************************************************/

     std::map<plint, Cell3D<T,Descriptor>* > cellIdToCell3D = cellField3D.getCellIdToCell3D();
     long int Nc = cellField3D.getNumberOfCells();
     std::vector<long int> cellIds;
     for (unsigned int i = 0; i < cellField3D.getCellIds().size(); ++i)
     {
        cellIds.push_back(cellField3D.getCellIds()[i]);
     }

     if (Nc == 0) { return; }
     long int firstCell = cellIds[0];
     std::vector<plint> const& scalarCcrIds = cellIdToCell3D[firstCell]->getScalarCcrIds();
     std::vector<plint> const& vectorCcrIds = cellIdToCell3D[firstCell]->getVectorCcrIds();
     std::vector<plint> const& tensorCcrIds = cellIdToCell3D[firstCell]->getTensorCcrIds();

     /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

     hsize_t dimVertices[2]; 
     dimVertices[0] = Nc; dimVertices[1] = 3;

    hsize_t chunk[2];  
    chunk[0] = (Nc/p) > 1 ? (Nc/p) : 1 ; chunk[1] = 3;

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     long int iterHDF5=iter;
     H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

     H5LTset_attribute_long (file_id, "/", "numberOfCells", &Nc, 1);
     H5LTmake_dataset_long(file_id, "cellIds", 1, dimVertices, &cellIds[0]);

     /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
     /*  Take care of Vectors    */
     plint vN = vectorCcrIds.size();
     double * matrixTensor = new double [3 * Nc];
     Array<T,3> array;
     for (plint ivN = 0; ivN < vN; ++ivN) {
         plint ccrId=vectorCcrIds[ivN];
         plint itr=0;
         for (plint iC = 0; iC < Nc; ++iC) {
            plint cellId = cellIds[iC];
            array = cellIdToCell3D[cellId]->get3D(ccrId);
            // TODO: Change in XDMF file.
            matrixTensor[itr++] = array[0];
            matrixTensor[itr++] = array[1];
            matrixTensor[itr++] = array[2];
         }
#ifdef NO_COMPRESSION
         H5LTmake_dataset_double(file_id, ccrNames[ccrId].c_str(), 2, dimVertices, matrixTensor);
#else               
        int sid = H5Screate_simple(2,dimVertices,NULL);
        int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,ccrNames[ccrId].c_str(),H5T_NATIVE_DOUBLE,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,matrixTensor);
        H5Dclose(did);
        H5Sclose(sid);
#endif
     }
     delete [] matrixTensor;

     /*  Take care of Scalars    */
     plint sN = scalarCcrIds.size();
     double * scalarTensor = new double [ Nc ];
     for (plint isN = 0; isN < sN; ++isN) {
         plint ccrId=scalarCcrIds[isN];
         plint itr=0;
         for (plint iC = 0; iC < Nc; ++iC) {
            plint cellId = cellIds[iC];
             scalarTensor[itr++] = cellIdToCell3D[cellId]->get1D(ccrId);
         }
#ifdef NO_COMPRESSION
         H5LTmake_dataset_double(file_id, ccrNames[ccrId].c_str(), 1, dimVertices, scalarTensor);
#else
            
        int sid = H5Screate_simple(1,dimVertices,NULL);
        int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,ccrNames[ccrId].c_str(),H5T_NATIVE_DOUBLE,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,scalarTensor);
        H5Dclose(did);
        H5Sclose(sid);
#endif
     }
     delete [] scalarTensor;

     /*  Take care of Tensors    */
     plint tN = tensorCcrIds.size();
     Array<T,3> vector;
     for (plint itN = 0; itN < tN; ++itN) {
         plint ccrId=tensorCcrIds[itN];
         plint numT = getReductionDimension(ccrId);
         double * tensorTensor = new double [numT * Nc];
         plint itr=0;
         for (plint iC = 0; iC < Nc; ++iC) {
            plint cellId = cellIds[iC];
             std::vector<T> const& vector = cellIdToCell3D[cellId]->getND(ccrId);
             for (plint iT = 0; iT < numT; ++iT) {
                 tensorTensor[itr++] = vector[iT];
             }
         }
         hsize_t dimTensor[2]; dimTensor[0] = Nc; dimTensor[1] = numT;
#ifdef NO_COMPRESSION
        H5LTmake_dataset_double(file_id, ccrNames[ccrId].c_str(), 2, dimTensor, tensorTensor);
#else            
        int sid = H5Screate_simple(2,dimTensor,NULL);
        int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,ccrNames[ccrId].c_str(),H5T_NATIVE_DOUBLE,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tensorTensor);
        H5Dclose(did);
        H5Sclose(sid);
#endif
         delete [] tensorTensor;
     }
     H5Fclose(file_id);

}

template<typename T, template<typename U> class Descriptor>
WriteCell3DInMultipleHDF5Files<T,Descriptor>* WriteCell3DInMultipleHDF5Files<T,Descriptor>::clone() const {
    return new WriteCell3DInMultipleHDF5Files<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void WriteCell3DInMultipleHDF5Files<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT WriteCell3DInMultipleHDF5Files<T,Descriptor>::appliesTo () const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor>
void writeCell3D_HDF5(CellField3D<T, Descriptor>& cellField3D, T dx, T dt, plint iter, std::string preString)
{
	std::string identifier = preString + cellField3D.getIdentifier() + "_Cell3D";
    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            new WriteCell3DInMultipleHDF5Files<T,Descriptor> (cellField3D, iter, identifier, dx, dt),
            cellField3D.getBoundingBox(), cellField3D.getParticleArg() );

}





#endif  // FICSION_CELL_HDF5IO_CPP

//  Get the number of processes.

