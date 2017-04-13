#ifndef FICSION_PARTICLE_HDF5IO_HH
#define FICSION_PARTICLE_HDF5IO_HH

#include "ParticleHdf5IO.h"

/* ******** WriteCellField3DInMultipleHDF5Files *********************************** */
WriteCellField3DInMultipleHDF5Files::WriteCellField3DInMultipleHDF5Files (
        HemoCellField & cellField3D_,
        plint iter_, std::string identifier_,
        double dx_, double dt_, int ctype_) :
        cellField3D(cellField3D_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_), ctype(ctype_) {};

void WriteCellField3DInMultipleHDF5Files::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{

    PLB_PRECONDITION( blocks.size() > 0 );
    int id = global::mpi().getRank();
    long int size = global::mpi().getSize();

      HemoParticleField3D& particleField =
        *dynamic_cast<HemoParticleField3D*>(blocks[0]);
   
    /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.", particleField.atomicBlockId,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     long int iterHDF5=iter;
     H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);
     H5LTset_attribute_long (file_id, "/", "numberOfProcessors", &size, 1);
          
     hsize_t dimVertices[2];
     hsize_t chunk[2];  
 
     vector<vector<double>> positions;
     
    /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
         
    for (pluint i = 0; i < cellField3D.desiredOutputVariables.size(); i++) {
        vector<vector<double>> * output = new vector<vector<double>>();
        std::string vectorname;
        particleField.passthroughpass(cellField3D.desiredOutputVariables[i],domain,*output,cellField3D.ctype,vectorname);
        dimVertices[0] = (*output).size();
        dimVertices[1] = dimVertices[0] == 0 ? 0 :(*output)[0].size();
        chunk[0] = 1000 < (*output).size() ? 1000 : (*output).size();
        chunk[1] = dimVertices[1];
        chunk[0] = chunk[0] > 1 ? chunk[0] : 1;
        chunk[1] = chunk[1] > 1 ? chunk[1] : 1;

        double* output_formatted = new double[dimVertices[0] * dimVertices[1]];
        
        int fmt_cnt = 0;
        for (pluint x=0; x< dimVertices[0];x++) {
            for (pluint y=0; y < dimVertices[1] ; y++) {
                output_formatted[fmt_cnt] = (*output)[x][y];
                fmt_cnt++;
            }
        }
        int sid = H5Screate_simple(2,dimVertices,NULL);
        int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,vectorname.c_str(),H5T_NATIVE_DOUBLE,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,output_formatted);
        H5Dclose(did);
        H5Sclose(sid);
            
        if (i == 0) {
            long int nP = (*output).size();
            H5LTset_attribute_long (file_id, "/", "numberOfParticles", &nP, 1);
        }
        if (cellField3D.desiredOutputVariables[i] == OUTPUT_POSITION) {
            positions = (*output);
        }
        delete output;
        delete[] output_formatted;
    }
     
    if (cellField3D.outputTriangles) { //Treat triangles seperately because of double/int issues
        vector<vector<plint>> * output = new vector<vector<plint>>();
        std::string vectorname;
        particleField.outputTriangles(domain,*output,positions, cellField3D.ctype,vectorname);
        
        dimVertices[0] = output->size();
        dimVertices[1] = dimVertices[0] == 0 ? 0 :(*output)[0].size();
        chunk[0] = 1000 < output->size() ? 1000 : output->size();
        chunk[1] = dimVertices[1];
        chunk[0] = chunk[0] > 1 ? chunk[0] : 1;
        chunk[1] = chunk[1] > 1 ? chunk[1] : 1;
           
        int* output_formatted = new int[dimVertices[0] * dimVertices[1]];
        int fmt_cnt = 0;
        for (pluint x=0; x< dimVertices[0];x++) {
            for (pluint y=0; y < dimVertices[1] ; y++) {
                output_formatted[fmt_cnt] = (*output)[x][y];
                fmt_cnt++;
            }
        }
        
        int sid = H5Screate_simple(2,dimVertices,NULL);
        int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,vectorname.c_str(),H5T_NATIVE_INT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,output_formatted);
        H5Dclose(did);
        H5Sclose(sid);
        
        long int nT = output->size();
        H5LTset_attribute_long (file_id, "/", "numberOfTriangles", &nT, 1);
        delete output;
        delete[] output_formatted;
     }

     H5Fclose(file_id);

}

WriteCellField3DInMultipleHDF5Files* WriteCellField3DInMultipleHDF5Files::clone() const {
    return new WriteCellField3DInMultipleHDF5Files(*this);
}

void WriteCellField3DInMultipleHDF5Files::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    for (pluint i = 0; i < modified.size(); ++i) {
        modified[i] = modif::nothing;
    }
}

BlockDomain::DomainT WriteCellField3DInMultipleHDF5Files::appliesTo () const {
    return BlockDomain::bulk;
}



void writeCellField3D_HDF5(CellFields3D& cellFields, double dx, double dt, plint iter, std::string preString)
{
    for (pluint i = 0; i < cellFields.size(); i++) {
	std::string identifier = preString + cellFields[i]->getIdentifier();
        WriteCellField3DInMultipleHDF5Files * bprf = new WriteCellField3DInMultipleHDF5Files(*cellFields[i], iter, identifier, dx, dt, i);
        vector<MultiBlock3D*> wrapper;
        wrapper.push_back(cellFields[i]->getParticleArg());
        applyProcessingFunctional (bprf,cellFields[i]->getParticleField3D()->getBoundingBox(), wrapper );
    }
}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

