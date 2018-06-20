/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "ParticleHdf5IO.h"

#include <hdf5.h>
#include <hdf5_hl.h>

/* ******** WriteCellField3DInMultipleHDF5Files *********************************** */
WriteCellField3DInMultipleHDF5Files::WriteCellField3DInMultipleHDF5Files (
        HemoCellField & cellField3D_,
        plint iter_, std::string identifier_,
        T dx_, T dt_, int ctype_) :
        cellField3D(cellField3D_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_), ctype(ctype_) {};

void WriteCellField3DInMultipleHDF5Files::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{

    PLB_PRECONDITION( blocks.size() > 0 );
    int id = global::mpi().getRank();
    long int size = global::mpi().getSize();

      HemoCellParticleField& particleField =
        *dynamic_cast<HemoCellParticleField*>(blocks[0]);
   
      if (cellField3D.desiredOutputVariables.size() == 0) {
          return; //No output desired, no problem
      }
    /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/
     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + to_string(iter) + '/' + identifier + "."  + to_string(iter) + ".p." + to_string(particleField.atomicBlockId) + ".h5";
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
 
     vector<vector<T>> positions;
     
    /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
         
    for (pluint i = 0; i < cellField3D.desiredOutputVariables.size(); i++) {
        vector<vector<T>> * output = new vector<vector<T>>();
        std::string vectorname;
        particleField.passthroughpass(cellField3D.desiredOutputVariables[i],domain,*output,cellField3D.ctype,vectorname);
        dimVertices[0] = (*output).size();
        dimVertices[1] = dimVertices[0] == 0 ? 0 :(*output)[0].size();
        chunk[0] = 1000 < (*output).size() ? 1000 : (*output).size();
        chunk[1] = dimVertices[1];
        chunk[0] = chunk[0] > 1 ? chunk[0] : 1;
        chunk[1] = chunk[1] > 1 ? chunk[1] : 1;

        float* output_formatted = new float[dimVertices[0] * dimVertices[1]];
        
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
        int did = H5Dcreate2(file_id,vectorname.c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,output_formatted);
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
     
    if (cellField3D.outputTriangles) { //Treat triangles seperately because of T/int issues
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

     if (cellField3D.outputLines) { //Treat triangles seperately because of T/int issues
        vector<vector<plint>> * output = new vector<vector<plint>>();
        std::string vectorname;
        particleField.outputLines(domain,*output,positions, cellField3D.ctype,vectorname);
        
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
        H5LTset_attribute_long (file_id, "/", "numberOfLines", &nT, 1);
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



void writeCellField3D_HDF5(HemoCellFields& cellFields, T dx, T dt, plint iter, std::string preString)
{
    for (pluint i = 0; i < cellFields.size(); i++) {
	std::string identifier = preString + cellFields[i]->getIdentifier();
        WriteCellField3DInMultipleHDF5Files * bprf = new WriteCellField3DInMultipleHDF5Files(*cellFields[i], iter, identifier, dx, dt, i);
        vector<MultiBlock3D*> wrapper;
        wrapper.push_back(cellFields[i]->getParticleArg());
        applyProcessingFunctional (bprf,cellFields[i]->getParticleField3D()->getBoundingBox(), wrapper );
    }
}