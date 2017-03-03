#ifndef FICSION_PARTICLE_HDF5IO_HH
#define FICSION_PARTICLE_HDF5IO_HH

#include "ParticleHdf5IO.h"
#include "ficsionInit.h"



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

      HemoParticleField3D<double,DESCRIPTOR>& particleField =
        *dynamic_cast<HemoParticleField3D<double,DESCRIPTOR>*>(blocks[0]);
    /************************************************************/
   /**            Fill triangle and particle lists            **/
  /************************************************************/
/*
     std::map<plint, Cell3D<double,DESCRIPTOR>* > cellIdToCell3D = cellField3D.getCellIdToCell3D();
     std::vector< long int > triangles;
     std::vector<SurfaceParticle3D<double,DESCRIPTOR>* > particles;
     plint sumLocalVertices=0;
     plint numCells=cellIdToCell3D.size();
     long int NpBulk = 0;

     std::map<plint, Array<double,3> > correctPBPosition;
     typename std::map<plint, Cell3D<double,DESCRIPTOR>* >::iterator itrtr;
     for (itrtr  = cellIdToCell3D.begin(); itrtr != cellIdToCell3D.end(); ++itrtr) {
       Cell3D<double,DESCRIPTOR> * cell3d = (itrtr->second);
      //   Array<double,3> cellPosition = cell3d->getPosition();
      //   NpBulk += cell3d->getNumVertices_LocalBulk();
      //   correctPBPosition[itrtr->first].resetToZero();

      //   if (cellPosition[0] > Nx) { correctPBPosition[itrtr->first][0] = -int(cellPosition[0]/Nx)*Nx;}
      //   if (cellPosition[0] <  0) { correctPBPosition[itrtr->first][0] =  int(cellPosition[0]/Nx)*Nx;}
      //   if (cellPosition[1] > Ny) { correctPBPosition[itrtr->first][1] = -int(cellPosition[1]/Ny)*Ny;}
      //   if (cellPosition[1] <  0) { correctPBPosition[itrtr->first][1] =  int(cellPosition[1]/Ny)*Ny;}
      //   if (cellPosition[2] > Nz) { correctPBPosition[itrtr->first][2] = -int(cellPosition[2]/Nz)*Nz;}
      //   if (cellPosition[2] <  0) { correctPBPosition[itrtr->first][2] =  int(cellPosition[2]/Nz)*Nz;}

        std::vector<plint> const& cellVertices = cell3d->getVertices();
        std::vector<plint> const& cellTriangles = cell3d->getTriangles();
        for (std::vector<plint>::const_iterator iVP = cellVertices.begin(); iVP != cellVertices.end(); ++iVP)
        {
             particles.push_back( dynamic_cast<SurfaceParticle3D<double,DESCRIPTOR>*>(cell3d->getParticle3D(*iVP)) );
         }
         std::map<plint, plint> iv = cell3d->getInvertVertices();
         for (pluint iT=0; iT < cellTriangles.size(); iT++) {
             plint iTriangle=cellTriangles[iT];
             plint v0 = cell3d->getVertexId(iTriangle, 0);
             plint v1 = cell3d->getVertexId(iTriangle, 1);
             plint v2 = cell3d->getVertexId(iTriangle, 2);
             plint t0 = iv[v0]+sumLocalVertices;
             plint t1 = iv[v1]+sumLocalVertices;
             plint t2 = iv[v2]+sumLocalVertices;
             triangles.push_back(t0);
             triangles.push_back(t1);
             triangles.push_back(t2);
         }
         sumLocalVertices += cellVertices.size();
     }

     long int Np = particles.size();
     long int Nt = triangles.size() / 3;
*/
  
   
    /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.",id,3) + ".h5";
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
        dimVertices[1] = (*output)[0].size();
        chunk[0] = 1000 < (*output).size() ? 1000 : (*output).size();
        chunk[1] = dimVertices[1];
        
        double output_formatted[dimVertices[0] *dimVertices[1]];
        
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
        free(output);
    }
     
    if (cellField3D.outputTriangles) { //Treat triangles seperately because of double/int issues
        vector<vector<plint>> * output = new vector<vector<plint>>();
        std::string vectorname;
        particleField.outputTriangles(domain,*output,positions, cellField3D.ctype,vectorname);
        
        dimVertices[0] = output->size();
        dimVertices[1] = (*output)[0].size();
        chunk[0] = 1000 < output->size() ? 1000 : output->size();
        chunk[1] = dimVertices[1];
        
        int output_formatted[dimVertices[0] * dimVertices[1]];
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
        free(output);

     }
#if 0
     /*  Take care of Vectors    */
     if (numCells > 0 and particles.size() > 0) {
         SurfaceParticle3D<double,DESCRIPTOR> * icParticle = particles[0];
         plint vN = icParticle->getVectorsNumber();
         float * matrixTensor = new float [3 * Np];

         Array<double,3> vector;
         for (plint ivN = 0; ivN < vN; ++ivN) {
             plint itr=0;
             for (plint iP = 0; iP < Np; ++iP) {
                icParticle = particles[iP];
                icParticle->getVector(ivN, vector);
                // If is pbcPosition
                if (ivN == 0) { vector = vector + correctPBPosition[icParticle->get_cellId()] ; }
                // TODO: Change in XDMF file.
                matrixTensor[itr++] = vector[0];
                matrixTensor[itr++] = vector[1];
                matrixTensor[itr++] = vector[2];
             }_
#ifdef NO_COMPRESSION
            H5LTmake_dataset_float(file_id, icParticle->getVectorName(ivN).c_str(), 2, dimVertices, matrixTensor);
#else            
           
#endif
         }
         
         delete [] matrixTensor;

         /*  Take care of Scalars    */http://stackoverflow.com/questions/594089/does-stdvector-clear-do-delete-free-memory-on-each-element
         plint sN = icParticle->getScalarsNumber();
         float * scalarTensor = new float [Np];
         for (plint isN = 0; isN < sN; ++isN) {
             double scalarValue;
             for (plint iP = 0; iP < Np; ++iP) {
                 particles[iP]->getScalar(isN, scalarValue);http://stackoverflow.com/questions/594089/does-stdvector-clear-do-delete-free-memory-on-each-element
                 scalarTensor[iP] = scalarValue;
             }
#ifdef NO_COMPRESSION
             H5LTmake_dataset_float(file_id, icParticle->getScalarName(isN).c_str(), 1, dimVertices, scalarTensor);
#else            
            int sid = H5Screate_simple(1,dimVertices,NULL);
            int plist_id = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, chunk); 
            H5Pset_deflate(plist_id, 7);
            int did = H5Dcreate2(file_id,icParticle->getScalarName(isN).c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
            H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,scalarTensor);
            H5Dclose(did);
            H5Sclose(sid);
#endif
         }
         delete [] scalarTensor;
     }
#endif

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
        applyProcessingFunctional (bprf,cellFields[i]->getBoundingBox(), wrapper );
    }
}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

