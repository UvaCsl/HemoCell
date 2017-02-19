#ifndef FICSION_PARTICLE_HDF5IO_HH
#define FICSION_PARTICLE_HDF5IO_HH

#include "ParticleHdf5IO.h"
#include "ficsionInit.h"



/* ******** WriteCellField3DInMultipleHDF5Files *********************************** */
WriteCellField3DInMultipleHDF5Files::WriteCellField3DInMultipleHDF5Files (
        HemoCellField & cellField3D_,
        plint iter_, std::string identifier_,
        double dx_, double dt_) :
        cellField3D(cellField3D_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_) {};

void WriteCellField3DInMultipleHDF5Files::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size() > 0 );
    int p = global::mpi().getSize();
    int id = global::mpi().getRank();
    plint Nx = cellField3D.getParticleArg()->getNx();
    plint Ny = cellField3D.getParticleArg()->getNy();
    plint Nz = cellField3D.getParticleArg()->getNz();

    /************************************************************/
   /**            Fill triangle and particle lists            **/
  /************************************************************/

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
         Array<double,3> cellPosition = cell3d->getPosition();
         NpBulk += cell3d->getNumVertices_LocalBulk();
         correctPBPosition[itrtr->first].resetToZero();

         if (cellPosition[0] > Nx) { correctPBPosition[itrtr->first][0] = -int(cellPosition[0]/Nx)*Nx;}
         if (cellPosition[0] <  0) { correctPBPosition[itrtr->first][0] =  int(cellPosition[0]/Nx)*Nx;}
         if (cellPosition[1] > Ny) { correctPBPosition[itrtr->first][1] = -int(cellPosition[1]/Ny)*Ny;}
         if (cellPosition[1] <  0) { correctPBPosition[itrtr->first][1] =  int(cellPosition[1]/Ny)*Ny;}
         if (cellPosition[2] > Nz) { correctPBPosition[itrtr->first][2] = -int(cellPosition[2]/Nz)*Nz;}
         if (cellPosition[2] <  0) { correctPBPosition[itrtr->first][2] =  int(cellPosition[2]/Nz)*Nz;}

//         std::cout << MPI::COMM_WORLD.Get_rank() << " Cell Volume " << cell3d->getVolume()
//                         << " surface " << cell3d->getSurface()
//                         << " CCR_ANGLE_MEAN " << cell3d->getMeanAngle()
//                         << " getMeanEdgeLength " << cell3d->getMeanEdgeLength()
//                         << " getEnergy " << cell3d->getEnergy()
//                         << " CCR_FORCE (" << cell3d->getForce()[0]
//                         << ", " << cell3d->getForce()[1]
//                         << ", " << cell3d->getForce()[2] << ") "
//                 << std::endl;


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

     /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

     hsize_t dimTriangles[2]; 
     dimTriangles[0] = Nt; dimTriangles[1] = 3;
     
     hsize_t dimVertices[2]; 
     dimVertices[0] = Np; dimVertices[1] = 3;

     hsize_t chunk[2];  
     chunk[0] = (Np/p) > 1 ? (Np/p) : 1; chunk[1] = 3;

     H5LTmake_dataset_long(file_id, "triangles", 2, dimTriangles, &triangles[0]);

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     long int iterHDF5=iter;
     H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

     H5LTset_attribute_long (file_id, "/", "numberOfParticles", &Np, 1);
     H5LTset_attribute_long (file_id, "/", "numberOfParticlesBulk", &NpBulk, 1);
     H5LTset_attribute_long (file_id, "/", "numberOfTriangles", &Nt, 1);

     /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
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
             }
#ifdef NO_COMPRESSION
            H5LTmake_dataset_float(file_id, icParticle->getVectorName(ivN).c_str(), 2, dimVertices, matrixTensor);
#else            
            int sid = H5Screate_simple(2,dimVertices,NULL);
            int plist_id = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 2, chunk); 
            H5Pset_deflate(plist_id, 7);
            int did = H5Dcreate2(file_id,icParticle->getVectorName(ivN).c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
            H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,matrixTensor);
            H5Dclose(did);
            H5Sclose(sid);
#endif
         }
         // TODO: This part is commented out, since it causes some MPI problems!
         /*  Take care of Vertex Normals */
         /*
         plint itr=0;
         for (plint iP = 0; iP < Np; ++iP) {
            icParticle = particles[iP];

            //TODO: This fails when a particle moves to an other core!
            Cell3D<T,Descriptor> * cell3d = cellIdToCell3D[icParticle->get_cellId()];
            vector = cell3d->computeVertexNormal(icParticle->getVertexId());

            // TODO: Change in XDMF file.
            matrixTensor[itr++] = vector[0];
            matrixTensor[itr++] = vector[1];
            matrixTensor[itr++] = vector[2];
         }
         H5LTmake_dataset_float(file_id, "normal", 2, dimVertices, matrixTensor);
        */

         delete [] matrixTensor;

         /*  Take care of Scalars    */
         plint sN = icParticle->getScalarsNumber();
         float * scalarTensor = new float [Np];
         for (plint isN = 0; isN < sN; ++isN) {
             double scalarValue;
             for (plint iP = 0; iP < Np; ++iP) {
                 particles[iP]->getScalar(isN, scalarValue);
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



void writeCellField3D_HDF5(HemoCellField& cellField3D, double dx, double dt, plint iter, std::string preString)
{
	std::string identifier = preString + cellField3D.getIdentifier();
    applyProcessingFunctional ( // compute force applied on the fluid by the particles
            (BoxProcessingFunctional3D*)new WriteCellField3DInMultipleHDF5Files(cellField3D, iter, identifier, dx, dt),
            cellField3D.getBoundingBox(), cellField3D.getParticleArg() );

}





#endif  // FICSION_HDF5IO_HH

//  Get the number of processes.

