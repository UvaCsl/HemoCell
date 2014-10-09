#ifndef FICSION_PARTICLE_HDF5IO_HH
#define FICSION_PARTICLE_HDF5IO_HH

#include "ParticleHdf5IO.h"
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

     std::map<plint, Cell3D<T,Descriptor>* > cellIdToCell3D = cellField3D.getCellIdToCell3D();
     std::vector< plint > triangles;
     std::vector<Particle3D<T,Descriptor>* > particles;
     plint sumLocalVertices=0;
     plint numCells=cellIdToCell3D.size();

     typename std::map<plint, Cell3D<T,Descriptor>* >::iterator itrtr;
     for (itrtr  = cellIdToCell3D.begin(); itrtr != cellIdToCell3D.end(); ++itrtr) {
         Cell3D<T,Descriptor> * cell3d = (itrtr->second);

         std::cout << MPI::COMM_WORLD.Get_rank() << " Cell Volume " << cell3d->getVolume()
                         << " surface " << cell3d->getSurface()
                         << " CCR_ANGLE_MEAN " << cell3d->getMeanAngle()
                         << " getMeanEdgeLength " << cell3d->getMeanEdgeLength()
                         << " getEnergy " << cell3d->getEnergy()
                         << " CCR_FORCE (" << cell3d->getForce()[0]
                         << ", " << cell3d->getForce()[1]
                         << ", " << cell3d->getForce()[2] << ") "
                 << std::endl;


         std::vector<plint> const& cellVertices = cell3d->getVertices();
         std::vector<plint> const& cellTriangles = cell3d->getTriangles();
         for (std::vector<plint>::const_iterator iVP = cellVertices.begin(); iVP != cellVertices.end(); ++iVP)
         {
             particles.push_back( cell3d->getParticle3D(*iVP) );
         }
         std::map<plint, plint> iv = cell3d->getInvertVertices();
         for (pluint iT=0; iT < cellTriangles.size(); iT++) {
             plint iTriangle=cellTriangles[iT];
             plint v0 = cell3d->getVertexId(iTriangle, 0);
             plint v1 = cell3d->getVertexId(iTriangle, 1);
             plint v2 = cell3d->getVertexId(iTriangle, 2);
//             if (castParticleToICP3D(cell3d->getParticle3D(v0))->get_processor() == id &&
//                     castParticleToICP3D(cell3d->getParticle3D(v1))->get_processor() == id &&
//                     castParticleToICP3D(cell3d->getParticle3D(v2))->get_processor() == id) {
                 plint t0 = iv[v0]+sumLocalVertices;
                 plint t1 = iv[v1]+sumLocalVertices;
                 plint t2 = iv[v2]+sumLocalVertices;
                 triangles.push_back(t0);
                 triangles.push_back(t1);
                 triangles.push_back(t2);
//             }
         }
         sumLocalVertices += cellVertices.size();
     }

     plint Np = particles.size();
     plint Nt = triangles.size() / 3;

     /************************************************************/
    /**            Initialise HDF5 file                        **/
   /************************************************************/

     std::string fileName = global::directories().getOutputDir() + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.",id,3) + ".h5";
     hid_t file_id;
     file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

     hsize_t dimTriangles[2]; dimTriangles[0] = Nt; dimTriangles[1] = 3;
     hsize_t dimVertices[2]; dimVertices[0] = Np; dimVertices[1] = 3;
     H5LTmake_dataset_long(file_id, "triangles", 2, dimTriangles, &triangles[0]);

     H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
     H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
     H5LTset_attribute_long (file_id, "/", "iteration", &iter, 1);
     H5LTset_attribute_int (file_id, "/", "numberOfProcessors", &p, 1);
     H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

     H5LTset_attribute_long (file_id, "/", "numberOfParticles", &Np, 1);
     H5LTset_attribute_long (file_id, "/", "numberOfTriangles", &Nt, 1);

     /************************************************************/
    /**            Write output to HDF5 file                   **/
   /************************************************************/
     /*  Take care of Vectors    */
     if (numCells > 0) {
         ImmersedCellParticle3D<T,Descriptor> * icParticle = castParticleToICP3D(particles[0]);
         plint vN = icParticle->getVectorsNumber();
         double * matrixTensor = new double [3 * Np];

         Array<T,3> vector;
         for (plint ivN = 0; ivN < vN; ++ivN) {
             plint itr=0;
             for (plint iP = 0; iP < Np; ++iP) {
                castParticleToICP3D(particles[iP])->getVector(ivN, vector);
                matrixTensor[itr++] = vector[0];
                matrixTensor[itr++] = vector[1];
                matrixTensor[itr++] = vector[2];
             }
             H5LTmake_dataset_double(file_id, icParticle->getVectorName(ivN).c_str(), 2, dimVertices, matrixTensor);
         }
         delete [] matrixTensor;

         /*  Take care of Scalars    */
         plint sN = icParticle->getScalarsNumber();
         double * scalarTensor = new double [Np];
         for (plint isN = 0; isN < sN; ++isN) {
             for (plint iP = 0; iP < Np; ++iP) {
                castParticleToICP3D(particles[iP])->getScalar(isN, scalarTensor[iP]);
             }
             H5LTmake_dataset_double(file_id, icParticle->getScalarName(isN).c_str(), 1, dimVertices, scalarTensor);
         }
         delete [] scalarTensor;
     }
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

