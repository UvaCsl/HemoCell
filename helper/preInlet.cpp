#include "preInlet.h"
#include <hdf5.h>
#include <hdf5_hl.h>

#ifndef H5_HAVE_PARALLEL
herr_t H5Pset_fapl_mpio( hid_t fapl_id, MPI_Comm comm, MPI_Info info ) {
  if (global::mpi().getSize() > 1) {
    pcerr << "Not compiled with HDF5 OpenMPI version, cowardly refusing to generate corrupted hdf5 files" << endl; 
    exit(0); 
  }
  return 0;
}
herr_t H5Pset_dxpl_mpio( hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode ) {
  if (global::mpi().getSize() > 1) {
    pcerr << "Not compiled with HDF5 OpenMPI version, cowardly refusing to generate corrupted hdf5 files" << endl; 
    exit(0); 
  }
  return 0;
}
#endif

PreInlet::PreInlet(Box3D domain_, string sourceFileName_, int particlePositionTimestep_, Direction flowDir_, HemoCell& hemocell_) 
  : hemocell(hemocell_) {
  this->domain = domain_;
  this->sourceFileName = global::directories().getOutputDir() + "/" + sourceFileName_  + ".hdf5";
  this->particlePositionTimeStep = particlePositionTimestep_;
  this->flowDir = flowDir_;

  fluidDomain = domain;
  switch(flowDir) {
    case Xpos:
      fluidDomain.x1 = fluidDomain.x0;
      break;
    case Xneg:
      fluidDomain.x0 = fluidDomain.x1;
      break;
    case Ypos:
      fluidDomain.y1 = fluidDomain.y0;
      break;
    case Yneg:
      fluidDomain.y0 = fluidDomain.y1;
      break;
    case Zpos:
      fluidDomain.z1 = fluidDomain.z0;
      break;
    case Zneg:
      fluidDomain.z0 = fluidDomain.z1;
      break;
  }

  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
  boundary->setVelocityConditionOnBlockBoundaries(*hemocell.lattice,fluidDomain);
  setBoundaryVelocity(*hemocell.lattice, fluidDomain, plb::Array<T,3>(0.,0.,0.));
  
  hid_t plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_file_id, global::mpi().getGlobalCommunicator(), info);  
      
    file_id = H5Fopen(sourceFileName.c_str(), H5F_ACC_RDONLY, plist_file_id);
    if(file_id < 0) { cout << "Error opening Preinlet hdf5 file." << endl; exit(0); }
  H5Pclose(plist_file_id);
  
  Direction stored_flow;
  if(H5LTget_attribute_int(file_id,"/","flowDirection",(int *)&stored_flow) < 0) { 
    cout << "Error reading flowDirection attribute" << endl; exit(0); 
  }
  if (stored_flow != flowDir) {
    cout << "Flow direction does not match, exiting ..." << endl;
    exit(0);
  }
  int stored_size[3];
  if(H5LTget_attribute_int(file_id,"/","boxSize",stored_size) < 0) { 
    cout << "Error reading boxSize attribute" << endl; exit(0); 
  }
  if (domain.getNx() != stored_size[0] ||
      domain.getNy() != stored_size[1] ||
      domain.getNz() != stored_size[2]) {
    cout << "Box dimensions do not match, exiting ..." << endl;
    exit(0);
  }
    
  dataset_velocity_id = H5Dopen2( file_id, "Velocity Boundary", H5P_DEFAULT);  
  dataspace_velocity_id = H5Dget_space(dataset_velocity_id);
}

void PreInlet::PreInletFunctional::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  BlockLattice3D<double,DESCRIPTOR> * fluidfield = dynamic_cast<BlockLattice3D<double,DESCRIPTOR>*>(blocks[0]);
  HemoCellParticleField * particlefield = dynamic_cast<HemoCellParticleField*>(blocks[1]);
   
  //When we are also saving particles
  if (parent.hemocell.iter%parent.particlePositionTimeStep==0) {
    
  }
  
  Box3D fluidBox; 
  Box3D shiftedBox = domain.shift(fluidfield->getLocation().x,fluidfield->getLocation().y,fluidfield->getLocation().z);
  
  if (!intersect(shiftedBox,parent.fluidDomain,fluidBox)) {return;}
  
  hsize_t offset[5]={parent.hemocell.iter,(hsize_t)fluidBox.x0-parent.fluidDomain.x0,(hsize_t)fluidBox.y0-parent.fluidDomain.y0,(hsize_t)fluidBox.z0-parent.fluidDomain.z0,0};
  hsize_t count[5]= {1,(hsize_t)fluidBox.getNx(),(hsize_t)fluidBox.getNy(),(hsize_t)fluidBox.getNz(),3};
  if(H5Sselect_hyperslab(parent.dataspace_velocity_id, H5S_SELECT_SET, offset, NULL, count, NULL ) < 0 ) {
    cerr << "Error selecting hyperslab, exiting.." << endl; 
    exit(0);
  }
  hid_t memspace_id = H5Screate_simple (5, count, NULL); 

  double *data = new double[count[1]*count[2]*count[3]*count[4]];

  H5Dread(parent.dataset_velocity_id,H5T_NATIVE_DOUBLE,memspace_id,parent.dataspace_velocity_id,H5P_DEFAULT,data);
  
  plb::Array<double,3> vel;
  Box3D unshiftedfluidBox = fluidBox.shift(-fluidfield->getLocation().x,-fluidfield->getLocation().y,-fluidfield->getLocation().z);
  plint xx = unshiftedfluidBox.x0;
  plint yy = unshiftedfluidBox.y0;
  plint zz = unshiftedfluidBox.z0;
  plint index = 0;
  for (plint z = 0 ; z < fluidBox.getNz() ; z++) {
    for (plint y = 0 ; y < fluidBox.getNy() ; y++) {
      for (plint x = 0; x < fluidBox.getNx() ; x++) {
        //Skip boundaries
        /*if (fluidfield->get(xx,yy,zz).getDynamics().isBoundary()) {
          index+=3;
          xx++;
          continue;
        };*/

        vel[0] = data[index];
        vel[1] = data[index+1];
        vel[2] = data[index+2];

        fluidfield->get(xx,yy,zz).defineVelocity(vel);

        index+=3;
        xx++;
      }
      xx = unshiftedfluidBox.x0;
      yy++;
    }
    yy = unshiftedfluidBox.y0;
    zz++;
  }
  H5Sclose(memspace_id);
  free(data);
}
void PreInlet::update() {
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell.cellfields->lattice);
  
  wrapper.push_back(hemocell.cellfields->immersedParticles);
  applyProcessingFunctional(new PreInletFunctional(*this),domain,wrapper);
}

PreInlet::~PreInlet() {
  H5Dclose(dataset_velocity_id);
  H5Sclose(dataspace_velocity_id);
  H5Fclose(file_id);
}

PreInlet::PreInletFunctional * PreInlet::PreInletFunctional::clone() const {
  return new PreInletFunctional(*this);
}