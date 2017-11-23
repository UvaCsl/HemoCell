#include "preInlet.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <hdf5/openmpi/H5LTpublic.h>

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

createPreInlet::createPreInlet(Box3D domain_, string outputFileName_,int particlePositionTimestep_, Direction flowDir_, HemoCell & hemocell_) 
: hemocell(hemocell_) {
  this->domain = domain_;
  this->outputFileName = global::directories().getOutputDir() + "/" + outputFileName_  + ".hdf5";
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
  
  
  hid_t plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_file_id, global::mpi().getGlobalCommunicator(), info);  
    file_id = H5Fcreate(outputFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);
  H5Pclose(plist_file_id);   

  //Attributes:
  H5LTset_attribute_int (file_id, "/", "flowDirection", (int *)&flowDir, 1);
  H5LTset_attribute_int(file_id, "/", "Number Of Cells", &hemocell.cellfields->number_of_cells,1);
  
  //Velocity, 1st dim == timestep, 2-4th == space of fluid nodes, one dim will be 1, 5 == all three velocity directions 
  hsize_t dims[5] = {0,(hsize_t)fluidDomain.getNx(),(hsize_t)fluidDomain.getNy(),(hsize_t)fluidDomain.getNz(),3};
  hsize_t maximum_size[5] = {H5S_UNLIMITED,(hsize_t)fluidDomain.getNx(),(hsize_t)fluidDomain.getNy(),(hsize_t)fluidDomain.getNz(),3};
  hsize_t chunks[5] = {1,dims[1],dims[2],dims[3],dims[4]};
  int sizes[3] = {(int)domain.getNx(),(int)domain.getNy(),(int)domain.getNz()};
  
  H5LTset_attribute_int (file_id, "/", "boxSize", sizes, 3);  
  
  hid_t plist_dataset_id = H5Pcreate(H5P_DATASET_CREATE);
    if (H5Pset_chunk(plist_dataset_id, 5, chunks) < 0) { 
      cerr << "Error creating PreInlet HDF5 File, exiting.." << endl; 
      exit(0);
    }
  //  H5Pset_deflate(plist_dataset_id, 7);

    dataspace_velocity_id = H5Screate_simple(5, dims, maximum_size);
    dataset_velocity_id = H5Dcreate2(file_id,"Velocity Boundary",H5T_NATIVE_DOUBLE,dataspace_velocity_id,H5P_DEFAULT,plist_dataset_id,H5P_DEFAULT);
  H5Pclose(plist_dataset_id);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);

  //Create Particle Type
  size_t particle_type_size = sizeof(particle_hdf5_t);
  particle_type_mem = H5Tcreate(H5T_COMPOUND,particle_type_size);
  H5Tinsert (particle_type_mem, "location_X", HOFFSET (particle_hdf5_t, location[0]), H5T_NATIVE_FLOAT);
  H5Tinsert (particle_type_mem, "location_Y", HOFFSET (particle_hdf5_t, location[1]), H5T_NATIVE_FLOAT);
  H5Tinsert (particle_type_mem, "location_Z", HOFFSET (particle_hdf5_t, location[2]), H5T_NATIVE_FLOAT);
  H5Tinsert (particle_type_mem, "velocity_X", HOFFSET (particle_hdf5_t, velocity[0]), H5T_NATIVE_FLOAT);
  H5Tinsert (particle_type_mem, "velocity_Y", HOFFSET (particle_hdf5_t, velocity[1]), H5T_NATIVE_FLOAT);
  H5Tinsert (particle_type_mem, "velocity_Z", HOFFSET (particle_hdf5_t, velocity[2]), H5T_NATIVE_FLOAT);
  H5Tinsert (particle_type_mem, "cell_id", HOFFSET (particle_hdf5_t, cell_id), H5T_NATIVE_INT);
  H5Tinsert (particle_type_mem, "vertex_id", HOFFSET (particle_hdf5_t, vertex_id), H5T_NATIVE_INT);
  H5Tinsert (particle_type_mem, "particle_type", HOFFSET (particle_hdf5_t, particle_type), H5T_NATIVE_UINT);

  
  //TODO immutable particle type, but its at least consistent over ubuntu installations
  particle_type_size = H5Tget_size(H5T_IEEE_F32LE)*6+H5Tget_size(H5T_STD_I32LE)*2+H5Tget_size(H5T_STD_U32LE);
  particle_type_h5 = H5Tcreate(H5T_COMPOUND,particle_type_size);
  if (particle_type_h5  < 0) {
    cerr << "Error creating particle type hdf5, exiting.." << endl;
    exit(0);
  };
  H5Tinsert (particle_type_h5, "location_X", 0, H5T_IEEE_F32LE); 
  H5Tinsert (particle_type_h5, "location_Y", H5Tget_size(H5T_IEEE_F32LE), H5T_IEEE_F32LE);
  H5Tinsert (particle_type_h5, "location_Z", H5Tget_size(H5T_IEEE_F32LE)*2, H5T_IEEE_F32LE);
  H5Tinsert (particle_type_h5, "velocity_X", H5Tget_size(H5T_IEEE_F32LE)*3, H5T_IEEE_F32LE);
  H5Tinsert (particle_type_h5, "velocity_Y", H5Tget_size(H5T_IEEE_F32LE)*4, H5T_IEEE_F32LE);
  H5Tinsert (particle_type_h5, "velocity_Z", H5Tget_size(H5T_IEEE_F32LE)*5, H5T_IEEE_F32LE);
  H5Tinsert (particle_type_h5, "cell_id", H5Tget_size(H5T_IEEE_F32LE)*6, H5T_STD_I32LE);
  H5Tinsert (particle_type_h5, "vertex_id", H5Tget_size(H5T_IEEE_F32LE)*6+H5Tget_size(H5T_STD_I32LE), H5T_STD_I32LE);
  H5Tinsert (particle_type_h5, "particle_type", H5Tget_size(H5T_IEEE_F32LE)*6+H5Tget_size(H5T_STD_I32LE)*2, H5T_STD_U32LE);

      
  plist_dataset_collective_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_dataset_collective_id, H5FD_MPIO_INDEPENDENT);  
}

createPreInlet::~createPreInlet() {
  H5Pclose(plist_dataset_collective_id);
  H5Dclose(dataset_velocity_id);
  H5Sclose(dataspace_velocity_id);
  H5Fclose(file_id);
  H5Tclose(particle_type_h5);
  H5Tclose(particle_type_mem);
}

void createPreInlet::createPreInletFunctional::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  BlockLattice3D<double,DESCRIPTOR> * fluidfield = dynamic_cast<BlockLattice3D<double,DESCRIPTOR>*>(blocks[0]);
   
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
  plb::Array<double,3> vel;
  Box3D unshiftedfluidBox = fluidBox.shift(-fluidfield->getLocation().x,-fluidfield->getLocation().y,-fluidfield->getLocation().z);
  plint xx = unshiftedfluidBox.x0;
  plint yy = unshiftedfluidBox.y0;
  plint zz = unshiftedfluidBox.z0;
  plint index = 0;
  for (plint z = 0 ; z < fluidBox.getNz() ; z++) {
    for (plint y = 0 ; y < fluidBox.getNy() ; y++) {
      for (plint x = 0; x < fluidBox.getNx() ; x++) {
        fluidfield->grid[xx][yy][zz].computeVelocity(vel);
        data[index] = vel[0];
        data[index+1] = vel[1];
        data[index+2] = vel[2];
        index+=3;
        
        xx++;
      }
      xx = unshiftedfluidBox.x0;
      yy++;
    }
    yy = unshiftedfluidBox.y0;
    zz++;
  }

  if (H5Dwrite(parent.dataset_velocity_id,H5T_NATIVE_DOUBLE,memspace_id,parent.dataspace_velocity_id,parent.plist_dataset_collective_id,data) < 0) {
    cerr << "Error writing to hyperslab, exiting.." << endl; 
    exit(0);
  };
  H5Sclose(memspace_id);
  delete[] data;

}

void createPreInlet::saveCurrent() {    
   //When we are also saving particles
  if (hemocell.iter%particlePositionTimeStep==0) {
    map<int,int> particles_per_proc;
    vector<HemoCellParticle *> particles;
    hemocell.cellfields->getParticles(particles,domain);
    particles_per_proc[global::mpi().getRank()] = particles.size();
    
    HemoCellGatheringFunctional<int>::gather(particles_per_proc,global::mpi().getSize());
    hsize_t offset[1] = {0}, count[1] = {(hsize_t)particles.size()}, total[1] = {0};
    for (auto & pair : particles_per_proc)  {
      if (pair.first == global::mpi().getRank()) {
        break;
      }
      offset[0] += pair.second;
    }
    for (auto & pair : particles_per_proc) {
      total[0] += pair.second;
    }
    particle_hdf5_t * particles_hdf5 = new particle_hdf5_t[particles.size()];
    
    int i = 0;
    for (HemoCellParticle * particle : particles) {
      particles_hdf5[i].location[0] = particle->position[0]-domain.x0;
      particles_hdf5[i].location[1] = particle->position[1]-domain.y0;
      particles_hdf5[i].location[2] = particle->position[2]-domain.z0;
      particles_hdf5[i].velocity[0] = particle->v[0];
      particles_hdf5[i].velocity[1] = particle->v[1];
      particles_hdf5[i].velocity[2] = particle->v[2];
      particles_hdf5[i].cell_id = particle->cellId;
      particles_hdf5[i].vertex_id = particle->vertexId;
      particles_hdf5[i].particle_type = particle->celltype;
      i++;
    }
    hid_t plist_dataset_id = H5Pcreate(H5P_DATASET_CREATE);
      hid_t dataspace_particles_id = H5Screate_simple(1,total,NULL);
        hid_t memspace_id = H5Screate_simple (1, count, NULL); 

          std::string dataset_name = "particles_" + std::to_string(hemocell.iter);
          hid_t dataset_particles_id = H5Dcreate2(file_id, dataset_name.c_str() ,particle_type_h5,dataspace_particles_id,H5P_DEFAULT,plist_dataset_id,H5P_DEFAULT);    
            if (dataset_particles_id < 0) {
              cerr << "Error creating particle dataset, exiting.." <<endl;
              exit(0);
            }
            if (H5Sselect_hyperslab(dataspace_particles_id,H5S_SELECT_SET,offset,NULL,count,NULL) < 0 ) {
              cerr << "Error selecting particle hyperslab, exiting.." <<endl;
              exit(0);
            }
            if (H5Dwrite(dataset_particles_id,particle_type_mem,memspace_id,dataspace_particles_id,plist_dataset_collective_id,particles_hdf5) < 0) {
              cerr << "Error writing particles, exiting..." <<endl;
              exit(0);
            }
          H5Dclose(dataset_particles_id);
          H5Sclose(memspace_id);
      H5Sclose(dataspace_particles_id); 
    H5Pclose(plist_dataset_id);
    delete[] particles_hdf5;
  }
    
  hsize_t size[5] = {hemocell.iter+1,(hsize_t)fluidDomain.getNx(),(hsize_t)fluidDomain.getNy(),(hsize_t)fluidDomain.getNz(),3};
  if (H5Dset_extent( dataset_velocity_id, size ) < 0) {
    cerr << "Error enlarging extent, exiting.." << endl; 
    exit(0);
  }
  H5Sclose(dataspace_velocity_id);
  dataspace_velocity_id = H5Dget_space(dataset_velocity_id);

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell.cellfields->lattice);
  
  applyProcessingFunctional(new createPreInletFunctional(*this),domain,wrapper);
  //H5Fflush(file_id,H5F_SCOPE_GLOBAL);
}

createPreInlet::createPreInletFunctional * createPreInlet::createPreInletFunctional::clone() const {
  return new createPreInletFunctional(*this);
}