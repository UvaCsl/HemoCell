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
#include "preInlet.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5LTpublic.h>

#ifndef H5_HAVE_PARALLEL
herr_t H5Pset_fapl_mpio( hid_t fapl_id, MPI_Comm comm, MPI_Info info ) {
  if (global::mpi().getSize() > 1) {
    pcerr << "Not compiled with HDF5 OpenMPI version, cowardly refusing to generate corrupted hdf5 files" << endl; 
    exit(1);
  }
  return 0;
}
herr_t H5Pset_dxpl_mpio( hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode ) {
  if (global::mpi().getSize() > 1) {
    pcerr << "Not compiled with HDF5 OpenMPI version, cowardly refusing to generate corrupted hdf5 files" << endl; 
    exit(1);
  }
  return 0;
}
#endif

createPreInlet::createPreInlet(Box3D domain_, string outputFileName_,int particlePositionTimestep_, Direction flowDir_, HemoCell & hemocell_, int desired_iterations_, bool reducedPrecision_) 
: hemocell(hemocell_) {
  this->desired_iterations = desired_iterations_;
  this->domain = domain_;
  this->outputFileName = global::directories().getOutputDir() + "/" + outputFileName_  + ".hdf5";
  this->particlePositionTimeStep = particlePositionTimestep_;
  this->flowDir = flowDir_;
  this->reducedPrecision = reducedPrecision_;

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
    
    if(access(outputFileName.c_str(),F_OK) == 0) {
      pcout << "WARNING: Appending to existing preinlet file, sanity and dimension checks are not performed" << endl;
      file_id = H5Fopen(outputFileName.c_str(),H5F_ACC_RDWR,plist_file_id);
      if (file_id < 0) {
        pcout << "Error opening existing preinlet file, exiting" << endl;
        exit(1);
      }
    } else {
      file_id = H5Fcreate(outputFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);      
      if (file_id < 0) {
        pcout << "Error creating preinlet file, exiting" << endl;
        exit(1);
      }
    }
  H5Pclose(plist_file_id);   

  //Attributes:
  H5LTset_attribute_int (file_id, "/", "flowDirection", (int *)&flowDir, 1);
  H5LTset_attribute_int(file_id, "/", "Number Of Cells", &hemocell.cellfields->number_of_cells,1);
  int rP;
  if (reducedPrecision) {
    rP = 1;
    reducedPrecisionDirection = (int)flowDir/2;
  } else {
    rP = 0;
  }
  H5LTset_attribute_int (file_id, "/", "reduced precision", &rP, 1);

  int sizes[3] = {(int)domain.getNx(),(int)domain.getNy(),(int)domain.getNz()};
  
  H5LTset_attribute_int (file_id, "/", "boxSize", sizes, 3);  
  
  


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
    exit(1);
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
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Pclose(plist_dataset_collective_id);
  H5Dclose(dataset_velocity_id);
  H5Fclose(file_id);
  H5Tclose(particle_type_h5);
  H5Tclose(particle_type_mem);
}

void createPreInlet::createPreInletFunctional::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  BlockLattice3D<T,DESCRIPTOR> * fluidfield = dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
   
  Box3D fluidBox; 
  Box3D shiftedBox = domain.shift(fluidfield->getLocation().x,fluidfield->getLocation().y,fluidfield->getLocation().z);
  
  if (!intersect(shiftedBox,parent.fluidDomain,fluidBox)) {return;}
  
  hsize_t offset[5]={parent.hemocell.iter%DSET_SLICE,(hsize_t)fluidBox.x0-parent.fluidDomain.x0,(hsize_t)fluidBox.y0-parent.fluidDomain.y0,(hsize_t)fluidBox.z0-parent.fluidDomain.z0,0};
  hsize_t count[5]= {1,(hsize_t)fluidBox.getNx(),(hsize_t)fluidBox.getNy(),(hsize_t)fluidBox.getNz(),3};
  if(parent.reducedPrecision) {
    count[4]=1;
  }
  if(H5Sselect_hyperslab(parent.dataspace_velocity_id, H5S_SELECT_SET, offset, NULL, count, NULL ) < 0 ) {
    cerr << "Error selecting hyperslab, exiting.." << endl; 
    exit(1);
  }
  hid_t memspace_id = H5Screate_simple (5, count, NULL); 
  
  plb::Array<T,3> vel;
  Box3D unshiftedfluidBox = fluidBox.shift(-fluidfield->getLocation().x,-fluidfield->getLocation().y,-fluidfield->getLocation().z);
  plint xx = unshiftedfluidBox.x0;
  plint yy = unshiftedfluidBox.y0;
  plint zz = unshiftedfluidBox.z0;
  plint index = 0;
  herr_t err = 0;
  if (parent.reducedPrecision) {
    static_assert(sizeof(float)==4, "Plattform float is not 32 bit, preinlet will not work");
    float * data = new float[count[1]*count[2]*count[3]*count[4]];
    
    for (plint x = 0 ; x < fluidBox.getNx() ; x++) {
      for (plint y = 0 ; y < fluidBox.getNy() ; y++) {
        for (plint z = 0; z < fluidBox.getNz() ; z++) {
          fluidfield->get(xx,yy,zz).computeVelocity(vel);
          data[index] = vel[parent.reducedPrecisionDirection];
          index++;
    
          zz++;
        }
        zz = unshiftedfluidBox.z0;
        yy++;
      }
      yy = unshiftedfluidBox.y0;
      xx++;
    }
    err = H5Dwrite(parent.dataset_velocity_id,H5T_IEEE_F32LE,memspace_id,parent.dataspace_velocity_id,parent.plist_dataset_collective_id,data);
    delete[] data;
  } else {
    static_assert(sizeof(double)==8, "Plattform double is not 64 bit, preinlet will not work");
    double * data = new double[count[1]*count[2]*count[3]*count[4]];
    for (plint x = 0 ; x < fluidBox.getNx() ; x++) {
      for (plint y = 0 ; y < fluidBox.getNy() ; y++) {
        for (plint z = 0; z < fluidBox.getNz() ; z++) {
          fluidfield->get(xx,yy,zz).computeVelocity(vel);
          data[index] = vel[0];
          data[index+1] = vel[1];
          data[index+2] = vel[2];
          index+=3;
        
          zz++;
        }
        zz = unshiftedfluidBox.z0;
        yy++;
      }
      yy = unshiftedfluidBox.y0;
      xx++;
    }
    err = H5Dwrite(parent.dataset_velocity_id,H5T_IEEE_F64LE,memspace_id,parent.dataspace_velocity_id,parent.plist_dataset_collective_id,data);
    delete[] data;
  }
 
  if ( err < 0) {
    cerr << "Error writing to hyperslab, exiting.." << endl; 
    exit(1);
  };
  H5Sclose(memspace_id);

}

void createPreInlet::saveCurrent() {    
   //When we are also saving particles
  if (hemocell.iter%particlePositionTimeStep==0) {
    map<int,int> particles_per_proc;
    vector<HemoCellParticle *> particles;
    hemocell.cellfields->getParticles(particles,domain);
    particles_per_proc[global::mpi().getRank()] = particles.size();
    
    HemoCellGatheringFunctional<int>::gather(particles_per_proc);
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
      particles_hdf5[i].location[0] = particle->sv.position[0]-domain.x0;
      particles_hdf5[i].location[1] = particle->sv.position[1]-domain.y0;
      particles_hdf5[i].location[2] = particle->sv.position[2]-domain.z0;
      particles_hdf5[i].velocity[0] = particle->sv.v[0];
      particles_hdf5[i].velocity[1] = particle->sv.v[1];
      particles_hdf5[i].velocity[2] = particle->sv.v[2];
      particles_hdf5[i].cell_id = particle->sv.cellId;
      particles_hdf5[i].vertex_id = particle->sv.vertexId;
      particles_hdf5[i].particle_type = particle->sv.celltype;
      i++;
    }
    hid_t plist_dataset_mpi_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_dataset_mpi_id, H5FD_MPIO_COLLECTIVE);  
      hid_t plist_dataset_id = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dataspace_particles_id = H5Screate_simple(1,total,NULL);
          hid_t memspace_id = H5Screate_simple (1, count, NULL); 

            std::string dataset_name = "particles_" + std::to_string(hemocell.iter);
            if(H5Lexists(file_id,dataset_name.c_str(),H5P_DEFAULT)) {
              H5Ldelete(file_id,dataset_name.c_str(),H5P_DEFAULT);
            }
            hid_t dataset_particles_id = H5Dcreate2(file_id, dataset_name.c_str() ,particle_type_h5,dataspace_particles_id,H5P_DEFAULT,plist_dataset_id,H5P_DEFAULT);    
              if (dataset_particles_id < 0) {
                cerr << "Error creating particle dataset, exiting.." <<endl;
                exit(1);
              }
              if (H5Sselect_hyperslab(dataspace_particles_id,H5S_SELECT_SET,offset,NULL,count,NULL) < 0 ) {
                cerr << "Error selecting particle hyperslab, exiting.." <<endl;
                exit(1);
              }
              if (H5Dwrite(dataset_particles_id,particle_type_mem,memspace_id,dataspace_particles_id,plist_dataset_mpi_id,particles_hdf5) < 0) {
                cerr << "Error writing particles, exiting..." <<endl;
                exit(1);
              }
            H5Dclose(dataset_particles_id);
            H5Sclose(memspace_id);
        H5Sclose(dataspace_particles_id); 
      H5Pclose(plist_dataset_id);
    H5Pclose(plist_dataset_mpi_id);
    delete[] particles_hdf5;
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  }
  
  //Velocity, 1st dim == timestep, 2-4th == space of fluid nodes, one dim will be 1, 5 == all three velocity directions 
  hsize_t dims[5] = {(hsize_t)DSET_SLICE,(hsize_t)fluidDomain.getNx(),(hsize_t)fluidDomain.getNy(),(hsize_t)fluidDomain.getNz(),3};
  if (reducedPrecision) {
    dims[4] = 1;
  }
  dataspace_velocity_id = H5Screate_simple(5, dims, NULL);
 
    //Open a new dataset for every slice of DSET_SLICE timesteps
    if ((int)(hemocell.iter/DSET_SLICE) != current_velocity_field) {
      current_velocity_field = hemocell.iter/DSET_SLICE;

      if (dataset_velocity_id != -1) {
        H5Dclose(dataset_velocity_id);
      }
      string dataset_name = "Velocity Boundary_" + to_string(current_velocity_field * DSET_SLICE);
  
      dataset_velocity_id = -1;
      if(H5Lexists(file_id,dataset_name.c_str(),H5P_DEFAULT)) {
        dataset_velocity_id = H5Dopen(file_id,dataset_name.c_str(),H5P_DEFAULT);
      }
      if(dataset_velocity_id < 0) {
        if (reducedPrecision) {
          dataset_velocity_id = H5Dcreate2(file_id,dataset_name.c_str(),H5T_IEEE_F32LE,dataspace_velocity_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        } else {
          dataset_velocity_id = H5Dcreate2(file_id,dataset_name.c_str(),H5T_IEEE_F64LE,dataspace_velocity_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        }
      }
    }

    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(hemocell.cellfields->lattice);

    applyProcessingFunctional(new createPreInletFunctional(*this),domain,wrapper);

  H5Sclose(dataspace_velocity_id);
  
  H5LTset_attribute_uint (file_id, "/", "Iterations", &hemocell.iter, 1);
}

createPreInlet::createPreInletFunctional * createPreInlet::createPreInletFunctional::clone() const {
  return new createPreInletFunctional(*this);
}
