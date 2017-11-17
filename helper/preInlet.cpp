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
  if(H5LTget_attribute_int(file_id,"/","Number Of Cells",&nCellsSelf) < 0) { 
    cout << "Error reading number of cells attribute" << endl; exit(0); 
  }
  nCellsOffset = hemocell.cellfields->number_of_cells;
  hemocell.cellfields->number_of_cells += nCellsSelf;
     
  dataset_velocity_id = H5Dopen2( file_id, "Velocity Boundary", H5P_DEFAULT);  
  dataspace_velocity_id = H5Dget_space(dataset_velocity_id);
  
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
}

void PreInlet::ImmersePreInletParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HemoCellParticleField * particlefield = dynamic_cast<HemoCellParticleField*>(blocks[0]);
  const map<int,vector<int>> & ppc = particlefield->get_preinlet_particles_per_cell();
  
  for (auto & pair : ppc) {
    bool complete = true;
    for (int entry : pair.second) {
      if (entry == -1) {
        complete = false;
        break;
      }
    }
    if (complete) {
      for (int entry:pair.second) {
        particlefield->particles[entry].fromPreInlet = false;
      }
    }
  }
}

void PreInlet::DeletePreInletParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HemoCellParticleField * particlefield = dynamic_cast<HemoCellParticleField*>(blocks[0]);
  
  for (HemoCellParticle & particle : particlefield->particles) {
    if (particle.fromPreInlet) {
      particle.tag = 1;
    }
  }
  particlefield->removeParticles(1);
}

void PreInlet::PreInletFunctional::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  BlockLattice3D<double,DESCRIPTOR> * fluidfield = dynamic_cast<BlockLattice3D<double,DESCRIPTOR>*>(blocks[0]);
  HemoCellParticleField * particlefield = dynamic_cast<HemoCellParticleField*>(blocks[1]);
  
  //When we are also saving particles
  if (parent.hemocell.iter%parent.particlePositionTimeStep==0) {
    const map<int,vector<int>> & ppc = particlefield->get_particles_per_cell();
    for (int i = 0 ; i < parent.particles_size ; i++) {
      particle_hdf5_t * particle = &parent.particles[i];
      hemo::Array<double,3> loc; 
      loc[0] = particle->location[0];
      loc[1] = particle->location[1];
      loc[2] = particle->location[2];
      
      
      if (ppc.find(particle->cell_id) != ppc.end()) {
        if (ppc.at(particle->cell_id)[particle->vertex_id] != -1) {
          if (! particlefield->particles[ppc.at(particle->cell_id)[particle->vertex_id]].fromPreInlet) {
            continue;
          }
        }
      }

      HemoCellParticle particle_h = HemoCellParticle(loc,particle->cell_id,particle->vertex_id,particle->particle_type);
      particle_h.fromPreInlet = true;
      hemo::Array<double,3> vel; 
      vel[0] = particle->velocity[0];
      vel[1] = particle->velocity[1];
      vel[2] = particle->velocity[2];
      particle_h.v = vel;
      particlefield->addParticle(domain,&particle_h);
    }
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
  delete[] data;
}
void PreInlet::update() {  
  //When we are also saving particles
  if (hemocell.iter%particlePositionTimeStep==0) {
    std::string filename = "particles_" + std::to_string(hemocell.iter);
    dataset_particles_id = H5Dopen2( file_id, filename.c_str() , H5P_DEFAULT);  
    dataspace_particles_id = H5Dget_space(dataset_particles_id);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace_particles_id,dims,NULL);
    particles_size = dims[0];
    particles = new particle_hdf5_t[particles_size];
    H5Dread(dataset_particles_id,particle_type_mem,H5S_ALL,dataspace_particles_id,H5P_DEFAULT,particles);
    for(int i = 0 ; i < particles_size ; i++) {
      particles[i].location[0] += domain.x0;
      particles[i].location[1] += domain.y0;
      particles[i].location[2] += domain.z0;
      
      //Generate unique id, we know everything to calculate it, quite complex
      int wraps = 0;
      if (particles[i].cell_id < 0) {
        wraps = ((particles[i].cell_id*-1)/nCellsSelf)*-1; //Negative rounding is implementation defined, therefore avoid it;
      } else {
        wraps = particles[i].cell_id/nCellsSelf;
      }
      int particleBaseId = (((particles[i].cell_id % nCellsSelf) + particles[i].cell_id) % nCellsSelf) + nCellsOffset; //Modulo from remainder operation
      particles[i].cell_id = particleBaseId+(wraps*hemocell.cellfields->number_of_cells);
    }
  }
  
  
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell.cellfields->immersedParticles);
  applyProcessingFunctional(new DeletePreInletParticles(),hemocell.cellfields->immersedParticles->getBoundingBox(),wrapper);
  
  wrapper.clear();
  wrapper.push_back(hemocell.cellfields->lattice);
  wrapper.push_back(hemocell.cellfields->immersedParticles);
  applyProcessingFunctional(new PreInletFunctional(*this),domain,wrapper);
  
  hemocell.cellfields->syncEnvelopes();
  
  wrapper.clear();
  wrapper.push_back(hemocell.cellfields->immersedParticles);
  applyProcessingFunctional(new ImmersePreInletParticles(),domain,wrapper);
  
  if (hemocell.iter%particlePositionTimeStep==0) {
   H5Dclose(dataset_particles_id); 
   H5Sclose(dataspace_particles_id);
   delete[] particles;
  }
}

PreInlet::~PreInlet() {
  H5Dclose(dataset_velocity_id);
  H5Sclose(dataspace_velocity_id);
  H5Fclose(file_id);
}

PreInlet::PreInletFunctional * PreInlet::PreInletFunctional::clone() const { return new PreInletFunctional(*this); }
PreInlet::DeletePreInletParticles * PreInlet::DeletePreInletParticles::clone() const { return new DeletePreInletParticles(*this); }
PreInlet::ImmersePreInletParticles * PreInlet::ImmersePreInletParticles::clone() const { return new ImmersePreInletParticles(*this); }
