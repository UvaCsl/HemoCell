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

#include "palabos3D.h"
#include "palabos3D.hh"
#include "preInlet.h"

#include "boundaryCondition/boundaryInstantiator3D.h"

namespace hemo {
  struct Box3D_simple {
    plint x0,x1,y0,y1,z0,z1;
    inline Box3D_simple & operator=(const Box3D & box) {
      x0 = box.x0;
      x1 = box.x1;
      y0 = box.y0;
      y1 = box.y1;
      z0 = box.z0;
      z1 = box.z1;
      return *this;
    }
    inline Box3D getBox3D() const {
      return Box3D(x0,x1,y0,y1,z0,z1);
    }
  };
}

#include <hdf5.h>
#include <hdf5_hl.h>


namespace hemo {

#ifndef H5_HAVE_PARALLEL
herr_t H5Pset_fapl_mpio( hid_t fapl_id, MPI_Comm comm, MPI_Info info ) {
  if (plb::global::mpi().getSize() > 1) {
    hlog << "Not compiled with HDF5 OpenMPI version, cowardly refusing to generate corrupted hdf5 files" << std::endl; 
    exit(1);
  }
  return 0;
}
#endif

void applyPreInletVelocityBoundary(HemoCell & hemocell) {
  MPI_Barrier(MPI_COMM_WORLD);
  vector<MPI_Request> reqs;
  Box3D & domain = hemocell.preInlet.fluidInlet;
  Box3D result;
  int z_o = domain.getNz();
  int y_o = domain.getNy();
  int tag;
  plb::Array<T,3> vel;
  for (int bId : hemocell.lattice->getLocalInfo().getBlocks()) {
    Box3D bulk = hemocell.lattice->getMultiBlockManagement().getBulk(bId);
    if (!intersect(domain,bulk,result)) { continue; }
  
    const Dot3D & loc = hemocell.lattice->getComponent(bId).getLocation();
    result = result.shift(-loc.x,-loc.y,-loc.z);
    
    for (int x  = result.x0 ; x <= result.x1 ; x++) {
     for (int y  = result.y0 ; y <= result.y1 ; y++) {
      for (int z  = result.z0 ; z <= result.z1 ; z++) {
        tag = ((z_o * ((z + loc.z) - domain.z0)) + ( y + loc.y ) - domain.y0) * y_o + x + loc.x - domain.x0;
        if (!hemocell.lattice->get(x,y,z).getDynamics().isBoundary()) {
          if (hemocell.partOfpreInlet) {
            hemocell.lattice->getComponent(bId).get(x,y,z).computeVelocity(vel);
            int dest = hemocell.lattice_management->getThreadAttribution().getMpiProcess(hemocell.lattice_management->getSparseBlockStructure().locate(x+loc.x,y+loc.y,z+loc.z));
            reqs.emplace_back();
            MPI_Isend(&vel[0],3*sizeof(T),MPI_CHAR,dest,tag,MPI_COMM_WORLD,&reqs.back());
          } else {
            Box3D point(x,x,y,y,z,z);
            int source = hemocell.preinlet_management->getThreadAttribution().getMpiProcess(hemocell.preinlet_management->getSparseBlockStructure().locate(x+loc.x,y+loc.y,z+loc.z));
            MPI_Recv(&vel[0],3*sizeof(T),MPI_CHAR,source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            setBoundaryVelocity(hemocell.lattice->getComponent(bId),point,vel);
          }
        }
      }
     }
    }
  }
  MPI_Waitall(reqs.size(),&reqs[0],MPI_STATUS_IGNORE);
}
void createPreInletVelocityBoundary(plb::MultiBlockLattice3D<T,DESCRIPTOR> * fluid, plb::MultiScalarField3D<int> * flagMatrix,plb::Array<double,3> speed, HemoCell & hemocell) {
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* bc =
        createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
  Box3D domain;
  intersect(flagMatrix->getBoundingBox(),hemocell.preInlet.location,domain);
  domain.z0 = domain.z1; 
  hemocell.preInlet.fluidInlet = domain;
  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
   for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
    for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
      if (flagMatrix->get(x,y,z) == 1 && !hemocell.partOfpreInlet) {
        Box3D point(x,x,y,y,z,z);
        bc->addVelocityBoundary2N(point,*hemocell.lattice);
        setBoundaryVelocity(*hemocell.lattice, point, speed );
      }
    }
   }
  }
}

void PreInlet::CreatePreInletBoundingBox::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  ScalarField3D<int> * sf = dynamic_cast<ScalarField3D<int>*>(blocks[0]);
  
  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
    for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
      for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
        if (sf->get(x,y,z)) {
          int xx = x + sf->getLocation().x;
          int yy = y + sf->getLocation().y;
          int zz = z + sf->getLocation().z;

          fluidslice[xx][yy] = true;
          if(!foundPreInlet) {
            boundingBox.x0 = boundingBox.x1 = xx;
            boundingBox.y0 = boundingBox.y1 = yy;
            boundingBox.z0 = boundingBox.z1 = zz;
            foundPreInlet = true;
          } else {
            if (xx < boundingBox.x0) { boundingBox.x0 = xx; }
            if (xx > boundingBox.x1) { boundingBox.x1 = xx; }
            if (yy < boundingBox.y0) { boundingBox.y0 = yy; }
            if (yy > boundingBox.y1) { boundingBox.y1 = yy; }
            if (zz < boundingBox.z0) { boundingBox.z0 = zz; }
            if (zz > boundingBox.z1) { boundingBox.z1 = zz; }
          }
        }
      }
    }
  }
}

PreInlet::PreInlet(MultiScalarField3D<int>& flagMatrix) {
  initialized = true;
  Box3D fluidDomain = flagMatrix.getMultiBlockManagement().getBoundingBox();
 
  fluidDomain.z1 = fluidDomain.z0+1;
  fluidDomain.z0 = fluidDomain.z1;
  fluidslice.resize((fluidDomain.x1-fluidDomain.x0)+1);
  for (auto & slice : fluidslice) {
    slice.resize((fluidDomain.y1-fluidDomain.y0)+1,false);
  }
  Box3D preInletDomain;
  bool foundPreInlet = false;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(&flagMatrix);
  
  applyProcessingFunctional(new CreatePreInletBoundingBox(preInletDomain,foundPreInlet,fluidslice),fluidDomain,wrapper);
  
  map<int,Box3D_simple> values;
  if (foundPreInlet) {
    values[global::mpi().getRank()] = preInletDomain;
  }
  HemoCellGatheringFunctional<Box3D_simple>::gather(values);
  
  if (values.size() == 0) {
    hlog << "(PreInlet) no preinlet found, does the stl go up to the Z-Negative wall?" << endl;
    exit(1);
  }
  preInletDomain =  values.begin()->second.getBox3D();
  for (auto & pair : values ) {
    Box3D value  = pair.second.getBox3D();
    if (value.x0 < preInletDomain.x0) { preInletDomain.x0 = value.x0; }
    if (value.x1 > preInletDomain.x1) { preInletDomain.x1 = value.x1; }
    if (value.y0 < preInletDomain.y0) { preInletDomain.y0 = value.y0; }
    if (value.y1 > preInletDomain.y1) { preInletDomain.y1 = value.y1; }
    if (value.z0 < preInletDomain.z0) { preInletDomain.z0 = value.z0; }
    if (value.z1 > preInletDomain.z1) { preInletDomain.z1 = value.z1; }
  }  
  
  location = preInletDomain.enlarge(1,0);
  location = location.enlarge(1,1);

}


void PreInlet::createBoundary(plb::MultiBlockLattice3D<T,DESCRIPTOR> * fluid, plb::MultiScalarField3D<int> * flagMatrix) {
  plb::Box3D domain = location; 
  plb::Box3D flagBox = flagMatrix->getBoundingBox();
  domain.x0 = domain.x0 < flagBox.x0 ? flagBox.x0 : domain.x0;
  domain.x1 = domain.x1 > flagBox.x1 ? flagBox.x1 : domain.x1;
  domain.y0 = domain.y0 < flagBox.y0 ? flagBox.y0 : domain.y0;
  domain.y1 = domain.y1 > flagBox.y1 ? flagBox.y1 : domain.y1;

  for (int x  = location.x0-1 ; x <= location.x1 ; x++) {
    for (int y  = location.y0-1 ; y <= location.y1 ; y++) {
      for (int z  = location.z0-1 ; z <= location.z1 ; z++) {
        if (partOfpreInlet) {
          if (x == location.x1 || y == location.y1 || x == location.x0-1 || y == location.y0-1) {
            defineDynamics(*fluid,x,y,z, new BounceBack<T,DESCRIPTOR>(1.));
          }
        }
      }
    }
  }

  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
    for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
      for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
        if(flagMatrix->get(x,y,flagMatrix->getBoundingBox().z0+1) == 0) {
          if(partOfpreInlet) {
            defineDynamics(*fluid,x,y,z, new BounceBack<T,DESCRIPTOR>(1.));
          }
        }
      }
    }
  }
  if (partOfpreInlet) {
    fluid->periodicity().toggle(2,true);  
  }
}

PreInlet_old::PreInlet_old(Box3D domain_, string sourceFileName_, int particlePositionTimestep_, Direction flowDir_, HemoCell& hemocell_, bool reducedPrecision_) 
  : hemocell(hemocell_) {
  this->domain = domain_;
  this->sourceFileName = global::directories().getOutputDir() + "/" + sourceFileName_  + ".hdf5";
  this->particlePositionTimeStep = particlePositionTimestep_;
  this->flowDir = flowDir_;
  this->reducedPrecision = reducedPrecision_;
  if (reducedPrecision) {
    this->reducedPrecisionDirection = flowDir/2;
  }
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

  plb::OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
  boundary->setVelocityConditionOnBlockBoundaries(*hemocell.lattice,fluidDomain);
  setBoundaryVelocity(*hemocell.lattice, fluidDomain, plb::Array<T,3>(0.,0.,0.));
  
  hid_t plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_file_id, global::mpi().getGlobalCommunicator(), info);  
      
    file_id = H5Fopen(sourceFileName.c_str(), H5F_ACC_RDONLY, plist_file_id);
    if(file_id < 0) { cout << "Error opening Preinlet hdf5 file." << endl; exit(1); }
  H5Pclose(plist_file_id);
  
  
  //If strange errors here, its the file_id datatype, should be hid_t, but is plint
  Direction stored_flow;
  if(H5LTget_attribute_int(file_id,"/","flowDirection",(int *)&stored_flow) < 0) { 
    cout << "Error reading flowDirection attribute" << endl; exit(1); 
  }
  if (stored_flow != flowDir) {
    cout << "Flow direction does not match, exiting ..." << endl;
    exit(1);
  }
  int stored_size[3];
  if(H5LTget_attribute_int(file_id,"/","boxSize",stored_size) < 0) { 
    cout << "Error reading boxSize attribute" << endl; exit(1); 
  }
  if (domain.getNx() != stored_size[0] ||
      domain.getNy() != stored_size[1] ||
      domain.getNz() != stored_size[2]) {
    cout << "Box dimensions do not match, exiting ..." << endl;
    exit(1);
  }
  if(H5LTget_attribute_int(file_id,"/","Number Of Cells",&nCellsSelf) < 0) { 
    cout << "Error reading number of cells attribute" << endl; exit(1); 
  }
  
  int rp;
  if(H5LTget_attribute_int(file_id,"/","reduced precision",&rp) < 0) { 
    cout << "Error reading reduced precision attribute" << endl; exit(1); 
  }
  if ((rp == 0 && reducedPrecision) || (rp == 1 && !reducedPrecision)) {
    cout << "Dataset is reduced precision, none requested or the other way around" << endl; exit(1);
  }
  
  nCellsOffset = hemocell.cellfields->number_of_cells;
  hemocell.cellfields->number_of_cells += nCellsSelf;
  
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
}

void PreInlet_old::ImmersePreInletParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
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
      particlefield->invalidate_preinlet_ppc();
      particlefield->invalidate_ppc();
      for (int entry:pair.second) {
        particlefield->particles[entry].sv.fromPreInlet = false;
      }
    }
  }
}

void PreInlet_old::DeletePreInletParticles::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HemoCellParticleField * particlefield = dynamic_cast<HemoCellParticleField*>(blocks[0]);
  
  for (HemoCellParticle & particle : particlefield->particles) {
    if (particle.sv.fromPreInlet) {
      particle.tag = 1;
    }
  }
  particlefield->removeParticles(1);
}

void PreInlet_old::PreInletFunctional::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  BlockLattice3D<T,DESCRIPTOR> * fluidfield = dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
  HemoCellParticleField * particlefield = dynamic_cast<HemoCellParticleField*>(blocks[1]);
  
  //When we are also saving particles
  if (parent.hemocell.iter%parent.particlePositionTimeStep==0) {
    const map<int,vector<int>> & ppc = particlefield->get_particles_per_cell();
    for (int i = 0 ; i < parent.particles_size ; i++) {
      particle_hdf5_t * particle = &parent.particles[i];
      hemo::Array<T,3> loc; 
      loc[0] = particle->location[0];
      loc[1] = particle->location[1];
      loc[2] = particle->location[2];
      
      
      if (ppc.find(particle->cell_id) != ppc.end()) {
        if (ppc.at(particle->cell_id)[particle->vertex_id] != -1) {
          if (! particlefield->particles[ppc.at(particle->cell_id)[particle->vertex_id]].sv.fromPreInlet) {
            continue;
          }
        }
      }

      HemoCellParticle particle_h = HemoCellParticle(loc,particle->cell_id,particle->vertex_id,particle->particle_type);
      particle_h.sv.fromPreInlet = true;
      hemo::Array<T,3> vel; 
      vel[0] = particle->velocity[0];
      vel[1] = particle->velocity[1];
      vel[2] = particle->velocity[2];
      particle_h.sv.v = vel;
      particlefield->addParticle(domain,&particle_h);
    }
  }
  
  Box3D fluidBox; 
  Box3D shiftedBox = domain.shift(fluidfield->getLocation().x,fluidfield->getLocation().y,fluidfield->getLocation().z);
  
  if (!intersect(shiftedBox,parent.fluidDomain,fluidBox)) {return;}
  
  hsize_t offset[5]={parent.hemocell.iter%DSET_SLICE,(hsize_t)fluidBox.x0-parent.fluidDomain.x0,(hsize_t)fluidBox.y0-parent.fluidDomain.y0,(hsize_t)fluidBox.z0-parent.fluidDomain.z0,0};
  hsize_t count[5]= {1,(hsize_t)fluidBox.getNx(),(hsize_t)fluidBox.getNy(),(hsize_t)fluidBox.getNz(),3};
  if (parent.reducedPrecision) {
    count[4] = 1;
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

  if(parent.reducedPrecision) {
    static_assert(sizeof(float) == 4, "float size not 32 bits, exiting");
    float *data = new float[count[1]*count[2]*count[3]*count[4]];
    H5Dread(parent.dataset_velocity_id,H5T_IEEE_F32LE,memspace_id,parent.dataspace_velocity_id,H5P_DEFAULT,data);
  
  
    for (plint x = 0 ; x < fluidBox.getNx() ; x++) {
      for (plint y = 0 ; y < fluidBox.getNy() ; y++) {
        for (plint z = 0; z < fluidBox.getNz() ; z++) {
          //Skip boundaries
          /*if (fluidfield->get(xx,yy,zz).getDynamics().isBoundary()) {
            index+=3;
            xx++;
            continue;
          };*/

          vel[0] = 0;
          vel[1] = 0;
          vel[2] = 0;
          vel[parent.reducedPrecisionDirection] = data[index];
          fluidfield->get(xx,yy,zz).defineVelocity(vel);

          index++;
          zz++;
        }
        zz = unshiftedfluidBox.z0;
        yy++;
      }
      yy = unshiftedfluidBox.y0;
      xx++;
    }
    delete[] data;
  } else {
    static_assert(sizeof(double) == 8, "double size not 64 bits, exiting");
    double *data = new double[count[1]*count[2]*count[3]*count[4]];
    H5Dread(parent.dataset_velocity_id,H5T_IEEE_F64LE,memspace_id,parent.dataspace_velocity_id,H5P_DEFAULT,data);
  
  
    for (plint x = 0 ; x < fluidBox.getNx() ; x++) {
      for (plint y = 0 ; y < fluidBox.getNy() ; y++) {
        for (plint z = 0; z < fluidBox.getNz() ; z++) {
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
          zz++;
        }
        zz = unshiftedfluidBox.z0;
        yy++;
      }
      yy = unshiftedfluidBox.y0;
      xx++;
    }
    delete[] data;
  }
  
  H5Sclose(memspace_id);
}
void PreInlet_old::update() {  
 
  //When we are also reading particles
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
  
  //Load in correct dataset for the fluid field
  if ((int)(hemocell.iter/DSET_SLICE) != current_velocity_field) {
    current_velocity_field = hemocell.iter/DSET_SLICE;
    
    if(dataset_velocity_id != -1){
      H5Dclose(dataset_velocity_id);
    }
    string dataset_name = "Velocity Boundary_" + to_string(current_velocity_field * DSET_SLICE);
    
    dataset_velocity_id = H5Dopen(file_id,dataset_name.c_str(),H5P_DEFAULT);
    if (dataset_velocity_id < 0) {
      pcout << "Error opening dataset " << dataset_name << " from preinlet file. Exiting ..." << endl;
      exit(1);
    }
  }
    
  
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell.cellfields->immersedParticles);
  applyProcessingFunctional(new DeletePreInletParticles(),hemocell.cellfields->immersedParticles->getBoundingBox(),wrapper);
  
  
  dataspace_velocity_id = H5Dget_space(dataset_velocity_id);
  
    wrapper.clear();
    wrapper.push_back(hemocell.cellfields->lattice);
    wrapper.push_back(hemocell.cellfields->immersedParticles);
    applyProcessingFunctional(new PreInletFunctional(*this),domain,wrapper);

  H5Sclose(dataspace_velocity_id);
  
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

PreInlet_old::~PreInlet_old() {
  H5Dclose(dataset_velocity_id);
  H5Fclose(file_id);
}

PreInlet_old::PreInletFunctional * PreInlet_old::PreInletFunctional::clone() const { return new PreInletFunctional(*this); }
PreInlet_old::DeletePreInletParticles * PreInlet_old::DeletePreInletParticles::clone() const { return new DeletePreInletParticles(*this); }
PreInlet_old::ImmersePreInletParticles * PreInlet_old::ImmersePreInletParticles::clone() const { return new ImmersePreInletParticles(*this); }
PreInlet::CreatePreInletBoundingBox *        PreInlet::CreatePreInletBoundingBox::clone() const { return new PreInlet::CreatePreInletBoundingBox(*this);}
}
