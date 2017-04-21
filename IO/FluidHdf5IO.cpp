#ifndef FLUID_HDF5_IO_CPP
#define FLUID_HDF5_IO_CPP

#include "FluidHdf5IO.h"

void writeFluidField_HDF5(CellFields3D& cellfields, double dx, double dt, plint iter, string preString) {
  WriteFluidField * wff = new WriteFluidField(cellfields, cellfields.lattice,iter,"Fluid",dx,dt);
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(&cellfields.lattice);
  wrapper.push_back(cellfields.immersedParticles); //Needed for the atomicblock id, nothing else
  applyProcessingFunctional(wff,cellfields.lattice.getBoundingBox(),wrapper);
}

BlockDomain::DomainT WriteFluidField::appliesTo () const {
  return BlockDomain::bulk;
}

WriteFluidField* WriteFluidField::clone() const {
  return new WriteFluidField(*this);
}

void WriteFluidField::getTypeOfModification ( vector<modif::ModifT>& modified ) const {
  for (pluint i = 0; i < modified.size(); i++) {
    modified[i] = modif::nothing;
  }
}

WriteFluidField::WriteFluidField(CellFields3D& cellfields_,
                                 MultiBlockLattice3D<double,DESCRIPTOR>& fluid_,
                                 plint iter_, string identifier_, 
                                 double dx_, double dt_) :
cellfields(cellfields_), fluid(fluid_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_)
{ }

void WriteFluidField::processGenericBlocks( Box3D domain, vector<AtomicBlock3D*> blocks ) {

  int id = global::mpi().getRank();
	ablock = dynamic_cast<BlockLattice3D<double,DESCRIPTOR>*>(blocks[0]);
  int blockid = dynamic_cast<HemoParticleField3D*>(blocks[1])->atomicBlockId; //Nasty trick to prevent us from having to overload the fluid field ( palabos domain)

  this->odomain = &domain; //Access for output functions
  if (cellfields.desiredFluidOutputVariables.size() == 0 ) {
    return; //No output needed? ok
  }

  std::string fileName = global::directories().getOutputDir() + "/hdf5/" + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.", blockid,3) + ".h5";
  hid_t file_id;
  file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
	H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
	long int iterHDF5=iter;
	H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
	H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

  hsize_t Nx = domain.x1 - domain.x0;
  hsize_t Ny = domain.y1 - domain.y0;
  hsize_t Nz = domain.z1 - domain.z0;
  int subdomainSize[]  = {int(Nx), int(Ny), int(Nz)};
  hsize_t nCells = Nx*Ny*Nz;
  this->nCells = &nCells;
  Dot3D rp_temp = blocks[0]->getLocation();
  long int relativePosition[3] = {rp_temp.x,rp_temp.y,rp_temp.z};

	H5LTset_attribute_int (file_id, "/", "subdomainSize", subdomainSize, 3);
  H5LTset_attribute_long(file_id, "/", "relativePosition", relativePosition, 3);

  //Also compute chunking here
  hsize_t chunk[4];
  chunk[0] = 1000 < Nx ? 1000 : Nx;
  chunk[1] = 1000 < Ny ? 1000 : Ny;
  chunk[2] = 1000 < Nz ? 1000 : Nz;

  //I could do fancy schmancy function pointer lookup like with the particles, but it
  //takes time, just do it here
  for (int outputVariable : cellfields.desiredFluidOutputVariables) {
    //These variables should be set by the functions:
    float * output = 0;
    hsize_t dim[] = {Nx,Ny,Nz,0};
    string name;
    
    switch(outputVariable) {
      case OUTPUT_VELOCITY:
        outputVelocity(output);
        name = "Velocity";
        dim[3] = 3;
        break;
      default:
        continue;
    }

    //We can calulate nvalues through the dims
    chunk[3] = dim[3];
    unsigned int nvalues = dim[0]*dim[1]*dim[2]*dim[3];

    int sid = H5Screate_simple(4,dim,NULL);
      int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,name.c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,output);
      H5Dclose(did);
    H5Sclose(sid);

    delete[] output;

  }
  H5Fclose(file_id);
}

void WriteFluidField::outputVelocity(float * output) {
  output = new float [(*nCells)*3];
  unsigned int n = 0;
  Array<double,3> vel;
  for (plint iX=odomain->x0; iX<odomain->x1; ++iX) {
    for (plint iY=odomain->y0; iY<odomain->y1; ++iY) {
      for (plint iZ=odomain->z0; iZ<odomain->z1; ++iZ) {

        ablock->grid[iX][iY][iZ].computeVelocity(vel);
        output[n] = vel[0];
        output[n+1] = vel[1];
        output[n+2] = vel[2];
        n += 3;
      }
    }
  }
}
#endif
