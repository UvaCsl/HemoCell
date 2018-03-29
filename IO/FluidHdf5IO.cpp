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
#include "FluidHdf5IO.h"
#include "hemocell.h"

#include <hdf5.h>
#include <hdf5_hl.h>

void outputHDF5(hsize_t* dim, hsize_t* chunk, hid_t& file_id, string& name, float* output) {
    //We can calulate nvalues through the dims
    chunk[3] = dim[3];
    //unsigned int nvalues = dim[0]*dim[1]*dim[2]*dim[3];

    int sid = H5Screate_simple(4,dim,NULL);
      int plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk); 
        H5Pset_deflate(plist_id, 7);
        int did = H5Dcreate2(file_id,name.c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,output);
      H5Dclose(did);
    H5Sclose(sid);
}

void writeFluidField_HDF5(HemoCellFields& cellfields, T dx, T dt, plint iter, string preString) {
  if(std::find(cellfields.desiredFluidOutputVariables.begin(), cellfields.desiredFluidOutputVariables.end(), OUTPUT_FORCE) != cellfields.desiredFluidOutputVariables.end()) {
    if(LOG_LEVEL >= 2)
      pcout << "(FluidOutput) (OutputForce) The force on the fluid field is reset to zero, If there is a bodyforce, reset it after this function" << endl; 
    cellfields.spreadParticleForce();
  }
  WriteFluidField * wff = new WriteFluidField(cellfields, *cellfields.lattice,iter,"Fluid",dx,dt);
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfields.lattice);
  wrapper.push_back(cellfields.immersedParticles); //Needed for the atomicblock id, nothing else
  applyProcessingFunctional(wff,cellfields.lattice->getBoundingBox(),wrapper);
  if(std::find(cellfields.desiredFluidOutputVariables.begin(), cellfields.desiredFluidOutputVariables.end(), OUTPUT_FORCE) != cellfields.desiredFluidOutputVariables.end()) {
    // Reset Forces on the lattice, TODO do own efficient implementation
    setExternalVector(*cellfields.hemocell.lattice, (*cellfields.hemocell.lattice).getBoundingBox(),
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
  }
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

WriteFluidField::WriteFluidField(HemoCellFields& cellfields_,
                                 MultiBlockLattice3D<T,DESCRIPTOR>& fluid_,
                                 plint iter_, string identifier_, 
                                 T dx_, T dt_) :
cellfields(cellfields_), fluid(fluid_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_)
{ }

void WriteFluidField::processGenericBlocks( Box3D domain, vector<AtomicBlock3D*> blocks ) {

  int id = global::mpi().getRank();
  ablock = dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
  particlefield = dynamic_cast<HemoCellParticleField*>(blocks[1]);
  blockid = particlefield->atomicBlockId; //Nasty trick to prevent us from having to overload the fluid field ( palabos domain)

  this->odomain = &domain; //Access for output functions
  if (cellfields.desiredFluidOutputVariables.size() == 0 ) {
    return; //No output needed? ok
  }

  std::string fileName = global::directories().getOutputDir() + "/hdf5/" + zeroPadNumber(iter) + '/'  + createFileName((identifier+".").c_str(),iter,8) + createFileName(".p.", blockid,3) + ".h5";
  hid_t file_id;
  file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	H5LTset_attribute_double (file_id, "/", "dx", &dx, 1);
	H5LTset_attribute_double (file_id, "/", "dt", &dt, 1);
	long int iterHDF5=iter;
	H5LTset_attribute_long (file_id, "/", "iteration", &iterHDF5, 1);
	H5LTset_attribute_int (file_id, "/", "processorId", &id, 1);

  hsize_t Nx = domain.x1 - domain.x0+1 +2;//+1 for = and <=, +2 for an envelope of 1 on each side for paraview
  hsize_t Ny = domain.y1 - domain.y0+1 +2;
  hsize_t Nz = domain.z1 - domain.z0+1 +2;
  hsize_t nCells = Nx*Ny*Nz;
  this->nCells = &nCells;
  int ncells = Nx*Ny*Nz;
  Dot3D rp_temp = blocks[0]->getLocation();
  int subdomainSize[]  = {int(Nz), int(Ny), int(Nx)}; //Reverse for paraview
  float dxdydz[3] = {1.,1.,1.};
  float relativePosition[3] = {float(rp_temp.z+domain.z0-1.5),
                               float(rp_temp.y+domain.y0-1.5),
                               float(rp_temp.x+domain.x0-1.5)}; //Reverse for paraview

  if (cellfields.hemocell.outputInSiUnits) {
    relativePosition[0] *= param::dx;
    relativePosition[1] *= param::dx;
    relativePosition[2] *= param::dx;
    dxdydz[0] = param::dx;
    dxdydz[1] = param::dx;
    dxdydz[2] = param::dx;
  }
  
  H5LTset_attribute_int (file_id, "/", "numberOfCells", &ncells, 1);
  H5LTset_attribute_int (file_id, "/", "subdomainSize", subdomainSize, 3);
  H5LTset_attribute_float(file_id, "/", "relativePosition", relativePosition, 3);
  H5LTset_attribute_float(file_id,"/","dxdydz",dxdydz,3);

  //Also compute chunking here
  hsize_t chunk[4];
  chunk[2] = 1000 < Nx ? 1000 : Nx;
  chunk[1] = 1000 < Ny ? 1000 : Ny;
  chunk[0] = 1000 < Nz ? 1000 : Nz;

  //I could do fancy schmancy function pointer lookup like with the particles, but it
  //takes time, just do it here
  for (int outputVariable : cellfields.desiredFluidOutputVariables) {
    //These variables should be set by the functions:
    float * output = 0;
    hsize_t dim[4] = {Nz,Ny,Nx,0};
    string name;
    
    switch(outputVariable) {
      case OUTPUT_VELOCITY:
        output = outputVelocity();
        name = "Velocity";
        dim[3] = 3;
        break;
      case OUTPUT_FORCE:
        output = outputForce();
        name = "Force";
        dim[3] = 3;
        break;
      case OUTPUT_DENSITY:
        output = outputDensity();
        name = "Density";
        dim[3] = 1;
        break;
      case OUTPUT_CELL_DENSITY:
        for (unsigned int i = 0 ; i < cellfields.size() ; i++) {
          output = outputCellDensity(cellfields[i]->name);
          name = "CellDensity_" + cellfields[i]->name;
          dim[3] = 1;
          outputHDF5(dim,chunk,file_id,name,output);
          delete[] output;
 	}
	continue;
      case OUTPUT_SHEAR_STRESS:
	output = outputShearStress();
	name = "ShearStress";
        dim[3] = 6;
        break;
      default:
        continue;
    }


    outputHDF5(dim,chunk,file_id,name,output);
    delete[] output;

  }
  H5Fclose(file_id);
}

float * WriteFluidField::outputVelocity() {
  float * output = new float [(*nCells)*3];
  unsigned int n = 0;
  plb::Array<T,3> vel;
  for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
    for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
      for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

        ablock->get(iX,iY,iZ).computeVelocity(vel);
        output[n] = vel[0];
        output[n+1] = vel[1];
        output[n+2] = vel[2];
        n += 3;
      }
    }
  }
  
  if (cellfields.hemocell.outputInSiUnits) {
    for (unsigned int i = 0 ; i < (*nCells)*3 ; i++) {
      output[i] = output[i]*param::dx/param::dt;
    }
  }
  
  return output;
}

float * WriteFluidField::outputForce() {
  float * output = new float [(*nCells)*3];
  unsigned int n = 0;
  hemo::Array<T,3> vel;
  for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
    for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
      for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

        output[n] = ablock->get(iX,iY,iZ).external.data[0];
        output[n+1] = ablock->get(iX,iY,iZ).external.data[1];
        output[n+2] = ablock->get(iX,iY,iZ).external.data[2];
        n += 3;
      }
    }
  }
  
  if (cellfields.hemocell.outputInSiUnits) {
    for (unsigned int i = 0 ; i < (*nCells)*3 ; i++) {
      output[i] = output[i]*param::df;
    }
  }
  
  return output;
}

float * WriteFluidField::outputDensity() {
  float * output = new float [(*nCells)];
  unsigned int n = 0;
  for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
    for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
      for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

        output[n] = ablock->get(iX,iY,iZ).computeDensity();
        n++;
      }
    }
  }
  
  if (cellfields.hemocell.outputInSiUnits) {
    for (unsigned int i = 0 ; i < (*nCells) ; i++) {
      output[i] = output[i]*param::dm;
    }
  }
  
  return output;
}

float * WriteFluidField::outputCellDensity(string name) {
  float * output = new float [(*nCells)];
  memset(output, 0, sizeof(float)*(*nCells));
  
  vector<HemoCellParticle*> found;
  particlefield->findParticles(particlefield->localDomain,found,cellfields[name]->ctype);

  int Ystride = ((odomain->x1-odomain->x0)+3);
  int Zstride = Ystride*((odomain->y1-odomain->y0)+3);

  for (HemoCellParticle * particle : found) {
    plint iX,iY,iZ;
    //Coordinates are relative
    const Dot3D tmpDot = ablock->getLocation(); 
    iX = plint((particle->sv.position[0]-tmpDot.x)+0.5);
    iY = plint((particle->sv.position[1]-tmpDot.y)+0.5);
    iZ = plint((particle->sv.position[2]-tmpDot.z)+0.5);

    output[(iX)+(iY)*Ystride+(iZ)*Zstride] += 1;
  }
    
  if (cellfields.hemocell.outputInSiUnits) {
    for (unsigned int i = 0 ; i < (*nCells) ; i++) {
      output[i] *= cellfields[name]->volumeFractionOfLspPerNode;
    }
  }
  
  return output;
}

float * WriteFluidField::outputShearStress() {
  float * output = new float [(*nCells)*6];
  unsigned int n = 0;
  plb::Array<T,6> stress;
  for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
    for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
      for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

        ablock->get(iX,iY,iZ).computeShearStress(stress);
	//Array<T,6> stress_vec = stress.get(iX,iY,iZ);
	output[n] = stress[0];
	output[n+1] = stress[1];
	output[n+2] = stress[2];
	output[n+3] = stress[3];
	output[n+4] = stress[4];
	output[n+5] = stress[5];
	n += 6;
      }
    }
  }

  if (cellfields.hemocell.outputInSiUnits) {
    for (unsigned int i = 0 ; i < (*nCells)*6 ; i++) {
      output[i] = output[i]*param::dm*((param::dx*param::dx)/(param::dt*param::dt));
    }
  }

  return output;
}

