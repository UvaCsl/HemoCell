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
#ifndef FLUID_HDF5_IO_HH
#define FLUID_HDF5_IO_HH

#include "FluidHdf5IO.hh"
#include "palabos3D.h"
#include "palabos3D.hh"

#include <hdf5.h>
#include <hdf5_hl.h>

namespace hemo {
  
void outputHDF5(hsize_t* dim, hsize_t* chunk, hid_t& file_id, string& name, float* output) {
    //We can calulate nvalues through the dims
    chunk[3] = dim[3];
    //unsigned int nvalues = dim[0]*dim[1]*dim[2]*dim[3];

    hid_t sid = H5Screate_simple(4,dim,NULL);
      hid_t plist_id = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk); 
        H5Pset_deflate(plist_id, 7);
        hid_t did = H5Dcreate2(file_id,name.c_str(),H5T_NATIVE_FLOAT,sid,H5P_DEFAULT,plist_id,H5P_DEFAULT);
        H5Dwrite(did,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,output);
      H5Dclose(did);
    H5Sclose(sid);
}

template<template<class U> class DD>
class WriteFluidField : public BoxProcessingFunctional3D
{
public:
 WriteFluidField(HemoCellFields& cellfields_, MultiBlock3D & fluid_, plint iter_, string identifier_, T dx_, T dt_, vector<int> & outputVariables_) :
    cellfields(cellfields_), fluid(fluid_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_),outputVariables(outputVariables_) { }
 
 ~WriteFluidField(){};

  BlockDomain::DomainT appliesTo () const {
    return BlockDomain::bulk;
  }

  WriteFluidField* clone() const {
    return new WriteFluidField(*this);
  }

  void getTypeOfModification ( vector<modif::ModifT>& modified ) const {
    for (pluint i = 0; i < modified.size(); i++) {
      modified[i] = modif::nothing;
    }
  }

  void processGenericBlocks( Box3D domain, vector<AtomicBlock3D*> blocks ) {

    int id = global::mpi().getRank();
    ablock = dynamic_cast<BlockLattice3D<T,DD>*>(blocks[0]);
    particlefield = dynamic_cast<HemoCellParticleField*>(blocks[1]);
    blockid = particlefield->atomicBlockId; //Nasty trick to prevent us from having to overload the fluid field ( palabos domain)

    this->odomain = &domain; //Access for output functions
    if (outputVariables.size() == 0 ) {
      return; //No output needed? ok
    }
    
    if (cellfields.hemocell.partOfpreInlet) {
      identifier += "_PRE";
    }
    std::string fileName = global::directories().getOutputDir() + "/hdf5/" + zeroPadNumber(iter) + '/' + identifier + "."  + zeroPadNumber(iter) + ".p." + to_string(blockid) + ".h5";
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
    for (int outputVariable : outputVariables) {
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
        case OUTPUT_BOUNDARY:
          output = outputBoundary();
          name = "Boundary";
          dim[3] = 1;
        break;
        case OUTPUT_OMEGA:
          output = outputOmega();
          name = "Omega";
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
        case OUTPUT_SHEAR_RATE:
          output = outputShearRate();
          name = "ShearRate";
          dim[3] = 6;
          break;
        case OUTPUT_STRAIN_RATE:
          output = outputStrainRate();
          name = "StrainRate";
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

private:

  float * outputVelocity() {
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

  float * outputForce() {
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

  float * outputDensity() {
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
        output[i] = output[i]*(param::df/(param::dx*param::dx));
      }
    }

    return output;
  }
  
  float * outputBoundary() {
    float * output = new float [(*nCells)];
    unsigned int n = 0;
    for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
      for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
        for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

          if (ablock->get(iX,iY,iZ).getDynamics().isBoundary()) {
            output[n] = 1;
          } else {
            output[n] = 0;
          }
          n++;
        }
      }
    }
    return output;
  }
  
  float * outputOmega() {
    float * output = new float [(*nCells)];
    unsigned int n = 0;
    for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
      for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
        for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

          output[n] = ablock->get(iX,iY,iZ).getDynamics().getOmega();
          n++;
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0 ; i < (*nCells) ; i++) {
        output[i] = output[i]*(param::df/(param::dx*param::dx));
      }
    }

    return output;
  }
    
  float * outputCellDensity(string name) {
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

  float * outputShearStress() {
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
        output[i] = output[i]*(param::df/(param::dx*param::dx));
      }
    }

    return output;
  }
     
  float * outputShearRate() {
    float * output = new float [(*nCells)*6];
    unsigned int n = 0;
    plb::Array<T,6> shearrate;
        
    for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
      for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
        for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {

          ablock->get(iX,iY,iZ).computeShearRate(*ablock,shearrate,iX,iY,iZ);
          //Array<T,6> stress_vec = stress.get(iX,iY,iZ);
          output[n] = shearrate[0];
          output[n+1] = shearrate[1];
          output[n+2] = shearrate[2];
          output[n+3] = shearrate[3];
          output[n+4] = shearrate[4];
          output[n+5] = shearrate[5];
          n += 6;
        }
      }
    }   
     
    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0 ; i < (*nCells)*6 ; i++) {
        output[i] = output[i]*(1/(param::dt));
      }
    }

    return output;
  }
  
  float * outputStrainRate() {
    float * output = new float [(*nCells)*6];
    unsigned int n = 0;
    // calculate tensorfield strain rate
    std::auto_ptr<TensorField3D<T,6> > strainrate (computeStrainRateFromStress(*ablock));
    // calculate norm of tensorfield
    std::auto_ptr<ScalarField3D<T> > shearrate (computeSymmetricTensorNorm(*strainrate));
    
    //strainrate.get()
    
    for (plint iZ=odomain->z0-1; iZ<=odomain->z1+1; ++iZ) {
      for (plint iY=odomain->y0-1; iY<=odomain->y1+1; ++iY) {
        for (plint iX=odomain->x0-1; iX<=odomain->x1+1; ++iX) {
            
            
            if ((iX < fluid.getBoundingBox().x0-2) || (iX > fluid.getBoundingBox().x1+2) ||
	    (iY < fluid.getBoundingBox().y0-2) || (iY > fluid.getBoundingBox().y1+2) ||
	    (iZ < fluid.getBoundingBox().z0-2) || (iZ > fluid.getBoundingBox().z1+2) )
            {   output[n] = 0;
                output[n+1] = 0;
                output[n+2] = 0;
                output[n+3] = 0;
                output[n+4] = 0;
                output[n+5] = 0;
                n += 6;
           }
            else {
    
            output[n] = strainrate->get(iX,iY,iZ)[0];
            output[n+1] = strainrate->get(iX,iY,iZ)[1];
            output[n+2] = strainrate->get(iX,iY,iZ)[2];
            output[n+3] = strainrate->get(iX,iY,iZ)[3];
            output[n+4] = strainrate->get(iX,iY,iZ)[4];
            output[n+5] = strainrate->get(iX,iY,iZ)[5];
            n += 6;
            }
        }
      }
    }   
     
    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0 ; i < (*nCells) ; i++) {
        output[i] = output[i]*(1/(param::dt));
      }
    }

    return output;
  }
   
 
    HemoCellFields& cellfields;
    MultiBlock3D& fluid;
    plint iter;
    string identifier;
    double dx;
    double dt;
    Box3D * odomain;
    BlockLattice3D<T,DD> * ablock;
    HemoCellParticleField * particlefield;
    int blockid;
    hsize_t * nCells;
    vector<int> & outputVariables;
};
}
#endif
