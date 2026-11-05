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

#include "hdf5Output.h"
#include "xdmfOutput.h"
#include "palabos3D.h"
#include "palabos3D.hh"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <cmath>

namespace hemo {

// Write a 4D dataset (z, y, x, component), one chunk per block (up to 1000 nodes per direction).
inline void outputHDF5(hsize_t * dim, hsize_t * chunk, hid_t & file_id, string & name, float * output) {
  chunk[3] = dim[3];
  writeCompressedDataset(file_id, name, H5T_NATIVE_FLOAT, 4, dim, chunk, output);
}

/**
 * Writes the selected fields of one lattice (fluid or CEPAC) to one HDF5 file
 * per atomic block: the values of the block's bulk nodes, `(z, y, x,
 * component)`. Arguments: [0] lattice, [1] particle field (for the block ID).
 */
template <template <class U> class DD>
class WriteFluidField : public BoxProcessingFunctional3D {
public:
  /// `densityToSI_` converts the lattice density to SI output (1 for a
  /// concentration). `written_` receives one record per block file written;
  /// it is owned by the caller and shared by all clones of the functional.
  WriteFluidField(HemoCellFields & cellfields_, MultiBlock3D & fluid_, plint iter_, string identifier_, T dx_, T dt_,
                  vector<int> & outputVariables_, T densityToSI_, std::vector<XdmfBlockFile> * written_)
      : cellfields(cellfields_), fluid(fluid_), iter(iter_), identifier(identifier_), dx(dx_), dt(dt_),
        densityToSI(densityToSI_), written(written_), outputVariables(outputVariables_) {}

  ~WriteFluidField(){};

  BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

  WriteFluidField * clone() const { return new WriteFluidField(*this); }

  void getTypeOfModification(vector<modif::ModifT> & modified) const {
    for (pluint i = 0; i < modified.size(); i++) { modified[i] = modif::nothing; }
  }

  void processGenericBlocks(Box3D domain, vector<AtomicBlock3D *> blocks) {
    int id = global::mpi().getRank();
    ablock = dynamic_cast<BlockLattice3D<T, DD> *>(blocks[0]);
    particlefield = dynamic_cast<HemoCellParticleField *>(blocks[1]);
    blockid = particlefield
                  ->atomicBlockId; //Nasty trick to prevent us from having to overload the fluid field ( palabos domain)

    // Output functions access the domain and size through these pointers; they
    // point to members, not to this function's locals, so they stay valid.
    this->odomainStorage = domain;
    this->odomain = &this->odomainStorage;
    if (outputVariables.size() == 0) {
      return; //No output needed? ok
    }

    XdmfBlockFile record;
    record.identifier = identifier;
    if (cellfields.hemocell.partOfpreInlet) { record.identifier += "_PRE"; }
    record.blockId = blockid;
    record.structured = true;
    record.file = "hdf5/" + zeroPadNumber(iter) + '/' + record.identifier + "." + zeroPadNumber(iter) + ".p." +
                  to_string(blockid) + ".h5";
    std::string fileName = global::directories().getOutputDir() + "/" + record.file;
    hid_t file_id;
    file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5LTset_attribute_double(file_id, "/", "dx", &dx, 1);
    H5LTset_attribute_double(file_id, "/", "dt", &dt, 1);
    long int iterHDF5 = iter;
    H5LTset_attribute_long(file_id, "/", "iteration", &iterHDF5, 1);
    H5LTset_attribute_int(file_id, "/", "processorId", &id, 1);

    // The bulk nodes only: neighbouring blocks write disjoint nodes. Each node
    // is the centre of a cell of size dx, so the blocks tile the domain.
    hsize_t Nx = domain.x1 - domain.x0 + 1;
    hsize_t Ny = domain.y1 - domain.y0 + 1;
    hsize_t Nz = domain.z1 - domain.z0 + 1;
    hsize_t nCells = Nx * Ny * Nz;
    this->nCellsStorage = nCells;
    this->nCells = &this->nCellsStorage;
    int ncells = Nx * Ny * Nz;
    Dot3D rp_temp = blocks[0]->getLocation();
    int subdomainSize[] = {int(Nz), int(Ny), int(Nx)}; //Reverse for paraview
    // Lower corner of the first node's cell, z, y, x.
    const double scale = cellfields.hemocell.outputInSiUnits ? double(param::dx) : 1.;
    const double corner[3] = {(rp_temp.z + domain.z0 - 0.5) * scale, (rp_temp.y + domain.y0 - 0.5) * scale,
                              (rp_temp.x + domain.x0 - 0.5) * scale};
    float dxdydz[3] = {float(scale), float(scale), float(scale)};
    float relativePosition[3] = {float(corner[0]), float(corner[1]), float(corner[2])};
    record.nodes[0] = Nz;
    record.nodes[1] = Ny;
    record.nodes[2] = Nx;
    for (int i = 0; i < 3; i++) {
      record.origin[i] = corner[i];
      record.spacing[i] = scale;
    }

    H5LTset_attribute_int(file_id, "/", "numberOfCells", &ncells, 1);
    H5LTset_attribute_int(file_id, "/", "subdomainSize", subdomainSize, 3);
    H5LTset_attribute_float(file_id, "/", "relativePosition", relativePosition, 3);
    H5LTset_attribute_float(file_id, "/", "dxdydz", dxdydz, 3);

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
      hsize_t dim[4] = {Nz, Ny, Nx, 0};
      string name;

      switch (outputVariable) {
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
        case OUTPUT_BINDING_SITES:
          output = outputBindingSites();
          name = "BindingSites";
          dim[3] = 1;
          break;
        case OUTPUT_INTERIOR_POINTS:
          output = outputInteriorPoints();
          name = "InteriorPoints";
          dim[3] = 1;
          break;
        case OUTPUT_OMEGA:
          output = outputOmega();
          name = "Omega";
          dim[3] = 1;
          break;
        case OUTPUT_CELL_DENSITY:
          for (unsigned int i = 0; i < cellfields.size(); i++) {
            output = outputCellDensity(cellfields[i]->name);
            name = "CellDensity_" + cellfields[i]->getIdentifier();
            dim[3] = 1;
            outputHDF5(dim, chunk, file_id, name, output);
            record.attributes.push_back(attribute(name, 1));
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
          dim[3] = 9;
          break;
        case OUTPUT_STRAIN_RATE:
          output = outputStrainRate();
          name = "StrainRate";
          dim[3] = 6;
          break;

        default:

          continue;
      }


      outputHDF5(dim, chunk, file_id, name, output);
      record.attributes.push_back(attribute(name, dim[3]));
      delete[] output;
    }
    H5Fclose(file_id);
    written->push_back(record);
  }

private:
  static XdmfAttribute attribute(const string & name, hsize_t components) {
    XdmfAttribute a;
    a.name = name;
    a.components = int(components);
    return a;
  }

  float * outputVelocity() {
    float * output = new float[(*nCells) * 3];
    unsigned int n = 0;
    plb::Array<T, 3> vel;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          ablock->get(iX, iY, iZ).computeVelocity(vel);
          output[n] = vel[0];
          output[n + 1] = vel[1];
          output[n + 2] = vel[2];
          n += 3;
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells) * 3; i++) { output[i] = output[i] * param::dx / param::dt; }
    }

    return output;
  }

  float * outputForce() {
    float * output = new float[(*nCells) * 3];
    unsigned int n = 0;
    hemo::Array<T, 3> vel;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          output[n] = ablock->get(iX, iY, iZ).external.data[0];
          output[n + 1] = ablock->get(iX, iY, iZ).external.data[1];
          output[n + 2] = ablock->get(iX, iY, iZ).external.data[2];
          n += 3;
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells) * 3; i++) { output[i] = output[i] * param::df; }
    }

    return output;
  }

  float * outputDensity() {
    float * output = new float[(*nCells)];
    unsigned int n = 0;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          output[n] = ablock->get(iX, iY, iZ).computeDensity();
          n++;
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells); i++) { output[i] = output[i] * densityToSI; }
    }

    return output;
  }

  float * outputBoundary() {
    float * output = new float[(*nCells)];
    unsigned int n = 0;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          if (ablock->get(iX, iY, iZ).getDynamics().isBoundary()) {
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

  float * outputBindingSites() {
    float * output = new float[(*nCells)];
    if (!particlefield->bindingField) {
      pcout << "(FluidHdf5) (Error) OUTPUT_BINDING_SITES requested, but binding sites not used, outputting a zero field"
            << endl;
      for (hsize_t i = 0; i < *nCells; i++) { output[i] = 0; }
      return output;
    }
    unsigned int n = 0;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          if (particlefield->bindingField->get(iX, iY, iZ)) {
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

  float * outputInteriorPoints() {
    float * output = new float[(*nCells)];
    if (!particlefield->interiorViscosityField) {
      pcout
          << "(FluidHdf5) (Error) OUTPUT_INTERIOR_POINTS requested, but interior viscosity not used, outputting a zero field"
          << endl;
      for (hsize_t i = 0; i < *nCells; i++) { output[i] = 0; }
      return output;
    }
    unsigned int n = 0;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          output[n] = (particlefield->interiorViscosityField->get(iX, iY, iZ));
          n++;
        }
      }
    }
    return output;
  }

  float * outputOmega() {
    float * output = new float[(*nCells)];
    unsigned int n = 0;
    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          output[n] = ablock->get(iX, iY, iZ).getDynamics().getOmega();
          n++;
        }
      }
    }
    // omega = 1/tau is dimensionless: the same in SI and lattice units.

    return output;
  }

  float * outputCellDensity(string name) {
    float * output = new float[(*nCells)];
    memset(output, 0, sizeof(float) * (*nCells));

    vector<HemoCellParticle *> found;
    particlefield->findParticles(particlefield->localDomain, found, cellfields[name]->ctype);

    const plint nx = odomain->x1 - odomain->x0 + 1;
    const plint ny = odomain->y1 - odomain->y0 + 1;
    const plint nz = odomain->z1 - odomain->z0 + 1;
    const Dot3D location = ablock->getLocation();

    for (HemoCellParticle * particle : found) {
      // Node whose cell (half a node in each direction, upper face included,
      // as in isContainedABS()) holds the particle, relative to the first
      // bulk node of this block.
      const plint iX = plint(std::ceil(particle->sv.position[0] - location.x - 0.5)) - odomain->x0;
      const plint iY = plint(std::ceil(particle->sv.position[1] - location.y - 0.5)) - odomain->y0;
      const plint iZ = plint(std::ceil(particle->sv.position[2] - location.z - 0.5)) - odomain->z0;
      if (iX < 0 || iX >= nx || iY < 0 || iY >= ny || iZ < 0 || iZ >= nz) { continue; }
      output[iX + nx * (iY + ny * iZ)] += 1;
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells); i++) { output[i] *= cellfields[name]->volumeFractionOfLspPerNode; }
    }

    return output;
  }

  float * outputShearStress() {
    float * output = new float[(*nCells) * 6];
    unsigned int n = 0;
    plb::Array<T, 6> stress;

    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          ablock->get(iX, iY, iZ).computeShearStress(stress);
          //Array<T,6> stress_vec = stress.get(iX,iY,iZ);
          output[n] = stress[0];
          output[n + 1] = stress[1];
          output[n + 2] = stress[2];
          output[n + 3] = stress[3];
          output[n + 4] = stress[4];
          output[n + 5] = stress[5];
          n += 6;
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells) * 6; i++) {
        output[i] = output[i] * (param::df / (param::dx * param::dx));
      }
    }

    return output;
  }

  float * outputShearRate() {
    float * output = new float[(*nCells) * 9];
    unsigned int n = 0;
    plb::Array<T, 9> shearrate;
    plb::Array<T, 3> velp1;
    plb::Array<T, 3> veln1;
    plb::Array<T, 3> velp2;
    plb::Array<T, 3> veln2;
    plb::Array<T, 3> velp3;
    plb::Array<T, 3> veln3;

    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          // velocity gradients in y-direction
          ablock->get(iX + 1, iY, iZ).computeVelocity(velp1); //vel at cell+1
          ablock->get(iX - 1, iY, iZ).computeVelocity(veln1); //vel at cell+1

          // velocity gradients in y-direction
          ablock->get(iX, iY + 1, iZ).computeVelocity(velp2); //vel at cell+1
          ablock->get(iX, iY - 1, iZ).computeVelocity(veln2); //vel at cell+1

          // velocity gradients in z-direction
          ablock->get(iX, iY, iZ + 1).computeVelocity(velp3); //vel at cell+1
          ablock->get(iX, iY, iZ - 1).computeVelocity(veln3); //vel at cell+1

          shearrate[0] = (velp1[0] - veln1[0]) / 2; // dv_x /2*dx
          shearrate[3] = (velp1[1] - veln1[1]) / 2; // dv_y/ 2*dx
          shearrate[6] = (velp1[2] - veln1[2]) / 2; // dv_z /2*dx
          shearrate[1] = (velp2[0] - veln2[0]) / 2; // dv_x / 2*dy
          shearrate[4] = (velp2[1] - veln2[1]) / 2; // dv_y / 2*dy
          shearrate[7] = (velp2[2] - veln2[2]) / 2; // dv_z / 2*dy
          shearrate[2] = (velp3[0] - veln3[0]) / 2; // dv_x / 2*dz
          shearrate[5] = (velp3[1] - veln3[1]) / 2; // dv_y / 2*dz
          shearrate[8] = (velp3[2] - veln3[2]) / 2; // dv_z / 2*dz

          output[n] = shearrate[0];
          output[n + 1] = shearrate[1];
          output[n + 2] = shearrate[2];
          output[n + 3] = shearrate[3];
          output[n + 4] = shearrate[4];
          output[n + 5] = shearrate[5];
          output[n + 6] = shearrate[6];
          output[n + 7] = shearrate[7];
          output[n + 8] = shearrate[8];
          n += 9;
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells) * 9; i++) { output[i] = output[i] * (1 / (param::dt)); }
    }

    return output;
  }

  float * outputStrainRate() {
    float * output = new float[(*nCells) * 6];
    unsigned int n = 0;
    // calculate tensorfield strain rate
    std::unique_ptr<TensorField3D<T, 6>> strainrate(computeStrainRateFromStress(*ablock));
    // calculate norm of tensorfield
    std::unique_ptr<ScalarField3D<T>> shearrate(computeSymmetricTensorNorm(*strainrate));

    //strainrate.get()

    for (plint iZ = odomain->z0; iZ <= odomain->z1; ++iZ) {
      for (plint iY = odomain->y0; iY <= odomain->y1; ++iY) {
        for (plint iX = odomain->x0; iX <= odomain->x1; ++iX) {
          if ((iX < fluid.getBoundingBox().x0 - 2) || (iX > fluid.getBoundingBox().x1 + 2) ||
              (iY < fluid.getBoundingBox().y0 - 2) || (iY > fluid.getBoundingBox().y1 + 2) ||
              (iZ < fluid.getBoundingBox().z0 - 2) || (iZ > fluid.getBoundingBox().z1 + 2)) {
            output[n] = 0;
            output[n + 1] = 0;
            output[n + 2] = 0;
            output[n + 3] = 0;
            output[n + 4] = 0;
            output[n + 5] = 0;
            n += 6;
          } else {
            output[n] = strainrate->get(iX, iY, iZ)[0];
            output[n + 1] = strainrate->get(iX, iY, iZ)[1];
            output[n + 2] = strainrate->get(iX, iY, iZ)[2];
            output[n + 3] = strainrate->get(iX, iY, iZ)[3];
            output[n + 4] = strainrate->get(iX, iY, iZ)[4];
            output[n + 5] = strainrate->get(iX, iY, iZ)[5];
            n += 6;
          }
        }
      }
    }

    if (cellfields.hemocell.outputInSiUnits) {
      for (unsigned int i = 0; i < (*nCells) * 6; i++) { output[i] = output[i] * (1 / (param::dt)); }
    }

    return output;
  }


  HemoCellFields & cellfields;
  MultiBlock3D & fluid;
  plint iter;
  string identifier;
  double dx;
  double dt;
  T densityToSI;
  std::vector<XdmfBlockFile> * written;
  Box3D * odomain;
  Box3D odomainStorage;
  BlockLattice3D<T, DD> * ablock;
  HemoCellParticleField * particlefield;
  int blockid;
  hsize_t * nCells;
  hsize_t nCellsStorage = 0;
  vector<int> & outputVariables;
};
}
#endif
