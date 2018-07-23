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
#include "hemocell.h"
#include "FluidHdf5IO.hh"

#include "dataProcessors/dataInitializerWrapper3D.hh"
#include "dataProcessors/dataInitializerFunctional3D.hh"

namespace hemo {

void writeCEPACField_HDF5(HemoCellFields& cellfields, T dx, T dt, plint iter, string preString) {
  global.statistics.getCurrent()["writeCEPACField"].start();

  WriteFluidField<plb::descriptors::AdvectionDiffusionD3Q19Descriptor> * wff = new WriteFluidField<plb::descriptors::AdvectionDiffusionD3Q19Descriptor>(cellfields, *cellfields.CEPACfield,iter,"CEPAC",dx,dt,cellfields.desiredCEPACfieldOutputVariables);
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfields.CEPACfield);
  wrapper.push_back(cellfields.immersedParticles); //Needed for the atomicblock id, nothing else
  applyProcessingFunctional(wff,cellfields.CEPACfield->getBoundingBox(),wrapper);
  
  global.statistics.getCurrent().stop();
}
void writeFluidField_HDF5(HemoCellFields& cellfields, T dx, T dt, plint iter, string preString) {
  global.statistics.getCurrent()["writeFluidField"].start();

  if(std::find(cellfields.desiredFluidOutputVariables.begin(), cellfields.desiredFluidOutputVariables.end(), OUTPUT_FORCE) != cellfields.desiredFluidOutputVariables.end()) {
    hlogfile << "(FluidOutput) (OutputForce) The force on the fluid field is reset to zero, If there is a bodyforce, reset it after this output function (FluidField write force, OUTPUT_FORCE)" << endl; 
    cellfields.spreadParticleForce();
  }
  WriteFluidField<DESCRIPTOR> * wff = new WriteFluidField<DESCRIPTOR>(cellfields, *cellfields.lattice,iter,"Fluid",dx,dt,cellfields.desiredFluidOutputVariables);
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(cellfields.lattice);
  wrapper.push_back(cellfields.immersedParticles); //Needed for the atomicblock id, nothing else
  applyProcessingFunctional(wff,cellfields.lattice->getBoundingBox(),wrapper);
  if(std::find(cellfields.desiredFluidOutputVariables.begin(), cellfields.desiredFluidOutputVariables.end(), OUTPUT_FORCE) != cellfields.desiredFluidOutputVariables.end()) {
    // Reset Forces on the lattice, TODO do own efficient implementation
    plb::setExternalVector(*cellfields.hemocell.lattice, (*cellfields.hemocell.lattice).getBoundingBox(),
          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
          plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
  }
  
  global.statistics.getCurrent().stop();
}

}
