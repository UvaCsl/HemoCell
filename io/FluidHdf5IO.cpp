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

using namespace std;
using namespace plb;

namespace hemo {

namespace {
// The force output spreads the particle forces onto the lattice to show the
// total force on the fluid. The lattice's external force must be identical
// afterwards (e.g. a body force set by the case before writeOutput()), so the
// external field of every node of every local atomic block (bulk and
// envelope) is copied out before and back in after the output. This exact
// copy only costs memory and time while output is written.
typedef std::vector<std::vector<T>> ExternalForceBackup;

void saveExternalForce(plb::MultiBlockLattice3D<T, DESCRIPTOR> & lattice, ExternalForceBackup & backup) {
  const int nScalars = DESCRIPTOR<T>::ExternalField::numScalars;
  const std::vector<plint> & blocks = lattice.getLocalInfo().getBlocks();
  backup.resize(blocks.size());
  for (pluint iBlock = 0; iBlock < blocks.size(); iBlock++) {
    plb::BlockLattice3D<T, DESCRIPTOR> & block = lattice.getComponent(blocks[iBlock]);
    std::vector<T> & saved = backup[iBlock];
    saved.resize(block.getNx() * block.getNy() * block.getNz() * nScalars);
    pluint i = 0;
    for (plint x = 0; x < block.getNx(); x++) {
      for (plint y = 0; y < block.getNy(); y++) {
        for (plint z = 0; z < block.getNz(); z++) {
          const T * external = block.get(x, y, z).getExternal(0);
          for (int s = 0; s < nScalars; s++) { saved[i++] = external[s]; }
        }
      }
    }
  }
}

void restoreExternalForce(plb::MultiBlockLattice3D<T, DESCRIPTOR> & lattice, const ExternalForceBackup & backup) {
  const int nScalars = DESCRIPTOR<T>::ExternalField::numScalars;
  const std::vector<plint> & blocks = lattice.getLocalInfo().getBlocks();
  for (pluint iBlock = 0; iBlock < blocks.size(); iBlock++) {
    plb::BlockLattice3D<T, DESCRIPTOR> & block = lattice.getComponent(blocks[iBlock]);
    const std::vector<T> & saved = backup[iBlock];
    pluint i = 0;
    for (plint x = 0; x < block.getNx(); x++) {
      for (plint y = 0; y < block.getNy(); y++) {
        for (plint z = 0; z < block.getNz(); z++) {
          T * external = block.get(x, y, z).getExternal(0);
          for (int s = 0; s < nScalars; s++) { external[s] = saved[i++]; }
        }
      }
    }
  }
}
// Time of an output step in the XDMF files: seconds, or iterations for lattice-unit output.
double outputTime(HemoCellFields & cellfields, plint iter, T dt) {
  return cellfields.hemocell.outputInSiUnits ? iter * dt : double(iter);
}
} // namespace

void writeCEPACField_HDF5(HemoCellFields & cellfields, T dx, T dt, plint iter, string preString) {
  global.statistics.getCurrent()["writeCEPACField"].start();

  // The CEPAC "density" is a concentration, written as it is.
  std::vector<XdmfBlockFile> written;
  WriteFluidField<plb::descriptors::AdvectionDiffusionD3Q19Descriptor> * wff =
      new WriteFluidField<plb::descriptors::AdvectionDiffusionD3Q19Descriptor>(
          cellfields, *cellfields.CEPACfield, iter, "CEPAC", dx, dt, cellfields.desiredCEPACfieldOutputVariables, 1.,
          &written);
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(cellfields.CEPACfield);
  wrapper.push_back(cellfields.immersedParticles); //Needed for the atomicblock id, nothing else
  applyProcessingFunctional(wff, cellfields.CEPACfield->getBoundingBox(), wrapper);
  writeXdmf(written, iter, outputTime(cellfields, iter, dt));

  global.statistics.getCurrent().stop();
}
void writeFluidField_HDF5(HemoCellFields & cellfields, T dx, T dt, plint iter, string preString) {
  global.statistics.getCurrent()["writeFluidField"].start();

  const bool outputForce =
      std::find(cellfields.desiredFluidOutputVariables.begin(), cellfields.desiredFluidOutputVariables.end(),
                OUTPUT_FORCE) != cellfields.desiredFluidOutputVariables.end();
  ExternalForceBackup externalForce;
  if (outputForce) {
    // Show the total force (external/body force + particle forces). The
    // lattice force is restored afterwards, so a body force set by the case
    // before writeOutput() survives and output does not alter the simulation.
    saveExternalForce(*cellfields.lattice, externalForce);
    cellfields.spreadParticleForce();
  }
  // Density in kg/m^3: the lattice density times the mass per lattice volume (rhoP).
  const T densityToSI = param::dm / (param::dx * param::dx * param::dx);
  std::vector<XdmfBlockFile> written;
  WriteFluidField<DESCRIPTOR> * wff =
      new WriteFluidField<DESCRIPTOR>(cellfields, *cellfields.lattice, iter, "Fluid", dx, dt,
                                      cellfields.desiredFluidOutputVariables, densityToSI, &written);
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(cellfields.lattice);
  wrapper.push_back(cellfields.immersedParticles); //Needed for the atomicblock id, nothing else
  applyProcessingFunctional(wff, cellfields.lattice->getBoundingBox(), wrapper);
  if (outputForce) { restoreExternalForce(*cellfields.lattice, externalForce); }
  writeXdmf(written, iter, outputTime(cellfields, iter, dt));

  global.statistics.getCurrent().stop();
}

}
