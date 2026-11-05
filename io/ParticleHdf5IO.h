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
#ifndef PARTICLE_HDF5IO_H
#define PARTICLE_HDF5IO_H

#include "hemoCellParticle.h"
#include "hemoCellFields.h"
#include "hemoCellField.h"
#include "xdmfOutput.h"

namespace hemo {
/**
 * @brief Write the selected outputs of every cell type to HDF5 files (one per
 * atomic block) and describe them in XDMF files (see writeXdmf()). Collective.
 *
 * Every cell is written once, by the block that holds its vertex 0 (see
 * HemoCellParticleField::getOutputCells()).
 *
 * @param[in] cellFields All cell types.
 * @param[in] dx Lattice spacing in metres (stored as attribute).
 * @param[in] dt Time step in seconds (stored as attribute; gives the time of the step).
 * @param[in] iter Iteration of the output step.
 * @param[in] preString Prefix of the output names.
 */
void writeCellField3D_HDF5(HemoCellFields & cellFields, T dx, T dt, plint iter, std::string preString = "");


class WriteCellField3DInMultipleHDF5Files : public BoxProcessingFunctional3D {
public:
  /// `written_` receives one record per block file written; it is owned by the
  /// caller and shared by all clones of the functional.
  WriteCellField3DInMultipleHDF5Files(HemoCellField & cellField3D_, plint iter_, std::string identifier_, T dx_, T dt_,
                                      std::vector<XdmfBlockFile> * written_);
  /// Arguments: [0] Particle-field.
  virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
  virtual WriteCellField3DInMultipleHDF5Files * clone() const;
  virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const;
  virtual BlockDomain::DomainT appliesTo() const;

private:
  HemoCellField & cellField3D;
  plint iter;
  std::string identifier;
  double dx;
  double dt;
  std::vector<XdmfBlockFile> * written;
};
}
#endif // PARTICLE_HDF5IO_H
