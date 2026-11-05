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
#ifndef HEMOCELL_NOOP_H
#define HEMOCELL_NOOP_H
#include "cellMechanics.h"

namespace hemo {
/**
 * @brief Null cell model: no membrane forces.
 *
 * ParticleMechanics() adds nothing, so the vertices only carry the forces of
 * other sources (repulsion) and otherwise follow the fluid like tracers. The
 * cell does not keep its shape. The material file `<name>.xml` must still
 * exist, as for every cell type (the mesh and the common constants are read
 * from it).
 */
class NoOp : public CellMechanics {
public:
  NoOp(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_, modelCfg_), cellField(cellField_) {}

  void ParticleMechanics(std::map<plint, std::vector<HemoCellParticle *>> &, const std::map<plint, bool> &,
                         pluint) override {}
  void statistics() override {
    hlog << "(Cell-mechanics model) NoOp for " << cellField.name << " cellfield: no membrane forces" << std::endl;
  }

private:
  HemoCellField & cellField;
};
}

#endif
