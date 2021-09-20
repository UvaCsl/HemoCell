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
#ifndef HEMOCELL_PLTSIMPLEMODEL_H
#define HEMOCELL_PLTSIMPLEMODEL_H

#include "config.h"
#include "cellMechanics.h"
#include "constant_defaults.h"
#include "hemoCellField.h"

namespace hemo {
class PltSimpleModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const T k_volume;
  const T k_area;
  const T k_link;
  const T k_bend;
  const T eta_m;
  
  //Constructor
  public:
  PltSimpleModel(Config & modelCfg_, HemoCellField & cellField_);

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> &particles_per_cell, const map<int,bool> &lpc, pluint ctype);
#ifdef SOLIDIFY_MECHANICS
  void solidifyMechanics(const std::map<int,std::vector<int>>&,std::vector<HemoCellParticle>&,plb::BlockLattice3D<T,DESCRIPTOR> *,plb::BlockLattice3D<T,CEPAC_DESCRIPTOR> *, pluint ctype, HemoCellParticleField&);
#endif
  void statistics();

};
}
#endif
