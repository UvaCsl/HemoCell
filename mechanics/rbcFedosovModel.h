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
#ifndef HEMOCELL_RBCFEDOSOVMODEL_H
#define HEMOCELL_RBCFEDOSOVMODEL_H

#include "config.h"
#include "cellMechanics.h"
#include "hemoCellField.h"
#include "constant_defaults.h"

namespace hemo {

class RbcFedosovModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;

  const T eta_m;
  const T kd;
  const T ka;
  const T kv;
  const T kb;
  const T x0;
  const T m;

  public:
  RbcFedosovModel(Config & modelCfg_, HemoCellField & cellField_);

  private:
  void ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> &lpc, size_t ctype) ;
  void statistics();

};
}
#endif
