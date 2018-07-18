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
#ifndef HEMOCELL_WBCHIGHORDERMODEL_H
#define HEMOCELL_WBCHIGHORDERMODEL_H

#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"
namespace hemo {
class WbcHighOrderModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const T k_volume;
  const T k_area;
  const T k_link;
  const T k_bend;
  const T eta_m;
  const T k_inner_rigid;
  const T k_cytoskeleton;
  const T core_radius;
  const T radius;

  public:
  WbcHighOrderModel(Config & modelCfg_, HemoCellField & cellField_) ;

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> &lpc, size_t ctype) ;

  void statistics();

  static T calculate_coreRadius(Config & cfg );
  static T calculate_radius(Config & cfg );
  static T calculate_kInnerRigid(Config & cfg );
  static T calculate_kCytoskeleton(Config & cfg );
};
}
#endif
