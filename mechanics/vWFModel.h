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
#ifndef HEMOCELL_VWFMODEL_H
#define HEMOCELL_VFWMODEL_H

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "hemoCellFields.h"

class vWFModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const double k_link;
  const double k_bend;
  const double link_eq;
  const double bend_eq;
  const double eta_m;
  const double eta_v;

  public:
  vWFModel(Config & modelCfg_, HemoCellField & cellField_) ;

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> &lpc, pluint ctype) ;

  void statistics();

  static double calculate_kBend(Config & cfg);
  static double calculate_kLink(Config & cfg);
  static double calculate_bend_eq(Config & cfg);
  static double calculate_link_eq(Config & cfg);
  static double calculate_etaM(Config & cfg );
  static double calculate_etaV(Config & cfg );

};

#endif
