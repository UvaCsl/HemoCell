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
class NoOp : public CellMechanics {
  public:
  NoOp() :CellMechanics() {};
  NoOp(Config & cfg, HemoCellField & cellfield) :CellMechanics() {};


  inline void ParticleMechanics(map<int,vector<HemoCellParticle *>>,map<int,bool>, pluint ctype) {} ;
  inline void statistics () {
    cerr << "Mechanical model is NoOp";
  }
};

#endif
