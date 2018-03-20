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
#ifndef HEMOCELLSTRETCH_H
#define HEMOCELLSTRETCH_H

#include "hemocell_internal.h"
#include "hemoCellFields.h"
#include "hemoCellFunctional.h"

class HemoCellStretch {
  public:

  HemoCellStretch(HemoCellField & cellfield_, unsigned int n_forced_lsps_, T external_force_);
  void applyForce();
  
  class FindForcedLsps: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   FindForcedLsps * clone() const;
  };
  
  class ForceForcedLsps: public HemoCellFunctional {
   void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
   ForceForcedLsps * clone() const;
  };
  
  //Vertex id of lsps
  static vector<plint> lower_lsps;
  static vector<plint> upper_lsps;
  
  HemoCellField & cellfield;
  static unsigned int n_forced_lsps;
  static T external_force;
};

#endif
