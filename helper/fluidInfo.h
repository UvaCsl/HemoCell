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
#ifndef FLUIDINFO_H
#define FLUIDINFO_H

#include "hemocell.h"
#include "hemoCellFunctional.h"
#include <vector>

namespace hemo {

struct FluidStatistics {
  double min;
  double max;
  double avg;
  pluint ncells;
};

class GatherFluidVelocity : public HemoCellGatheringFunctional<FluidStatistics> {
  public:  
    using HemoCellGatheringFunctional<FluidStatistics>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    GatherFluidVelocity * clone() const;
};
class GatherFluidForce : public HemoCellGatheringFunctional<FluidStatistics> {
  public:  
    using HemoCellGatheringFunctional<FluidStatistics>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    GatherFluidForce * clone() const;
};
class FluidInfo {
public:
  static FluidStatistics calculateVelocityStatistics(HemoCell * hemocell_);
  static FluidStatistics calculateForceStatistics(HemoCell * hemocell_);
};

}
#endif /* FLUIDINFO_H */

