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
#ifndef PARTICLEINFO_H
#define PARTICLEINFO_H

#include "hemocell_internal.h"
#include "hemoCellFunctional.h"

struct ParticleStatistics {
  double min;
  double max;
  double avg;
  pluint ncells;
};

class GatherParticleVelocity : public HemoCellGatheringFunctional<ParticleStatistics> {
  public:  
    using HemoCellGatheringFunctional<ParticleStatistics>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherParticleVelocity * clone() const;
};
class GatherParticleForce : public HemoCellGatheringFunctional<ParticleStatistics> {
  public:  
    using HemoCellGatheringFunctional<ParticleStatistics>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherParticleForce * clone() const;
};

class ParticleInfo {
public:
  static ParticleStatistics calculateVelocityStatistics(HemoCell * hemocell_);
  static ParticleStatistics calculateForceStatistics(HemoCell * hemocell_);
};
#endif /* PARTICLEINFO_H */

