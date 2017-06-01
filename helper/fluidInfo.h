/* 
 * File:   FluidInfo.h
 * Author: vikko
 *
 * Created on May 31, 2017, 12:13 PM
 */

#ifndef FLUIDINFO_H
#define FLUIDINFO_H

#include "hemocell_internal.h"
#include "hemoCellFunctional.h"

struct FluidStatistics {
  double min;
  double max;
  double avg;
  pluint ncells;
};

class GatherFluidVelocity : public HemoCellGatheringFunctional<FluidStatistics> {
  public:  
    using HemoCellGatheringFunctional<FluidStatistics>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherFluidVelocity * clone() const;
};
class GatherFluidForce : public HemoCellGatheringFunctional<FluidStatistics> {
  public:  
    using HemoCellGatheringFunctional<FluidStatistics>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherFluidForce * clone() const;
};
class FluidInfo {
public:
  static FluidStatistics calculateVelocityStatistics(HemoCell * hemocell_);
  static FluidStatistics calculateForceStatistics(HemoCell * hemocell_);
};

#endif /* FLUIDINFO_H */

