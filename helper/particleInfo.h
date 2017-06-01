/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particleInfo.h
 * Author: vikko
 *
 * Created on 01 June 2017, 13:53
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

