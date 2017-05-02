#ifndef LOADBALANCER_H
#define LOADBALANCER_H

class LoadBalancer;

#include "hemocell_internal.h"
#include "hemocell.h"

class LoadBalancer {
  public:
  LoadBalancer(HemoCell & hemocell);
  double calculateFractionalLoadImbalance();
  void doLoadBalance();

  //Functionals for gathering data
  class GatherNumberOfLSPs : public HemoCellGatheringFunctional<int> {
  public:  
    using HemoCellGatheringFunctional<int>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherNumberOfLSPs * clone() const;
  };
  private:
  HemoCell & hemocell;
};

#endif
