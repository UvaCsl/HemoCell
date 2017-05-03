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

  struct TOAB_t{
    double fluid_time;
    double particle_time;
    int mpi_proc;
  };
  class GatherTimeOfAtomicBlocks : public HemoCellGatheringFunctional<TOAB_t> {
  public:  
    using HemoCellGatheringFunctional<TOAB_t>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherTimeOfAtomicBlocks * clone() const;
  };
  private:
  HemoCell & hemocell;
};

#endif
