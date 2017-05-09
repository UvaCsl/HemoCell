#ifndef LOADBALANCER_H
#define LOADBALANCER_H

class LoadBalancer;

#include "hemocell_internal.h"
#include "hemocell.h"
#include <parmetis.h>

class LoadBalancer {
  public:
  LoadBalancer(HemoCell & hemocell);
  double calculateFractionalLoadImbalance();
  void doLoadBalance();

  //Functionals for gathering data
  struct TOAB_t{
    double fluid_time;
    double particle_time;
    int n_lsp;
    int mpi_proc;
    int n_neighbours;
    int neighbours[HEMOCELL_MAX_NEIGHBOURS];
  };
  class GatherTimeOfAtomicBlocks : public HemoCellGatheringFunctional<TOAB_t> {
  public:  
    using HemoCellGatheringFunctional<TOAB_t>::HemoCellGatheringFunctional; //Inherit Constructor
    void processGenericBlocks(Box3D, vector<AtomicBlock3D*>);
    GatherTimeOfAtomicBlocks * clone() const;
  };
  private:
  bool FLI_iscalled = false;
  map<int,TOAB_t> gatherValues;
  HemoCell & hemocell;
};

#endif
