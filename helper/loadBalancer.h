#ifndef LOADBALANCER_H
#define LOADBALANCER_H
#define HEMO_PARMETIS

class LoadBalancer;

#include "hemocell_internal.h"
#include "hemocell.h"

class LoadBalancer {  
  public:
  LoadBalancer(HemoCell & hemocell);
  double calculateFractionalLoadImbalance();
  void doLoadBalance();
  
  /**
   * Restructure blocks to reduce communication on one processor
   * Set checkpoint_available to false if not called in the same iteration right after doLoadBalance()
   */
  void restructureBlocks(bool checkpoint_available=true);

  //Functionals for gathering data
  struct TOAB_t{
    double fluid_time;
    double particle_time;
    int n_lsp;
    int mpi_proc;
    int n_neighbours;
    int location[3];
    int neighbours[HEMOCELL_MAX_NEIGHBOURS];
  };
  struct Box3D_simple {
    plint x0,x1,y0,y1,z0,z1;
    inline Box3D_simple & operator=(const Box3D & box) {
      x0 = box.x0;
      x1 = box.x1;
      y0 = box.y0;
      y1 = box.y1;
      z0 = box.z0;
      z1 = box.z1;
      return *this;
    }
    inline Box3D getBox3D() const {
      return Box3D(x0,x1,y0,y1,z0,z1);
    }
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
  bool original_block_stored = false;
  SparseBlockStructure3D original_block_structure;
  ThreadAttribution* original_thread_attribution;
};

#endif
