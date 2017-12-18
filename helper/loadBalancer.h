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
#ifndef LOADBALANCER_H
#define LOADBALANCER_H

class LoadBalancer;

#include "hemocell_internal.h"
#include "hemocell.h"

class LoadBalancer {  
  public:
  LoadBalancer(HemoCell & hemocell_) : hemocell(hemocell_), original_block_structure(0,0,0) { }
#ifdef HEMO_PARMETIS
  double calculateFractionalLoadImbalance();
  /**
   * Restructure blocks to reduce communication on one processor
   * Set checkpoint_available to false if not called in the same iteration right after doLoadBalance()
   */
  void restructureBlocks(bool checkpoint_available=true);
#else
  double calculateFractionalLoadImbalance() {return 0.0;}
  void restructureBlocks(bool checkpoint_available=true) {}
#endif

  void doLoadBalance();

  /**
   * used to reload a checkpoint, but first reload the config file
   */
  void reloadCheckpoint();
  
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
