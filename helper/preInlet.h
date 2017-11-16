/* 
 * File:   preInlet.h
 * Author: vikko
 *
 * Created on 14 September 2017, 11:36
 */
#include "hemocell.h"


#ifndef PREINLET_H
#define PREINLET_H

enum Direction : int {
  Xpos, Xneg, Ypos, Yneg, Zpos, Zneg
};

struct particle_hdf5_t {
  float location[3];
  float velocity[3];
  int cell_id;
  int vertex_id;
  unsigned int particle_type;
};

class createPreInlet {
  class createPreInletFunctional: public HemoCellFunctional {
    createPreInlet & parent;
    void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    createPreInletFunctional * clone() const;
  public:
    createPreInletFunctional(createPreInlet & parent_) : parent(parent_) {}
  };
  
  HemoCell & hemocell;
  Box3D domain;
  Box3D fluidDomain;
  string outputFileName;
  int particlePositionTimeStep;
  Direction flowDir;
  int file_id, dataspace_velocity_id,dataset_velocity_id,plist_dataset_collective_id;
  int particle_type_mem, particle_type_h5;
  
public:
  createPreInlet(Box3D domain, string outputFileName, int particlePositionTimestep, Direction flow_direction, HemoCell & hemocell);
  void saveCurrent();
  ~createPreInlet();
};



class PreInlet {
  class PreInletFunctional: public HemoCellFunctional {
    PreInlet & parent;
    void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    PreInletFunctional * clone() const;
  public:
    PreInletFunctional(PreInlet & parent_) : parent(parent_) {}
  };
  HemoCell & hemocell;
  Box3D domain;
  Box3D fluidDomain;
  string sourceFileName;
  int particlePositionTimeStep;
  Direction flowDir;
  int file_id,dataspace_velocity_id,dataset_velocity_id,dataset_particles_id,dataspace_particles_id;
  int particle_type_mem, particle_type_h5, particles_size;
  particle_hdf5_t * particles;


public:  
  PreInlet(Box3D domain, string sourceFileName, int particlePositionTimestep, Direction flow_direction, HemoCell & hemocell);
  void update();
  ~PreInlet();
};

#endif /* PREINLET_H */

