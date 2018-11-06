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
#ifndef PREINLET_H
#define PREINLET_H

class PreInlet;
#include "config/constant_defaults.h"
#include "core/hemoCellFunctional.h"

#include "core/geometry3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "config.h"

#ifndef HEMOCELL_H
namespace hemo {
  class HemoCell;
}
#endif

enum Direction : int {
  Xpos, Xneg, Ypos, Yneg, Zpos, Zneg
};
  

#define DSET_SLICE 1000

namespace hemo {


inline plint cellsInBoundingBox(plb::Box3D const & box) {
  return abs((box.x1 - box.x0)*(box.y1-box.y0)*(box.z1-box.z0));
}


class PreInlet {
public:
    class CreateDrivingForceFunctional: public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CreateDrivingForceFunctional * clone() const;
    plint & count;
    public:
      CreateDrivingForceFunctional(plint & count_) :
                                count(count_) {count = 0;}
  };
  class CreatePreInletBoundingBox: public HemoCellFunctional {
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CreatePreInletBoundingBox * clone() const;
    plb::Box3D & boundingBox;
    bool & foundPreInlet;
    public:
      CreatePreInletBoundingBox(plb::Box3D & b_, bool & fp_) :
                                boundingBox(b_), foundPreInlet(fp_) {}
  };
  
  PreInlet(plb::MultiScalarField3D<int> & flagMatrix, hemo::HemoCell * hemocell_);
  inline plint getNumberOfNodes() { return cellsInBoundingBox(location);}
  void createBoundary(plb::MultiBlockLattice3D<T,DESCRIPTOR> *,plb::MultiScalarField3D<int> * flagMatrix);
  void setDrivingForce();
  void calculateDrivingForce();
  plb::Box3D location;
  plb::Box3D fluidInlet;
  int nProcs = 0;
  bool initialized = false;
  double drivingForce = 0.0;
  std::map<plint,plint> BlockToMpi;
  std::map<int,plint> particleSendMpi; // the value equals the number of atomic blocks that are sent.
  std::map<int,bool> particleReceiveMpi;
  std::vector<plint> communicationBlocks;
  int sendingBlocks = 0;
  bool partOfpreInlet = false;
  int inflow_length = 0;
  int preinlet_length = 0;
  bool communications_mapped = false;
  std::vector<int> particle_receivers;
  std::vector<int> particle_senders;
  HemoCell * hemocell;
};

}

  #include "hemocell.h"

namespace hemo {
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
  plint file_id, dataspace_velocity_id,dataset_velocity_id = -1,plist_dataset_collective_id;
  plint particle_type_mem, particle_type_h5;
  bool reducedPrecision;
  int reducedPrecisionDirection;
  int desired_iterations;
  //ofstream counter;
  int current_velocity_field = -1;  
public:
  createPreInlet(Box3D domain, string outputFileName, int particlePositionTimestep, Direction flow_direction, HemoCell & hemocell, int desired_iterations, bool reducedPrecision_ = false);
  void saveCurrent();
  ~createPreInlet();
};

void mapPreInletParticleBoundary(HemoCell & hemocell);
void createPreInletVelocityBoundary(plb::MultiBlockLattice3D<T,DESCRIPTOR> * fluid, plb::MultiScalarField3D<int> * flagmatrix,plb::Array<double,3> speed, HemoCell & hemocell);
void applyPreInletVelocityBoundary(HemoCell & hemocell);
void applyPreInletParticleBoundary(HemoCell & hemocell);

}
#endif /* PREINLET_H */

