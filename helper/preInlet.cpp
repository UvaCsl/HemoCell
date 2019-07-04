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

#include "palabos3D.h"
#include "palabos3D.hh"
#include "preInlet.h"

#include "boundaryCondition/boundaryInstantiator3D.h"
#include "hemoCellFields.h"

namespace hemo {
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
}

#include <hdf5.h>
#include <hdf5_hl.h>


namespace hemo {

void PreInlet::initializePreInletParticleBoundary() {
  MPI_Barrier(MPI_COMM_WORLD);
  Box3D domain = fluidInlet;
  if (hemocell->partOfpreInlet) {
    if (direction == Direction::Xneg) {
      domain.x0 = location.x0;
    }
    if (direction == Direction::Yneg) {
      domain.y0 = location.y0;
    }
    if (direction == Direction::Zneg) {
      domain.z0 = location.z0;
    }
    if (direction == Direction::Xpos) {
      domain.x1 = location.x1;
    }
    if (direction == Direction::Ypos) {
      domain.y1 = location.y1;
    }
    if (direction == Direction::Zpos) {
      domain.z1 = location.z1;
    }
  }
  if (direction == Direction::Xneg) {
    domain.x0 += 1;
    domain.x1 = domain.x0 + inflow_length;
  }
  if (direction == Direction::Yneg) {
    domain.y0 += 1;
    domain.y1 = domain.y0 + inflow_length;
  }
  if (direction == Direction::Zneg) {
    domain.z0 += 1;
    domain.z1 = domain.z0 + inflow_length;
  }
  if (direction == Direction::Xpos) {
    domain.x1 -= 1;
    domain.x0 = domain.x1 - inflow_length;
  }
  if (direction == Direction::Ypos) {
    domain.y1 -= 1;
    domain.y0 = domain.y1 - inflow_length;
  }
  if (direction == Direction::Zpos) {
    domain.z1 -= 1;
    domain.z0 = domain.z1 - inflow_length;
  }
  
 
  
  for (int bId : hemocell->lattice->getLocalInfo().getBlocks()) {
    Box3D bulk = hemocell->lattice->getMultiBlockManagement().getBulk(bId);
    Box3D result;
    if (!intersect(domain,bulk,result)) { continue; }
    if (hemocell->partOfpreInlet) {
      particleSendMpi[global::mpi().getRank()] += 1;
    } else {
      particleReceiveMpi[global::mpi().getRank()] = true;
    }
    communicationBlocks.push_back(bId);
  }
  
  HemoCellGatheringFunctional<plint>::gather(particleSendMpi);
  HemoCellGatheringFunctional<bool>::gather(particleReceiveMpi);
  for ( auto & pair : particleSendMpi) {
    sendingBlocks += pair.second;
  }
}

void PreInlet::applyPreInletParticleBoundary() {
  global.statistics.getCurrent()["applyPreInletParticleBoundary"].start();
  vector<MPI_Request> requests;
  vector<vector<char>> buffers;
  MPI_Barrier(MPI_COMM_WORLD);
  if (partOfpreInlet) {
    if(particleSendMpi.find(global::mpi().getRank()) != particleSendMpi.end()) {
      for (plint bid : communicationBlocks) {
        buffers.push_back(vector<char>());
        Box3D domain = fluidInlet;
        if (direction == Direction::Xneg) {
          domain.x0 = domain.x0 - preinlet_length;
          domain.x1 = domain.x0 + inflow_length;
        }
        if (direction == Direction::Yneg) {
          domain.y0 = domain.y0 - preinlet_length;
          domain.y1 = domain.y0 + inflow_length;       
        }
        if (direction == Direction::Zneg) {
          domain.z0 = domain.z0 - preinlet_length;
          domain.z1 = domain.z0 + inflow_length;
        }
        if (direction == Direction::Xpos) {
          domain.x1 = domain.x1 + preinlet_length;
          domain.x0 = domain.x1 - inflow_length;       
        }
        if (direction == Direction::Ypos) {
          domain.y1 = domain.y1 + preinlet_length;
          domain.y0 = domain.y1 - inflow_length;        
        }
        if (direction == Direction::Zpos) {
          domain.z1 = domain.z1 + preinlet_length;
          domain.z0 = domain.z1 - inflow_length;       
        }

        Dot3D shift = hemocell->cellfields->immersedParticles->getComponent(bid).getLocation();
        domain = domain.shift(-shift.x,-shift.y,-shift.z);
        hemocell->cellfields->immersedParticles->getComponent(bid).particleDataTransfer.send(domain,buffers.back(),modif::hemocell);
        for (auto & pair : particleReceiveMpi) {
          const int & pid = pair.first;
          requests.push_back(MPI_Request());
          MPI_Isend(&buffers.back()[0],buffers.back().size(),MPI_CHAR,pid,0,MPI_COMM_WORLD,&requests.back());
        }
      }
    }
    MPI_Waitall(requests.size(),&requests[0],MPI_STATUSES_IGNORE);
  } else {
    if(particleReceiveMpi.find(global::mpi().getRank()) != particleReceiveMpi.end()) {
      for (int counter = 0 ; counter < sendingBlocks ; counter++) {
        MPI_Status status;
        int count;
        MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status,MPI_CHAR,&count);
        buffers.push_back(vector<char>(count));
        MPI_Recv(&buffers.back()[0],count,MPI_CHAR,status.MPI_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        Dot3D offset(0,0,0);
        if (direction == Direction::Xneg) {
          offset.x = preinlet_length;
        }
        if (direction == Direction::Yneg) {
          offset.y = preinlet_length;
        }
        if (direction == Direction::Zneg) {
          offset.z = preinlet_length;
        }
        if (direction == Direction::Xpos) {
          offset.x = -preinlet_length;
        }
        if (direction == Direction::Ypos) {
          offset.y = -preinlet_length;
        }
        if (direction == Direction::Zpos) {
          offset.z = -preinlet_length;
        }
        for (int bId : hemocell->cellfields->immersedParticles->getLocalInfo().getBlocks()) {
  
          hemocell->cellfields->immersedParticles->getComponent(bId).particleDataTransfer.receive(&buffers.back()[0],buffers.back().size(),modif::hemocell,offset);
          hemocell->cellfields->immersedParticles->getComponent(bId).invalidate_ppc();
          hemocell->cellfields->immersedParticles->getComponent(bId).invalidate_lpc();
          hemocell->cellfields->immersedParticles->getComponent(bId).invalidate_pg();
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  global.statistics.getCurrent().stop();
}
void PreInlet::applyPreInletVelocityBoundary() {
  global.statistics.getCurrent()["applyPreInletVelocityBoundary"].start();
  MPI_Barrier(MPI_COMM_WORLD);
  vector<MPI_Request> reqs;
  Box3D & domain = fluidInlet;
  Box3D result;
  int z_o = domain.getNz();
  int y_o = domain.getNy();
  int tag;
  plb::Array<T,3> vel;
  for (int bId : hemocell->lattice->getLocalInfo().getBlocks()) {
    Box3D bulk = hemocell->lattice->getMultiBlockManagement().getBulk(bId);
    if (!intersect(domain,bulk,result)) { continue; }
  
    const Dot3D & loc = hemocell->lattice->getComponent(bId).getLocation();
    result = result.shift(-loc.x,-loc.y,-loc.z);
    
    for (int x  = result.x0 ; x <= result.x1 ; x++) {
     for (int y  = result.y0 ; y <= result.y1 ; y++) {
      for (int z  = result.z0 ; z <= result.z1 ; z++) {
        tag = ((z_o * ((z + loc.z) - domain.z0)) + ( y + loc.y ) - domain.y0) * y_o + x + loc.x - domain.x0;
        if (!hemocell->lattice->get(x,y,z).getDynamics().isBoundary()) {
          if (hemocell->partOfpreInlet) {
            hemocell->lattice->getComponent(bId).get(x,y,z).computeVelocity(vel);
            int dest = hemocell->domain_lattice_management->getThreadAttribution().getMpiProcess(hemocell->domain_lattice_management->getSparseBlockStructure().locate(x+loc.x,y+loc.y,z+loc.z));
            reqs.emplace_back();
            MPI_Isend(&vel[0],3*sizeof(T),MPI_CHAR,dest,tag,MPI_COMM_WORLD,&reqs.back());
          } else {
            Box3D point(x,x,y,y,z,z);
            int source = hemocell->preinlet_lattice_management->getThreadAttribution().getMpiProcess(hemocell->preinlet_lattice_management->getSparseBlockStructure().locate(x+loc.x,y+loc.y,z+loc.z));
            MPI_Recv(&vel[0],3*sizeof(T),MPI_CHAR,source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            setBoundaryVelocity(hemocell->lattice->getComponent(bId),point,vel);
          }
        }
      }
     }
    }
  }
  MPI_Waitall(reqs.size(),&reqs[0],MPI_STATUS_IGNORE);
  global.statistics.getCurrent().stop();
}
void PreInlet::initializePreInletVelocityBoundary() {
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* bc =
        createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
  Box3D domain;
  intersect(flagMatrix->getBoundingBox(),location,domain);
  if (direction == Direction::Xneg) {
    domain.x0 = domain.x1;
  }
  if (direction == Direction::Yneg) {
    domain.y0 = domain.y1;
  }
  if (direction == Direction::Zneg) {
    domain.z0 = domain.z1; 
  }
  if (direction == Direction::Xpos) {
    domain.x1 = domain.x0;
  }
  if (direction == Direction::Ypos) {
    domain.y1 = domain.y0;
  }
  if (direction == Direction::Zpos) {
    domain.z1 = domain.z0; 
  }
  fluidInlet = domain;
  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
   for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
    for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
      if (flagMatrix->get(x,y,z) == 1 && !hemocell->partOfpreInlet) {
        Box3D point(x,x,y,y,z,z);
        if (direction == Direction::Xneg) {
          bc->addVelocityBoundary0N(point,*hemocell->lattice);
        }
        if (direction == Direction::Yneg) {
          bc->addVelocityBoundary1N(point,*hemocell->lattice);
        }
        if (direction == Direction::Zneg) {
          bc->addVelocityBoundary2N(point,*hemocell->lattice);
        }
        if (direction == Direction::Xpos) {
          bc->addVelocityBoundary0P(point,*hemocell->lattice);
        }
        if (direction == Direction::Ypos) {
          bc->addVelocityBoundary1P(point,*hemocell->lattice);
        }
        if (direction == Direction::Zpos) {
          bc->addVelocityBoundary2P(point,*hemocell->lattice);
        }
        setBoundaryVelocity(*hemocell->lattice, point, {0.,0.,0.} );
      }
    }
   }
  }
}

void PreInlet::preInletFromSlice(Direction direction_, Box3D boundary) {
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* bc =
        createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
  direction = direction_;
  initialized = true;

  switch(direction) {
    case Direction::Xneg:
    case Direction::Xpos:
      if (boundary.x0 != boundary.x1) {
        hlog << "Not a flat slice, refusing to create preInlet" << endl;
        exit(1);
      }
      break;
    case Direction::Yneg:
    case Direction::Ypos:
      if (boundary.y0 != boundary.y1) {
        hlog << "Not a flat slice, refusing to create preInlet" << endl;
        exit(1);
      }
      break;
    case Direction::Zneg:
    case Direction::Zpos:
      if (boundary.z0 != boundary.z1) {
        hlog << "Not a flat slice, refusing to create preInlet" << endl;
        exit(1);
      }
      break;
    default:
      hlog << "This should not happen" << endl;
      exit(1);
      break;
  }
  
  inflow_length = (*hemocell->cfg)["domain"]["particleEnvelope"].read<int>();
  Box3D preInletDomain;
  bool foundPreInlet = false;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(flagMatrix);
  
  applyProcessingFunctional(new CreatePreInletBoundingBox(preInletDomain,foundPreInlet),boundary,wrapper);
  
  map<int,Box3D_simple> values;
  if (foundPreInlet) {
    values[global::mpi().getRank()] = preInletDomain;
  }
  HemoCellGatheringFunctional<Box3D_simple>::gather(values);
  
  if (values.size() == 0) {
    hlog << "(PreInlet) no preinlet found, is it in the correct location?" << endl;
    exit(1);
  }
  preInletDomain =  values.begin()->second.getBox3D();
  for (auto & pair : values ) {
    Box3D value  = pair.second.getBox3D();
    if (value.x0 < preInletDomain.x0) { preInletDomain.x0 = value.x0; }
    if (value.x1 > preInletDomain.x1) { preInletDomain.x1 = value.x1; }
    if (value.y0 < preInletDomain.y0) { preInletDomain.y0 = value.y0; }
    if (value.y1 > preInletDomain.y1) { preInletDomain.y1 = value.y1; }
    if (value.z0 < preInletDomain.z0) { preInletDomain.z0 = value.z0; }
    if (value.z1 > preInletDomain.z1) { preInletDomain.z1 = value.z1; }
  }  
  
  location = preInletDomain.enlarge(1);
  
  if (direction == Direction::Xneg) {
    location.x0 -= preinlet_length;
  } 
  if (direction == Direction::Yneg) {
    location.y0 -= preinlet_length;
  } 
  if (direction == Direction::Zneg) {
    location.z0 -= preinlet_length;
  }
  if (direction == Direction::Xpos) {
    location.x1 += preinlet_length;
  } 
  if (direction == Direction::Ypos) {
    location.y1 += preinlet_length;
  } 
  if (direction == Direction::Zpos) {
    location.z1 += preinlet_length;
  }
  Box3D system = hemocell->lattice->getBoundingBox();
  if (direction == Direction::Xneg || direction == Direction::Xpos) {
    if (location.y0 < system.y0 ) { location.y0 =system.y0; }
    if (location.y1 > system.y1 ) { location.y1 =system.y1; }
    if (location.z0 < system.z0 ) { location.z0 =system.z0; }
    if (location.z1 > system.z1 ) { location.z1 =system.z1; }
  }
  if (direction == Direction::Yneg || direction == Direction::Ypos) {
    if (location.x0 < system.x0 ) { location.x0 =system.x0; }
    if (location.x1 > system.x1 ) { location.x1 =system.x1; }
    if (location.z0 < system.z0 ) { location.z0 =system.z0; }
    if (location.z1 > system.z1 ) { location.z1 =system.z1; }
  }
  if (direction == Direction::Zneg || direction == Direction::Zpos) {
    if (location.y0 < system.y0 ) { location.y0 =system.y0; }
    if (location.y1 > system.y1 ) { location.y1 =system.y1; }
    if (location.x0 < system.x0 ) { location.x0 =system.x0; }
    if (location.x1 > system.x1 ) { location.x1 =system.x1; }
  }
}

void PreInlet::CreatePreInletBoundingBox::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  ScalarField3D<int> * sf = dynamic_cast<ScalarField3D<int>*>(blocks[0]);
  
  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
    for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
      for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
        if (sf->get(x,y,z)) {
          int xx = x + sf->getLocation().x;
          int yy = y + sf->getLocation().y;
          int zz = z + sf->getLocation().z;

          if(!foundPreInlet) {
            boundingBox.x0 = boundingBox.x1 = xx;
            boundingBox.y0 = boundingBox.y1 = yy;
            boundingBox.z0 = boundingBox.z1 = zz;
            foundPreInlet = true;
          } else {
            if (xx < boundingBox.x0) { boundingBox.x0 = xx; }
            if (xx > boundingBox.x1) { boundingBox.x1 = xx; }
            if (yy < boundingBox.y0) { boundingBox.y0 = yy; }
            if (yy > boundingBox.y1) { boundingBox.y1 = yy; }
            if (zz < boundingBox.z0) { boundingBox.z0 = zz; }
            if (zz > boundingBox.z1) { boundingBox.z1 = zz; }
          }
        }
      }
    }
  }
}

void PreInlet::autoPreinletFromBoundary(Direction dir_) {
  direction = dir_;
  initialized = true;
  Box3D fluidDomain = flagMatrix->getMultiBlockManagement().getBoundingBox();
 
  if (direction == Direction::Xneg) {
    fluidDomain.x1 = fluidDomain.x0+1;
    fluidDomain.x0 = fluidDomain.x1;
  } 
  if (direction == Direction::Yneg) {
    fluidDomain.y1 = fluidDomain.y0+1;
    fluidDomain.y0 = fluidDomain.y1;
  } 
  if (direction == Direction::Zneg) {
    fluidDomain.z1 = fluidDomain.z0+1;
    fluidDomain.z0 = fluidDomain.z1;
  }
  if (direction == Direction::Xpos) {
    fluidDomain.x0 = fluidDomain.x1-1;
    fluidDomain.x1 = fluidDomain.x0;
  } 
  if (direction == Direction::Ypos) {
    fluidDomain.y0 = fluidDomain.y1-1;
    fluidDomain.y1 = fluidDomain.y0;
  } 
  if (direction == Direction::Zpos) {
    fluidDomain.z0 = fluidDomain.z1-1;
    fluidDomain.z1 = fluidDomain.z0;
  }
  
  inflow_length = (*hemocell->cfg)["domain"]["particleEnvelope"].read<int>();
  Box3D preInletDomain;
  bool foundPreInlet = false;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(flagMatrix);
  
  applyProcessingFunctional(new CreatePreInletBoundingBox(preInletDomain,foundPreInlet),fluidDomain,wrapper);
  
  map<int,Box3D_simple> values;
  if (foundPreInlet) {
    values[global::mpi().getRank()] = preInletDomain;
  }
  HemoCellGatheringFunctional<Box3D_simple>::gather(values);
  
  if (values.size() == 0) {
    hlog << "(PreInlet) no preinlet found, does fluid domain extend to the wall?" << endl;
    exit(1);
  }
  preInletDomain =  values.begin()->second.getBox3D();
  for (auto & pair : values ) {
    Box3D value  = pair.second.getBox3D();
    if (value.x0 < preInletDomain.x0) { preInletDomain.x0 = value.x0; }
    if (value.x1 > preInletDomain.x1) { preInletDomain.x1 = value.x1; }
    if (value.y0 < preInletDomain.y0) { preInletDomain.y0 = value.y0; }
    if (value.y1 > preInletDomain.y1) { preInletDomain.y1 = value.y1; }
    if (value.z0 < preInletDomain.z0) { preInletDomain.z0 = value.z0; }
    if (value.z1 > preInletDomain.z1) { preInletDomain.z1 = value.z1; }
  }  
  
  location = preInletDomain.enlarge(1);
  
  if (direction == Direction::Xneg) {
    location.x0 -= preinlet_length;
  } 
  if (direction == Direction::Yneg) {
    location.y0 -= preinlet_length;
  } 
  if (direction == Direction::Zneg) {
    location.z0 -= preinlet_length;
  }
  if (direction == Direction::Xpos) {
    location.x1 += preinlet_length;
  } 
  if (direction == Direction::Ypos) {
    location.y1 += preinlet_length;
  } 
  if (direction == Direction::Zpos) {
    location.z1 += preinlet_length;
  }
    Box3D system = hemocell->lattice->getBoundingBox();
  if (direction == Direction::Xneg || direction == Direction::Xpos) {
    if (location.y0 < system.y0 ) { location.y0 =system.y0; }
    if (location.y1 > system.y1 ) { location.y1 =system.y1; }
    if (location.z0 < system.z0 ) { location.z0 =system.z0; }
    if (location.z1 > system.z1 ) { location.z1 =system.z1; }
  }
  if (direction == Direction::Yneg || direction == Direction::Ypos) {
    if (location.x0 < system.x0 ) { location.x0 =system.x0; }
    if (location.x1 > system.x1 ) { location.x1 =system.x1; }
    if (location.z0 < system.z0 ) { location.z0 =system.z0; }
    if (location.z1 > system.z1 ) { location.z1 =system.z1; }
  }
  if (direction == Direction::Zneg || direction == Direction::Zpos) {
    if (location.y0 < system.y0 ) { location.y0 =system.y0; }
    if (location.y1 > system.y1 ) { location.y1 =system.y1; }
    if (location.x0 < system.x0 ) { location.x0 =system.x0; }
    if (location.x1 > system.x1 ) { location.x1 =system.x1; }
  }
}

PreInlet::PreInlet(HemoCell * hemocell_, plb::MultiScalarField3D<int> * flagMatrix_) {
  hemocell = hemocell_;
  flagMatrix = flagMatrix_;
  preinlet_length = (*hemocell->cfg)["preInlet"]["parameters"]["lengthN"].read<int>();
}

PreInlet::PreInlet(HemoCell * hemocell_, plb::MultiBlockManagement3D & management) {
  hemocell = hemocell_;
  flagMatrix = new plb::MultiScalarField3D<int>(management,
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<int>(),
            1);
  vector<MultiBlock3D *> wrapper;
  wrapper.push_back(hemocell->lattice);
  wrapper.push_back(flagMatrix);
  applyProcessingFunctional(new FillFlagMatrix(),hemocell->lattice->getBoundingBox(),wrapper);
  preinlet_length = (*hemocell->cfg)["preInlet"]["parameters"]["lengthN"].read<int>();
}

void PreInlet::CreateDrivingForceFunctional::processGenericBlocks(plb::Box3D domain, std::vector<plb::AtomicBlock3D*> blocks) {
  BlockLattice3D<T,DESCRIPTOR> * ff = dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
  for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
    for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
      for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
        if (!ff->get(iX,iY,iZ).getDynamics().isBoundary()) {
          count++;
        }
      }
    }
  }
}

void PreInlet::FillFlagMatrix::processGenericBlocks(plb::Box3D domain, std::vector<plb::AtomicBlock3D*> blocks) {
  BlockLattice3D<T,DESCRIPTOR> * ff = dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
  ScalarField3D<int> * sf = dynamic_cast<ScalarField3D<int>*>(blocks[1]);
  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
    for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
      for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
        if (ff->get(x,y,z).getDynamics().isBoundary()) {
          sf->get(x,y,z) = 0;
        }
      }
    }
  }
}

void PreInlet::calculateDrivingForce() {
  map<int,plint> areas;
  plint fluidArea = 0;
  double re = (*hemocell->cfg)["preInlet"]["parameters"]["Re"].read<T>();

  if (partOfpreInlet) {
    
    plb::Box3D domain = hemocell->lattice->getBoundingBox();
    if (direction == Direction::Xneg) {
      domain.x1 = domain.x0 = domain.x0 + 2;
    } 
    if (direction == Direction::Yneg) {
      domain.y1 = domain.y0 = domain.y0 + 2;
    } 
    if (direction == Direction::Zneg) {
      domain.z1 = domain.z0 = domain.z0 + 2;
    }
    if (direction == Direction::Xpos) {
      domain.x0 = domain.x1 = domain.x1 - 2;
    } 
    if (direction == Direction::Ypos) {
      domain.y0 = domain.y1 = domain.y1 - 2;
    } 
    if (direction == Direction::Zpos) {
      domain.z0 = domain.z1 = domain.z1 - 2;
    } 
    vector<MultiBlock3D*> wrapper;
    wrapper.push_back(hemocell->lattice);
    applyProcessingFunctional(new CreateDrivingForceFunctional(fluidArea),domain,wrapper);
    areas[global::mpi().getRank()] = fluidArea;
  }
  
  HemoCellGatheringFunctional<plint>::gather(areas);
  
    fluidArea = 0;
    for (auto & pair: areas) {
      fluidArea += pair.second;
    }
    
    T pipe_radius = sqrt(fluidArea/PI);
    hlog << "(Parameters) Your preInlet pipe has a calculated radius of " << pipe_radius << " LU, assuming a perfect circle" << std::endl;
  if(partOfpreInlet) {     
    T u_lbm_max = re * param::nu_lbm / (pipe_radius*2);
    drivingForce =  8 * param::nu_lbm * (u_lbm_max * 0.5) / pipe_radius / pipe_radius;
  }
}

void PreInlet::setDrivingForce() {
  if (partOfpreInlet) {
    plb::Array<T,3> force(0.,0.,0.);
    if (direction == Direction::Xneg) {
      force[0] = drivingForce;
    } 
    if (direction == Direction::Yneg) {
      force[1] = drivingForce;
    } 
    if (direction == Direction::Zneg) {
      force[2] = drivingForce;
    }
    if (direction == Direction::Xpos) {
      force[0] = -drivingForce;
    } 
    if (direction == Direction::Ypos) {
      force[1] = -drivingForce;
    } 
    if (direction == Direction::Zpos) {
      force[2] = -drivingForce;
    }
    setExternalVector(*hemocell->lattice, (*hemocell->lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    force);
  }
}

void PreInlet::createBoundary() {
  plb::Box3D domain = location; 

  for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
    for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
      for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
        bool fluid = false;
        if (direction == Direction::Xneg || direction == Direction::Xpos) {
          if (y == location.y0-1 || z == location.z0-1) {
            defineDynamics(*hemocell->lattice,x,y,z, new BounceBack<T,DESCRIPTOR>(1.));
            continue;
          }
          fluid = flagMatrix->get(fluidInlet.x0,y,z);
        }
        if (direction == Direction::Yneg || direction == Direction::Ypos) {
          if (x == location.x0-1 || z == location.z0-1) {
            defineDynamics(*hemocell->lattice,x,y,z, new BounceBack<T,DESCRIPTOR>(1.));
            continue;
          }
          fluid = flagMatrix->get(x,fluidInlet.y0,z);
        }
        if (direction == Direction::Zneg || direction == Direction::Zpos) {
          if (x == location.x0-1 || y == location.y0-1) {
            defineDynamics(*hemocell->lattice,x,y,z, new BounceBack<T,DESCRIPTOR>(1.));
            continue;
          }
          fluid = flagMatrix->get(x,y,fluidInlet.z0);
        }
        if(partOfpreInlet && !fluid) {
          defineDynamics(*hemocell->lattice,x,y,z, new BounceBack<T,DESCRIPTOR>(1.));
        }
      }
    }
  }
     
  if (partOfpreInlet) {
    if (direction == Direction::Xneg || direction == Direction::Xpos) {
      hemocell->lattice->periodicity().toggle(0,true);  
    }
    if (direction == Direction::Yneg || direction == Direction::Ypos) {
      hemocell->lattice->periodicity().toggle(1,true);  
    }
    if (direction == Direction::Zneg || direction == Direction::Zpos) {
      hemocell->lattice->periodicity().toggle(2,true);  
    }
  }
}

PreInlet::CreatePreInletBoundingBox *        PreInlet::CreatePreInletBoundingBox::clone() const { return new PreInlet::CreatePreInletBoundingBox(*this);}
PreInlet::CreateDrivingForceFunctional *        PreInlet::CreateDrivingForceFunctional::clone() const { return new PreInlet::CreateDrivingForceFunctional(*this);}
PreInlet::FillFlagMatrix * PreInlet::FillFlagMatrix::clone() const { return new PreInlet::FillFlagMatrix(*this);}
}
