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

#include "core/geometry3D.h"
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


/*
 * Initialise the communication structure between ranks near the inlets of the
 * pre-inlet and main domain. These ranks are responsible to transmit (copy)
 * the cells (RBC, PLT, etc.) from the periodic pre-inlet into the main domain.
 * To minimise communication, send/receive pairs are only defined for ranks for
 * which their (bulk) domains overlaps, i.e. their bounding boxes intersect.
 */
void PreInlet::initializePreInletParticleBoundary() {
  MPI_Barrier(MPI_COMM_WORLD);

  Box3D preinlet_boundary = fluidInlet;

  if (hemocell->partOfpreInlet) {
    switch (direction) {
      case Direction::Xneg:
        preinlet_boundary.x0 = location.x0;
        break;
      case Direction::Yneg:
        preinlet_boundary.y0 = location.y0;
        break;
      case Direction::Zneg:
        preinlet_boundary.z0 = location.z0;
        break;
      case Direction::Xpos:
        preinlet_boundary.x1 = location.x1;
        break;
      case Direction::Ypos:
        preinlet_boundary.y1 = location.y1;
        break;
      case Direction::Zpos:
        preinlet_boundary.z1 = location.z1;
        break;
    }
  }

  switch (direction) {
    case Direction::Xneg:
      preinlet_boundary.x0 += 1;
      preinlet_boundary.x1 = preinlet_boundary.x0 + inflow_length;
      break;
    case Direction::Yneg:
      preinlet_boundary.y0 += 1;
      preinlet_boundary.y1 = preinlet_boundary.y0 + inflow_length;
      break;
    case Direction::Zneg:
      preinlet_boundary.z0 += 1;
      preinlet_boundary.z1 = preinlet_boundary.z0 + inflow_length;
      break;
    case Direction::Xpos:
      preinlet_boundary.x1 -= 1;
      preinlet_boundary.x0 = preinlet_boundary.x1 - inflow_length;
      break;
    case Direction::Ypos:
      preinlet_boundary.y1 -= 1;
      preinlet_boundary.y0 = preinlet_boundary.y1 - inflow_length;
      break;
    case Direction::Zpos:
      preinlet_boundary.z1 -= 1;
      preinlet_boundary.z0 = preinlet_boundary.z1 - inflow_length;
      break;
  }

  auto rank = global::mpi().getRank();

  // Extract all ranks and their domains that intersect the pre-inlet interface.
  for (int block_id : hemocell->lattice->getLocalInfo().getBlocks()) {
    Box3D bulk = hemocell->lattice->getMultiBlockManagement().getBulk(block_id);

    domain_at_rank[rank] = bulk;

    /*
    FIXME: This assumes only a single bulk per rank is considered, which
    seems to be the cases most of the time. An improviement would be to
    implementt the subsequent if-statement: if the bulk is already set for
    the current rank -> merge the bulks, otherwise inser the current bulk.
    If situations do exists where such merges would fail, we might need to
    opt to store an vector of bulks and loop over all present entries.

    if (domain_at_rank.find(rank) == domain_at_rank.end()) {
      domain_at_rank[rank] = bulk;
    } else {
      auto merged = merge(domain_at_rank[rank], bulk);
      if (not merged) {
        hlog << "Pre-inlet cannot handle not-aligned bulks efficiently." << endl;
        exit(1);
      }
    }
    */

    if (!doesIntersect(preinlet_boundary, bulk)) {
      continue;
    }

    if (hemocell->partOfpreInlet) {
      particleSendMpi[rank] += 1;
    } else {
      particleReceiveMpi[rank] = true;
    }

    communicating_blocks.push_back(block_id);
  }

  // Exchange the maps: make them available on all other ranks.
  HemoCellGatheringFunctional<plint>::gather(particleSendMpi);
  HemoCellGatheringFunctional<bool>::gather(particleReceiveMpi);
  HemoCellGatheringFunctional<Box3D>::gather(domain_at_rank);

  // FIXME: The implementation is currently separated between the "preinlet"
  // and "!preinlet" sides of the domain. It is, most-likely, possible to
  // collapse this if-statement, as both operations are identical except for
  // which maps they considered. Alternatively only the maps are swapped, while
  // the remainder of the routine is kept the same.
  if (hemocell->partOfpreInlet) {
    if (particleSendMpi[rank] > 0) {
      Box3D domain = domain_at_rank[rank];

      for (auto & pair : particleReceiveMpi) {
        auto other = domain_at_rank[pair.first];
        bool intersects = false;

        // As the considered bulks are already adjacent with the interface of the
        // pre-inlet, only intersection in the plane orthogonal to its
        // configuration need to be tested for.
        switch (direction) {
          case Direction::Xneg:
          case Direction::Xpos:
            intersects = doesIntersect(
              Box2D(domain.y0, domain.y1, domain.z0, domain.z1),
              Box2D(other.y0, other.y1, other.z0, other.z1));
            break;
          case Direction::Yneg:
          case Direction::Ypos:
            intersects = doesIntersect(
              Box2D(domain.x0, domain.x1, domain.z0, domain.z1),
              Box2D(other.x0, other.x1, other.z0, other.z1));
            break;
          case Direction::Zneg:
          case Direction::Zpos:
            intersects = doesIntersect(
              Box2D(domain.x0, domain.x1, domain.y0, domain.y1),
              Box2D(other.x0, other.x1, other.y0, other.y1));
            break;
        }

        if (pair.second && intersects) {
          my_send_blocks.push_back(pair.first);
        }
      }
    }
  } else {
    if (particleReceiveMpi[rank] == true) {
      auto domain = domain_at_rank[rank];

      for (auto & pair : particleSendMpi) {
        auto other = domain_at_rank[pair.first];
        bool intersects = false;

        switch (direction) {
          case Direction::Xneg:
          case Direction::Xpos:
            intersects = doesIntersect(
              Box2D(domain.y0, domain.y1, domain.z0, domain.z1),
              Box2D(other.y0, other.y1, other.z0, other.z1));
            break;
          case Direction::Yneg:
          case Direction::Ypos:
            intersects = doesIntersect(
              Box2D(domain.x0, domain.x1, domain.z0, domain.z1),
              Box2D(other.x0, other.x1, other.z0, other.z1));
            break;
          case Direction::Zneg:
          case Direction::Zpos:
            intersects = doesIntersect(
              Box2D(domain.x0, domain.x1, domain.y0, domain.y1),
              Box2D(other.x0, other.x1, other.y0, other.y1));
            break;
        }

        if (pair.second && intersects) {
          my_recv_blocks.push_back(pair.first);
        }
      }
    }
  }
}


/*
 * Transmit cells (RBC, PLT, ...) located within the periodic pre-inlet to
 * their corresponding location in the main simulation lattice. Cells are only
 * transmitted between ranks of which their domains overlap near the interface
 * between the pre-inlet and the main domain. All other ranks return early.
 *
 * The communication only happens in a single direction, where cells are only
 * send from the pre-inlet towards the main simulation domain.
 */
void PreInlet::applyPreInletParticleBoundary() {
  global.statistics.getCurrent()["applyPreInletParticleBoundary"].start();
  if (!partOfpreInlet) {
    hemocell->cellfields->syncEnvelopes();
    hemocell->cellfields->deleteIncompleteCells(false);
  }
  vector<MPI_Request> requests;
  vector<vector<char>> buffers;
  MPI_Barrier(MPI_COMM_WORLD);
  if (partOfpreInlet) {
    if (particleSendMpi.find(global::mpi().getRank()) != particleSendMpi.end()) {
      for (plint bid : communicating_blocks) {
        buffers.push_back(vector<char>());
        Box3D domain = fluidInlet;

        switch (direction) {
          case Direction::Xneg:
            domain.x0 = domain.x0 - preinlet_length;
            domain.x1 = domain.x0 + inflow_length;
            break;
          case Direction::Yneg:
            domain.y0 = domain.y0 - preinlet_length;
            domain.y1 = domain.y0 + inflow_length;
            break;
          case Direction::Zneg:
            domain.z0 = domain.z0 - preinlet_length;
            domain.z1 = domain.z0 + inflow_length;
            break;
          case Direction::Xpos:
            domain.x1 = domain.x1 + preinlet_length;
            domain.x0 = domain.x1 - inflow_length;
            break;
          case Direction::Ypos:
            domain.y1 = domain.y1 + preinlet_length;
            domain.y0 = domain.y1 - inflow_length;
            break;
          case Direction::Zpos:
            domain.z1 = domain.z1 + preinlet_length;
            domain.z0 = domain.z1 - inflow_length;
            break;
        }

        Dot3D shift = hemocell->cellfields->immersedParticles->getComponent(bid).getLocation();
        domain = domain.shift(-shift.x,-shift.y,-shift.z);
        hemocell->cellfields->immersedParticles->getComponent(bid).particleDataTransfer.send_preinlet(domain,buffers.back(),modif::hemocell);

        for (auto & pid : my_send_blocks) {
          requests.push_back(MPI_Request());
          MPI_Isend(&buffers.back()[0],buffers.back().size(),MPI_CHAR,pid,0,MPI_COMM_WORLD,&requests.back());
        }
      }
    }
    MPI_Waitall(requests.size(),&requests[0],MPI_STATUSES_IGNORE);
  } else {
    for (size_t i = 0 ; i < my_recv_blocks.size() ; i++) {
      MPI_Status status;
      int count;
      MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
      MPI_Get_count(&status,MPI_CHAR,&count);
      buffers.push_back(vector<char>(count));
      MPI_Recv(&buffers.back()[0],count,MPI_CHAR,status.MPI_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      Dot3D offset(0,0,0);

      switch (direction) {
        case Direction::Xneg:
          offset.x = preinlet_length;
          break;
        case Direction::Yneg:
          offset.y = preinlet_length;
          break;
        case Direction::Zneg:
          offset.z = preinlet_length;
          break;
        case Direction::Xpos:
          offset.x = -preinlet_length;
          break;
        case Direction::Ypos:
          offset.y = -preinlet_length;
          break;
        case Direction::Zpos:
          offset.z = -preinlet_length;
          break;
      }
      // NOTE: Reading into all blocks should be OK as long as a single block
      // or adjacent blocks are considered. For blocks that are sparse and/or
      // spatially "far away", additional checks should be added to avoid
      // excess extracting of cells that are not present in the current block.
      for (int bId : hemocell->cellfields->immersedParticles->getLocalInfo().getBlocks()) {
        hemocell->cellfields->immersedParticles->getComponent(bId).particleDataTransfer.receivePreInlet(&buffers.back()[0],buffers.back().size(),modif::hemocell,offset);
        hemocell->cellfields->immersedParticles->getComponent(bId).invalidate_ppc();
        hemocell->cellfields->immersedParticles->getComponent(bId).invalidate_lpc();
        hemocell->cellfields->immersedParticles->getComponent(bId).invalidate_pg();
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  global.statistics.getCurrent().stop();
}

void PreInlet::applyPreInletVelocityBoundary() {
  global.statistics.getCurrent()["applyPreInletVelocityBoundary"].start();
  MPI_Barrier(MPI_COMM_WORLD);

  Box3D & domain = fluidInlet;
  Box3D result;
  int z_o = domain.getNz();
  int y_o = domain.getNy();
  int tag;

  plb::Array<T,3> receiver;
  plb::Array<T,3> vel;
  std::vector<plb::Array<T,3>> buffer;
  std::vector<MPI_Request> requests;

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

            // send velocity from preInlet boundary
            buffer.push_back(vel);
            requests.push_back(MPI_Request());
            MPI_Isend(&buffer.back()[0],3*sizeof(T),MPI_CHAR,dest,tag,MPI_COMM_WORLD,&requests.back());
          } else {
            Box3D point(x,x,y,y,z,z);
            int source = hemocell->preinlet_lattice_management->getThreadAttribution().getMpiProcess(hemocell->preinlet_lattice_management->getSparseBlockStructure().locate(x+loc.x,y+loc.y,z+loc.z));

            // receive velocity from preInlet
            MPI_Recv(&receiver,3*sizeof(T),MPI_CHAR,source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            // update velocity lattice main domain
            setBoundaryVelocity(hemocell->lattice->getComponent(bId),point,receiver);
          }
        }
      }
     }
    }
  }
  MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
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
  if (!hemocell->lattice) {
    hlog << "(PreInlet::preInletFromSlice) preInletFromSlice called without setting up a lattice first" << endl;
  } else {
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
  if (!hemocell->lattice) {
    hlog << "(PreInlet::autoPreinletFromBoundary) autoPreinletFromBoundary called without setting up a lattice first" << endl;
  } else {
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

// JON: average function...
double PreInlet::average(vector<double> values) {
  double sum = 0.0;
  int N = values.size();
  for (int i = 0; i < N; i++) {
    sum += values[i];
  }
  return sum / N;
}

// JON: Read pulsatility data into normalizedVelocityTimes / normalizedVelocityValues arrays which are declared in preInlet.h
// Also computes average which is a constant that is used every time that setDrivingForceTimeDependent() is called
bool PreInlet::readNormalizedVelocities() {

  // Open file
  std::string pulseFileName = (*hemocell->cfg)["preInlet"]["parameters"]["pulseFileName"].read<std::string>();
  std::ifstream pulseFile(pulseFileName);

  // Check if pulsatility file was found
  if (!pulseFile.is_open()) {
    cout << "*** WARNING! pulsatility data file " << pulseFileName << " does not exist!" << endl;
    return false;
  }

  // Read and parse data into arrays
  std::string line;
  while (std::getline(pulseFile, line))
  {
      std::string DELIM = " ";
      std::string t = line.substr(0, line.find(DELIM));
      std::string v = line.substr(line.find(DELIM), line.size());

      normalizedVelocityTimes .push_back(stod(t));
      normalizedVelocityValues.push_back(stod(v));
  }
  pulseFile.close();

  // Set average velocity
  average_vel = PreInlet::average(normalizedVelocityValues);

  // Length of heartbeat pulse is also important
  pulseEndTime = normalizedVelocityTimes[normalizedVelocityTimes.size() - 1];

  // Read in pulsatility frequency value. Just leave it at its value as per the input pulse data if it's not in the XML.
  try  { pFrequency = (*hemocell->cfg)["preInlet"]["parameters"]["pFrequency"].read<double>(); }
  catch (const std::invalid_argument& e)  { pFrequency = 1.0 / pulseEndTime; }

  return true;
}

// JON: got this interpolation from the internet somewhere
double PreInlet::interpolate(vector<double> &xData, vector<double> &yData, double x, bool extrapolate)
{
   int size = xData.size();

   // find left end of interval for interpolation
   int i = 0;
   // special case: beyond right end
   if ( x >= xData[size - 2] )
   {
      i = size - 2;
   }
   else
   {
      while ( x > xData[i+1] ) i++;
   }
   // points on either side (unless beyond ends)
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];
   // if beyond ends of array and not extrapolating
   if ( !extrapolate )
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   // gradient
   double dydx = ( yR - yL ) / ( xR - xL );

   // linear interpolation
   return yL + dydx * ( x - xL );
}

// JON: this function is similar to setDrivingForce() but it reads in a value
// from pulsatile velocity data and then interpolates driving force
void PreInlet::setDrivingForceTimeDependent(double t) {
  if (partOfpreInlet) {
    // Multiply time value by frequency and length of input pulse in order to get appropriate point
    // input beat curve
    t *= pFrequency * pulseEndTime;

    // Time value must wrap around normalized velocity data in order to generate periodic heart beat
    t = fmod(t, pulseEndTime);

    // JON: Compute current driving force based on ratio of data-derived "current velocity"
    // and average velocity
    double current_vel = PreInlet::interpolate(normalizedVelocityTimes, normalizedVelocityValues, t, false);
    double currentDrivingForce = (current_vel / average_vel) * drivingForce;

    plb::Array<T,3> force(0.,0.,0.);
    if (direction == Direction::Xneg) {
      force[0] = currentDrivingForce;
    }
    if (direction == Direction::Yneg) {
      force[1] = currentDrivingForce;
    }
    if (direction == Direction::Zneg) {
      force[2] = currentDrivingForce;
    }
    if (direction == Direction::Xpos) {
      force[0] = -currentDrivingForce;
    }
    if (direction == Direction::Ypos) {
      force[1] = -currentDrivingForce;
    }
    if (direction == Direction::Zpos) {
      force[2] = -currentDrivingForce;
    }
    setExternalVector(*hemocell->lattice, (*hemocell->lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    force);
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
