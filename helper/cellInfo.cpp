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
#include "cellInfo.h"
#include "hemocell.h"

void CellInformationFunctionals::clear_list() {
  info_per_cell.clear();
}
void CellInformationFunctionals::calculate_vol_pos_area(HemoCell* hemocell) {
  hemocell->cellfields->deleteIncompleteCells(false);
  calculateCellArea(hemocell);
  calculateCellVolume(hemocell);
  calculateCellPosition(hemocell);
}

void CellInformationFunctionals::CellVolume::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);

  for (const auto & pair : pf->get_lpc()) {
    T volume = 0.;
    const int & cid = pair.first;
    const vector<int> & cell = pf->get_particles_per_cell().at(cid);
    const pluint ctype = pf->particles[cell[0]].sv.celltype;
    for (hemo::Array<plint,3> triangle : (*hemocell->cellfields)[ctype]->mechanics->cellConstants.triangle_list) {
      const hemo::Array<T,3> & v0 = pf->particles[cell[triangle[0]]].sv.position;
      const hemo::Array<T,3> & v1 = pf->particles[cell[triangle[1]]].sv.position;
      const hemo::Array<T,3> & v2 = pf->particles[cell[triangle[2]]].sv.position;
      
      //Volume
      const T v210 = v2[0]*v1[1]*v0[2];
      const T v120 = v1[0]*v2[1]*v0[2];
      const T v201 = v2[0]*v0[1]*v1[2];
      const T v021 = v0[0]*v2[1]*v1[2];
      const T v102 = v1[0]*v0[1]*v2[2];
      const T v012 = v0[0]*v1[1]*v2[2];
      volume += (1.0/6.0)*(-v210+v120+v201-v021-v102+v012);
    }
    info_per_cell[cid].volume = volume;
  }
}
void CellInformationFunctionals::CellArea::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->get_lpc()) {
    T total_area = 0.;
    const int & cid = pair.first;
    const vector<int> & cell = pf->get_particles_per_cell().at(cid);
    const pluint ctype = pf->particles[cell[0]].sv.celltype;
    for (hemo::Array<plint,3> triangle : (*hemocell->cellfields)[ctype]->mechanics->cellConstants.triangle_list) {
      const hemo::Array<T,3> & v0 = pf->particles[cell[triangle[0]]].sv.position;
      const hemo::Array<T,3> & v1 = pf->particles[cell[triangle[1]]].sv.position;
      const hemo::Array<T,3> & v2 = pf->particles[cell[triangle[2]]].sv.position;

      total_area += computeTriangleArea(v0,v1,v2);  
    }
    info_per_cell[cid].area = total_area;
  }
}
void CellInformationFunctionals::CellPosition::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->get_lpc()) {
    hemo::Array<T,3> position = {0.,0.,0.};
    const int & cid = pair.first;
    const vector<int> & cell = pf->get_particles_per_cell().at(cid);
    unsigned int size = 0;
    for (const int pid : cell ) {
      if (pid == -1) { continue; }
      size++;
      position += pf->particles[pid].sv.position;
    }
    if ( info_per_cell.find(cid) == info_per_cell.end() || !info_per_cell[cid].centerLocal) {
      info_per_cell[cid].position = position/T(size);
      info_per_cell[cid].centerLocal = pf->isContainedABS(info_per_cell[cid].position,pf->localDomain);
    } else {
     
    }
  }
}
void CellInformationFunctionals::CellStretch::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->get_lpc()) {
    T max_stretch = 0.;
    const int & cid = pair.first;
    const vector<int> & cell = pf->get_particles_per_cell().at(cid);
    for (unsigned int i = 0 ; i < cell.size() - 1 ; i++ ) {
      for (unsigned int j = i + 1 ; j < cell.size() ; j ++) {
        if (cell[i] == -1 || cell[j] == -1) {continue;}
        T distance = sqrt( pow(pf->particles[cell[i]].sv.position[0]-pf->particles[cell[j]].sv.position[0],2) +
                                pow(pf->particles[cell[i]].sv.position[1]-pf->particles[cell[j]].sv.position[1],2) +
                                pow(pf->particles[cell[i]].sv.position[2]-pf->particles[cell[j]].sv.position[2],2));
        max_stretch = max_stretch < distance ? distance : max_stretch;
      }
    }
    info_per_cell[cid].stretch = max_stretch;
  }
}
void CellInformationFunctionals::CellBoundingBox::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->get_lpc()) {
    hemo::Array<T,6> bbox;
    const int & cid = pair.first;
    const vector<int> & cell = pf->get_particles_per_cell().at(cid);
    HemoCellParticle * particle = &pf->particles[cell[0]];
    
    bbox[0] = particle->sv.position[0];
    bbox[1] = particle->sv.position[0];
    bbox[2] = particle->sv.position[1];
    bbox[3] = particle->sv.position[1];
    bbox[4] = particle->sv.position[2];
    bbox[5] = particle->sv.position[2];
    
    for (const int pid : cell ) {
      particle = &pf->particles[pid];
      bbox[0] = bbox[0] > particle->sv.position[0] ? particle->sv.position[0] : bbox[0];
      bbox[1] = bbox[1] < particle->sv.position[0] ? particle->sv.position[0] : bbox[1];
      bbox[2] = bbox[2] > particle->sv.position[1] ? particle->sv.position[1] : bbox[2];
      bbox[3] = bbox[3] < particle->sv.position[1] ? particle->sv.position[1] : bbox[3];
      bbox[4] = bbox[4] > particle->sv.position[2] ? particle->sv.position[2] : bbox[4];
      bbox[5] = bbox[5] < particle->sv.position[2] ? particle->sv.position[2] : bbox[5];

    }
    info_per_cell[cid].bbox = bbox;
  }
}
void CellInformationFunctionals::CellAtomicBlock::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->get_lpc()) {
    const int & cid = pair.first;

    info_per_cell[cid].blockId = pf->atomicBlockId;
  }
}
void CellInformationFunctionals::CellType::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->get_lpc()) {
    const int & cid = pair.first;

    info_per_cell[cid].cellType = pf->particles[pf->get_particles_per_cell().at(cid)[0]].sv.celltype;
  }
}

void CellInformationFunctionals::calculateCellVolume(HemoCell * hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellVolume(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::calculateCellArea(HemoCell * hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellArea(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::calculateCellPosition(HemoCell * hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellPosition(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::calculateCellStretch(HemoCell * hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellStretch(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::calculateCellBoundingBox(HemoCell * hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellBoundingBox(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::calculateCellAtomicBlock(HemoCell* hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellAtomicBlock(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::calculateCellType(HemoCell* hemocell_) {
  hemocell = hemocell_;
  hemocell->cellfields->syncEnvelopes();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellType(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
pluint CellInformationFunctionals::getTotalNumberOfCells(HemoCell* hemocell) {
  info_per_cell.clear(); //TODO thread safe n such
  calculateCellPosition(hemocell); //Any functional will do
  pluint localCells = 0;
  for (const auto & pair : info_per_cell) {
    const CellInformation & cinfo = pair.second;
    if (cinfo.centerLocal) {
      localCells++;
    }
  }
  map<int,pluint> cells_per_proc;
  cells_per_proc[global::mpi().getRank()] = localCells;
  HemoCellGatheringFunctional<pluint>::gather(cells_per_proc);
  pluint total = 0;
  for (const auto & pair : cells_per_proc) {
    total += pair.second;
  }
  info_per_cell.clear();
  return total;
}
pluint CellInformationFunctionals::getNumberOfCellsFromType(HemoCell* hemocell, string type) {
  info_per_cell.clear(); //TODO thread safe n such
  calculateCellPosition(hemocell); 
  calculateCellType(hemocell);
  pluint localCells = 0;
  for (const auto & pair : info_per_cell) {
    const CellInformation & cinfo = pair.second;
    if (cinfo.centerLocal && (cinfo.cellType == (*hemocell->cellfields)[type]->ctype)) {
      localCells++;
    }
  }
  
  map<int,pluint> cells_per_proc;
  cells_per_proc[global::mpi().getRank()] = localCells;
  HemoCellGatheringFunctional<pluint>::gather(cells_per_proc);
  pluint total = 0;
  for (const auto & pair : cells_per_proc) {
    total += pair.second;
  }
  info_per_cell.clear();
  return total;
}


CellInformationFunctionals::CellVolume * CellInformationFunctionals::CellVolume::clone() const { return new CellInformationFunctionals::CellVolume(*this);}
CellInformationFunctionals::CellArea * CellInformationFunctionals::CellArea::clone() const { return new CellInformationFunctionals::CellArea(*this);}
CellInformationFunctionals::CellPosition * CellInformationFunctionals::CellPosition::clone() const { return new CellInformationFunctionals::CellPosition(*this);}
CellInformationFunctionals::CellStretch * CellInformationFunctionals::CellStretch::clone() const { return new CellInformationFunctionals::CellStretch(*this);}
CellInformationFunctionals::CellBoundingBox * CellInformationFunctionals::CellBoundingBox::clone() const { return new CellInformationFunctionals::CellBoundingBox(*this);}
CellInformationFunctionals::CellAtomicBlock * CellInformationFunctionals::CellAtomicBlock::clone() const { return new CellInformationFunctionals::CellAtomicBlock(*this);}
CellInformationFunctionals::CellType * CellInformationFunctionals::CellType::clone() const { return new CellInformationFunctionals::CellType(*this);}


map<int,CellInformation> CellInformationFunctionals::info_per_cell = map<int,CellInformation>();
HemoCell * CellInformationFunctionals::hemocell = 0;