#include "cellInformation.h"

void CellInformationFunctionals::clear_list() {
  info_per_cell.clear();
}
void CellInformationFunctionals::calculate_vol_pos_area(HemoCell* hemocell) {
  getCellArea(hemocell);
  getCellVolume(hemocell);
  getCellPosition(hemocell);
}

void CellInformationFunctionals::CellVolume::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    double volume = 0.;
    const int & cid = pair.first;
    const vector<HemoCellParticle*> & cell = pf->particles_per_cell[cid];
    const pluint ctype = pf->particles_per_cell[cid][0]->celltype;
    for (Array<plint,3> triangle : (*hemocell->cellfields)[ctype]->mechanics->cellConstants.triangle_list) {
      const Array<double,3> & v0 = cell[triangle[0]]->position;
      const Array<double,3> & v1 = cell[triangle[1]]->position;
      const Array<double,3> & v2 = cell[triangle[2]]->position;
      
      //Volume
      const double v210 = v2[0]*v1[1]*v0[2];
      const double v120 = v1[0]*v2[1]*v0[2];
      const double v201 = v2[0]*v0[1]*v1[2];
      const double v021 = v0[0]*v2[1]*v1[2];
      const double v102 = v1[0]*v0[1]*v2[2];
      const double v012 = v0[0]*v1[1]*v2[2];
      volume += (1.0/6.0)*(-v210+v120+v201-v021-v102+v012);
    }
    info_per_cell[cid].volume = volume;
  }
}
void CellInformationFunctionals::CellArea::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    double total_area = 0.;
    const int & cid = pair.first;
    const vector<HemoCellParticle*> & cell = pf->particles_per_cell[cid];
    const pluint ctype = pf->particles_per_cell[cid][0]->celltype;
    for (Array<plint,3> triangle : (*hemocell->cellfields)[ctype]->mechanics->cellConstants.triangle_list) {
      const Array<double,3> & v0 = cell[triangle[0]]->position;
      const Array<double,3> & v1 = cell[triangle[1]]->position;
      const Array<double,3> & v2 = cell[triangle[2]]->position;
      
      total_area += computeTriangleArea(v0,v1,v2);  
    }
    info_per_cell[cid].area = total_area;
  }
}
void CellInformationFunctionals::CellPosition::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    Array<double,3> position = {0.,0.,0.};
    const int & cid = pair.first;
    const vector<HemoCellParticle*> & cell = pf->particles_per_cell[cid];
    for (const HemoCellParticle * particle : cell ) {
      position += particle->position;
    }
    info_per_cell[cid].position = position/double(cell.size());
    info_per_cell[cid].centerLocal = pf->isContainedABS(info_per_cell[cid].position,pf->localDomain);
  }
}
void CellInformationFunctionals::CellStretch::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    double max_stretch = 0.;
    const int & cid = pair.first;
    const vector<HemoCellParticle*> & cell = pf->particles_per_cell[cid];
    for (unsigned int i = 0 ; i < cell.size() - 1 ; i++ ) {
      for (unsigned int j = i + 1 ; j < cell.size() ; j ++) {
        double distance = sqrt( pow(cell[i]->position[0]-cell[j]->position[0],2) +
                                pow(cell[i]->position[1]-cell[j]->position[1],2) +
                                pow(cell[i]->position[2]-cell[j]->position[2],2));
        max_stretch = max_stretch < distance ? distance : max_stretch;
      }
    }
    info_per_cell[cid].stretch = max_stretch;
  }
}
void CellInformationFunctionals::CellBoundingBox::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    Array<double,6> bbox;
    const int & cid = pair.first;
    const vector<HemoCellParticle*> & cell = pf->particles_per_cell[cid];
    
    bbox[0] = cell[0]->position[0];
    bbox[1] = cell[0]->position[0];
    bbox[2] = cell[0]->position[1];
    bbox[3] = cell[0]->position[1];
    bbox[4] = cell[0]->position[2];
    bbox[5] = cell[0]->position[2];
    
    for (const HemoCellParticle * particle : cell ) {
      bbox[0] = bbox[0] > particle->position[0] ? particle->position[0] : bbox[0];
      bbox[1] = bbox[1] < particle->position[0] ? particle->position[0] : bbox[1];
      bbox[2] = bbox[2] > particle->position[1] ? particle->position[1] : bbox[2];
      bbox[3] = bbox[3] < particle->position[1] ? particle->position[1] : bbox[3];
      bbox[4] = bbox[4] > particle->position[2] ? particle->position[2] : bbox[4];
      bbox[5] = bbox[5] < particle->position[2] ? particle->position[2] : bbox[5];

    }
    info_per_cell[cid].bbox = bbox;
  }
}
void CellInformationFunctionals::CellAtomicBlock::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    const int & cid = pair.first;

    info_per_cell[cid].blockId = pf->atomicBlockId;
  }
}
void CellInformationFunctionals::CellType::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
  HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
  
  for (const auto & pair : pf->lpc) {
    const int & cid = pair.first;

    info_per_cell[cid].cellType = pf->particles_per_cell[cid][0]->celltype;
  }
}

void CellInformationFunctionals::getCellVolume(HemoCell * hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellVolume(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::getCellArea(HemoCell * hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellArea(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::getCellPosition(HemoCell * hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellPosition(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::getCellStretch(HemoCell * hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellStretch(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::getCellBoundingBox(HemoCell * hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellBoundingBox(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::getCellAtomicBlock(HemoCell* hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellAtomicBlock(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
void CellInformationFunctionals::getCellType(HemoCell* hemocell_) {
  hemocell = hemocell_;
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new CellType(),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
}
pluint CellInformationFunctionals::getTotalNumberOfCells(HemoCell* hemocell) {
  info_per_cell.clear(); //TODO thread safe n such
  getCellPosition(hemocell); //Any functional will do
  pluint localCells = 0;
  for (const auto & pair : info_per_cell) {
    const CellInformation & cinfo = pair.second;
    if (cinfo.centerLocal) {
      localCells++;
    }
  }
  
  map<int,pluint> cells_per_proc;
  cells_per_proc[global::mpi().getRank()] = localCells;
  HemoCellGatheringFunctional<pluint>::gather(cells_per_proc,global::mpi().getSize());
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