#include "cellInformation.h"

void CellInformationFunctionals::clear_list() {
  info_per_cell.clear();
}
void CellInformationFunctionals::calculate_all(HemoCell* hemocell) {
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

CellInformationFunctionals::CellVolume * CellInformationFunctionals::CellVolume::clone() const { return new CellInformationFunctionals::CellVolume(*this);}
CellInformationFunctionals::CellArea * CellInformationFunctionals::CellArea::clone() const { return new CellInformationFunctionals::CellArea(*this);}
CellInformationFunctionals::CellPosition * CellInformationFunctionals::CellPosition::clone() const { return new CellInformationFunctionals::CellPosition(*this);}

map<int,CellInformation> CellInformationFunctionals::info_per_cell = map<int,CellInformation>();
HemoCell * CellInformationFunctionals::hemocell = 0;