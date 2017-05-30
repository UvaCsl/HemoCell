/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   writeCellInfoCSV.cpp
 * Author: vikko
 * 
 * Created on 30 May 2017, 14:27
 */

#include "writeCellInfoCSV.h"

void writeCellInfo_CSV(HemoCell * hemocell) {
  WriteCellInfoCSV * wcic = new WriteCellInfoCSV();
  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyProcessingFunctional(wcic,hemocell->lattice->getBoundingBox(),wrapper);
}

void WriteCellInfoCSV::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
}
WriteCellInfoCSV * WriteCellInfoCSV::clone() const { return new WriteCellInfoCSV(*this);}
