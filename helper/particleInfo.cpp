/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particleInfo.cpp
 * Author: vikko
 * 
 * Created on 01 June 2017, 13:53
 */

#include "particleInfo.h"
#include "hemoCellParticleField.h"
#include "hemocell.h"

void GatherParticleVelocity::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
    vector<HemoCellParticle *> localParticles;
    pf->findParticles(pf->localDomain,localParticles);
    
    if (localParticles.size() > 0) {
      //initial value
      Array<double,3> vel_vec = localParticles[0]->v;
      double vel = sqrt(vel_vec[0]*vel_vec[0]+vel_vec[1]*vel_vec[1]+vel_vec[2]*vel_vec[2]);
      double min=vel,max=vel,avg=0.;


      for (const HemoCellParticle * particle : localParticles) {
        vel_vec = particle->v;
        vel = sqrt(vel_vec[0]*vel_vec[0]+vel_vec[1]*vel_vec[1]+vel_vec[2]*vel_vec[2]);
        min = min > vel ? vel : min;
        max = max < vel ? vel : max;
        avg += vel;   
      }


      gatherValues[pf->atomicBlockId].min = min;
      gatherValues[pf->atomicBlockId].max = max;
      gatherValues[pf->atomicBlockId].avg = avg/localParticles.size();
      gatherValues[pf->atomicBlockId].ncells = localParticles.size();
    }
}
void GatherParticleForce::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
    vector<HemoCellParticle *> localParticles;
    pf->findParticles(pf->localDomain,localParticles);
    
    if (localParticles.size() > 0) {
      //initial value
      Array<double,3> force_vec = localParticles[0]->force +localParticles[0]->force_repulsion;
      double force = sqrt(force_vec[0]*force_vec[0]+force_vec[1]*force_vec[1]+force_vec[2]*force_vec[2]);
      double min=force,max=force,avg=0.;


      for (const HemoCellParticle * particle : localParticles) {
        force_vec = particle->force + particle->force_repulsion;
        force = sqrt(force_vec[0]*force_vec[0]+force_vec[1]*force_vec[1]+force_vec[2]*force_vec[2]);
        min = min > force ? force : min;
        max = max < force ? force : max;
        avg += force;   
      }


      gatherValues[pf->atomicBlockId].min = min;
      gatherValues[pf->atomicBlockId].max = max;
      gatherValues[pf->atomicBlockId].avg = avg/localParticles.size();
      gatherValues[pf->atomicBlockId].ncells = localParticles.size();
    }
}
ParticleStatistics ParticleInfo::calculateVelocityStatistics(HemoCell* hemocell) {
  map<int,ParticleStatistics> gatherValues;
  
  //hemocell->cellfields->interpolateFluidVelocity();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new GatherParticleVelocity(gatherValues),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
  int numAtomicBlock = hemocell->lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
  HemoCellGatheringFunctional<ParticleStatistics>::gather(gatherValues,numAtomicBlock);
  
  ParticleStatistics result = gatherValues.begin()->second;
  result.avg = 0.;
  result.ncells = 0;
  
  for (const auto & pair : gatherValues) {
    const ParticleStatistics & cur = pair.second;
    result.avg += (cur.avg * cur.ncells);
    result.ncells += cur.ncells;
    result.min = result.min > cur.min ? cur.min : result.min;
    result.max = result.max < cur.max ? cur.max : result.max;
  }
  
  result.avg /= result.ncells;

  return result;
}

ParticleStatistics ParticleInfo::calculateForceStatistics(HemoCell* hemocell) {
  map<int,ParticleStatistics> gatherValues;
  
  //hemocell->cellfields->interpolateFluidVelocity();

  vector<MultiBlock3D*> wrapper;
  wrapper.push_back(hemocell->cellfields->immersedParticles);
  applyTimedProcessingFunctional(new GatherParticleForce(gatherValues),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
  int numAtomicBlock = hemocell->lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
  HemoCellGatheringFunctional<ParticleStatistics>::gather(gatherValues,numAtomicBlock);
  
  ParticleStatistics result = gatherValues.begin()->second;
  result.avg = 0.;
  result.ncells = 0;
  
  for (const auto & pair : gatherValues) {
    const ParticleStatistics & cur = pair.second;
    result.avg += (cur.avg * cur.ncells);
    result.ncells += cur.ncells;
    result.min = result.min > cur.min ? cur.min : result.min;
    result.max = result.max < cur.max ? cur.max : result.max;
  }
  
  result.avg /= result.ncells;

  return result;
}

GatherParticleVelocity * GatherParticleVelocity::clone() const { return new GatherParticleVelocity(*this); }
GatherParticleForce * GatherParticleForce::clone() const { return new GatherParticleForce(*this); }
