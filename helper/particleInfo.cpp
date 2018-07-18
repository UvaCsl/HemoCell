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
#include "particleInfo.h"
#include "hemoCellParticleField.h"
#include "hemocell.h"

namespace hemo {
  
void GatherParticleVelocity::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
    HEMOCELL_PARTICLE_FIELD* pf = dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[0]);
    vector<HemoCellParticle *> localParticles;
    pf->findParticles(pf->localDomain,localParticles);
    
    if (localParticles.size() > 0) {
      //initial value
      hemo::Array<T,3> vel_vec = localParticles[0]->sv.v;
      T vel = sqrt(vel_vec[0]*vel_vec[0]+vel_vec[1]*vel_vec[1]+vel_vec[2]*vel_vec[2]);
      T min=vel,max=vel,avg=0.;


      for (const HemoCellParticle * particle : localParticles) {
        vel_vec = particle->sv.v;
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
      hemo::Array<T,3> force_vec = localParticles[0]->sv.force +localParticles[0]->sv.force_repulsion;
      T force = sqrt(force_vec[0]*force_vec[0]+force_vec[1]*force_vec[1]+force_vec[2]*force_vec[2]);
      T min=force,max=force,avg=0.;


      for (const HemoCellParticle * particle : localParticles) {
        force_vec = particle->sv.force + particle->sv.force_repulsion;
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
  applyProcessingFunctional(new GatherParticleVelocity(gatherValues),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
  //int numAtomicBlock = hemocell->lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
  HemoCellGatheringFunctional<ParticleStatistics>::gather(gatherValues);
  
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
  applyProcessingFunctional(new GatherParticleForce(gatherValues),hemocell->cellfields->immersedParticles->getBoundingBox(),wrapper);
  //int numAtomicBlock = hemocell->lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
  HemoCellGatheringFunctional<ParticleStatistics>::gather(gatherValues);
  
  if (gatherValues.size() > 0) {
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
  } else {
    ParticleStatistics result;
    result.avg = 0;
    result.max = 0;
    result.min = 0;
    result.ncells = 0;
    return result;
  }
}

GatherParticleVelocity * GatherParticleVelocity::clone() const { return new GatherParticleVelocity(*this); }
GatherParticleForce * GatherParticleForce::clone() const { return new GatherParticleForce(*this); }

}