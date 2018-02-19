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
#include "hemoCellParticleField.h"
#include "hemocell.h"

/* *************** class HemoParticleField3D ********************** */

HemoCellParticleField::HemoCellParticleField(plint nx, plint ny, plint nz)
    : AtomicBlock3D(nx,ny,nz, &this->dataTransfer)
{ 
    dataTransfer.setBlock(*this);
    AddOutputMap(); 
}

HemoCellParticleField::HemoCellParticleField(HemoCellParticleField const& rhs)
    : AtomicBlock3D(rhs)
{
    HemoCellParticle tmp;
    dataTransfer.setBlock(*this);
    for (const HemoCellParticle & particle : rhs.particles) {
      tmp = particle;
      addParticle(this->getBoundingBox(),&tmp);
    }
    ppc_up_to_date = false;
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    AddOutputMap();
}

HemoCellParticleField::~HemoCellParticleField()
{
  AtomicBlock3D::dataTransfer = new HemoCellParticleDataTransfer();
    //for (pluint i=0; i<particles.size(); ++i) {
    //    delete particles[i];
    //}
}

HemoCellParticleField& HemoCellParticleField::operator=(HemoCellParticleField const& rhs){
 HemoCellParticleField *copy = new HemoCellParticleField(rhs);
 delete this;
  return *copy;
}

HemoCellParticleField* HemoCellParticleField::clone() const
{
    return new HemoCellParticleField(*this);
}

const vector<vector<unsigned int>> & HemoCellParticleField::get_particles_per_type() { 
    if (!ppt_up_to_date) { update_ppt(); }
    return _particles_per_type;
  }
const map<int,vector<int>> & HemoCellParticleField::get_particles_per_cell() { 
    if (!ppc_up_to_date) { update_ppc(); }
    return _particles_per_cell;
  }
const map<int,vector<int>> & HemoCellParticleField::get_preinlet_particles_per_cell() { 
    update_preinlet_ppc();
    return _preinlet_particles_per_cell;
  }
const map<int,bool> & HemoCellParticleField::get_lpc() { 
    if (!lpc_up_to_date) { update_lpc(); }
    return _lpc;
  }
void HemoCellParticleField::update_lpc() {
  _lpc.clear();
  for (const HemoCellParticle & particle : particles) {
     if (isContainedABS(particle.position, localDomain) && ! particle.fromPreInlet) {
       _lpc[particle.cellId] = true;
     }
  }
  lpc_up_to_date = true;
}
void HemoCellParticleField::update_ppt() {
  _particles_per_type.clear();
  _particles_per_type.resize(cellFields->size());
  
  for (unsigned int i = 0 ; i <  particles.size() ; i++) { 
    _particles_per_type[particles[i].celltype].push_back(i);
  }
  ppt_up_to_date = true;
}
void HemoCellParticleField::update_ppc() {
  _particles_per_cell.clear();
  
  for (unsigned int i = 0 ; i <  particles.size() ; i++) { 
    if (!particles[i].fromPreInlet) {
     insert_ppc(&particles[i],i);
    }
  }
  ppc_up_to_date = true;
}
void HemoCellParticleField::update_preinlet_ppc() {
  _preinlet_particles_per_cell.clear();
  
  for (unsigned int i = 0 ; i <  particles.size() ; i++) {
    if (particles[i].fromPreInlet) {
      insert_preinlet_ppc(&particles[i],i);
    }
  }
  preinlet_ppc_up_to_date = true;
  ppc_up_to_date = false;
  lpc_up_to_date = false;
}

void HemoCellParticleField::addParticle(Box3D domain, HemoCellParticle* particle) {
  //Box3D finalDomain;
  //plint x,y,z;
  HemoCellParticle * local_sparticle;
  hemo::Array<double,3> pos = particle->position;
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();

  cellFields->celltype_per_cell[particle->cellId] = particle->celltype; 
  
  if( this->isContainedABS(pos, this->getBoundingBox()) )
  {
    //check if we have particle already, if so, we must overwrite but not
    //forget to delete the old entry
    if ((!(particles_per_cell.find(particle->cellId) == 
      particles_per_cell.end())) && particles_per_cell.at(particle->cellId)[particle->vertexId] != -1) {
      local_sparticle =  &particles[particles_per_cell.at(particle->cellId)[particle->vertexId]];

      //If our particle is local, do not replace it, envelopes are less important
      if (isContainedABS(local_sparticle->position, localDomain)) {

      } else {
        //We have the particle already, replace it
        *local_sparticle = *particle;
        particle = local_sparticle;

        //Invalidate lpc hemo::Array
        lpc_up_to_date = false;

      }
    } else {
      //new entry
      particles.push_back(*particle);
      particle = &particles.back();
      
      //invalidate ppt
      ppt_up_to_date=false;
      //_particles_per_type[particle->celltype].push_back(particles.size()-1); //last entry
      if(!particle->fromPreInlet) {
        if(this->isContainedABS(pos, localDomain)) {
          _lpc[particle->cellId] = true;
        }
        insert_ppc(particle, particles.size()-1);
      }
    }
    particle->setTag(-1);
  }
}

void inline HemoCellParticleField::insert_ppc(HemoCellParticle* sparticle, unsigned int index) {
  if (_particles_per_cell.find(sparticle->cellId) == _particles_per_cell.end()) {
    _particles_per_cell[sparticle->cellId].resize((*cellFields)[sparticle->celltype]->numVertex);
    for (unsigned int i = 0; i < _particles_per_cell[sparticle->cellId].size(); i++) {
      _particles_per_cell[sparticle->cellId][i] = -1;
    }
  }
  _particles_per_cell.at(sparticle->cellId)[sparticle->vertexId] = index;

}
void inline HemoCellParticleField::insert_preinlet_ppc(HemoCellParticle* sparticle, unsigned int index) {
  if (_preinlet_particles_per_cell.find(sparticle->cellId) == _preinlet_particles_per_cell.end()) {
    _preinlet_particles_per_cell[sparticle->cellId].resize((*cellFields)[sparticle->celltype]->numVertex);
    for (unsigned int i = 0; i < _preinlet_particles_per_cell[sparticle->cellId].size(); i++) {
      _preinlet_particles_per_cell[sparticle->cellId][i] = -1;
    }
  }
  _preinlet_particles_per_cell.at(sparticle->cellId)[sparticle->vertexId] = index;

}

void HemoCellParticleField::removeParticles(plint tag) {
//Almost the same, but we save a lot of branching by making a seperate function

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (particles[i].getTag() == tag) {
      particles[i] = particles.back();
      particles.pop_back();
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
  } 
}

void HemoCellParticleField::removeParticles(Box3D domain, plint tag) {
//Almost the same, but we save a lot of branching by making a seperate function
  Box3D finalDomain;
  
  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (particles[i].getTag() == tag && this->isContainedABS(particles[i].position,finalDomain)) {
      particles[i] = particles.back();
      particles.pop_back();
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
  } 
}

void HemoCellParticleField::removeParticles(Box3D domain) {
//Almost the same, but we save a lot of branching by making a seperate function

  Box3D finalDomain;
  
  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (this->isContainedABS(particles[i].position,finalDomain)) {
      particles[i] = particles.back();
      particles.pop_back();
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
  } 
}

//remove everything outside this domain
void HemoCellParticleField::removeParticles_inverse(Box3D domain) {
//Almost the same, but we save a lot of branching by making a seperate function

  Box3D finalDomain;
  
  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (!this->isContainedABS(particles[i].position,finalDomain)) {
      particles[i] = particles.back();
      particles.pop_back();
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
  } 
}

void HemoCellParticleField::syncEnvelopes() {
  removeParticles_inverse(localDomain);
}


void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<HemoCellParticle*>& found ) 
{
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    for (HemoCellParticle & particle : particles) {
        if (this->isContainedABS(particle.position,domain)) {
            found.push_back(&particle);
        }
    }
}
void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<HemoCellParticle*>& found, pluint type)
{
    
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    //hemo::Array<double,3> pos; 
    const vector<vector<unsigned int>> & particles_per_type = get_particles_per_type();
    if (!(particles_per_type.size() > type)) 
      {return;} 
    else {
      for (const unsigned int i : particles_per_type[type]) {
          if (this->isContainedABS(particles[i].position,domain)) {
              found.push_back(&(particles[i]));
          }
      }
    }
    
}

inline plint HemoCellParticleField::nearestCell(double const pos) const {
  return int(pos + 0.5);
}

inline void HemoCellParticleField::computeGridPosition (
            hemo::Array<double,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const
{
      Dot3D const& location = this->getLocation();
      *iX = nearestCell(position[0]) - location.x;
      *iY = nearestCell(position[1]) - location.y;
      *iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::computeGridPosition (
            hemo::Array<double,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const
{
      Dot3D const& location = this->getLocation();
      iX = nearestCell(position[0]) - location.x;
      iY = nearestCell(position[1]) - location.y;
      iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::issueWarning(HemoCellParticle & p){
	cout << "(HemoCell) (Delete Cells) WARNING! Particle deleted from local domain. This means the whole cell will be deleted!" << endl;
        cout << "\t Particle ID:" << p.cellId << endl;
    cout << "\t Position: " << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << "; vel.: " << p.v[0] << ", " <<  p.v[1] << ", " << p.v[2] << "; force: " << p.force[0] << ", " << p.force[1] << ", " << p.force[2] << endl;
}

int HemoCellParticleField::deleteIncompleteCells(pluint ctype, bool verbose) {
  int deleted = 0;

  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    if (!(*cellFields)[cellFields->celltype_per_cell[cellid]]->deleteIncomplete) { continue; }
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {
        broken = true;
        break;
      }
    }
    if (!broken) {continue;}

    //actually add to tobedeleted list
    bool warningIssued = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {continue;}
      if (particles[particles_per_cell.at(cellid)[i]].celltype != ctype) {break;} //certainly a entry, therefore we check here if it is the right type, if not, exit
      if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].position,localDomain) && verbose && !warningIssued) {
      	issueWarning(particles[particles_per_cell.at(cellid)[i]]);
        warningIssued = true;
      }
      particles[particles_per_cell.at(cellid)[i]].setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function
  removeParticles(1);

  return deleted; 
}

int HemoCellParticleField::deleteIncompleteCells(bool verbose) {
  int deleted = 0;

  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    if (!(*cellFields)[cellFields->celltype_per_cell[cellid]]->deleteIncomplete) { continue; }
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {
        broken = true;
        break;
      }
    }
    if (!broken) {continue;}

    //actually add to tobedeleted list
    bool warningIssued = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {continue;}
      if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].position,localDomain) && verbose && !warningIssued) {
		issueWarning(particles[particles_per_cell.at(cellid)[i]]);
        warningIssued = true;
      }
      particles[particles_per_cell.at(cellid)[i]].setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function

  removeParticles(1);

  return deleted; 
}

void HemoCellParticleField::setlocalDomain(Box3D & localDomain_) {
  localDomain = localDomain_;
  localDomain.x0 -= this->getLocation().x;
  localDomain.x1 -= this->getLocation().x;
  localDomain.y0 -= this->getLocation().y;
  localDomain.y1 -= this->getLocation().y;
  localDomain.z0 -= this->getLocation().z;
  localDomain.z1 -= this->getLocation().z;
}


void HemoCellParticleField::advanceParticles() {
  for(HemoCellParticle & particle:particles){
    particle.advance();
  }
  removeParticles(1);
}

void HemoCellParticleField::separateForceVectors() {
  //Also save the total force, therfore recalculate in advance
  applyConstitutiveModel();

  for (HemoCellParticle & sparticle : particles) {
    //Save Total Force
    sparticle.force_total = sparticle.force + sparticle.force_repulsion;

    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want
    ////TODO this can leak if particle is deleted between seperate and unify,
    //rewrite to reference
    sparticle.force_volume = new hemo::Array<double,3>({0.0,0.0,0.0});
    sparticle.force_link = new hemo::Array<double,3>({0.0,0.0,0.0});
    sparticle.force_area = new hemo::Array<double,3>({0.0,0.0,0.0});
    sparticle.force_bending = new hemo::Array<double,3>({0.0,0.0,0.0});
    sparticle.force_visc = new hemo::Array<double,3>({0.0,0.0,0.0});
    sparticle.force_inner_link = new hemo::Array<double,3>({0.0,0.0,0.0});
  }
}

void HemoCellParticleField::unifyForceVectors() {
  for (HemoCellParticle& sparticle : particles) {
    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want
    delete sparticle.force_volume;
    delete sparticle.force_link;
    delete sparticle.force_area;
    delete sparticle.force_bending;
    delete sparticle.force_visc;
    delete sparticle.force_inner_link;
    sparticle.repoint_force_vectors();
  }
}

void HemoCellParticleField::applyConstitutiveModel(bool forced) {

  deleteIncompleteCells();

  map<int,vector<HemoCellParticle*>> * ppc_new = new map<int,vector<HemoCellParticle*>>();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  const map<int,bool> & lpc = get_lpc();
  
  //Fill it here, probably needs optimization, ah well ...
  for (const auto & pair : particles_per_cell) {
    const int & cid = pair.first;
    const vector<int> & cell = pair.second; 
    (*ppc_new)[cid].resize(cell.size());
    for (unsigned int i = 0 ; i < cell.size() ; i++) {
      if (cell[i] == -1) {
        (*ppc_new)[cid][i] = NULL;
      } else {
        (*ppc_new)[cid][i] = &particles[cell[i]];
      }
    }
  }
  
  for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++) {
    if ((*cellFields).hemocell.iter % (*cellFields)[ctype]->timescale == 0 || forced) {
      vector<HemoCellParticle*> found;
      findParticles(getBoundingBox(),found,ctype);
      for (HemoCellParticle* particle : found) {
        particle->force = {0.,0.,0.};        
      }
      (*cellFields)[ctype]->mechanics->ParticleMechanics(*ppc_new,lpc,ctype);
    }
  }
  
  delete ppc_new;
  
}

void HemoCellParticleField::applyRepulsionForce(bool forced) {
  const double r_const = cellFields->repulsionConstant;
  const double r_cutoff = cellFields->repulsionCutoff;
  
  vector<HemoCellParticle *> found;
  Box3D localDomainplusOne = localDomain;
  localDomainplusOne.enlarge(1);
  findParticles(localDomainplusOne,found);
  
  if (cellFields->hemocell.iter % cellFields->repulsionTimescale == 0 || forced) {  
    for (HemoCellParticle * particle : found) {
      particle->force_repulsion = {0.,0.,0.};
    }
    
    if (found.size() > 1) {
      for (unsigned int i = 0 ; i < found.size() - 1 ; i ++) {
        for (unsigned int j = i + 1 ; j < found.size() ; j++ ) {
          if (found[i]->cellId == found[j]->cellId) { continue; } //No incell repulsion
          const hemo::Array<double,3> dv = found[i]->position - found[j]->position;
          const double distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
          if (distance < r_cutoff) {
            const hemo::Array<double, 3> rfm = r_const * (1/(distance/r_cutoff))  * (dv/distance);
            found[i]->force_repulsion = found[i]->force_repulsion + rfm;
            found[j]->force_repulsion = found[j]->force_repulsion - rfm;
          } 
        }
      }
    }
  }
}

void HemoCellParticleField::interpolateFluidVelocity(Box3D domain) {
  vector<HemoCellParticle*> localParticles;
  findParticles(domain,localParticles);
  //TODO, remove casting
  HemoCellParticle * sparticle;
  //Prealloc is nice
  hemo::Array<double,3> velocity;
  plb::Array<double,3> velocity_comp;

  for (pluint i = 0; i < localParticles.size(); i++ ) {
    sparticle = localParticles[i];
    if (sparticle->fromPreInlet) { continue; }
    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[sparticle->celltype]->kernelMethod(*atomicLattice,sparticle);

    //We have the kernels, now calculate the velocity of the particles.
    //Palabos developers, sorry for not using a functional...
    velocity = {0.0,0.0,0.0};
    for (pluint j = 0; j < sparticle->kernelLocations.size(); j++) {
      //Yay for direct access
      sparticle->kernelLocations[j]->computeVelocity(velocity_comp);
      velocity += (velocity_comp * sparticle->kernelWeights[j]);
    }
    sparticle->v = velocity;
  }

}

void HemoCellParticleField::spreadParticleForce(Box3D domain) {
  HemoCellParticle * sparticle;
  for (pluint i = 0; i < particles.size(); i++ ) {
    sparticle = &particles[i];
    if (sparticle->fromPreInlet) { continue; }

    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[sparticle->celltype]->kernelMethod(*atomicLattice,sparticle);

    // Capping force to ensure stability -> NOTE: this introduces error!
#ifdef FORCE_LIMIT
    const double force_mag = norm(sparticle->force);
    if(force_mag > param::f_limit)
      sparticle->force *= param::f_limit/force_mag;
#endif

    //Directly change the force on a node , Palabos developers hate this one
    //quick non-functional trick.
    for (pluint j = 0; j < sparticle->kernelLocations.size(); j++) {
      //Yay for direct access
      sparticle->kernelLocations[j]->external.data[0] += ((sparticle->force_repulsion[0] + sparticle->force[0]) * sparticle->kernelWeights[j]);
      sparticle->kernelLocations[j]->external.data[1] += ((sparticle->force_repulsion[1] + sparticle->force[1]) * sparticle->kernelWeights[j]);
      sparticle->kernelLocations[j]->external.data[2] += ((sparticle->force_repulsion[2] + sparticle->force[2]) * sparticle->kernelWeights[j]);
    }

  }
}

HemoCellParticleDataTransfer& HemoCellParticleField::getDataTransfer() {
    return dataTransfer;
}
HemoCellParticleDataTransfer const& HemoCellParticleField::getDataTransfer() const {
    return dataTransfer;
}

std::string HemoCellParticleField::getBlockName() {
    return std::string("HemoParticleField3D");
}

HemoCellFields* HemoCellParticleField::cellFields=0;