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
  #include "octree.h"
  #include "mollerTrumbore.h"
  #include "bindingField.h"
  #include "interiorViscosity.h"
  #include <Eigen3/Eigenvalues>

  namespace hemo { 
  /* *************** class HemoParticleField3D ********************** */

  HemoCellParticleField::HemoCellParticleField(plint nx, plint ny, plint nz)
      : AtomicBlock3D(nx,ny,nz,0), particleDataTransfer(*(new HemoCellParticleDataTransfer()))
  { 
      boundingBox = Box3D(0,this->getNx()-1, 0, this->getNy()-1, 0, this->getNz()-1);
      dataTransfer = &particleDataTransfer;
      particleDataTransfer.setBlock(*this);
      AddOutputMap(); 
  }

  HemoCellParticleField::HemoCellParticleField(HemoCellParticleField const& rhs)
      : AtomicBlock3D(rhs), particleDataTransfer(*(new HemoCellParticleDataTransfer()))
  {
      boundingBox = Box3D(0,this->getNx()-1, 0, this->getNy()-1, 0, this->getNz()-1);
      dataTransfer = &particleDataTransfer;
      particleDataTransfer.setBlock(*this);
      for (const HemoCellParticle & particle : rhs.particles) {
        addParticle(particle.sv);
      }
      ppc_up_to_date = false;
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      pg_up_to_date = false;
      AddOutputMap();
  }

  HemoCellParticleField::~HemoCellParticleField()
  {
    //AtomicBlock3D::dataTransfer = new HemoCellParticleDataTransfer();
    if (particle_grid) {
      delete[] particle_grid;
      particle_grid = 0;
    }
    if(particle_grid_size) {
      delete[] particle_grid_size;
      particle_grid_size = 0;
    }
    
    // Sanitize for MultiBlockLattice destructor (releasememory). It can't handle releasing non-background dynamics that are not singular
    if (global.enableInteriorViscosity) {
      for (const Dot3D & internalPoint : internalPoints ) {
        atomicLattice->get(internalPoint.x,internalPoint.y,internalPoint.z);
        atomicLattice->get(internalPoint.x,internalPoint.y,internalPoint.z).attributeDynamics(&atomicLattice->getBackgroundDynamics());
      }
      InteriorViscosityHelper::get(*cellFields).empty(*this);
    }
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

  const map<int,bool> & HemoCellParticleField::get_lpc() { 
      if (!lpc_up_to_date) { update_lpc(); }
      return _lpc;
    }
  void HemoCellParticleField::update_lpc() {
    _lpc.clear();
    for (const HemoCellParticle & particle : particles) {
       if (isContainedABS(particle.sv.position, localDomain)) {
         _lpc[particle.sv.cellId] = true;
       }
    }
    lpc_up_to_date = true;
  }
  void HemoCellParticleField::update_ppt() {
    _particles_per_type.clear();
    _particles_per_type.resize(cellFields->size());
    
    for (unsigned int i = 0 ; i <  particles.size() ; i++) { 
      _particles_per_type[particles[i].sv.celltype].push_back(i);
    }
    ppt_up_to_date = true;
  }
  void HemoCellParticleField::update_ppc() {
    _particles_per_cell.clear();
    
    for (unsigned int i = 0 ; i <  particles.size() ; i++) { 
       insert_ppc(&particles[i],i);
    }
    ppc_up_to_date = true;
  }

  void HemoCellParticleField::update_pg() {
    //Check if map exists, otherwise create
    if (!particle_grid) {
      particle_grid = new hemo::Array<unsigned int, 10>[this->atomicLattice->getNx()*this->atomicLattice->getNy()*this->atomicLattice->getNz()];
    }
    if (!particle_grid_size) {
      particle_grid_size = new unsigned int[this->atomicLattice->getNx()*this->atomicLattice->getNy()*this->atomicLattice->getNz()];
    }
    if (!this->atomicLattice) {
      return;
    }
    
    memset(particle_grid_size,0,sizeof(unsigned int)*this->atomicLattice->getNx()*this->atomicLattice->getNy()*this->atomicLattice->getNz());
    Dot3D const& location = this->atomicLattice->getLocation();
    hemo::Array<T,3> * pos;
    
    for (unsigned int i = 0 ; i <  particles.size() ; i++) {
      pos = &particles[i].sv.position;
      int x = pos->operator[](0)-location.x+0.5;
      int y = pos->operator[](1)-location.y+0.5;
      int z = pos->operator[](2)-location.z+0.5;
      if ((x >= 0) && (x < this->atomicLattice->getNx()) &&
  	(y >= 0) && (y < this->atomicLattice->getNy()) &&
  	(z >= 0) && (z < this->atomicLattice->getNz()) ) 
      {
        unsigned int index = grid_index(x,y,z);
        particle_grid[index][particle_grid_size[index]] = i;
        particle_grid_size[index]++;
      }
    }
    pg_up_to_date = true;
  }

  void HemoCellParticleField::addParticle(HemoCellParticle* particle) {
    addParticle(particle->sv);
  }  
  void HemoCellParticleField::addParticle(const HemoCellParticle::serializeValues_t & sv) {
    HemoCellParticle * local_sparticle, * particle;
    const hemo::Array<T,3> & pos = sv.position;
    const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();

    if( this->isContainedABS(pos, this->getBoundingBox()) )
    {
      //check if we have particle already, if so, we must overwrite but not
      //forget to delete the old entry
      if ((!(particles_per_cell.find(sv.cellId) == particles_per_cell.end()))) { 
        if (particles_per_cell.at(sv.cellId)[sv.vertexId] != -1) {
          local_sparticle =  &particles[particles_per_cell.at(sv.cellId)[sv.vertexId]];

          //If our particle is local, do not replace it, envelopes are less important
          if (isContainedABS(local_sparticle->sv.position, localDomain)) {
            return;
          } else {
            //We have the particle already, replace it
            local_sparticle->sv = sv;
            particle = local_sparticle;
            particle->setTag(-1);

            //Invalidate lpc hemo::Array
            lpc_up_to_date = false;
            pg_up_to_date = false;

          }
        } else {
          goto outer_else;
        }
      } else {
  outer_else:
        //new entry
        particles.emplace_back(sv);
        particle = &particles.back();
        
        //invalidate ppt
        ppt_up_to_date=false;
          if(this->isContainedABS(pos, localDomain)) {
            _lpc[particle->sv.cellId] = true;
          }
          if (ppc_up_to_date) { //Otherwise its rebuild anyway
           insert_ppc(particle, particles.size()-1);
          }
        
        if (pg_up_to_date) {
          Dot3D const& location = this->atomicLattice->getLocation();
          hemo::Array<T,3>  & pos = particle->sv.position;
          int x = pos[0]-location.x+0.5;
          int y = pos[1]-location.y+0.5;
          int z = pos[2]-location.z+0.5;
          if ((x >= 0) && (x <= this->atomicLattice->getNx()) &&
              (y >= 0) && (y <= this->atomicLattice->getNy()) &&
              (z >= 0) && (z <= this->atomicLattice->getNz()) ) 
          {
            unsigned int index = grid_index(x,y,z);
            particle_grid[index][particle_grid_size[index]] = particles.size()-1;
            particle_grid_size[index]++;
          }
        }
      }
    }
  }

  void inline HemoCellParticleField::insert_ppc(HemoCellParticle* sparticle, unsigned int index) {
    if (_particles_per_cell.find(sparticle->sv.cellId) == _particles_per_cell.end()) {
      _particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex,-1);
    }
    _particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;

  }
  void inline HemoCellParticleField::insert_preinlet_ppc(HemoCellParticle* sparticle, unsigned int index) {
    if (_preinlet_particles_per_cell.find(sparticle->sv.cellId) == _preinlet_particles_per_cell.end()) {
      _preinlet_particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex);
      for (unsigned int i = 0; i < _preinlet_particles_per_cell[sparticle->sv.cellId].size(); i++) {
        _preinlet_particles_per_cell[sparticle->sv.cellId][i] = -1;
      }
    }
    _preinlet_particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;

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
      pg_up_to_date = false;
    } 
  }

  void HemoCellParticleField::removeParticles(Box3D domain, plint tag) {
  //Almost the same, but we save a lot of branching by making a seperate function
    Box3D finalDomain;
    
    intersect(domain, this->getBoundingBox(), finalDomain);

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0 ; i < particles.size() ; i++) {
      if (particles[i].getTag() == tag && this->isContainedABS(particles[i].sv.position,finalDomain)) {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size) {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
    } 
  }

  void HemoCellParticleField::removeParticles(Box3D domain) {
  //Almost the same, but we save a lot of branching by making a seperate function

    Box3D finalDomain;
    
    intersect(domain, this->getBoundingBox(), finalDomain);

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0 ; i < particles.size() ; i++) {
      if (this->isContainedABS(particles[i].sv.position,finalDomain)) {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size) {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
    } 
  }

  //remove everything outside this domain
  void HemoCellParticleField::removeParticles_inverse(Box3D domain) {
  //Almost the same, but we save a lot of branching by making a seperate function

    Box3D finalDomain;
    
    intersect(domain, this->getBoundingBox(), finalDomain);

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0 ; i < particles.size() ; i++) {
      if (!this->isContainedABS(particles[i].sv.position,finalDomain)) {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size) {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
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
          if (this->isContainedABS(particle.sv.position,domain)) {
              found.push_back(&particle);
          }
      }
  }

  void HemoCellParticleField::findParticles (
          Box3D domain, std::vector<const HemoCellParticle*>& found ) const
  {
      found.clear();
      PLB_ASSERT( contained(domain, this->getBoundingBox()) );
      for (const HemoCellParticle & particle : particles) {
          if (this->isContainedABS(particle.sv.position,domain)) {
              found.push_back(&particle);
          }
      }
  }
  void HemoCellParticleField::findParticles (
          Box3D domain, std::vector<HemoCellParticle*>& found, pluint type)
  {
      
      found.clear();
      PLB_ASSERT( contained(domain, this->getBoundingBox()) );
      //hemo::Array<T,3> pos; 
      const vector<vector<unsigned int>> & particles_per_type = get_particles_per_type();
      if (!(particles_per_type.size() > type)) 
        {return;} 
      else {
        for (const unsigned int i : particles_per_type[type]) {
            if (this->isContainedABS(particles[i].sv.position,domain)) {
                found.push_back(&(particles[i]));
            }
        }
      }
      
  }

  inline plint HemoCellParticleField::nearestCell(T const pos) const {
    return int(pos + 0.5);
  }

  inline void HemoCellParticleField::computeGridPosition (
              hemo::Array<T,3> const& position,
                      plint* iX, plint* iY, plint* iZ ) const
  {
        Dot3D const& location = this->getLocation();
        *iX = nearestCell(position[0]) - location.x;
        *iY = nearestCell(position[1]) - location.y;
        *iZ = nearestCell(position[2]) - location.z;
  }

  void HemoCellParticleField::computeGridPosition (
              hemo::Array<T,3> const& position,
                      plint& iX, plint& iY, plint& iZ ) const
  {
        Dot3D const& location = this->getLocation();
        iX = nearestCell(position[0]) - location.x;
        iY = nearestCell(position[1]) - location.y;
        iZ = nearestCell(position[2]) - location.z;
  }

  void HemoCellParticleField::issueWarning(HemoCellParticle & p){
  	cout << "(HemoCell) (Delete Cells) WARNING! Particle deleted from local domain. This means the whole cell will be deleted!" << endl;
          cout << "\t Particle ID:" << p.sv.cellId << endl;
      cout << "\t Position: " << p.sv.position[0] << ", " << p.sv.position[1] << ", " << p.sv.position[2] << "; vel.: " << p.sv.v[0] << ", " <<  p.sv.v[1] << ", " << p.sv.v[2] << "; force: " << p.sv.force[0] << ", " << p.sv.force[1] << ", " << p.sv.force[2] << endl;
  }

  int HemoCellParticleField::deleteIncompleteCells(pluint ctype, bool verbose) {
    int deleted = 0;

    const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
    //Warning, TODO, high complexity, should be rewritten 
    //For now abuse tagging and the remove function
    for ( const auto &lpc_it : particles_per_cell ) {
      int cellid = lpc_it.first;
      bool broken = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
        if (particles_per_cell.at(cellid)[i] == -1) {
          broken = true;
          break;
        }
      }
      if (!broken) {continue;}

      bool warningIssued = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
        if (particles_per_cell.at(cellid)[i] == -1) {continue;}

        //issue warning
        if (verbose) {
          if (!warningIssued) {
            if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].sv.position,localDomain)) {
                    issueWarning(particles[particles_per_cell.at(cellid)[i]]);
              warningIssued = true;
            }
          }
        }
        
        //actually add to tobedeleted list
        particles[particles_per_cell.at(cellid)[i]].setTag(1);
        deleted++;
      }
    } 

    //We have our list, now abuse the removeall function
    removeParticles(1);

    return deleted; 
  }

  int HemoCellParticleField::deleteIncompleteCells(const bool verbose) {
    int deleted = 0;
    const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
    //Warning, TODO, high complexity, should be rewritten 
    //For now abuse tagging and the remove function
    for ( const auto &lpc_it : particles_per_cell ) {
      int cellid = lpc_it.first;
      bool broken = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
        if (particles_per_cell.at(cellid)[i] == -1) {
          broken = true;
          break;
        }
      }
      if (!broken) {continue;}

      bool warningIssued = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
        if (particles_per_cell.at(cellid)[i] == -1) {continue;}

        //issue warning
        if (verbose) {
          if (!warningIssued) {
            if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].sv.position,localDomain)) {
                    issueWarning(particles[particles_per_cell.at(cellid)[i]]);
              warningIssued = true;
            }
          }
        }
        
        //actually add to tobedeleted list
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
      //By lack of better place, check if it is on a boundary, if so, delete it
      plb::Box3D const box = atomicLattice->getBoundingBox();
      plb::Dot3D const& location = atomicLattice->getLocation();
      plint x = (particle.sv.position[0]-location.x)+0.5;
      plint y = (particle.sv.position[1]-location.y)+0.5;
      plint z = (particle.sv.position[2]-location.z)+0.5;

      if ((x >= box.x0) && (x <= box.x1) &&
  	(y >= box.y0) && (y <= box.y1) &&
  	(z >= box.z0) && (z <= box.z1)) {
        if (atomicLattice->get(x,y,z).getDynamics().isBoundary()) {
          particle.tag = 1;
        }
      }
    }
    removeParticles(1);
    
    lpc_up_to_date = false;
    pg_up_to_date = false;
  }

  void HemoCellParticleField::separateForceVectors() {
    //Also save the total force, therfore recalculate in advance
    applyConstitutiveModel();

    for (HemoCellParticle & sparticle : particles) {
      //Save Total Force
      sparticle.force_total = sparticle.sv.force + sparticle.sv.force_repulsion;

      //Just repoint all possible outputs for now //TODO only repoint the ones we
      //want

      sparticle.force_volume = new hemo::Array<T,3>({0.0,0.0,0.0});
      allocated_for_output.push_back(sparticle.force_volume);
      sparticle.force_link = new hemo::Array<T,3>({0.0,0.0,0.0});
      allocated_for_output.push_back(sparticle.force_link);
      sparticle.force_area = new hemo::Array<T,3>({0.0,0.0,0.0});
      allocated_for_output.push_back(sparticle.force_area);
      sparticle.force_bending = new hemo::Array<T,3>({0.0,0.0,0.0});
      allocated_for_output.push_back(sparticle.force_bending);
      sparticle.force_visc = new hemo::Array<T,3>({0.0,0.0,0.0});
      allocated_for_output.push_back(sparticle.force_visc);
      sparticle.force_inner_link = new hemo::Array<T,3>({0.0,0.0,0.0});
      allocated_for_output.push_back(sparticle.force_inner_link);
    }
  }

  void HemoCellParticleField::updateResidenceTime(unsigned int rtime) {
    for (HemoCellParticle & sparticle : particles) {
      sparticle.sv.restime += rtime;
    }
  }


  void HemoCellParticleField::unifyForceVectors() {
    for (const hemo::Array<T,3>* mem : allocated_for_output) {
      delete mem;
    }
    allocated_for_output.clear();
    for (HemoCellParticle& sparticle : particles) {
      sparticle.repoint_force_vectors();
    }
  }

  void HemoCellParticleField::applyConstitutiveModel(bool forced) {
    map<int,vector<HemoCellParticle*>> * ppc_new = new map<int,vector<HemoCellParticle*>>();
    const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
    map<int,bool> lpc;
    //Fill it here, probably needs optimization, ah well ...
    for (const auto & pair : particles_per_cell) {
      const int & cid = pair.first;
      const vector<int> & cell = pair.second; 
      (*ppc_new)[cid].resize(cell.size());
      for (unsigned int i = 0 ; i < cell.size() ; i++) {
        if (cell[i] == -1) {
          (*ppc_new).erase(cid); //not complete, remove entry
          goto no_add_lpc;
        } else {
          (*ppc_new)[cid][i] = &particles[cell[i]];
        }
      }
      lpc[cid]=true;
      no_add_lpc:;
    }
    
    for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++) {
      if ((*cellFields).hemocell.iter % (*cellFields)[ctype]->timescale == 0 || forced) {
        vector<HemoCellParticle*> found;
        findParticles(getBoundingBox(),found,ctype);
        if (found.size() > 0) {
          //only reset forces when the forces actually point at it.
          if (found[0]->force_area == &found[0]->sv.force) {
            for (HemoCellParticle* particle : found) {
              particle->sv.force = {0.,0.,0.};
  #ifdef INTERIOR_VISCOSITY
              particle->normalDirection = {0., 0., 0.};
  #endif
            }
          }
        }
        (*cellFields)[ctype]->mechanics->ParticleMechanics(*ppc_new,lpc,ctype);
      }
    }
    
    delete ppc_new;
    
  }

  #define inner_loop \
    const int & l_index = grid_index(x,y,z); \
    const int & n_index = grid_index(xx,yy,zz); \
    for (unsigned int i = 0; i < particle_grid_size[l_index];i++){ \
      for (unsigned int j = 0; j < particle_grid_size[n_index];j++){ \
        HemoCellParticle & lParticle = particles[particle_grid[l_index][i]]; \
        HemoCellParticle & nParticle = particles[particle_grid[n_index][j]]; \
        if (&nParticle == &lParticle) { continue; } \
        if (lParticle.sv.cellId == nParticle.sv.cellId) { continue; } \
        const hemo::Array<T,3> dv = lParticle.sv.position - nParticle.sv.position; \
        const T distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]); \
        if (distance < r_cutoff) { \
          const hemo::Array<T, 3> rfm = r_const * (1/(distance/r_cutoff))  * (dv/distance); \
          lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm; \
          nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - rfm; \
        } \
      } \
    }

  void HemoCellParticleField::applyRepulsionForce(bool forced) {
    const T r_const = cellFields->repulsionConstant;
    const T r_cutoff = cellFields->repulsionCutoff;
    if(!pg_up_to_date) {
      update_pg();
    }
    for (HemoCellParticle & particle : particles) {
      particle.sv.force_repulsion = {0.,0.,0.};
    }
    
    for (int x = 0; x < atomicLattice->getNx()-1; x++) {
      for (int y = 0; y < atomicLattice->getNy()-1; y++) {
        for (int z = 0; z < atomicLattice->getNz()-1; z++) {
          //Manual finding, we could make a map, but for now this should be fast enough
          //0, 0, 0 //0, 1, 0 //0, 0, 1 //0, 1, 1
          int xx = x;
          for (int yy = y; yy <= y+1; yy++) {
            for (int zz = z; zz <= z+1; zz++) {
              inner_loop
            }
          }
          //1, 0, 0
          //1, 1, 0
          //1, 0, 1
          //1, 1, 1
          //1, 0,-1
          //1,-1, 0
          //1,-1,-1
          //1, 1,-1
          //1,-1, 1
          xx = x+1;
          for (int yy = y-1; yy <= y+1; yy++) {
            for(int zz = z-1; zz <= z+1; zz++) {
              if (yy < 0) {continue;}
              if (zz < 0) {continue;}
              inner_loop
            }
          }
          xx = x;
          int yy = y+1, zz = z-1;
          if (zz<0) { continue; }
          //0, 1,-1
          inner_loop
          
        }
      }
    }
  }

  #ifdef INTERIOR_VISCOSITY
  void HemoCellParticleField::internalGridPointsMembrane(Box3D domain) {
    // This could be done less complex I guess?
    for (const HemoCellParticle & particle : particles) { // Go over each particle
       if (!(*cellFields)[particle.sv.celltype]->doInteriorViscosity) { continue; }

      for (unsigned int i = 0; i < particle.kernelCoordinates.size(); i++) {
        const hemo::Array<T, 3> latPos = particle.kernelCoordinates[i]-(particle.sv.position-atomicLattice->getLocation());
        const hemo::Array<T, 3> & normalP = particle.normalDirection;

        if (computeLength(latPos) > (*cellFields)[particle.sv.celltype]->mechanics->cellConstants.edge_mean_eq) {continue;}
        
        T dot1 = hemo::dot(latPos, normalP);

        if (dot1 < 0.) {  // Node is inside
          InteriorViscosityHelper::get(*cellFields).add(*this, {particle.kernelCoordinates[i][0],
                  particle.kernelCoordinates[i][1],
                  particle.kernelCoordinates[i][2]},
                  (*cellFields)[particle.sv.celltype]->interiorViscosityTau);
          particle.kernelLocations[i]->attributeDynamics((*cellFields)[particle.sv.celltype]->innerViscosityDynamics);
        } else {  // Node is outside
          InteriorViscosityHelper::get(*cellFields).remove(*this, {particle.kernelCoordinates[i][0],
                                                                  particle.kernelCoordinates[i][1],
                                                                  particle.kernelCoordinates[i][2]});
          particle.kernelLocations[i]->attributeDynamics(&atomicLattice->getBackgroundDynamics());
        }
      }
    }
  }

  // For performance reason, this is only executed once every n iterations to make
  // sure that there are no higher viscosity grid points left after substantial movement
  void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain) {
    // Reset all the lattice points to the orignal relaxation parameter
    for (const Dot3D & internalPoint : internalPoints ) {
      atomicLattice->get(internalPoint.x,internalPoint.y,internalPoint.z);
      atomicLattice->get(internalPoint.x,internalPoint.y,internalPoint.z).attributeDynamics(&atomicLattice->getBackgroundDynamics());
    }
    InteriorViscosityHelper::get(*cellFields).empty(*this);
    
    for (const auto & pair : get_lpc()) { // Go over each cell?
      const int & cid = pair.first;
      const vector<int> & cell = get_particles_per_cell().at(cid);
      const pluint ctype = particles[cell[0]].sv.celltype;

      // Plt and Wbc now have normal tau internal, so we don't have
      // to raycast these particles
      if (!(*cellFields)[ctype]->doInteriorViscosity) {
        continue;
      }
      
      hemo::OctreeStructCell octCell(3, 1, 30,
                                    (*cellFields)[ctype]->mechanics->cellConstants.triangle_list,
                                    particles, cell);

      std::set<Array<plint,3>> innerNodes;
      octCell.findInnerNodes(atomicLattice,particles,cell,innerNodes);
      for (const Array<plint,3> & node : innerNodes) {
        InteriorViscosityHelper::get(*cellFields).add(*this, {node[0],node[1],node[2]},(*cellFields)[ctype]->interiorViscosityTau );
        atomicLattice->get(node[0],node[1],node[2]).attributeDynamics((*cellFields)[ctype]->innerViscosityDynamics);
      }
    }
  }
  #else
  void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain) {
    pcout << "(HemoCellParticleField) (Error) findInternalParticleGridPoints called, but INTERIOR_VISCOSITY not defined, exiting..." << endl;
    exit(1);
  }
  void HemoCellParticleField::internalGridPointsMembrane(Box3D domain) {
    pcout << "(HemoCellParticleField) (Error) internalGridPointsMembrane called, but INTERIOR_VISCOSITY not defined, exiting..." << endl;
    exit(1);
  }
  #endif

  void HemoCellParticleField::interpolateFluidVelocity(Box3D domain) {
    //Prealloc is nice
    hemo::Array<T,3> velocity;
    plb::Array<T,3> velocity_comp;

    for (HemoCellParticle &particle:particles) {

      //Clever trick to allow for different kernels for different particle types.
      //(*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice,particle);

      //We have the kernels, now calculate the velocity of the particles.
      //Palabos developers, sorry for not using a functional...
      velocity = {0.0,0.0,0.0};
      for (pluint j = 0; j < particle.kernelLocations.size(); j++) {
        //Yay for direct access
        particle.kernelLocations[j]->computeVelocity(velocity_comp);
        velocity += (velocity_comp * particle.kernelWeights[j]);
      }
      particle.sv.v = velocity;
    }

  }

  void HemoCellParticleField::spreadParticleForce(Box3D domain) {
    for( HemoCellParticle &particle:particles) {

      //Clever trick to allow for different kernels for different particle types.
      (*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice,particle);

      // Capping force to ensure stability -> NOTE: this introduces error!
  #ifdef FORCE_LIMIT
      const T force_mag = norm(particle.sv.force);
      if(force_mag > param::f_limit)
        particle.sv.force *= param::f_limit/force_mag;
  #endif

      //Directly change the force on a node , Palabos developers hate this one
      //quick non-functional trick.
      for (pluint j = 0; j < particle.kernelLocations.size(); j++) {
        //Yay for direct access
        particle.kernelLocations[j]->external.data[0] += ((particle.sv.force_repulsion[0] + particle.sv.force[0]) * particle.kernelWeights[j]);
        particle.kernelLocations[j]->external.data[1] += ((particle.sv.force_repulsion[1] + particle.sv.force[1]) * particle.kernelWeights[j]);
        particle.kernelLocations[j]->external.data[2] += ((particle.sv.force_repulsion[2] + particle.sv.force[2]) * particle.kernelWeights[j]);
      }

    }
  }

  void HemoCellParticleField::populateBoundaryParticles() {
    //qqw2qswws32 <- Greatly appreciated input of Gábor

    for (int x = 0; x < this->atomicLattice->getNx()-1; x++) {
      for (int y = 0; y < this->atomicLattice->getNy()-1; y++) {
        for (int z = 0; z < this->atomicLattice->getNz()-1; z++) {
          if (this->atomicLattice->get(x,y,z).getDynamics().isBoundary()) {
            for (int xx = x-1; xx <= x+1; xx++) {
              if (xx < 0 || xx > this->atomicLattice->getNx()-1) {continue;}
              for (int yy = y-1; yy <= y+1; yy++) {
                if (yy < 0 || yy > this->atomicLattice->getNy()-1) {continue;}
                for (int zz = z-1; zz <= z+1; zz++) {
                  if (zz < 0 || zz > this->atomicLattice->getNz()-1) {continue;}
                  if (!this->atomicLattice->get(xx,yy,zz).getDynamics().isBoundary()) {
                    boundaryParticles.push_back({x,y,z});       
                    goto end_inner_loop;
                  }
                }
              }
            }
          }
  end_inner_loop:;
        }
      }
    }
  }

  void HemoCellParticleField::applyBoundaryRepulsionForce() {
      if(!pg_up_to_date) {
      update_pg();
    }
    const T & br_cutoff = cellFields->boundaryRepulsionCutoff;
    const T & br_const = cellFields->boundaryRepulsionConstant;
    for (Dot3D & b_particle : boundaryParticles) {
      for (int x = b_particle.x-1; x <= b_particle.x+1; x++) {
        if (x < 0 || x > this->atomicLattice->getNx()-1) {continue;}
        for (int y = b_particle.y-1; y <= b_particle.y+1; y++) {
          if (y < 0 || y > this->atomicLattice->getNy()-1) {continue;}
          for (int z = b_particle.z-1; z <= b_particle.z+1; z++) {
            if (z < 0 || z > this->atomicLattice->getNz()-1) {continue;}
            const int & index = grid_index(x,y,z);
            for (unsigned int i = 0 ; i < particle_grid_size[index] ; i++ ) {
              HemoCellParticle & lParticle = particles[particle_grid[index][i]];
              const hemo::Array<T,3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation()); 
              const T distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]); 
              if (distance < br_cutoff) { 
                const hemo::Array<T, 3> rfm = br_const * (1/(distance/br_cutoff))  * (dv/distance);
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm; 
              } 
            }
          }
        }
      }
    }
  }

  void HemoCellParticleField::populateBindingSites() {
    //qqw2qswws32 <- Greatly appreciated input of Gábor

    for (int x = 0; x < this->atomicLattice->getNx()-1; x++) {
      for (int y = 0; y < this->atomicLattice->getNy()-1; y++) {
        for (int z = 0; z < this->atomicLattice->getNz()-1; z++) {
          if (this->atomicLattice->get(x,y,z).getDynamics().isBoundary()) {
            for (int xx = x-1; xx <= x+1; xx++) {
              if (xx < 0 || xx > this->atomicLattice->getNx()-1) {continue;}
              for (int yy = y-1; yy <= y+1; yy++) {
                if (yy < 0 || yy > this->atomicLattice->getNy()-1) {continue;}
                for (int zz = z-1; zz <= z+1; zz++) {
                  if (zz < 0 || zz > this->atomicLattice->getNz()-1) {continue;}
                  if (!this->atomicLattice->get(xx,yy,zz).getDynamics().isBoundary()) {
                    bindingFieldHelper::get(*cellFields).add(*this, {x,y,z});
                    goto end_inner_loop;
                  }
                }
              }
            }
          }
  end_inner_loop:;
        }
      }
    }
  }


  T HemoCellParticleField::eigenValueFromCell(plb::Cell<T,DESCRIPTOR> & cell) {
      T n = 3;
      plb::Array<T,SymmetricTensor<T,DESCRIPTOR>::n> element;
      cell.computePiNeq(element);
      T omega     = cell.getDynamics().getOmega();
      T rhoBar    = cell.getDynamics().computeRhoBar(cell);
      T prefactor = - omega * DESCRIPTOR<T>::invCs2 *
                   DESCRIPTOR<T>::invRho(rhoBar) / (T)2;
          for (int iTensor=0; iTensor<SymmetricTensor<T,DESCRIPTOR>::n; ++iTensor) {
              element[iTensor] *= prefactor;
          }

      Array<Array<T,3>,3> S;  // Strain-rate tensor (symmetric).
      S[0][0] = element[0]; //s[xx];
      S[0][1] = element[1]; //s[xy];
      S[0][2] = element[2]; //s[xz];
                  
      S[1][0] = element[1]; //s[yx];
      S[1][1] = element[3]; //s[yy];
      S[1][2] = element[4]; //s[yz];

      S[2][0] = element[2]; //s[zx];
      S[2][1] = element[4]; //s[zy];
      S[2][2] = element[6]; //s[zz];
      
      Eigen::Matrix<T,3,3> A;
      for (plint i = 0; i < 3; i++) {
          for (plint j = 0; j < 3; j++) {
              A(i, j) = S[i][j];
          }
       }

      bool computeEigenvectors = false;
      Eigen::EigenSolver<Eigen::Matrix<T,3,3> > es(A, computeEigenvectors);
      std::vector<T> lambda(3);
      lambda[0] = std::real(es.eigenvalues()[0]);
      lambda[1] = std::real(es.eigenvalues()[1]);
      lambda[2] = std::real(es.eigenvalues()[2]);
      std::sort(lambda.begin(), lambda.end());
                  
  //    Array<Array<T,3>,3> x;  // Eigenvectors of S.
  //    Array<T,3> d;           // Eigenvalues of S.
  //    Eigen::eigenDecomposition(S, x, d);
  //    std::vector<T> lambda(3);
  //    lambda[0] = d[0];
  //    lambda[1] = d[1];
  //    lambda[2] = d[2];
   //   std::sort(lambda.begin(), lambda.end());
      T tresca = lambda[0]-lambda[2]/2;
      return tresca;
  }

  void HemoCellParticleField::solidifyCells() {
  #ifdef SOLIDIFY_MECHANICS
    for (HemoCellField * type : cellFields->cellFields) {
      ppc_up_to_date = false;
      if(type->doSolidifyMechanics) {
        type->mechanics->solidifyMechanics(get_particles_per_cell(),particles, this->atomicLattice, this->CEPAClattice, type->ctype, *this);
      }
    }
    removeParticles(1);
    if(!pg_up_to_date) {
      update_pg();
    }

    for (const Dot3D & b_particle : bindingSites) {
      for (int x = b_particle.x-1; x <= b_particle.x+1; x++) {
        if (x < 0 || x > this->atomicLattice->getNx()-1) {continue;}
        for (int y = b_particle.y-1; y <= b_particle.y+1; y++) {
          if (y < 0 || y > this->atomicLattice->getNy()-1) {continue;}
  	for (int z = b_particle.z-1; z <= b_particle.z+1; z++) {
  	  if (z < 0 || z > this->atomicLattice->getNz()-1) {continue;}
            const int & index = grid_index(x,y,z);
            for (unsigned int i = 0 ; i < particle_grid_size[index] ; i++ ) {
              HemoCellParticle & lParticle = particles[particle_grid[index][i]];
              const hemo::Array<T,3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation()); 
              const T distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]); 
              T tresca = eigenValueFromCell(this->atomicLattice->get(x,y,z));
   	    if ((distance <= (*cellFields)[lParticle.sv.celltype]->mechanics->cfg["MaterialModel"]["distanceThreshold"].read<T>())  
                      && (abs(tresca/1e-7) > (*cellFields)[lParticle.sv.celltype]->mechanics->cfg["MaterialModel"]["shearThreshold"].read<T>()) ) { 
  	      lParticle.sv.solidify = true;
              } 
            }
          }
        }
      }
    }
  #else
    hlog << "(HemoCellParticleField) SolidifyCells called but SOLIDIFY_MECHANICS not enabled" << endl;
    exit(1);
  #endif
  }


  HemoCellParticleDataTransfer& HemoCellParticleField::getDataTransfer() {
      return particleDataTransfer;
  }
  HemoCellParticleDataTransfer const& HemoCellParticleField::getDataTransfer() const {
      return particleDataTransfer;
  }

  std::string HemoCellParticleField::getBlockName() {
      return std::string("HemoParticleField3D");
  }

  HemoCellFields* HemoCellParticleField::cellFields=0;
  }
