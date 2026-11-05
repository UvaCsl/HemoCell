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

#include <cmath>
#include "hemocell.h"
#include "octree.h"
#include "mollerTrumbore.h"
#include "bindingField.h"
#include "interiorViscosity.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include <Eigen3/Eigenvalues>

using namespace std;
using namespace plb;
#pragma GCC diagnostic pop

namespace hemo {
/* *************** class HemoParticleField3D ********************** */

HemoCellParticleField::HemoCellParticleField(plint nx, plint ny, plint nz)
    : AtomicBlock3D(nx, ny, nz, 0), particleDataTransfer(*(new HemoCellParticleDataTransfer())) {
  boundingBox = Box3D(0, this->getNx() - 1, 0, this->getNy() - 1, 0, this->getNz() - 1);
  dataTransfer = &particleDataTransfer;
  particleDataTransfer.setBlock(*this);
  AddOutputMap();
}

HemoCellParticleField::HemoCellParticleField(HemoCellParticleField const & rhs)
    : AtomicBlock3D(rhs), particleDataTransfer(*(new HemoCellParticleDataTransfer())) {
  boundingBox = Box3D(0, this->getNx() - 1, 0, this->getNy() - 1, 0, this->getNz() - 1);
  delete dataTransfer; // clone made by AtomicBlock3D(rhs); replaced by our own (owned and deleted by AtomicBlock3D)
  dataTransfer = &particleDataTransfer;
  particleDataTransfer.setBlock(*this);
  for (const HemoCellParticle & particle : rhs.particles) { addParticle(particle.sv); }
  ppc_up_to_date = false;
  lpc_up_to_date = false;
  ppt_up_to_date = false;
  pg_up_to_date = false;
  AddOutputMap();
}

HemoCellParticleField::~HemoCellParticleField() {
  // Sanitize for MultiBlockLattice destructor (releasememory). It can't handle releasing non-background dynamics that are not singular
  if (global.enableInteriorViscosity) {
    for (const Dot3D & internalPoint : internalPoints) {
      atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z);
      atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z)
          .attributeDynamics(&atomicLattice->getBackgroundDynamics());
    }
    InteriorViscosityHelper::get(*cellFields).empty(*this);
  }
}

HemoCellParticleField & HemoCellParticleField::operator=(HemoCellParticleField const & rhs) {
  if (this == &rhs) { return *this; }
  particles.clear();
  lpc_up_to_date = false;
  ppt_up_to_date = false;
  ppc_up_to_date = false;
  preinlet_ppc_up_to_date = false;
  pg_up_to_date = false;
  for (const HemoCellParticle & particle : rhs.particles) { addParticle(particle.sv); }
  return *this;
}

HemoCellParticleField * HemoCellParticleField::clone() const {
  return new HemoCellParticleField(*this);
}

const vector<vector<unsigned int>> & HemoCellParticleField::get_particles_per_type() {
  if (!ppt_up_to_date) { update_ppt(); }
  return _particles_per_type;
}
const map<plint, vector<int>> & HemoCellParticleField::get_particles_per_cell() {
  if (!ppc_up_to_date) { update_ppc(); }
  return _particles_per_cell;
}

const map<plint, bool> & HemoCellParticleField::get_lpc() {
  if (!lpc_up_to_date) { update_lpc(); }
  return _lpc;
}
void HemoCellParticleField::update_lpc() {
  _lpc.clear();
  for (const HemoCellParticle & particle : particles) {
    if (isContainedABS(particle.sv.position, localDomain)) { _lpc[particle.sv.cellId] = true; }
  }
  lpc_up_to_date = true;
}
void HemoCellParticleField::update_ppt() {
  _particles_per_type.clear();
  _particles_per_type.resize(cellFields->size());

  for (unsigned int i = 0; i < particles.size(); i++) { _particles_per_type[particles[i].sv.celltype].push_back(i); }
  ppt_up_to_date = true;
}
void HemoCellParticleField::update_ppc() {
  _particles_per_cell.clear();

  for (unsigned int i = 0; i < particles.size(); i++) { insert_ppc(&particles[i], i); }
  ppc_up_to_date = true;
}

void HemoCellParticleField::update_pg() {
  if (!this->atomicLattice) { return; }
  particleGrid.build(particles, this->atomicLattice->getLocation(), this->atomicLattice->getNx(),
                     this->atomicLattice->getNy(), this->atomicLattice->getNz());
  pg_up_to_date = true;
}

void HemoCellParticleField::addParticle(HemoCellParticle * particle) {
  addParticle(particle->sv);
}
void HemoCellParticleField::addParticle(const HemoCellParticle::serializeValues_t & sv) {
  HemoCellParticle *local_sparticle, *particle;
  const hemo::Array<T, 3> & pos = sv.position;
  const map<plint, vector<int>> & particles_per_cell = get_particles_per_cell();

  if (this->isContainedABS(pos, this->getBoundingBox())) {
    //check if we have particle already, if so, we must overwrite but not
    //forget to delete the old entry
    if ((!(particles_per_cell.find(sv.cellId) == particles_per_cell.end()))) {
      if (particles_per_cell.at(sv.cellId)[sv.vertexId] != -1) {
        local_sparticle = &particles[particles_per_cell.at(sv.cellId)[sv.vertexId]];

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
      ppt_up_to_date = false;
      if (this->isContainedABS(pos, localDomain)) { _lpc[particle->sv.cellId] = true; }
      if (ppc_up_to_date) { //Otherwise its rebuild anyway
        insert_ppc(particle, particles.size() - 1);
      }

      // The particle grid is rebuilt lazily (update_pg) before its next use;
      // every step invalidates it in advanceParticles() anyway, so incremental
      // insertion would not save work.
      pg_up_to_date = false;
    }
  }
}

void HemoCellParticleField::addParticlePreinlet(const HemoCellParticle::serializeValues_t & sv) {
  HemoCellParticle * particle;
  const hemo::Array<T, 3> & pos = sv.position;
  const map<plint, vector<int>> & particles_per_cell = get_particles_per_cell();

  if (this->isContainedABS(pos, this->getBoundingBox())) {
    //check if we have particle already, if so, we must overwrite but not
    //forget to delete the old entry
    if ((!(particles_per_cell.find(sv.cellId) == particles_per_cell.end()))) {
      if (particles_per_cell.at(sv.cellId)[sv.vertexId] != -1) {
        return;
      } else {
        goto outer_else;
      }
    } else {
    outer_else:
      //new entry
      particles.emplace_back(sv);
      particle = &particles.back();

      //invalidate ppt
      ppt_up_to_date = false;
      if (this->isContainedABS(pos, localDomain)) { _lpc[particle->sv.cellId] = true; }
      if (ppc_up_to_date) { //Otherwise its rebuild anyway
        insert_ppc(particle, particles.size() - 1);
      }

      pg_up_to_date = false; // rebuilt lazily, see addParticle()
    }
  }
}

void inline HemoCellParticleField::insert_ppc(HemoCellParticle * sparticle, unsigned int index) {
  if (_particles_per_cell.find(sparticle->sv.cellId) == _particles_per_cell.end()) {
    _particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex, -1);
  }
  _particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;
}
void inline HemoCellParticleField::insert_preinlet_ppc(HemoCellParticle * sparticle, unsigned int index) {
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
  for (unsigned int i = 0; i < particles.size(); i++) {
    if (particles[i].getTag() == tag) {
      particles[i] = std::move(particles.back());
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

void HemoCellParticleField::removeParticles(plb::Box3D domain, plint tag) {
  //Almost the same, but we save a lot of branching by making a seperate function
  Box3D finalDomain;

  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0; i < particles.size(); i++) {
    if (particles[i].getTag() == tag && this->isContainedABS(particles[i].sv.position, finalDomain)) {
      particles[i] = std::move(particles.back());
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

void HemoCellParticleField::removeParticles(plb::Box3D domain) {
  //Almost the same, but we save a lot of branching by making a seperate function

  Box3D finalDomain;

  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0; i < particles.size(); i++) {
    if (this->isContainedABS(particles[i].sv.position, finalDomain)) {
      particles[i] = std::move(particles.back());
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
  for (unsigned int i = 0; i < particles.size(); i++) {
    if (!this->isContainedABS(particles[i].sv.position, finalDomain)) {
      particles[i] = std::move(particles.back());
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

void HemoCellParticleField::findParticles(plb::Box3D domain, std::vector<HemoCellParticle *> & found) {
  found.clear();
  PLB_ASSERT(contained(domain, this->getBoundingBox()));
  for (HemoCellParticle & particle : particles) {
    if (this->isContainedABS(particle.sv.position, domain)) { found.push_back(&particle); }
  }
}

void HemoCellParticleField::findParticles(plb::Box3D domain, std::vector<const HemoCellParticle *> & found) const {
  found.clear();
  PLB_ASSERT(contained(domain, this->getBoundingBox()));
  for (const HemoCellParticle & particle : particles) {
    if (this->isContainedABS(particle.sv.position, domain)) { found.push_back(&particle); }
  }
}
void HemoCellParticleField::findParticles(plb::Box3D domain, std::vector<HemoCellParticle *> & found, pluint type) {
  found.clear();
  PLB_ASSERT(contained(domain, this->getBoundingBox()));
  //hemo::Array<T,3> pos;
  const vector<vector<unsigned int>> & particles_per_type = get_particles_per_type();
  if (!(particles_per_type.size() > type)) {
    return;
  } else {
    for (const unsigned int i : particles_per_type[type]) {
      if (this->isContainedABS(particles[i].sv.position, domain)) { found.push_back(&(particles[i])); }
    }
  }
}

inline plint HemoCellParticleField::nearestCell(T const pos) const {
  return int(pos + 0.5);
}

inline void HemoCellParticleField::computeGridPosition(hemo::Array<T, 3> const & position, plint * iX, plint * iY,
                                                       plint * iZ) const {
  Dot3D const & location = this->getLocation();
  *iX = nearestCell(position[0]) - location.x;
  *iY = nearestCell(position[1]) - location.y;
  *iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::computeGridPosition(hemo::Array<T, 3> const & position, plint & iX, plint & iY,
                                                plint & iZ) const {
  Dot3D const & location = this->getLocation();
  iX = nearestCell(position[0]) - location.x;
  iY = nearestCell(position[1]) - location.y;
  iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::issueWarning(HemoCellParticle & p) {
  cout
      << "(HemoCell) (Delete Cells) WARNING! Particle deleted from local domain. This means the whole cell will be deleted!"
      << endl;
  cout << "\t Particle ID:" << p.sv.cellId << endl;
  cout << "\t Position: " << p.sv.position[0] << ", " << p.sv.position[1] << ", " << p.sv.position[2]
       << "; vel.: " << p.sv.v[0] << ", " << p.sv.v[1] << ", " << p.sv.v[2] << "; force: " << p.sv.force[0] << ", "
       << p.sv.force[1] << ", " << p.sv.force[2] << endl;
}

int HemoCellParticleField::deleteIncompleteCells(pluint ctype, bool verbose) {
  return deleteIncompleteCellsImpl(true, ctype, nullptr, verbose);
}

int HemoCellParticleField::deleteIncompleteCells(const bool verbose) {
  return deleteIncompleteCellsImpl(false, 0, nullptr, verbose);
}

int HemoCellParticleField::deleteBrokenCells(const plint envelope, const bool verbose) {
  const Box3D near = localDomain.enlarge(envelope);
  return deleteIncompleteCellsImpl(false, 0, &near, verbose);
}

void HemoCellParticleField::snapshotParticles() {
  // The originals are moved aside unchanged (the vector buffer is swapped, so
  // their addresses, and the force pointers aliasing them, stay valid).
  particlesSnapshot.swap(particles);
  particles = particlesSnapshot;
  lpc_up_to_date = false;
  ppt_up_to_date = false;
  ppc_up_to_date = false;
  preinlet_ppc_up_to_date = false;
  pg_up_to_date = false;
}

void HemoCellParticleField::restoreParticles() {
  particles.swap(particlesSnapshot);
  vector<HemoCellParticle>().swap(particlesSnapshot); // release the copy
  lpc_up_to_date = false;
  ppt_up_to_date = false;
  ppc_up_to_date = false;
  preinlet_ppc_up_to_date = false;
  pg_up_to_date = false;
}

int HemoCellParticleField::deleteIncompleteCellsImpl(const bool onlyType, const pluint ctype, const Box3D * near,
                                                     const bool verbose) {
  int deleted = 0;
  const map<plint, vector<int>> & particles_per_cell = get_particles_per_cell();
  //Warning, TODO, high complexity, should be rewritten
  //For now abuse tagging and the remove function
  for (const auto & lpc_it : particles_per_cell) {
    const vector<int> & cell = lpc_it.second;
    bool broken = false;
    int present = -1; // index of any particle of this cell that is present
    for (pluint i = 0; i < cell.size(); i++) {
      if (cell[i] == -1) {
        broken = true;
      } else if (present < 0) {
        present = cell[i];
      }
      if (broken && present >= 0) { break; }
    }
    if (!broken || present < 0) { continue; }
    if (onlyType && particles[present].sv.celltype != ctype) { continue; }
    if (near) {
      bool isNear = false;
      for (const int index : cell) {
        if (index >= 0 && isContainedABS(particles[index].sv.position, *near)) {
          isNear = true;
          break;
        }
      }
      if (!isNear) { continue; }
    }

    bool warningIssued = false;
    for (pluint i = 0; i < cell.size(); i++) {
      if (cell[i] == -1) { continue; }

      //issue warning
      if (verbose) {
        if (!warningIssued) {
          if (isContainedABS(particles[cell[i]].sv.position, localDomain)) {
            issueWarning(particles[cell[i]]);
            warningIssued = true;
          }
        }
      }

      //actually add to tobedeleted list
      particles[cell[i]].setTag(1);
      deleted++;
    }
  }

  //We have our list, now abuse the removeall function
  if (deleted > 0) { removeParticles(1); }

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
  plb::Box3D const box = atomicLattice->getBoundingBox();
  plb::Dot3D const location = atomicLattice->getLocation();
  bool anyOnBoundary = false;
  for (HemoCellParticle & particle : particles) {
    particle.advance();
    //By lack of better place, check if it is on a boundary, if so, delete it
    plint x = (particle.sv.position[0] - location.x) + 0.5;
    plint y = (particle.sv.position[1] - location.y) + 0.5;
    plint z = (particle.sv.position[2] - location.z) + 0.5;

    if ((x >= box.x0) && (x <= box.x1) && (y >= box.y0) && (y <= box.y1) && (z >= box.z0) && (z <= box.z1)) {
      if (atomicLattice->get(x, y, z).getDynamics().isBoundary()) {
        particle.tag = 1;
        anyOnBoundary = true;
      }
    }
  }
  // Tag 1 is only set above at this point (other users remove their tagged
  // particles immediately), so the removal pass can be skipped if none hit a wall.
  if (anyOnBoundary) { removeParticles(1); }

  lpc_up_to_date = false;
  pg_up_to_date = false;
}

void HemoCellParticleField::separateForceVectors() {
  // Save the total force as it will be applied in the next iteration. It is
  // deliberately *not* recomputed here: recomputing the constitutive model at
  // output time would change sv.force for cell types with a material timescale
  // separation and thereby make the trajectory depend on the output interval.
  for (HemoCellParticle & sparticle : particles) {
    //Save Total Force
    sparticle.force_total = sparticle.sv.force + sparticle.sv.force_repulsion;

    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want

    sparticle.force_volume = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
    allocated_for_output.push_back(sparticle.force_volume);
    sparticle.force_link = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
    allocated_for_output.push_back(sparticle.force_link);
    sparticle.force_area = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
    allocated_for_output.push_back(sparticle.force_area);
    sparticle.force_bending = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
    allocated_for_output.push_back(sparticle.force_bending);
    sparticle.force_visc = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
    allocated_for_output.push_back(sparticle.force_visc);
    sparticle.force_inner_link = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
    allocated_for_output.push_back(sparticle.force_inner_link);
  }
}

void HemoCellParticleField::updateResidenceTime(unsigned int rtime) {
  for (HemoCellParticle & sparticle : particles) { sparticle.sv.restime += rtime; }
}


void HemoCellParticleField::unifyForceVectors() {
  for (const hemo::Array<T, 3> * mem : allocated_for_output) { delete mem; }
  allocated_for_output.clear();
  for (HemoCellParticle & sparticle : particles) { sparticle.repoint_force_vectors(); }
}

void HemoCellParticleField::applyConstitutiveModel(bool forced) {
  map<plint, vector<HemoCellParticle *>> * ppc_new = new map<plint, vector<HemoCellParticle *>>();
  const map<plint, vector<int>> & particles_per_cell = get_particles_per_cell();
  map<plint, bool> lpc;
  //Fill it here, probably needs optimization, ah well ...
  for (const auto & pair : particles_per_cell) {
    const plint cid = pair.first;
    const vector<int> & cell = pair.second;
    (*ppc_new)[cid].resize(cell.size());
    for (unsigned int i = 0; i < cell.size(); i++) {
      if (cell[i] == -1) {
        (*ppc_new).erase(cid); //not complete, remove entry
        goto no_add_lpc;
      } else {
        (*ppc_new)[cid][i] = &particles[cell[i]];
      }
    }
    lpc[cid] = true;
  no_add_lpc:;
  }

  for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++) {
    if ((*cellFields).hemocell.iter % (*cellFields)[ctype]->timescale == 0 || forced) {
      vector<HemoCellParticle *> found;
      findParticles(getBoundingBox(), found, ctype);
      if (found.size() > 0) {
        //only reset forces when the forces actually point at it. During output
        //they point to separate arrays (separateForceVectors) and sv.force must
        //keep the force that is spread in the next iteration.
        if (found[0]->force_area == &found[0]->sv.force) {
          for (HemoCellParticle * particle : found) {
            particle->sv.force = {0., 0., 0.};
#ifdef INTERIOR_VISCOSITY
            particle->normalDirection = {0., 0., 0.};
#endif
          }
        }
      }
      (*cellFields)[ctype]->mechanics->ParticleMechanics(*ppc_new, lpc, ctype);
    }
  }

  delete ppc_new;
}

void HemoCellParticleField::applyRepulsionForce(bool forced) {
  const T r_const = cellFields->repulsionConstant;
  const T r_cutoff = cellFields->repulsionCutoff;
  if (!pg_up_to_date) { update_pg(); }
  for (HemoCellParticle & particle : particles) { particle.sv.force_repulsion = {0., 0., 0.}; }
  computePairRepulsion(particles, particleGrid, r_const, r_cutoff);
}

void HemoCellParticleField::capRepulsionForce(const T limit) {
  for (HemoCellParticle & particle : particles) {
    hemo::Array<T, 3> & force = particle.sv.force_repulsion;
    const T magnitude = sqrt(force[0] * force[0] + force[1] * force[1] + force[2] * force[2]);
    if (magnitude > limit) { force *= limit / magnitude; }
  }
}

#ifdef INTERIOR_VISCOSITY
void HemoCellParticleField::internalGridPointsMembrane(Box3D domain) {
  // This could be done less complex I guess?
  for (const HemoCellParticle & particle : particles) { // Go over each particle
    if (!(*cellFields)[particle.sv.celltype]->doInteriorViscosity) { continue; }

    for (unsigned int i = 0; i < particle.kernelCoordinates.size(); i++) {
      const hemo::Array<T, 3> latPos =
          particle.kernelCoordinates[i] - (particle.sv.position - atomicLattice->getLocation());
      const hemo::Array<T, 3> & normalP = particle.normalDirection;

      if (computeLength(latPos) > (*cellFields)[particle.sv.celltype]->mechanics->cellConstants.edge_mean_eq) {
        continue;
      }

      T dot1 = hemo::dot(latPos, normalP);

      if (dot1 < 0.) { // Node is inside
        InteriorViscosityHelper::get(*cellFields)
            .add(*this,
                 {particle.kernelCoordinates[i][0], particle.kernelCoordinates[i][1], particle.kernelCoordinates[i][2]},
                 (*cellFields)[particle.sv.celltype]->interiorViscosityTau);
        particle.kernelLocations[i]->attributeDynamics((*cellFields)[particle.sv.celltype]->innerViscosityDynamics);
      } else { // Node is outside
        InteriorViscosityHelper::get(*cellFields)
            .remove(*this, {particle.kernelCoordinates[i][0], particle.kernelCoordinates[i][1],
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
  for (const Dot3D & internalPoint : internalPoints) {
    atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z);
    atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z)
        .attributeDynamics(&atomicLattice->getBackgroundDynamics());
  }
  InteriorViscosityHelper::get(*cellFields).empty(*this);

  for (const auto & pair : get_lpc()) { // Go over each cell?
    const plint cid = pair.first;
    const vector<int> & cell = get_particles_per_cell().at(cid);
    const pluint ctype = particles[cell[0]].sv.celltype;

    // Plt and Wbc now have normal tau internal, so we don't have
    // to raycast these particles
    if (!(*cellFields)[ctype]->doInteriorViscosity) { continue; }

    hemo::OctreeStructCell octCell(3, 1, 30, (*cellFields)[ctype]->mechanics->cellConstants.triangle_list, particles,
                                   cell);

    std::set<Array<plint, 3>> innerNodes;
    octCell.findInnerNodes(atomicLattice, particles, cell, innerNodes);
    for (const Array<plint, 3> & node : innerNodes) {
      InteriorViscosityHelper::get(*cellFields)
          .add(*this, {node[0], node[1], node[2]}, (*cellFields)[ctype]->interiorViscosityTau);
      atomicLattice->get(node[0], node[1], node[2]).attributeDynamics((*cellFields)[ctype]->innerViscosityDynamics);
    }
  }
}
#else
void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain) {
  fatal("findInternalParticleGridPoints() needs a library built with INTERIOR_VISCOSITY.");
}
void HemoCellParticleField::internalGridPointsMembrane(Box3D domain) {
  fatal("internalGridPointsMembrane() needs a library built with INTERIOR_VISCOSITY.");
}
#endif

void HemoCellParticleField::interpolateFluidVelocity(Box3D domain) {
  // Preallocating
  hemo::Array<T, 3> velocity;
  plb::Array<T, 3> velocity_comp;

  for (HemoCellParticle & particle : particles) {
    // Only particles owned by this block: interpolation is always followed by
    // syncEnvelopes(), which drops the envelope copies and receives them again
    // with the velocity their owner interpolated.
    if (!isContainedABS(particle.sv.position, localDomain)) { continue; }

    // We have the kernels (computed in spreadParticleForce()), now calculate
    // the velocity of the particle.
    velocity = {0.0, 0.0, 0.0};
    for (pluint j = 0; j < particle.kernelLocations.size(); j++) {
      // Direct access
      particle.kernelLocations[j]->computeVelocity(velocity_comp);
      velocity += (velocity_comp * particle.kernelWeights[j]);
    }
    particle.sv.v = velocity;
  }
}

void HemoCellParticleField::spreadParticleForce(Box3D domain) {
#ifndef INTERIOR_VISCOSITY
  // Envelope copies more than two nodes outside the bulk only reach envelope
  // nodes from the second layer on. Those forces are never used: Palabos
  // collides the bulk plus the active fluid envelope, which HemoCell keeps at
  // one layer (initializeCellfield()), envelope populations are overwritten by
  // communication, and interpolation and output read at most the first
  // envelope layer. So their kernels and spreading are skipped. (Interior
  // viscosity uses the kernels of all particles, so it keeps computing them.)
  const plb::Box3D spreadDomain = localDomain.enlarge(2);
#endif
  // A non-finite force (e.g. from a degenerate triangle of a broken cell)
  // would make the whole lattice NaN; such a vertex spreads nothing and is
  // removed, which deletes its cell like any incomplete cell.
  plint nonFinite = 0;
  for (HemoCellParticle & particle : particles) {
    const T force_mag = norm(particle.sv.force);
#ifdef FORCE_LIMIT
    // Capping force to ensure stability -> NOTE: this can introduce an error if forces are large!
    // !(a <= b) is also true for NaN, so the finiteness check costs nothing
    // for forces within the limit.
    if (!(force_mag <= param::f_limit)) {
      if (std::isfinite(force_mag)) {
        particle.sv.force *= param::f_limit / force_mag;
      } else {
        dropNonFiniteParticle(particle, nonFinite);
        continue;
      }
    }
#else
    if (!std::isfinite(force_mag)) {
      dropNonFiniteParticle(particle, nonFinite);
      continue;
    }
#endif

#ifndef INTERIOR_VISCOSITY
    if (!isContainedABS(particle.sv.position, spreadDomain)) {
      particle.kernelLocations.clear();
      particle.kernelWeights.clear();
      continue;
    }
#endif

    //Trick to allow for different kernels for different particle types.
    (*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice, particle);

    // Directly change the force on a node, quick-and-dirty solution.
    const T force[3] = {particle.sv.force_repulsion[0] + particle.sv.force[0],
                        particle.sv.force_repulsion[1] + particle.sv.force[1],
                        particle.sv.force_repulsion[2] + particle.sv.force[2]};
    for (pluint j = 0; j < particle.kernelLocations.size(); j++) {
      // Direct access
      T * const external = particle.kernelLocations[j]->external.data;
      const T weight = particle.kernelWeights[j];
      external[0] += force[0] * weight;
      external[1] += force[1] * weight;
      external[2] += force[2] * weight;
    }
  }
  if (nonFinite > 0) { removeParticles(1); }
}

void HemoCellParticleField::dropNonFiniteParticle(HemoCellParticle & particle, plint & count) {
  if (count == 0) {
    warning("Non-finite force on vertex ", particle.sv.vertexId, " of cell ",
            cellFields->base_cell_id(particle.sv.cellId), " (type ", (*cellFields)[particle.sv.celltype]->name,
            ") at (", particle.sv.position[0], ", ", particle.sv.position[1], ", ", particle.sv.position[2],
            "); its cell is removed.");
  }
  count++;
  particle.sv.force = {0., 0., 0.};
  particle.kernelLocations.clear();
  particle.kernelWeights.clear();
  particle.setTag(1);
}

void HemoCellParticleField::populateBoundaryParticles() {
  for (int x = 0; x < this->atomicLattice->getNx() - 1; x++) {
    for (int y = 0; y < this->atomicLattice->getNy() - 1; y++) {
      for (int z = 0; z < this->atomicLattice->getNz() - 1; z++) {
        if (this->atomicLattice->get(x, y, z).getDynamics().isBoundary()) {
          for (int xx = x - 1; xx <= x + 1; xx++) {
            if (xx < 0 || xx > this->atomicLattice->getNx() - 1) { continue; }
            for (int yy = y - 1; yy <= y + 1; yy++) {
              if (yy < 0 || yy > this->atomicLattice->getNy() - 1) { continue; }
              for (int zz = z - 1; zz <= z + 1; zz++) {
                if (zz < 0 || zz > this->atomicLattice->getNz() - 1) { continue; }
                if (!this->atomicLattice->get(xx, yy, zz).getDynamics().isBoundary()) {
                  boundaryParticles.push_back({x, y, z});
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

void HemoCellParticleField::applyBoundaryRepulsionForce(bool resetFirst) {
  if (!pg_up_to_date) { update_pg(); }
  if (resetFirst) {
    for (HemoCellParticle & particle : particles) { particle.sv.force_repulsion = {0., 0., 0.}; }
  }
  const T & br_cutoff = cellFields->boundaryRepulsionCutoff;
  const T & br_const = cellFields->boundaryRepulsionConstant;
  for (Dot3D & b_particle : boundaryParticles) {
    for (int x = b_particle.x - 1; x <= b_particle.x + 1; x++) {
      if (x < 0 || x > this->atomicLattice->getNx() - 1) { continue; }
      for (int y = b_particle.y - 1; y <= b_particle.y + 1; y++) {
        if (y < 0 || y > this->atomicLattice->getNy() - 1) { continue; }
        for (int z = b_particle.z - 1; z <= b_particle.z + 1; z++) {
          if (z < 0 || z > this->atomicLattice->getNz() - 1) { continue; }
          const unsigned int index = particleGrid.index(x, y, z);
          const unsigned int * const index_end = particleGrid.end(index);
          for (const unsigned int * i = particleGrid.begin(index); i != index_end; ++i) {
            HemoCellParticle & lParticle = particles[*i];
            const hemo::Array<T, 3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation());
            const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
            if (distance < br_cutoff && distance > 0) { // see repulseParticlePair()
              const hemo::Array<T, 3> rfm = br_const * (1 / (distance / br_cutoff)) * (dv / distance);
              lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;
            }
          }
        }
      }
    }
  }
}

void HemoCellParticleField::populateBindingSites(plb::Box3D & domain_) {
  // Translate domain to make sense in the lattice domain.
  Dot3D shift = this->getLocation() - this->atomicLattice->getLocation();
  Box3D domain = domain_.shift(shift.x, shift.y, shift.z);

  for (int x = domain.x0; x <= domain.x1; x++) {
    for (int y = domain.y0; y <= domain.y1; y++) {
      for (int z = domain.z0; z <= domain.z1; z++) {
        if (this->atomicLattice->get(x, y, z).getDynamics().isBoundary()) {
          for (int xx = x - 1; xx <= x + 1; xx++) {
            if (xx < 0 || xx > this->atomicLattice->getNx() - 1) { continue; }
            for (int yy = y - 1; yy <= y + 1; yy++) {
              if (yy < 0 || yy > this->atomicLattice->getNy() - 1) { continue; }
              for (int zz = z - 1; zz <= z + 1; zz++) {
                if (zz < 0 || zz > this->atomicLattice->getNz() - 1) { continue; }
                if (!this->atomicLattice->get(xx, yy, zz).getDynamics().isBoundary()) {
                  bindingFieldHelper::get(*cellFields).add(*this, {x, y, z});
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


T HemoCellParticleField::eigenValueFromCell(plb::Cell<T, DESCRIPTOR> & cell) {
  plb::Array<T, SymmetricTensor<T, DESCRIPTOR>::n> element;
  cell.computePiNeq(element);
  T omega = cell.getDynamics().getOmega();
  T rhoBar = cell.getDynamics().computeRhoBar(cell);
  T prefactor = -omega * DESCRIPTOR<T>::invCs2 * DESCRIPTOR<T>::invRho(rhoBar) / (T)2;
  for (int iTensor = 0; iTensor < SymmetricTensor<T, DESCRIPTOR>::n; ++iTensor) { element[iTensor] *= prefactor; }

  Array<Array<T, 3>, 3> S; // Strain-rate tensor (symmetric).
  S[0][0] = element[0];    //s[xx];
  S[0][1] = element[1];    //s[xy];
  S[0][2] = element[2];    //s[xz];

  S[1][0] = element[1]; //s[yx];
  S[1][1] = element[3]; //s[yy];
  S[1][2] = element[4]; //s[yz];

  S[2][0] = element[2]; //s[zx];
  S[2][1] = element[4]; //s[zy];
  S[2][2] = element[5]; //s[zz];

  Eigen::Matrix<T, 3, 3> A;
  for (plint i = 0; i < 3; i++) {
    for (plint j = 0; j < 3; j++) { A(i, j) = S[i][j]; }
  }

  bool computeEigenvectors = false;
  Eigen::EigenSolver<Eigen::Matrix<T, 3, 3>> es(A, computeEigenvectors);
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
  T tresca = (lambda[2] - lambda[0]) / 2;
  return tresca;
}

void HemoCellParticleField::prepareSolidification() {
#ifdef SOLIDIFY_MECHANICS
  for (HemoCellField * type : cellFields->cellFields) {
    ppc_up_to_date = false;
    if (type->doSolidifyMechanics) {
      // Invoke cell-type specific solidification mechanics implementation.
      type->mechanics->solidifyMechanics(get_particles_per_cell(), particles, this->atomicLattice, this->CEPAClattice,
                                         type->ctype, *this);
    }
  }
#else
  fatal("prepareSolidification() needs a library built with SOLIDIFY_MECHANICS.");
#endif
}

void HemoCellParticleField::solidifyCells() {
#ifdef SOLIDIFY_MECHANICS

  // Remove any to be removed particles (tagged with `tag == 1`).
  removeParticles(1);
  if (!pg_up_to_date) { update_pg(); }

  // Detect particles to be solidified by looping over a block of 3x3 LBM cells
  // around each binding site. When a cell statisfies both:
  // - close enough in space to a binding site,
  // - shows a minimum tresca stress,
  // the particle is labelled to be solified.
  // The thresholds of a cell type are read from its material file once per
  // call, on first use (only the types found near a binding site need them).
  std::vector<T> distanceThreshold(cellFields->size(), -1), shearThreshold(cellFields->size(), -1);
  std::vector<bool> thresholdsRead(cellFields->size(), false);
  for (const Dot3D & b_particle : bindingSites) {
    for (int x = b_particle.x - 1; x <= b_particle.x + 1; x++) {
      if (x < 0 || x > this->atomicLattice->getNx() - 1) { continue; }

      for (int y = b_particle.y - 1; y <= b_particle.y + 1; y++) {
        if (y < 0 || y > this->atomicLattice->getNy() - 1) { continue; }

        for (int z = b_particle.z - 1; z <= b_particle.z + 1; z++) {
          if (z < 0 || z > this->atomicLattice->getNz() - 1) { continue; }

          const unsigned int index = particleGrid.index(x, y, z);
          const unsigned int * const index_end = particleGrid.end(index);

          for (const unsigned int * i = particleGrid.begin(index); i != index_end; ++i) {
            HemoCellParticle & lParticle = particles[*i];
            const hemo::Array<T, 3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation());
            const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
            T tresca = eigenValueFromCell(this->atomicLattice->get(x, y, z));

            const pluint ctype = lParticle.sv.celltype;
            if (!thresholdsRead[ctype]) {
              Config & material = (*cellFields)[ctype]->mechanics->cfg;
              distanceThreshold[ctype] = material["MaterialModel"]["distanceThreshold"].read<T>();
              shearThreshold[ctype] = material["MaterialModel"]["shearThreshold"].read<T>();
              thresholdsRead[ctype] = true;
            }
            if ((distance <= distanceThreshold[ctype]) && (abs(tresca / 1e-7) > shearThreshold[ctype])) {
              lParticle.sv.solidify = true;
            }
          }
        }
      }
    }
  }
#else
  fatal("solidifyCells() needs a library built with SOLIDIFY_MECHANICS.");
#endif
}


HemoCellParticleDataTransfer & HemoCellParticleField::getDataTransfer() {
  return particleDataTransfer;
}
HemoCellParticleDataTransfer const & HemoCellParticleField::getDataTransfer() const {
  return particleDataTransfer;
}

std::string HemoCellParticleField::getBlockName() {
  return std::string("HemoParticleField3D");
}

HemoCellFields * HemoCellParticleField::cellFields = 0;
}
