#ifndef HEMOCELL_PARTICLE_FIELD_CPP
#define HEMOCELL_PARTICLE_FIELD_CPP

#include "hemoCellParticleField.h"

/* *************** class HemoParticleDataTransfer3D ************************ */

HemoParticleDataTransfer3D::HemoParticleDataTransfer3D () {};

plint HemoParticleDataTransfer3D::staticCellSize() const {
    return 0;  // Particle containers have only dynamic data.
}

void HemoParticleDataTransfer3D::send (
        Box3D domain, std::vector<char>& buffer, modif::ModifT kind ) const
{
    buffer.clear();
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        std::vector<SurfaceParticle3D*> foundParticles;
        particleField->findParticles(domain, foundParticles);
        for (pluint iParticle=0; iParticle<foundParticles.size(); ++iParticle) {
        // The serialize function automatically reallocates memory for buffer.
        serialize(*foundParticles[iParticle], buffer);
        }
    }
}

void HemoParticleDataTransfer3D::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind )
{
    PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));

    // Clear the existing data before introducing the new data.
    particleField->removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            SurfaceParticle3D* newParticle = new SurfaceParticle3D();
            newParticle->unserialize(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField->addParticle(domain, newParticle);
        }
    }
}

void HemoParticleDataTransfer3D::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset )
{
        //Particle Locations should always be ABSOULUTE, an offset thus makes no
        //sense
        Array<T,3> realAbsoluteOffset((T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z);
    particleField->removeParticles(domain);
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            SurfaceParticle3D* newParticle = new SurfaceParticle3D();
            newParticle->unserialize(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle->position += realAbsoluteOffset;
            particleField->addParticle(domain, newParticle);
        }
    }
}

void HemoParticleDataTransfer3D::setBlock(AtomicBlock3D& block) {
   particleField = dynamic_cast<HemoParticleField3D*>(&block);
   PLB_ASSERT(particleField);
   constParticleField = particleField;
}

void HemoParticleDataTransfer3D::setConstBlock(AtomicBlock3D const& block) {
          constParticleField = dynamic_cast<HemoParticleField3D const*>(&block);
}

void HemoParticleDataTransfer3D::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind )
{
    Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
    std::vector<char> buffer;
    HemoParticleField3D const& fromParticleField =
        dynamic_cast<HemoParticleField3D const&>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

void HemoParticleDataTransfer3D::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset )
{
    Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
    std::vector<char> buffer;
    HemoParticleField3D const& fromParticleField =
        dynamic_cast<HemoParticleField3D const&>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}


/* *************** class HemoParticleField3D ********************** */
void HemoParticleField3D::AddOutputMap() {
  outputFunctionMap[OUTPUT_POSITION] = &HemoParticleField3D::outputPositions;
  outputFunctionMap[OUTPUT_FORCE] = &HemoParticleField3D::outputForces;
  outputFunctionMap[OUTPUT_FORCE_VOLUME] = &HemoParticleField3D::outputForceVolume;
  outputFunctionMap[OUTPUT_FORCE_AREA] = &HemoParticleField3D::outputForceArea;
  outputFunctionMap[OUTPUT_FORCE_INPLANE] = &HemoParticleField3D::outputForceInPlane;
  outputFunctionMap[OUTPUT_FORCE_BENDING] = &HemoParticleField3D::outputForceBending;
  
}

HemoParticleField3D::HemoParticleField3D(plint nx, plint ny, plint nz)
    : AtomicBlock3D(nx,ny,nz, &this->dataTransfer)
{ 
    dataTransfer.setBlock(*this);
    particle_grid.resize(getNx());
    for (auto & particlesnx : particle_grid) {
        particlesnx.resize(getNy());
        for (auto & particlesny : particlesnx) {
            particlesny.resize(getNz());
        }
    }
    AddOutputMap(); 

}

HemoParticleField3D::~HemoParticleField3D()
{
  AtomicBlock3D::dataTransfer = new HemoParticleDataTransfer3D();
    for (pluint i=0; i<particles.size(); ++i) {
        delete particles[i];
    }
}

HemoParticleField3D::HemoParticleField3D(HemoParticleField3D const& rhs)
    : AtomicBlock3D(rhs)
{
    dataTransfer.setBlock(*this);
    particle_grid.resize(getNx());
    for (auto & particlesnx : particle_grid) {
        particlesnx.resize(getNy());
        for (auto & particlesny : particlesnx) {
            particlesny.resize(getNz());
        }
    }

    for (pluint i=0; i<rhs.particles.size(); ++i) {
        addParticle(this->getBoundingBox(),rhs.particles[i]->clone());
    }
    AddOutputMap();
}

HemoParticleField3D& HemoParticleField3D::operator=(HemoParticleField3D const& rhs){
 HemoParticleField3D *copy = new HemoParticleField3D(rhs);
 delete this;
  return *copy;
}

HemoParticleField3D* HemoParticleField3D::clone() const
{
    return new HemoParticleField3D(*this);
}

void HemoParticleField3D::addParticle(Box3D domain, SurfaceParticle3D* particle) {
    Box3D finalDomain;
    plint x,y,z;
    SurfaceParticle3D * local_sparticle;
    Array<double,3> pos = particle->position;


    while (particles_per_type.size()<=particle->celltype) {
      particles_per_type.push_back(std::vector<SurfaceParticle3D*>());
    }

    if( this->isContainedABS(pos, this->getBoundingBox()) )
    {
      //check if we have particle already, if so, we must overwrite but not
      //forget to delete the old entry
      if ((!(particles_per_cell.find(particle->cellId) == 
           particles_per_cell.end())) && particles_per_cell[particle->cellId][particle->vertexId]) {
          local_sparticle =  particles_per_cell[particle->cellId][particle->vertexId];
          //Remove from bin as we finally add it back anyway
          x = local_sparticle->grid_pos[0];
          y = local_sparticle->grid_pos[1];
          z = local_sparticle->grid_pos[2];
          for (unsigned int i = 0; i < particle_grid[x][y][z].size(); i++) {
            if (particle_grid[x][y][z][i] == local_sparticle) {
                particle_grid[x][y][z][i] = particle_grid[x][y][z].back();
                particle_grid[x][y][z].pop_back();
            }
              
          }

          //If our particle is local, do not replace it, envelopes are less important
          if (isContainedABS(local_sparticle->position, getBoundingBox())) {

          } else {
            //We have the particle already, replace it
            *local_sparticle = *particle;
            local_sparticle->repoint_force_vectors();
          }
          //If the pointers are the same we shouldn't delete either
          if(!(particle==local_sparticle)) {
            delete particle;
            particle = local_sparticle;
          }
      } else {
        //new entry
        particles.push_back(particle);
          particles_per_type[particle->celltype].push_back(particle); //TODO, not accurate for already existing particles
          if(this->isContainedABS(pos, localDomain)) {
            lpc[particle->cellId] = true;
          }
          
        insert_ppc(particle);
      }

      computeGridPosition(pos,&x,&y,&z);          
      particle_grid[x][y][z].push_back(particle);
      particle->grid_pos = {x,y,z};
      particle->setTag(-1);
    }
    else {
        delete particle;
    }
}

void HemoParticleField3D::insert_ppc(SurfaceParticle3D* sparticle) {
  if (particles_per_cell.find(sparticle->cellId) == particles_per_cell.end()) {
    particles_per_cell[sparticle->cellId].resize((*cellFields)[sparticle->celltype]->numVertex);
    for (unsigned int i = 0; i < particles_per_cell[sparticle->cellId].size(); i++) {
      particles_per_cell[sparticle->cellId][i] = NULL;
    }
  }
  particles_per_cell[sparticle->cellId][sparticle->vertexId] = sparticle;

}


void HemoParticleField3D::removeParticles(plint tag) {
//Almost the same, but we save a lot of branching by manking a seperate function
    std::vector<SurfaceParticle3D*> remainingParticles = particles;
    SurfaceParticle3D * sparticle;
    plint x,y,z;
    for (pluint i=0; i < particles_per_type.size(); i++) {
      particles_per_type[i].clear();
    }
    particles_per_cell.clear();
    lpc.clear();
    particles.clear();

    for (pluint i=0; i<remainingParticles.size(); ++i) {
       if (remainingParticles[i]->getTag() == tag) {
         sparticle = remainingParticles[i];
         x = sparticle->grid_pos[0];
         y = sparticle->grid_pos[1];
         z = sparticle->grid_pos[2];
         for (unsigned int i = 0; i < particle_grid[x][y][z].size(); i++) {
           if (particle_grid[x][y][z][i] == remainingParticles[i]) {
              particle_grid[x][y][z][i] = particle_grid[x][y][z].back();
              particle_grid[x][y][z].pop_back();
           }
         }
         delete remainingParticles[i];
       }
       else {
           addParticle(this->getBoundingBox(), remainingParticles[i]);
       }
    }
}

void HemoParticleField3D::removeParticles(Box3D domain, plint tag) {
//Almost the same, but we save a lot of branching by manking a seperate function
    std::vector<SurfaceParticle3D*> remainingParticles = particles;
    SurfaceParticle3D * sparticle;
    plint x,y,z;
    Box3D finalDomain;
    Array<double,3> pos; 
    for (pluint i=0; i < particles_per_type.size(); i++) {
      particles_per_type[i].clear();
    }
    particles_per_cell.clear();
    lpc.clear();
    particles.clear();
    
    for (pluint i=0; i<remainingParticles.size(); ++i) {
       pos = remainingParticles[i]->position;
       intersect(domain, this->getBoundingBox(), finalDomain);
       if (this->isContainedABS(pos,finalDomain) && remainingParticles[i]->getTag() == tag) {
         sparticle = dynamic_cast<SurfaceParticle3D*>(remainingParticles[i]);
         x = sparticle->grid_pos[0];
         y = sparticle->grid_pos[1];
         z = sparticle->grid_pos[2];
         for (unsigned int i = 0; i < particle_grid[x][y][z].size(); i++) {
           if (particle_grid[x][y][z][i] == remainingParticles[i]) {
              particle_grid[x][y][z][i] = particle_grid[x][y][z].back();
              particle_grid[x][y][z].pop_back();
           }
         }
         delete remainingParticles[i];
       }
       else {
           addParticle(this->getBoundingBox(), remainingParticles[i]);
       }
    }
}

void HemoParticleField3D::removeParticles(Box3D domain) {
//Almost the same, but we save a lot of branching by making a seperate function
// Dont allow palabos to delete things
    //return;
    //if (!cellFields) return; //TODO shouldnt be necessary
    //if (!cellFields->hemocellfunction) return;
    std::vector<SurfaceParticle3D*> remainingParticles = particles;
    SurfaceParticle3D * sparticle;
    plint x,y,z;
    Box3D finalDomain;
    Array<double,3> pos; 
    for (pluint i=0; i < particles_per_type.size(); i++) {
      particles_per_type[i].clear();
    }
    particles_per_cell.clear();
    lpc.clear();
    particles.clear();
    
    for (pluint i=0; i<remainingParticles.size(); ++i) {
       pos = remainingParticles[i]->position;
       intersect(domain, this->getBoundingBox(), finalDomain);
       if (this->isContainedABS(pos,finalDomain)) {
         sparticle = dynamic_cast<SurfaceParticle3D*>(remainingParticles[i]);
         x = sparticle->grid_pos[0];
         y = sparticle->grid_pos[1];
         z = sparticle->grid_pos[2];
         for (unsigned int i = 0; i < particle_grid[x][y][z].size(); i++) {
           if (particle_grid[x][y][z][i] == remainingParticles[i]) {
              particle_grid[x][y][z][i] = particle_grid[x][y][z].back();
              particle_grid[x][y][z].pop_back();
           }
         }
         delete remainingParticles[i];
       }
       else {
           addParticle(this->getBoundingBox(), remainingParticles[i]);
       }
    }
}

void HemoParticleField3D::findParticles (
        Box3D domain, std::vector<SurfaceParticle3D*>& found ) 
{
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    Array<double,3> pos; 
    for (pluint i=0; i<particles.size(); ++i) {
        pos = particles[i]->position;
        if (this->isContainedABS(pos,domain)) {
            found.push_back(particles[i]);
        }
    }
}
void HemoParticleField3D::findParticles (
        Box3D domain, std::vector<SurfaceParticle3D*>& found, pluint type) const
{
    
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    Array<double,3> pos; 
    if (!(particles_per_type.size() > type)) {return;} else {
    for (pluint i=0; i<particles_per_type[type].size(); ++i) {
        pos = particles_per_type[type][i]->position;
        if (this->isContainedABS(pos,domain)) {
            found.push_back(particles_per_type[type][i]);
        }
    }
    }
    
}

bool HemoParticleField3D::isContainedABS(Array<double,3> pos, Box3D box) const {
		Dot3D const& location = this->getLocation();
    double x = pos[0]-location.x;
    double y = pos[1]-location.y;
    double z = pos[2]-location.z;
    if (box.z1 < -10000) {
      exit(0);
    } 
    if (location.x < -10000) {
      exit(0);
    } 

    return (x > box.x0-0.5) && (x <= box.x1+0.5) &&
           (y > box.y0-0.5) && (y <= box.y1+0.5) &&
           (z > box.z0-0.5) && (z <= box.z1+0.5);

}

inline plint HemoParticleField3D::nearestCell(double const pos) const {
  return int(pos + 0.5);
}

inline void HemoParticleField3D::computeGridPosition (
            Array<double,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const
{
      Dot3D const& location = this->getLocation();
      *iX = nearestCell(position[0]) - location.x;
      *iY = nearestCell(position[1]) - location.y;
      *iZ = nearestCell(position[2]) - location.z;
}
void HemoParticleField3D::computeGridPosition (
            Array<double,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const
{
      Dot3D const& location = this->getLocation();
      iX = nearestCell(position[0]) - location.x;
      iY = nearestCell(position[1]) - location.y;
      iZ = nearestCell(position[2]) - location.z;
}

int HemoParticleField3D::deleteIncompleteCells(pluint ctype, bool twice) {
  //Function must be called twice since addParticle can remove a particle
  //unintentionally, for now, catch it here; TODO, this can be done better
  int deleted = 0;
  if (!twice) { deleted = deleteIncompleteCells(true); }

  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (!particles_per_cell[cellid][i]) {
        broken = true;
        break;
      }
    }
    if (!broken) {continue;}

    //actually add to tobedeleted list
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (particles_per_cell[cellid][i] == NULL) {continue;}
      if (particles_per_cell[cellid][i]->celltype != ctype) {break;} //certainly a entry, therefore we check here if it is the right type, if not, exit

      particles_per_cell[cellid][i]->setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function
  removeParticles(1);

  return deleted; 
}

int HemoParticleField3D::deleteIncompleteCells(bool twice) {
  //Function must be called twice since addParticle can remove a particle
  //unintentionally, for now, catch it here; TODO, this can be done better
  int deleted = 0;
  if (!twice) {deleted = deleteIncompleteCells(true); }

  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (!particles_per_cell[cellid][i]) {
        broken = true;
        break;
      }
    }
    if (!broken) {continue;}

    //actually add to tobedeleted list
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (particles_per_cell[cellid][i] == NULL) {continue;}
      particles_per_cell[cellid][i]->setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function
  removeParticles(1);

  return deleted; 
}

void HemoParticleField3D::outputPositions(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  deleteIncompleteCells(ctype);
  name = "Position";
  output.clear();
  SurfaceParticle3D * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> pbv;
      pbv.push_back(sparticle->position[0]);
      pbv.push_back(sparticle->position[1]);
      pbv.push_back(sparticle->position[2]);
      output.push_back(pbv); //TODO, memory copy

    }
  }
}

void HemoParticleField3D::outputForceBending(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Bending Force";
  output.clear();
  SurfaceParticle3D * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_bending)[0]);
      tf.push_back((*sparticle->force_bending)[1]);
      tf.push_back((*sparticle->force_bending)[2]);
      output.push_back(tf);
    }
  }
}
void HemoParticleField3D::outputForceArea(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Area Force";
  output.clear();
  SurfaceParticle3D * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_area)[0]);
      tf.push_back((*sparticle->force_area)[1]);
      tf.push_back((*sparticle->force_area)[2]);
      output.push_back(tf);
    }
  }
}
void HemoParticleField3D::outputForceInPlane(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "In Plane Force";
  output.clear();
  SurfaceParticle3D * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_inplane)[0]);
      tf.push_back((*sparticle->force_inplane)[1]);
      tf.push_back((*sparticle->force_inplane)[2]);
      output.push_back(tf);
    }
  }
}
void HemoParticleField3D::outputForceVolume(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Volume Force";
  output.clear();
  SurfaceParticle3D * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_volume)[0]);
      tf.push_back((*sparticle->force_volume)[1]);
      tf.push_back((*sparticle->force_volume)[2]);
      output.push_back(tf);
    }
  }
}

void HemoParticleField3D::outputForces(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Total Force";
  output.clear();
  SurfaceParticle3D * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back(sparticle->force_total[0]);
      tf.push_back(sparticle->force_total[1]);
      tf.push_back(sparticle->force_total[2]);
      output.push_back(tf);
    }
  }
}

void HemoParticleField3D::outputTriangles(Box3D domain, vector<vector<plint>>& output, vector<vector<double>> & positions, pluint ctype, std::string & name) {
  name = "Triangles";
  output.clear();
  int counter = 0;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;

    for (pluint i = 0; i < (*cellFields)[ctype]->triangle_list.size(); i++) {
      vector<plint> triangle = {(*cellFields)[ctype]->triangle_list[i][0] + counter,
                          (*cellFields)[ctype]->triangle_list[i][1] + counter,
                          (*cellFields)[ctype]->triangle_list[i][2] + counter};

      //Do not add triangles over periodic boundaries
      bool toolarge = false;
      /*for (pluint x = 0; x < 3; x++) {
        for (pluint y = x +1; y < 3; y++) {
          if ((abs(positions[triangle[x]][0] - positions[triangle[y]][0]) > 2.0) ||
              (abs(positions[triangle[x]][1] - positions[triangle[y]][1]) > 2.0) ||
              (abs(positions[triangle[x]][2] - positions[triangle[y]][2]) > 2.0))
          {
            toolarge = true;
            break;
          }
        }
      }*/
      if (!toolarge) {
        output.push_back(triangle);
      }
    }
    counter += (*cellFields)[ctype]->numVertex;
  }
   
}

void HemoParticleField3D::setlocalDomain(Box3D & localDomain_) {
  localDomain = localDomain_;
  localDomain.x0 -= this->getLocation().x;
  localDomain.x1 -= this->getLocation().x;
  localDomain.y0 -= this->getLocation().y;
  localDomain.y1 -= this->getLocation().y;
  localDomain.z0 -= this->getLocation().z;
  localDomain.z1 -= this->getLocation().z;
}

void HemoParticleField3D::passthroughpass(int type, Box3D domain, vector<vector<double>>& output, pluint ctype, std::string & name) {
  //Too Much c++ will give you cancer like this function
  void (HemoParticleField3D::*cancerpointer)(Box3D,vector<vector<double>>&, pluint, std::string&) = outputFunctionMap[type];
  (this->*cancerpointer)(domain,output,ctype,name);
}

void HemoParticleField3D::advanceParticles() {
  for(SurfaceParticle3D* particle:particles){
    particle->advance();
  }
  removeParticles(1);
}

void HemoParticleField3D::separateForceVectors() {
  deleteIncompleteCells();
  //Also save the total force, therfore recalculate in advance
  applyConstitutiveModel();
  

  for (SurfaceParticle3D* sparticle : particles) {
    //Save Total Force
    sparticle->force_total = sparticle->force;

    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want
    ////TODO this can leak if particle is deleted between seperate and unify,
    //rewrite to reference
    sparticle->force_volume = new Array<double,3>(0.0,0.0,0.0);
    sparticle->force_inplane = new Array<double,3>(0.0,0.0,0.0);
    sparticle->force_area = new Array<double,3>(0.0,0.0,0.0);
    sparticle->force_bending = new Array<double,3>(0.0,0.0,0.0);

  }



}

void HemoParticleField3D::unifyForceVectors() {
  for (SurfaceParticle3D* sparticle : particles) {
    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want
    delete sparticle->force_volume;
    delete sparticle->force_inplane;
    delete sparticle->force_area;
    delete sparticle->force_bending;
    sparticle->repoint_force_vectors();
  }



}

void HemoParticleField3D::applyConstitutiveModel() {
  deleteIncompleteCells();
  for (SurfaceParticle3D* particle : particles) {
    particle->force = {0.0,0.0,0.0};
  }
  for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++) {
    (*cellFields)[ctype]->mechanics->ParticleMechanics(particles_per_cell,lpc,ctype);
  }
  
}

void HemoParticleField3D::applyRepulsionForce() {
    
}


void HemoParticleField3D::spreadParticleForce(Box3D domain) {
  vector<SurfaceParticle3D*> localParticles;
  findParticles(domain,localParticles);
  SurfaceParticle3D * sparticle;

  for (pluint i = 0; i < localParticles.size(); i++ ) {
    sparticle = localParticles[i];
    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[sparticle->celltype]->kernelMethod(*atomicLattice,sparticle);

    //Directly change the force on a node , Palabos developers hate this one
    //quick non-functional trick.
    for (pluint j = 0; j < sparticle->kernelLocations.size(); j++) {
      //Yay for direct access
      sparticle->kernelLocations[j]->external.data[0] += (sparticle->force[0] * sparticle->kernelWeights[j]);
      sparticle->kernelLocations[j]->external.data[1] += (sparticle->force[1] * sparticle->kernelWeights[j]);
      sparticle->kernelLocations[j]->external.data[2] += (sparticle->force[2] * sparticle->kernelWeights[j]);
    }

  }
}

void HemoParticleField3D::interpolateFluidVelocity(Box3D domain) {
  vector<SurfaceParticle3D*> localParticles;
  findParticles(domain,localParticles);
  //TODO, remove casting
  SurfaceParticle3D * sparticle;
  //Prealloc is nice
  Array<double,3> velocity;
  Array<double,3> velocity_comp;


  for (pluint i = 0; i < localParticles.size(); i++ ) {
    sparticle = localParticles[i];
    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[sparticle->celltype]->kernelMethod(*atomicLattice,sparticle);

    //We have the kernels, now calculate the velocity of the particles.
    //Palabos developers, sorry for not using a functional...
    velocity = {0.0,0.0,0.0};
    for (pluint j = 0; j < sparticle->kernelLocations.size(); j++) {
      //Yay for direct access
      sparticle->kernelLocations[j]->computeVelocity(velocity_comp);
      velocity += velocity_comp * sparticle->kernelWeights[j];
    }
    sparticle->v = velocity;
    /*
    if (sparticle->kernelLocations.size() == 0) {
      cerr << "Location: " << sparticle->position[0] << " " << sparticle->position[1] << " " << sparticle->position[2] << std::endl;
    Box3D temp = getBoundingBox();
    Dot3D tmp = getLocation();
    cerr << "Box: " << temp.x0 << " " << temp.x1 << " " << temp.y0  << " "<< temp.y1  << " "<< temp.z0 <<" "<< temp.z1 <<std::endl;
    cerr << "Loc: " << tmp.x << " " << tmp.y << " " << tmp.z << std::endl;
    }*/
  }

}

HemoParticleDataTransfer3D& HemoParticleField3D::getDataTransfer() {
    return dataTransfer;
}
HemoParticleDataTransfer3D const& HemoParticleField3D::getDataTransfer() const {
    return dataTransfer;
}

std::string HemoParticleField3D::getBlockName() {
    return std::string("HemoParticleField3D");
}

hemoCellFields* HemoParticleField3D::cellFields=0;

#endif