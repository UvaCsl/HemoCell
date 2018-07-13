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
#include "hemoCellParticleDataTransfer.h"
#include "hemoCellParticleField.h"
#include "hemocell.h"

///Used to circumvent buffer initialization of characters
struct NoInitChar
{
    char value;
    NoInitChar() {
        // do nothing
        static_assert(sizeof *this == sizeof value, "invalid size");
        static_assert(__alignof *this == __alignof value, "invalid alignment");
    }
};


/* *************** class HemoParticleDataTransfer3D ************************ */

inline plint HemoCellParticleDataTransfer::getOffset(Dot3D & absoluteOffset) {
      int offset = 0;
      //Offset only happens for wrapping around. Decide the new (or old) CellID to use
      if (absoluteOffset.x > 0) {
        offset ++;
      } else if (absoluteOffset.x < 0) {
        offset --;
      }
      if (absoluteOffset.y > 0) {
        offset += particleField->cellFields->periodicity_limit_offset_y;
      } else if (absoluteOffset.y < 0) {
        offset -= particleField->cellFields->periodicity_limit_offset_y;
      }
      if (absoluteOffset.z > 0) {
        offset += particleField->cellFields->periodicity_limit_offset_z;
      } else if (absoluteOffset.z < 0) {
        offset -= particleField->cellFields->periodicity_limit_offset_z;
      }

      offset*=particleField->cellFields->number_of_cells;

      return offset;
}

HemoCellParticleDataTransfer::HemoCellParticleDataTransfer () {};

plint HemoCellParticleDataTransfer::staticCellSize() const {
    return 0;  // Particle containers have only dynamic data.
}

void HemoCellParticleDataTransfer::send (
        Box3D domain, std::vector<char>& buffer, modif::ModifT kind ) const
{
  constParticleField->cellFields->hemocell.statistics.getCurrent()["MpiSend"].start();
    buffer.clear();
    std::vector<NoInitChar> * bufferNoInit = reinterpret_cast<std::vector<NoInitChar>*>(&buffer);
    
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::hemocell))
    {
        std::vector<HemoCellParticle*> foundParticles;
        particleField->findParticles(domain, foundParticles);
        bufferNoInit->resize(sizeof(HemoCellParticle::serializeValues_t)*foundParticles.size());
        pluint offset=0;
        for (HemoCellParticle * iParticle : foundParticles) {
          *((HemoCellParticle::serializeValues_t*)&(*bufferNoInit)[offset]) = iParticle->sv;
          offset += sizeof(HemoCellParticle::serializeValues_t);
        }
    }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::receive (Box3D const & domain, char * buffer, unsigned int size, modif::ModifT kind)
{
  PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));
  // Clear the existing data before introducing the new data.
  //particleField->removeParticles(domain);
  // Particles, by definition, are dynamic data, and they need to
  //   be reconstructed in any case. Therefore, the receive procedure
  //   is run whenever kind is one of the dynamic types.
 constParticleField->cellFields->hemocell.statistics.getCurrent()["MpiReceive"].start();

  if ( (kind==modif::hemocell) )
  {
      unsigned int posInBuffer = 0;

      HemoCellParticle::serializeValues_t * newParticle;;
      while (posInBuffer < size) {
        // 1. Generate dynamics object, and unserialize dynamic data.
        newParticle = (HemoCellParticle::serializeValues_t*)&buffer[posInBuffer];
        posInBuffer += sizeof(HemoCellParticle::serializeValues_t);       
        particleField->addParticle(domain, *newParticle);
      }
  }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::receive (Box3D const & domain, char * buffer, unsigned int size, modif::ModifT kind, Dot3D absoluteOffset )
{
  if (absoluteOffset.x == 0 && absoluteOffset.y == 0 && absoluteOffset.z == 0) {
    receive(domain,buffer,size,kind);
    return;
  }
  
  constParticleField->cellFields->hemocell.statistics.getCurrent()["MpiReceive"].start();

  if ( (kind==modif::hemocell) )
  {
    int offset = getOffset(absoluteOffset);
    hemo::Array<T,3> realAbsoluteOffset({(T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z});
    unsigned int posInBuffer = 0;
    
    HemoCellParticle::serializeValues_t * newParticle;;
    while (posInBuffer < size) {
        // 1. Generate dynamics object, and unserialize dynamic data.
      newParticle = (HemoCellParticle::serializeValues_t*)&buffer[posInBuffer];
      posInBuffer += sizeof(HemoCellParticle::serializeValues_t);
      //Edit in buffer, but it is not used again anyway
      newParticle->position += realAbsoluteOffset;
      newParticle->cellId += offset;
      particleField->addParticle(domain, *newParticle);
    }
  }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind )
{
  //cast to char* because we ARE going to edit
  receive(domain,(char*)&buffer[0],buffer.size(),kind);
}

void HemoCellParticleDataTransfer::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset )
{
  //cast to char* because we ARE going to edit
  receive(domain,(char*)&buffer[0],buffer.size(),kind,absoluteOffset);
}

void HemoCellParticleDataTransfer::setBlock(AtomicBlock3D& block) {
   particleField = dynamic_cast<HemoCellParticleField*>(&block);
   PLB_ASSERT(particleField);
   constParticleField = particleField;
}

void HemoCellParticleDataTransfer::setConstBlock(AtomicBlock3D const& block) {
          constParticleField = dynamic_cast<HemoCellParticleField const*>(&block);
}

void HemoCellParticleDataTransfer::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind )
{
  constParticleField->cellFields->hemocell.statistics.getCurrent()["LocalCommunication"].start();

    if ( kind==modif::hemocell )
    {
      Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
      HemoCellParticleField const& fromParticleField =
          dynamic_cast<HemoCellParticleField const&>(from);
      vector<const HemoCellParticle *> particles;
      fromParticleField.findParticles(fromDomain, particles);
      //Calling addParticle on self can invalidate particles pointer array on realloc from vector
      //Do for every local communication to accomodate overcoupling particle field in the future.
      vector<HemoCellParticle::serializeValues_t> sv_values;
      sv_values.reserve(particles.size());
      for (const HemoCellParticle * particle : particles) {
        sv_values.emplace_back(particle->sv);
      }     
      for (const HemoCellParticle::serializeValues_t & sv : sv_values) {
        particleField->addParticle(toDomain, sv);
      }
    }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset )
{
  if (absoluteOffset.x == 0 && absoluteOffset.y == 0 && absoluteOffset.z == 0) {
    attribute(toDomain,deltaX,deltaY,deltaZ,from,kind);
    return;
  }
  
  constParticleField->cellFields->hemocell.statistics.getCurrent()["LocalCommunication"].start();
  
  if ( kind==modif::hemocell )
  { 
  
    Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
    HemoCellParticleField const& fromParticleField =
        dynamic_cast<HemoCellParticleField const&>(from);
    vector<const HemoCellParticle *> particles;
    fromParticleField.findParticles(fromDomain, particles);
    int offset = getOffset(absoluteOffset);
    hemo::Array<T,3> realAbsoluteOffset({(T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z});

    //Calling addParticle on self can invalidate particles pointer array on realloc from vector
    //Do for every local communication to accomodate overcoupling particle field in the future.
    vector<HemoCellParticle::serializeValues_t> sv_values;
    sv_values.reserve(particles.size());
    for (const HemoCellParticle * particle : particles) {
      sv_values.emplace_back(particle->sv);
      sv_values.back().position += realAbsoluteOffset;
      sv_values.back().cellId += offset;
    }     
    for (const HemoCellParticle::serializeValues_t & sv : sv_values) {
      particleField->addParticle(toDomain, sv);
    }

  }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}