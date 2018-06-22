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
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        std::vector<HemoCellParticle*> foundParticles;
        particleField->findParticles(domain, foundParticles);
        for (pluint iParticle=0; iParticle<foundParticles.size(); ++iParticle) {
        // The serialize function automatically reallocates memory for buffer.
        serialize(*foundParticles[iParticle], buffer);
        }
    }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind )
{
    PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));
    // Clear the existing data before introducing the new data.
    //particleField->removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
   constParticleField->cellFields->hemocell.statistics.getCurrent()["MpiReceive"].start();

    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        pluint posInBuffer = 0;
        
        HemoCellParticle newParticle = HemoCellParticle();
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            newParticle.unserialize(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField->addParticle(domain, newParticle.sv);
        }
    }
    constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::receive (
        Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset )
{
  
  constParticleField->cellFields->hemocell.statistics.getCurrent()["MpiReceive"].start();

  if ( (kind==modif::dynamicVariables) ||
       (kind==modif::allVariables) ||
       (kind==modif::dataStructure) )
  {
    int offset = getOffset(absoluteOffset);
    hemo::Array<T,3> realAbsoluteOffset({(T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z});
    pluint posInBuffer = 0;
    HemoCellParticle newParticle = HemoCellParticle();
    while (posInBuffer < buffer.size()) {
      // 1. Generate dynamics object, and unserialize dynamic data.
      HierarchicUnserializer unserializer(buffer, posInBuffer);
      newParticle.unserialize(unserializer);
      posInBuffer = unserializer.getCurrentPos();
      newParticle.sv.position += realAbsoluteOffset;
      newParticle.sv.cellId += offset;
      particleField->addParticle(domain, newParticle.sv);
    }
  }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();

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

    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
      Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
      HemoCellParticleField const& fromParticleField =
          dynamic_cast<HemoCellParticleField const&>(from);
      vector<const HemoCellParticle *> particles;
      fromParticleField.findParticles(fromDomain, particles);
      for (const HemoCellParticle * particle : particles) {
        particleField->addParticle(toDomain,particle->sv);
      }
    }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}

void HemoCellParticleDataTransfer::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset )
{
  constParticleField->cellFields->hemocell.statistics.getCurrent()["LocalCommunication"].start();
  
  if ( (kind==modif::dynamicVariables) ||
       (kind==modif::allVariables) ||
       (kind==modif::dataStructure) )
  { 
  
    Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
    HemoCellParticleField const& fromParticleField =
        dynamic_cast<HemoCellParticleField const&>(from);
    vector<const HemoCellParticle *> particles;
    fromParticleField.findParticles(fromDomain, particles);
    int offset = getOffset(absoluteOffset);
    hemo::Array<T,3> realAbsoluteOffset({(T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z});

    HemoCellParticle::serializeValues_t sv;
    for (const HemoCellParticle * particle : particles) {
      sv = particle->sv;
      sv.position += realAbsoluteOffset;
      sv.cellId += offset;
      particleField->addParticle(toDomain, sv);
    }     
  }
  constParticleField->cellFields->hemocell.statistics.getCurrent().stop();
}