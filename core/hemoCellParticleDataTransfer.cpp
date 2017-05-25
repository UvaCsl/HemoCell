#ifndef HEMOCELL_PARTICLE_DATA_TRANSFER_CPP
#define HEMOCELL_PARTICLE_DATA_TRANSFER_CPP

#include "hemoCellParticleDataTransfer.h"
#include "hemoCellParticleField.h"

/* *************** class HemoParticleDataTransfer3D ************************ */

HemoCellParticleDataTransfer::HemoCellParticleDataTransfer () {};

plint HemoCellParticleDataTransfer::staticCellSize() const {
    return 0;  // Particle containers have only dynamic data.
}

void HemoCellParticleDataTransfer::send (
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
        std::vector<HemoCellParticle*> foundParticles;
        particleField->findParticles(domain, foundParticles);
        for (pluint iParticle=0; iParticle<foundParticles.size(); ++iParticle) {
        // The serialize function automatically reallocates memory for buffer.
        serialize(*foundParticles[iParticle], buffer);
        }
    }
}

void HemoCellParticleDataTransfer::receive (
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
            HemoCellParticle* newParticle = new HemoCellParticle();
            newParticle->unserialize(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField->addParticle(domain, newParticle);
        }
    }
}

void HemoCellParticleDataTransfer::receive (
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
            HemoCellParticle* newParticle = new HemoCellParticle();
            newParticle->unserialize(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle->position += realAbsoluteOffset;
            particleField->addParticle(domain, newParticle);
        }
    }
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
    Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
    std::vector<char> buffer;
    HemoCellParticleField const& fromParticleField =
        dynamic_cast<HemoCellParticleField const&>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

void HemoCellParticleDataTransfer::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset )
{
    Box3D fromDomain(toDomain.shift(deltaX,deltaY,deltaZ));
    std::vector<char> buffer;
    HemoCellParticleField const& fromParticleField =
        dynamic_cast<HemoCellParticleField const&>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}



#endif
