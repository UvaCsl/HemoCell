#ifndef HEMOCELL_PARTICLE_DATA_TRANSFER_H
#define HEMOCELL_PARTICLE_DATA_TRANSFER_H

#include "hemocell_internal.h"
class HemoCellParticleField;

class HemoCellParticleDataTransfer : public BlockDataTransfer3D {
public:
    HemoCellParticleDataTransfer();
    virtual plint staticCellSize() const;
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset);
    virtual void receive( Box3D domain, std::vector<char> const& buffer,
                          modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    virtual void setBlock(AtomicBlock3D& block);
    virtual void setConstBlock(AtomicBlock3D const& block) ;
    virtual HemoCellParticleDataTransfer* clone() const {
      return new HemoCellParticleDataTransfer(*this);
    }
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset);
private:
    HemoCellParticleField* particleField;
    HemoCellParticleField const * constParticleField;
};

#endif
