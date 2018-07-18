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
#ifndef HEMOCELL_PARTICLE_DATA_TRANSFER_H
#define HEMOCELL_PARTICLE_DATA_TRANSFER_H

namespace hemo {
  class HemoCellParticleDataTransfer;
}

#include "atomicBlock/atomicBlock3D.h"
#include "hemoCellParticleField.h"

namespace hemo {
using namespace plb;

class HemoCellParticleDataTransfer : public BlockDataTransfer3D {
public:
    HemoCellParticleDataTransfer();
    virtual plint staticCellSize() const;
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    //Much faster since we can circumvent memcpy (twice!)
    void receive(Box3D const & domain, char *, unsigned int size, modif::ModifT);
    void receive(Box3D const & domain, char *, unsigned int size, modif::ModifT, Dot3D absoluteOffset);
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
    inline plint getOffset(Dot3D &);
private:
    HemoCellParticleField* particleField;
    HemoCellParticleField const * constParticleField;
};
}
#endif
