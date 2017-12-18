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
#ifndef PARTICLE_HDF5IO_H
#define PARTICLE_HDF5IO_H

#include "hemocell_internal.h"
#include "hemoCellParticle.h"
#include "hemoCellFields.h"
#include "hemoCellField.h"

void writeCellField3D_HDF5(HemoCellFields& cellFields, double dx, double dt, plint iter, std::string preString="");


class WriteCellField3DInMultipleHDF5Files : public BoxProcessingFunctional3D
{
public:
    WriteCellField3DInMultipleHDF5Files (
            HemoCellField & cellField3D_,
            plint iter_, std::string identifier_,
            double dx_, double dt_, int i);
    /// Arguments: [0] Particle-field. [1] Lattice.
    ~WriteCellField3DInMultipleHDF5Files(){}; //Fuck C c++
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual WriteCellField3DInMultipleHDF5Files* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    HemoCellField& cellField3D;
    plint iter;
    std::string identifier;
    double dx;
    double dt;
    int ctype;
};
#endif  // FICSION_PARTICLE_HDF5IO_H

