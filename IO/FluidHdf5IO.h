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
#ifndef FLUID_HDF5_IO_H
#define FLUID_HDF5_IO_H

#include "hemocell_internal.h"
#include "hemoCellFields.h"

void writeFluidField_HDF5(HemoCellFields& cellFields, double dx, double dt, plint iter, string preString="");

#ifndef hsize_t
typedef long long unsigned int hsize_t;
#endif

class WriteFluidField : public BoxProcessingFunctional3D
{
public:
    WriteFluidField (
            HemoCellFields& cellfields,
						MultiBlockLattice3D<double,DESCRIPTOR>& fluid,
            plint iter_, string identifier_,
            double dx_, double dt_);
    ~WriteFluidField(){}; //Fuck C c++
    virtual void processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> fields);
    virtual WriteFluidField* clone() const;
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    HemoCellFields& cellfields;
    MultiBlockLattice3D<double,DESCRIPTOR>& fluid;
    plint iter;
    string identifier;
    double dx;
    double dt;
    Box3D * odomain;
    BlockLattice3D<double,DESCRIPTOR> * ablock;
    HemoCellParticleField * particlefield;
    int blockid;
    hsize_t * nCells;

    float * outputVelocity();
    float * outputForce();
    float * outputDensity();
    float * outputCellDensity(string name);
    
};

#endif
