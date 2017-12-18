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
#ifndef READ_POSISIONS_OF_BLOOD_CELLS_H
#define READ_POSISIONS_OF_BLOOD_CELLS_H

#include "hemocell_internal.h"
#include "hemoCellFields.h"
#include "config.h"
#include "constantConversion.h"
void readPositionsBloodCellField3D(HemoCellFields & cellFields, double dx, Config & cfg);
int getTotalNumberOfCells(HemoCellFields & cellFields);

void getReadPositionsBloodCellsVector(Box3D realDomain,
                                           std::vector<TriangularSurfaceMesh<double>* > & meshes,
                                           std::vector<plint> & Np,
                                           std::vector<std::vector<hemo::Array<double,3> > > & positions,
                                           std::vector<std::vector<plint> > & cellIds,
                                           std::vector<std::vector<hemo::Array<double,3> > > & randomAngles,
                                           Config & cfg, HemoCellFields & cellFields,
                                           HemoCellParticleField & particleField);

class ReadPositionsBloodCellField3D : public BoxProcessingFunctional3D
{
public:
    ReadPositionsBloodCellField3D (HemoCellFields & cellFields_, double dx_, Config & cfg_)
            : cellFields(cellFields_), cfg(cfg_) {dx = dx_;}
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ReadPositionsBloodCellField3D* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    double dx;
    HemoCellFields & cellFields;
    Config & cfg;
};

#endif
