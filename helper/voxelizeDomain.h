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
#ifndef VOXELIZEDOMAIN_H
#define VOXELIZEDOMAIN_H

#include "hemocell_internal.h"

class CopyFromNeighbor : public BoxProcessingFunctional3D_S<int> {
public:
    CopyFromNeighbor(hemo::Array<plint, 3> offset_) : offset(offset_) { };

    virtual void process(Box3D domain, ScalarField3D<int> &field1);

    virtual CopyFromNeighbor *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

    virtual BlockDomain::DomainT appliesTo() const;

private:
    hemo::Array<plint, 3> offset;
};

void getFlagMatrixFromSTL(std::string meshFileName, plint extendedEnvelopeWidth, plint refDirLength, plint refDir,
                          VoxelizedDomain3D<T> *&voxelizedDomain, MultiScalarField3D<int> *&flagMatrix, plint blockSize, int particleEnvelope = 0);
#endif
