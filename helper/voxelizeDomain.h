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
                          VoxelizedDomain3D<double> *&voxelizedDomain, MultiScalarField3D<int> *&flagMatrix, plint blockSize);
#endif
