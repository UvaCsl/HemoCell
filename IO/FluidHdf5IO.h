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
    hsize_t * nCells;

    float * outputVelocity();
    float * outputForce();

};

#endif
