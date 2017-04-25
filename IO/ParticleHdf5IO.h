#ifndef FICSION_PARTICLE_HDF5IO_H
#define FICSION_PARTICLE_HDF5IO_H

#include "hemocell_internal.h"
#include "hemoCellParticle.h"
#include "hemoCellFields.h"
#include "hemoCellParticleType.h"

#include <hdf5.h>
#include <hdf5_hl.h>

void writeCellField3D_HDF5(hemoCellFields& cellFields, double dx, double dt, plint iter, std::string preString="");


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

