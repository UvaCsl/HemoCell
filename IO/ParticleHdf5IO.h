/*#ifndef FICSION_PARTICLE_HDF5IO_H
#define FICSION_PARTICLE_HDF5IO_H

#include <limits>
// #include "ficsion.h"
#include "palabos3D.h"
#include "surfaceParticle3D.h"
#include "cellField3D.h"

#include <vector>
#include <string>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "cellFields3D.h"

using namespace plb;
using namespace std;

void writeCellField3D_HDF5(HemoCellField& cellField3D, double dx, double dt, plint iter, std::string preString="");


class WriteCellField3DInMultipleHDF5Files : public BoxProcessingFunctional3D
{
public:
    WriteCellField3DInMultipleHDF5Files (
            HemoCellField & cellField3D_,
            plint iter_, std::string identifier_,
            double dx_, double dt_);
    /// Arguments: [0] Particle-field. [1] Lattice.
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
};




#include "ParticleHdf5IO.cpp"
#endif  // FICSION_PARTICLE_HDF5IO_H
*/
