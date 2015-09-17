#ifndef FICSION_PARTICLE_HDF5IO_H
#define FICSION_PARTICLE_HDF5IO_H

#include <limits>
// #include "ficsion.h"
#include "palabos3D.h"
#include "immersedCellParticle3D.h"
#include "cellField3D.h"

#include <vector>
#include <string>
#include <hdf5.h>
#include <hdf5_hl.h>

using namespace plb;
using namespace std;

template<typename T, template<typename U> class Descriptor>
void writeCellField3D_HDF5(CellField3D<T, Descriptor>& cellField3D, T dx, T dt, plint iter, std::string preString="");


template<typename T, template<typename U> class Descriptor>
class WriteCellField3DInMultipleHDF5Files : public BoxProcessingFunctional3D
{
public:
    WriteCellField3DInMultipleHDF5Files (
            CellField3D<T, Descriptor>& cellField3D_,
            plint iter_, std::string identifier_,
            T dx_, T dt_);
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual WriteCellField3DInMultipleHDF5Files<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    CellField3D<T, Descriptor>& cellField3D;
    plint iter;
    std::string identifier;
    T dx;
    T dt;
};




#include "ParticleHdf5IO.cpp"
#endif  // FICSION_PARTICLE_HDF5IO_H
