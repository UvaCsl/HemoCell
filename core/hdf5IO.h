#ifndef FICSION_HDF5IO_H
#define FICSION_HDF5IO_H

#include <limits>
// #include "ficsion.h"
#include "palabos3D.h"

#include <vector>
#include <hdf5.h>
#include <hdf5_hl.h>

using namespace plb;
using namespace std;



template<typename T, template<typename U> class Descriptor>
class WriteInMultipleHDF5Files : public BoxProcessingFunctional3D
{
public:
    WriteInMultipleHDF5Files (
            std::vector<std::string> & hdf5ContainerNames_,
            std::vector<plint> & hdf5ContainerDimensions_,
            plint iter_);
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual WriteInMultipleHDF5Files<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<std::string> & hdf5ContainerNames;
    std::vector<plint> & hdf5ContainerDimensions;
    plint iter;
};

template<typename T, template<typename U> class Descriptor>
void writeHDF5(MultiBlockLattice3D<T, Descriptor>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter);

#include "hdf5IO.cpp"
#endif  // FICSION_HDF5IO_H
