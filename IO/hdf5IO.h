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

// Exporting files in XDMF switches X and Z axes.
template<typename T, template<typename U> class Descriptor>
void writeHDF5(MultiBlockLattice3D<T, Descriptor>& lattice,
              T dx, T dt, plint iter, bool invertXZ_for_XDMF=false);


template<typename T>
class WriteInMultipleHDF5Files : public BoxProcessingFunctional3D
{
public:
    WriteInMultipleHDF5Files (
            std::vector<std::string> & hdf5ContainerNames_,
            std::vector<plint> & hdf5ContainerDimensions_,
            plint iter_, T dx_, T dt_,
            plint envelopeWidth_=0, bool invertXZ_for_XDMF_=false);
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual WriteInMultipleHDF5Files<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<std::string> & hdf5ContainerNames;
    std::vector<plint> & hdf5ContainerDimensions;
    plint iter;
    T dx;
    T dt;
    bool invertXZ_for_XDMF;
    plint envelopeWidth;
};



#include "hdf5IO.cpp"
#endif  // FICSION_HDF5IO_H
