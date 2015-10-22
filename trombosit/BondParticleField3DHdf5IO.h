#ifndef FICSION_BONDPARTICLEFIELD3D_HDF5IO_H
#define FICSION_BONDPARTICLEFIELD3D_HDF5IO_H

#include <limits>
// #include "ficsion.h"
#include "palabos3D.h"
#include "BondParticle3D.h"
#include "immersedCellParticle3D.h"
#include "cellField3D.h"

#include <vector>
#include <string>
#include <hdf5.h>
#include <hdf5_hl.h>

using namespace plb;
using namespace std;

namespace trombocit {

template<typename T, template<typename U> class Descriptor>
void writeBondParticleField3D_HDF5(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField, T dx, T dt, plint iter, std::string identifier);


template<typename T, template<typename U> class Descriptor>
class WriteBondParticleField3D : public BoxProcessingFunctional3D
{
public:
	WriteBondParticleField3D (
            plint iter_, std::string identifier_,
            T dx_, T dt_);
    /// Arguments: [0] Particle-field. [1] Lattice.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual WriteBondParticleField3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iter;
    std::string identifier;
    T dx;
    T dt;
};


}

#include "BondParticleField3DHdf5IO.cpp"
#endif  // FICSION_BONDPARTICLEFIELD3D_HDF5IO_H
