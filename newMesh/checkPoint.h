#ifndef FCN_CHECKPOINT_H
#define FCN_CHECKPOINT_H
#include "palabos3D.h"
#include "palabos3D.hh"
#include "cellField3D.h"

#include <vector>
#include <string>


template<typename T, template<typename U> class Descriptor>
class FcnCheckpoint
{
public:
	FcnCheckpoint(XMLreader & documentXML);
	FcnCheckpoint(std::string paramXmlFileName);
	~FcnCheckpoint();


	void read(MultiBlockLattice3D<T, Descriptor> & lattice_);
	void read(MultiParticleField3D<DenseParticleField3D<T,Descriptor> > & particleField, std::string identifier);
	void save(MultiBlockLattice3D<T, Descriptor> & lattice_, std::vector<CellField3D<T, Descriptor>* > & cellFields, plint iter);
	/* data */
private:
    XMLwriter xmlw;

};



#include "checkPoint.hh"
#endif // FCN_CHECKPOINT_H
