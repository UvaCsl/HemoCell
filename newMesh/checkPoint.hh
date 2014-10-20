#ifndef FCN_CHECKPOINT_HH
#define FCN_CHECKPOINT_HH
#include "checkPoint.h"



template<typename T, template<typename U> class Descriptor>
void FcnCheckpoint<T, Descriptor>::save(MultiBlockLattice3D<T, Descriptor> & lattice_, std::vector<CellField3D<T, DESCRIPTOR>* > & cellFields, plint iter) {

}



#endif // FCN_CHECKPOINT_HH
