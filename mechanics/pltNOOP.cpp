#ifndef HEMO_pltNOOP_cpp
#define HEMO_pltNOOP_cpp
#include "cellMechanics.h"
class PltNOOP : public CellMechanics {
  public:
  PltNOOP() :CellMechanics() {};

  void ParticleMechanics(map<int,vector<SurfaceParticle3D *>>,map<int,bool>, pluint ctype) {} ;
  void statistics () {
    cerr << "PLT Mechinical model is NOOP";
  }
};

#endif
