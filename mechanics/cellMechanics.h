#ifndef HEMO_CELLMECHANICS
#define HEMO_CELLMECHANICS

class CellMechanics;
#include "hemocell_internal.h"
#include "hemoCellParticle.h"

class CellMechanics {
  public:
  virtual void ParticleMechanics(map<int,vector<SurfaceParticle3D *>>,map<int,bool>, pluint ctype) = 0 ;
  virtual void statistics() = 0;
};
#endif
