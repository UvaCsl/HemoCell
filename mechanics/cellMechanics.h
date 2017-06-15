#ifndef HEMO_CELLMECHANICS
#define HEMO_CELLMECHANICS

class CellMechanics;
#include "hemocell_internal.h"
#include "hemoCellParticle.h"
#include "hemoCellFields.h"
#include "commonCellConstants.h"

class CellMechanics {
  public:
  const CommonCellConstants cellConstants;
  
  CellMechanics(HemoCellField & cellfield);

  virtual void ParticleMechanics(map<int,vector<HemoCellParticle *>> &,const map<int,bool> &, pluint ctype) = 0 ;
  virtual void statistics() = 0;
};
#endif
