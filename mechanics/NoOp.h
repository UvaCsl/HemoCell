#ifndef HEMOCELL_NOOP_H
#define HEMOCELL_NOOP_H
#include "cellMechanics.h"
class NoOp : public CellMechanics {
  public:
  NoOp() :CellMechanics() {};

  inline void ParticleMechanics(map<int,vector<HemoCellParticle *>>,map<int,bool>, pluint ctype) {} ;
  inline void statistics () {
    cerr << "Mechanical model is NoOp";
  }
};

#endif
