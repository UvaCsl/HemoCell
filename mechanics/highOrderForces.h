#ifndef HEMO_highOrderForces
#define HEMO_highOrderForces

#include "cellMechanics.h"
#include "commonCellConstants.h"

class HighOrderForces : public CellMechanics {

  public:
  //Variables
  const CommonCellConstants cellConstants;
  HemoCellField & cellField;
  const double k_volume;
  const double k_area;
  const double k_inPlane;
  const double k_bend;


  //Constructor
  public:
  HighOrderForces(HemoCellField & cellField_, double k_volume_, double k_area_, double k_inPlane_, double k_bend_) ;
  void ParticleMechanics(map<int,vector<SurfaceParticle3D *>> particles_per_cell, map<int,bool> lpc, pluint ctype) ;

  void statistics();
};
#endif
