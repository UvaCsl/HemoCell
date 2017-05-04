#ifndef HEMO_highOrderForces
#define HEMO_highOrderForces

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"
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
  void ParticleMechanics(map<int,vector<HemoCellParticle *>> particles_per_cell, map<int,bool> lpc, pluint ctype) ;

  void statistics();
};

class HighOrderForcesXML : public HighOrderForces {
  public:
  HighOrderForcesXML(Config & cfg,HemoCellField & cellField_, string name) :
                        HighOrderForces(cellField_,
                                        param::calculate_kVolume(cfg,name,*cellField_.meshmetric),
                                        param::calculate_kArea(cfg,name,*cellField_.meshmetric),
                                        param::calculate_kInPlane(cfg,name,*cellField_.meshmetric),
                                        param::calculate_kBend(cfg,name)) {};
};
#endif
