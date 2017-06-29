#ifndef HEMOCELL_VWFMODEL_H
#define HEMOCELL_VFWMODEL_H

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "hemoCellFields.h"

class vWFModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const double k_link;
  const double k_bend;
  const double link_eq;
  const double bend_eq;
  const double eta_m;
  const double eta_v;

  public:
  vWFModel(Config & modelCfg_, HemoCellField & cellField_) ;

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> &lpc, pluint ctype) ;

  void statistics();

  static double calculate_kBend(Config & cfg);
  static double calculate_kLink(Config & cfg);
  static double calculate_bend_eq(Config & cfg);
  static double calculate_link_eq(Config & cfg);
  static double calculate_etaM(Config & cfg );
  static double calculate_etaV(Config & cfg );

};

#endif
