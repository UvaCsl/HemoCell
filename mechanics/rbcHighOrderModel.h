#ifndef HEMOCELL_RBCHIGHORDERMODEL_H
#define HEMOCELL_RBCHIGHORDERMODEL_H

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"

class RbcHighOrderModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const double k_volume;
  const double k_area;
  const double k_link;
  const double k_bend;
  const double eta_m;
  const double eta_v;
#ifdef FORCE_LIMIT
  const double areaRatio_limit;
  const double volumeRatio_limit; 
  const double linkRatio_limit; 
  const double bendingRatio_limit;
#endif
  public:
  RbcHighOrderModel(Config & modelCfg_, HemoCellField & cellField_) ;

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> &lpc, pluint ctype) ;

  void statistics();

  static double calculate_kBend(Config & cfg);
  static double calculate_kVolume(Config & cfg, MeshMetrics<double> &);
  static double calculate_kArea(Config & cfg, MeshMetrics<double> &);
  static double calculate_kLink(Config & cfg, MeshMetrics<double> &);
  static double calculate_etaM(Config & cfg );
  static double calculate_etaV(Config & cfg );
#ifdef FORCE_LIMIT
  static double calculate_areaRatio_limit(Config & cfg, MeshMetrics<double> &);
  static double calculate_volumeRatio_limit(Config & cfg, MeshMetrics<double> &);
  static double calculate_linkRatio_limit(Config & cfg, MeshMetrics<double> &);
  static double calculate_bendRatio_limit(Config & cfg);
#endif


};

#endif
