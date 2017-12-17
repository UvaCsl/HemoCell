#ifndef HEMOCELL_WBCHIGHORDERMODEL_H
#define HEMOCELL_WBCHIGHORDERMODEL_H

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"

class WbcHighOrderModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const double k_volume;
  const double k_area;
  const double k_link;
  const double k_bend;
  const double eta_m;
  const double eta_v;
  const double k_cytoskeleton;
  const double k_inner_rigid;
  const double core_radius;
  const double radius;

  public:
  WbcHighOrderModel(Config & modelCfg_, HemoCellField & cellField_) ;

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> &lpc, size_t ctype) ;

  void statistics();

  static double calculate_kBend(Config & cfg, MeshMetrics<double> &);
  static double calculate_kVolume(Config & cfg, MeshMetrics<double> &);
  static double calculate_kArea(Config & cfg, MeshMetrics<double> &);
  static double calculate_kLink(Config & cfg, MeshMetrics<double> &);
  static double calculate_etaM(Config & cfg );
  static double calculate_etaV(Config & cfg );
  static double calculate_coreRadius(Config & cfg );
  static double calculate_radius(Config & cfg );
  static double calculate_kInnerRigid(Config & cfg );
  static double calculate_kCytoskeleton(Config & cfg );
};

#endif
