#ifndef HEMOCELL_PLTSIMPLEMODEL_H
#define HEMOCELL_PLTSIMPLEMODEL_H

#include "hemocell_internal.h"
#include "constantConversion.h"
#include "config.h"
#include "cellMechanics.h"
#include "commonCellConstants.h"
#include "hemoCellFields.h"

class PltSimpleModel : public CellMechanics {

  public:
  //Variables
  HemoCellField & cellField;
  const double k_volume;
  const double k_area;
  const double k_link;
  const double k_bend;
  const double eta;


  //Constructor
  public:
  PltSimpleModel(Config & modelCfg_, HemoCellField & cellField_);

  void ParticleMechanics(map<int,vector<HemoCellParticle *>> &particles_per_cell, const map<int,bool> &lpc, pluint ctype);

  void statistics();

  static double calculate_kBend(Config & cfg);
  static double calculate_kVolume(Config & cfg, MeshMetrics<double> &);
  static double calculate_kArea(Config & cfg, MeshMetrics<double> &);
  static double calculate_kLink(Config & cfg, MeshMetrics<double> &);
  static double calculate_eta(Config & cfg );

};

#endif
