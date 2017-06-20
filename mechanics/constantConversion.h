#ifndef HEMO_constantConversion_h
#define HEMO_constantConversion_h
//This class converses SI things to LBM things for classes like the HO model.
#include "config.h"
#include "meshMetrics.h"
#include "palabos3D.h"

class Parameters {
  //Variables for the simulation
  public:
  static double dt,dx,dm,df;
  static double nu_p,rho_p,kBT_p;
  static double tau,re;
  static double nu_lbm, u_lbm_max, kBT_lbm, shearrate_lbm;
#ifdef FORCE_LIMIT
  static double f_limit;
#endif
  static double ef_lbm; //used in cellStretching


  static void lbm_pipe_parameters(Config & cfg, Box3D domainBox);
  static void lbm_shear_parameters(Config & cfg, double nx);
  static void lbm_base_parameters(Config & cfg);
  static void printParameters();
};
#endif
