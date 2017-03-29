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
  static double nu_lbm, u_lbm_max, kBT_lbm;


  static void lbm_parameters(Config & cfg, Box3D domainBox);
  static double calculate_kBend(Config & cfg, std::string );
  static double calculate_kVolume(Config & cfg, std::string , MeshMetrics<double> &);
  static double calculate_kArea(Config & cfg, std::string , MeshMetrics<double> &);
  static double calculate_kInPlane(Config & cfg, std::string , MeshMetrics<double> &);
};
#endif
