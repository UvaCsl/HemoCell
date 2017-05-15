#ifndef HEMO_constantConversion_cpp
#define HEMO_constantConversion_cpp
//This class converses SI things to LBM things for classes like the HO model.
#include "constantConversion.h"

void Parameters::lbm_pipe_parameters(Config & cfg, Box3D domainBox) {
    dt = cfg["domain"]["dt"].read<double>();
    dx = cfg["domain"]["dx"].read<double>();
    nu_p = cfg["domain"]["nuP"].read<double>();
    rho_p = cfg["domain"]["rhoP"].read<double>();
    re = cfg["domain"]["Re"].read<double>();
    kBT_p = cfg["domain"]["kBT"].read<double>();

    if (dt < 0.0 ) { //dt is not set, calculate it from nu_p and dx, tau is set to 1
      tau = 1.0;
      nu_lbm = 1.0/3.0 * (tau - 0.5);
      dt = nu_lbm /nu_p * (dx*dx);
      cerr << "(main) Tau is set to unity, and dt is derived from that!" << std::endl;
    } else { //dt is set, this means we must set tau
      nu_lbm = nu_p * dt/ (dx*dx);
      tau = 3.0 * nu_lbm + 0.5;
    }

    dm = rho_p * (dx * dx *dx);
    df = dm * dx / (dt * dt);

    // TODO: This is only true for a pipe set along the x-axis.
    //       Make some document note, that Re is calculated using Ny!
    u_lbm_max = re * nu_lbm / domainBox.getNy();

    kBT_lbm = kBT_p/(df*dx);
};

void Parameters::lbm_shear_parameters(Config & cfg,double nx) {
  double shearrate_p = cfg["domain"]["shearrate"].read<double>();
  nu_p = cfg["domain"]["nuP"].read<double>();
  dx = cfg["domain"]["dx"].read<double>();
  dt = cfg["domain"]["dt"].read<double>();
  rho_p = cfg["domain"]["rhoP"].read<double>();
  kBT_p = cfg["domain"]["kBT"].read<double>();


  re = (nx* (shearrate_p * (nx*0.5))) / nu_p;
  
  if (dt < 0.0 ) { //dt is not set, calculate it from nu_p and dx, tau is set to 1
    tau = 1.0;
    nu_lbm = 1.0/3.0 * (tau - 0.5);
    dt = nu_lbm /nu_p * (dx*dx);
    cerr << "(main) Tau is set to unity, and dt is derived from that!" << std::endl;
  } else {
    nu_lbm = nu_p * dt/ (dx*dx);
    tau = 3.0 * nu_lbm + 0.5;
  }
  
  dm = rho_p * (dx * dx *dx);
  df = dm * dx / (dt * dt);

  // TODO: This is only true for a pipe set along the x-axis.
  //       Make some document note, that Re is calculated using Ny!
  u_lbm_max = re * nu_lbm / nx;
  
  shearrate_lbm = shearrate_p*dt;
  kBT_lbm = kBT_p/(df*dx);

}

double Parameters::dt = 0.0;
double Parameters::dx = 0.0;
double Parameters::dm = 0.0;
double Parameters::df = 0.0;
double Parameters::nu_p = 0.0;
double Parameters::rho_p = 0.0;
double Parameters::tau = 0.0;
double Parameters::re = 0.0;
double Parameters::nu_lbm = 0.0;
double Parameters::u_lbm_max = 0.0;
double Parameters::shearrate_lbm = 0.0;

double Parameters::kBT_lbm = 0.0;
double Parameters::kBT_p = 0.0;
#endif
