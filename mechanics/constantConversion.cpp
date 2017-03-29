#ifndef HEMO_constantConversion_cpp
#define HEMO_constantConversion_cpp
//This class converses SI things to LBM things for classes like the HO model.
#include "constantConversion.h"

void Parameters::lbm_parameters(Config & cfg, Box3D domainBox) {
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

    u_lbm_max = re * nu_lbm / domainBox.getNy();

    kBT_lbm = kBT_p/(df*dx);
};

double Parameters::calculate_kBend(Config & cfg, std::string  cellname ){
  return cfg["MaterialModel"][cellname]["kBend"].read<double>() * kBT_lbm;
};

double Parameters::calculate_kVolume(Config & cfg, std::string  cellname, MeshMetrics<double> & meshmetric){
  double kVolume =  cfg["MaterialModel"][cellname]["kVolume"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kVolume *= kBT_lbm/(eqLength*eqLength*eqLength);
  return kVolume;
};

double Parameters::calculate_kArea(Config & cfg, std::string  cellname, MeshMetrics<double> & meshmetric){
  double kArea =  cfg["MaterialModel"][cellname]["kShear"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kArea *= kBT_lbm/(eqLength*eqLength);
  return kArea;
}

double Parameters::calculate_kInPlane(Config & cfg, std::string  cellname, MeshMetrics<double> & meshmetric){
  double kInPlane = cfg["MaterialModel"][cellname]["kWLC"].read<double>();
  double plf = cfg["MaterialModel"][cellname]["persistenceLengthFine"].read<double>();
  double plc = plf/dx * sqrt((meshmetric.getNumVertices()-2.0) / (23867-2.0)); //Kanadiakis magic
  kInPlane *= kBT_lbm/(4.0*plc);
  return kInPlane;
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

double Parameters::kBT_lbm = 0.0;
double Parameters::kBT_p = 0.0;
#endif
