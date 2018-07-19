/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//This class converses SI things to LBM things for classes like the HO model.
#include "constantConversion.h"
#include "meshMetrics.h"
#include "hemoCellFunctional.h"
#include "logfile.h"

#include "multiBlock/multiDataField3D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.hh"
#include "dataProcessors/dataAnalysisFunctional3D.hh"
#include "dataProcessors/dataAnalysisWrapper3D.hh"

namespace hemo {
void Parameters::lbm_base_parameters(Config & cfg) {
    dt = cfg["domain"]["dt"].read<T>();
    dx = cfg["domain"]["dx"].read<T>();
    nu_p = cfg["domain"]["nuP"].read<T>();
    rho_p = cfg["domain"]["rhoP"].read<T>();
    kBT_p = cfg["domain"]["kBT"].read<T>();

    if (dt < 0.0 ) { //dt is not set, calculate it from nu_p and dx, tau is set to 1
      tau = 1.0;
      nu_lbm = 1.0/3.0 * (tau - 0.5);
      dt = nu_lbm /nu_p * (dx*dx);
      hlog << "(HemoCell) dt is set to *auto*. Tau will be set to 1!" << std::endl;
    } else { //dt is set, this means we must set tau
      nu_lbm = nu_p * dt/ (dx*dx);
      tau = 3.0 * nu_lbm + 0.5;
    }

    dm = rho_p * (dx * dx *dx);
    df = dm * dx / (dt * dt);
#ifdef FORCE_LIMIT
    f_limit = FORCE_LIMIT / 1.0e12 / df; // Changing pN to lbm force
#endif
    kBT_lbm = kBT_p/(df*dx);
};

void Parameters::lbm_pipe_parameters(Config & cfg, plb::MultiScalarField3D<int> * sf) {
    Parameters::lbm_base_parameters(cfg);
    re = cfg["domain"]["Re"].read<T>();
    
    plb::Box3D domain = sf->getBoundingBox();
    domain.x1 = domain.x0;
    T fluidArea = plb::computeSum<int>(*sf,domain);
    
    pipe_radius = sqrt(fluidArea/PI);
    hlog << "(Parameters) Your pipe has a calculated radius of " << pipe_radius << " LU, assuming a perfect circle" << std::endl;
    u_lbm_max = re * nu_lbm / (pipe_radius*2);
};

void Parameters::lbm_pipe_parameters(Config & cfg, int nY) {
    Parameters::lbm_base_parameters(cfg);
    re = cfg["domain"]["Re"].read<T>();
    
    pipe_radius = nY;
    hlog << "(Parameters) Your pipe has a given radius of " << pipe_radius << " LU, given so might be wrong" << std::endl;
    u_lbm_max = re * nu_lbm / (pipe_radius*2);
};

void Parameters::lbm_shear_parameters(Config & cfg,T nx) {
  Parameters::lbm_base_parameters(cfg);
  T shearrate_p = cfg["domain"]["shearrate"].read<T>();
  re = (nx* (shearrate_p * (nx*0.5))) / nu_p;
  u_lbm_max = re * nu_lbm / nx;  
  shearrate_lbm = shearrate_p*dt;
}

void Parameters::printParameters() {
  hlog << "(HemoCell) System parameters:" << std::endl;
  hlog << "\t dx: \t" << dx << std::endl;
  hlog << "\t dt: \t" << dt << std::endl;
  hlog << "\t dm: \t" << dm << std::endl;
  hlog << "\t dN: \t" << df << std::endl;
  hlog << "\t tau: \t" << tau << std::endl;
  hlog << "\t nu_lbm: \t" << nu_lbm << std::endl;
  hlog << "\t u_lb_max: \t" << u_lbm_max << std::endl;
#ifdef FORCE_LIMIT
  hlog << "\t f_limit: \t" << f_limit << std::endl;
#endif
}

T Parameters::dt = 0.0;
T Parameters::dx = 0.0;
T Parameters::dm = 0.0;
T Parameters::df = 0.0;
T Parameters::nu_p = 0.0;
T Parameters::rho_p = 0.0;
T Parameters::tau = 0.0;
T Parameters::re = 0.0;
T Parameters::nu_lbm = 0.0;
T Parameters::u_lbm_max = 0.0;
T Parameters::shearrate_lbm = 0.0;

T Parameters::kBT_lbm = 0.0;
T Parameters::kBT_p = 0.0;
T Parameters::ef_lbm = 0.0;
#ifdef FORCE_LIMIT
T Parameters::f_limit = 0.0;
T Parameters::pipe_radius = 0.0;
}
#endif
