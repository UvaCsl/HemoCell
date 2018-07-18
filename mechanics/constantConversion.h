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
#ifndef HEMO_constantConversion_h
#define HEMO_constantConversion_h
//This class converses SI things to LBM things for classes like the HO model.
#include "config.h"
#include "constant_defaults.h"
#include "multiBlock/multiDataField3D.hh"

namespace hemo {
  class Parameters {

    //Variables for the simulation
    public:
    static T dt,dx,dm,df;
    static T nu_p,rho_p,kBT_p;
    static T tau,re;
    static T nu_lbm, u_lbm_max, kBT_lbm, shearrate_lbm;
  #ifdef FORCE_LIMIT
    static T f_limit;
  #endif
    static T ef_lbm; //used in cellStretching
    static T pipe_radius; //used for bodyforce calculation



    static void lbm_pipe_parameters(Config & cfg, plb::MultiScalarField3D<int> * sf);
    static void lbm_pipe_parameters(Config & cfg, int nY);
    static void lbm_shear_parameters(Config & cfg, T nx);
    static void lbm_base_parameters(Config & cfg);
    static void printParameters();
  };
}
#endif
