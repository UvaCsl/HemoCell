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
#include "vWFModel.h"
//TODO Make all inner hemo::Array variables constant as well

void vWFModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, pluint ctype) {

  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    if (cell[0]->celltype != ctype) continue; //only execute on correct particle

    //Since this is a very special class we only need the separate vertexes, no volume force etc
    for (unsigned int i = 0 ; i < cell.size() -1; i++) {
      if (cell[i] == NULL) { continue; }
      if (cell[i+1] == NULL) { i++; continue; }

      //Calculate stretching force here if applicable
      const hemo::Array<T,3> & p0 = cell[i]->position;
      const hemo::Array<T,3> & p1 = cell[i+1]->position;
      
      const hemo::Array<T,3> edge_vec = p1-p0;
      const T edge_length = norm(edge_vec);
      const hemo::Array<T,3> edge_uv = edge_vec/edge_length;
      
      const T edge_frac = (edge_length - link_eq) / link_eq;
      
#ifdef FORCE_LIMIT
      const T edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#else
      const T edge_force_scalar = k_link * ( edge_frac + edge_frac/(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#endif
      const hemo::Array<T,3> force = edge_uv*edge_force_scalar;
      *cell[i]->force_link += force;
      *cell[i+1]->force_link -= force;
    }
  }
};

void vWFModel::statistics() {
    pcout << "(Cell-mechanics model) High Order model parameters for " << cellField.name << " cellfield" << std::endl; 
    pcout << "\t k_link:   " << k_link << endl; 
    pcout << "\t k_bend:   " << k_bend << endl; 
    pcout << "\t bend_eq:  " << bend_eq << endl;
    pcout << "\t link_eq:  " << link_eq << endl;
    pcout << "\t eta_m:    " << eta_m << endl;
    pcout << "\t eta_v:    " << eta_v << endl;
};


// Provide methods to calculate and scale to coefficients from here

T vWFModel::calculate_etaV(Config & cfg ){
  return cfg["MaterialModel"]["eta_v"].read<T>() * param::dx * param::dt / param::dm; //== dx^2/dN/dt
};

T vWFModel::calculate_etaM(Config & cfg ){
  return cfg["MaterialModel"]["eta_m"].read<T>() * param::dx / param::dt / param::df;
};

T vWFModel::calculate_kBend(Config & cfg ){
  return cfg["MaterialModel"]["kBend"].read<T>() * param::kBT_lbm;
};

T vWFModel::calculate_bend_eq(Config & cfg ) {
  return cfg["MaterialModel"]["bend_eq"].read<T>();
}

T vWFModel::calculate_kLink(Config & cfg){
  T kLink = cfg["MaterialModel"]["kLink"].read<T>();
  T persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  T plc = persistenceLengthFine/param::dx;
  kLink *= param::kBT_lbm/plc;
  return kLink;
};

T vWFModel::calculate_link_eq(Config & cfg) {
  return cfg["MaterialModel"]["link_eq"].read<T>() * (1e-6/param::dx);
}
