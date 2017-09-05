#include "vWFModel.h"
//TODO Make all inner hemo::Array variables constant as well


vWFModel::vWFModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_, modelCfg_),
                  cellField(cellField_),
                  k_link( vWFModel::calculate_kLink(modelCfg_) ), 
                  k_bend( vWFModel::calculate_kBend(modelCfg_) ),
                  link_eq(vWFModel::calculate_link_eq(modelCfg_) ),
                  bend_eq(vWFModel::calculate_bend_eq(modelCfg_) ),
                  eta_m( vWFModel::calculate_etaM(modelCfg_) ),
                  eta_v( vWFModel::calculate_etaV(modelCfg_) )
{
  cellField.deleteIncomplete = false;
  cellField.numVertex = modelCfg_["MaxVertices"].read<unsigned int>();
};

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
      const hemo::Array<double,3> & p0 = cell[i]->position;
      const hemo::Array<double,3> & p1 = cell[i+1]->position;
      
      const hemo::Array<double,3> edge_vec = p1-p0;
      const double edge_length = norm(edge_vec);
      const hemo::Array<double,3> edge_uv = edge_vec/edge_length;
      
      const double edge_frac = (edge_length - link_eq) / link_eq;
      
#ifdef FORCE_LIMIT
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#else
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#endif
      const hemo::Array<double,3> force = edge_uv*edge_force_scalar;
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

double vWFModel::calculate_etaV(Config & cfg ){
  return cfg["MaterialModel"]["eta_v"].read<double>() * param::dx * param::dt / param::dm; //== dx^2/dN/dt
};

double vWFModel::calculate_etaM(Config & cfg ){
  return cfg["MaterialModel"]["eta_m"].read<double>() * param::dx / param::dt / param::df;
};

double vWFModel::calculate_kBend(Config & cfg ){
  return cfg["MaterialModel"]["kBend"].read<double>() * param::kBT_lbm;
};

double vWFModel::calculate_bend_eq(Config & cfg ) {
  return cfg["MaterialModel"]["bend_eq"].read<double>();
}

double vWFModel::calculate_kLink(Config & cfg){
  double kLink = cfg["MaterialModel"]["kLink"].read<double>();
  double persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  double plc = persistenceLengthFine/param::dx;
  kLink *= param::kBT_lbm/plc;
  return kLink;
};

double vWFModel::calculate_link_eq(Config & cfg) {
  return cfg["MaterialModel"]["link_eq"].read<double>() * (1e-6/param::dx);
}