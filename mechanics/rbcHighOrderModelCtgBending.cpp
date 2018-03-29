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
#include "rbcHighOrderModelCtgBending.h"
//TODO Make all inner hemo::Array variables constant as well


RbcHighOrderModelnewBending::RbcHighOrderModelnewBending(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_, modelCfg_),
                  cellField(cellField_),
                  k_volume( RbcHighOrderModelnewBending::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ),
                  k_area( RbcHighOrderModelnewBending::calculate_kArea(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( RbcHighOrderModelnewBending::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( RbcHighOrderModelnewBending::calculate_kBend(modelCfg_,*cellField_.meshmetric) ),
                  eta_m( RbcHighOrderModelnewBending::calculate_etaM(modelCfg_) ),
                  eta_v( RbcHighOrderModelnewBending::calculate_etaV(modelCfg_) )
    {};

void RbcHighOrderModelnewBending::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, size_t ctype) {

  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    if (cell[0]->sv.celltype != ctype) continue; //only execute on correct particle

    //Calculate Cell Values that need all particles (but do it most efficient
    //tailored to this class)
    T volume = 0.0;
    int triangle_n = 0;
    vector<T> triangle_areas;
    triangle_areas.reserve(cellConstants.triangle_list.size());
    vector<hemo::Array<T,3>> triangle_normals;
    triangle_normals.reserve(cellConstants.triangle_list.size());
    vector<T> edges_length;
    edges_length.reserve(cellConstants.edge_list.size());
    vector<hemo::Array<T,3>> edges_vector;
    edges_vector.reserve(cellConstants.edge_list.size());

    // Per-triangle calculations
    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      const hemo::Array<T,3> & v0 = cell[triangle[0]]->sv.position;
      const hemo::Array<T,3> & v1 = cell[triangle[1]]->sv.position;
      const hemo::Array<T,3> & v2 = cell[triangle[2]]->sv.position;
      
      //Volume
      const T v210 = v2[0]*v1[1]*v0[2];
      const T v120 = v1[0]*v2[1]*v0[2];
      const T v201 = v2[0]*v0[1]*v1[2];
      const T v021 = v0[0]*v2[1]*v1[2];
      const T v102 = v1[0]*v0[1]*v2[2];
      const T v012 = v0[0]*v1[1]*v2[2];
      volume += (-v210+v120+v201-v021-v102+v012); // the factor of 1/6 moved to after the summation -> saves a few flops
      
      //Area
      T area; 
      hemo::Array<T,3> t_normal;
      computeTriangleAreaAndUnitNormal(v0, v1, v2, area, t_normal);

      const T areaRatio = (area - /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n])
                               / /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n];      
       
      //area force magnitude
      const T afm = k_area * (areaRatio+areaRatio/std::fabs(0.09-areaRatio*areaRatio));

      hemo::Array<T,3> centroid;
      centroid[0] = (v0[0]+v1[0]+v2[0])/3.0;
      centroid[1] = (v0[1]+v1[1]+v2[1])/3.0;
      centroid[2] = (v0[2]+v1[2]+v2[2])/3.0;
      hemo::Array<T,3> av0 = centroid - v0;
      hemo::Array<T,3> av1 = centroid - v1;
      hemo::Array<T,3> av2 = centroid - v2;

      *cell[triangle[0]]->force_area += afm*av0;
      *cell[triangle[1]]->force_area += afm*av1;
      *cell[triangle[2]]->force_area += afm*av2;

      //Store values necessary later
      triangle_areas.push_back(area);
      triangle_normals.push_back(t_normal);

      triangle_n++;
    }
    
    volume *= (1.0/6.0);

    //Volume
    const T volume_frac = (volume-cellConstants.volume_eq)/cellConstants.volume_eq;
    const T volume_force = -k_volume * volume_frac/std::fabs(0.01-volume_frac*volume_frac);
    triangle_n = 0;

//Volume force loop
    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      // Scale volume force with local face area
      const hemo::Array<T, 3> local_volume_force = (volume_force*triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      *cell[triangle[0]]->force_volume += local_volume_force;
      *cell[triangle[1]]->force_volume += local_volume_force;
      *cell[triangle[2]]->force_volume += local_volume_force;

      triangle_n++;
    }
    
    // Per-edge calculations
    int edge_n=0;
    for (const hemo::Array<plint,2> & edge : cellConstants.edge_list) {
      const hemo::Array<T,3> & p0 = cell[edge[0]]->sv.position;
      const hemo::Array<T,3> & p1 = cell[edge[1]]->sv.position;

      // Link force
      const hemo::Array<T,3> edge_vec = p1-p0;
      const T edge_length = norm(edge_vec);
      edges_vector.push_back(edge_vec);
      edges_length.push_back(edge_length);
      
      const hemo::Array<T,3> edge_uv = edge_vec/edge_length;
      const T edge_frac = (edge_length - /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n])
                               / /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n];

      const T edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
      const hemo::Array<T,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_link += force;
      *cell[edge[1]]->force_link -= force;

      // Membrane viscosity of bilipid layer
      // F = eta * (dv/l) * l. 
      const hemo::Array<T,3> rel_vel = cell[edge[1]]->sv.v - cell[edge[0]]->sv.v;
      const hemo::Array<T,3> rel_vel_projection = dot(rel_vel, edge_uv) * edge_uv;
      const hemo::Array<T,3> Fvisc_memb = eta_m * rel_vel_projection;
      *cell[edge[0]]->force_visc += Fvisc_memb;
      *cell[edge[1]]->force_visc -= Fvisc_memb; 

      edge_n++;
    }
    
    //Bending force calculation (Method B from "On the bending algorithms for soft objects in flows" by Achim Guckenberger)

    for (unsigned int vertex = 0 ; vertex < cellConstants.vertex_vertexes.size() ; vertex++) {
      const hemo::Array<plint,6> & neighbour_edges = cellConstants.vertex_edges[vertex];
      const hemo::Array<signed int, 6> neighbour_edges_sign = cellConstants.vertex_edges_sign[vertex];
      const unsigned int & n_neighbours = cellConstants.vertex_n_vertexes[vertex];
      const hemo::Array<hemo::Array<plint,2>,6> & outer_edges_per_neighbour = cellConstants.vertex_outer_edges_per_vertex[vertex];
      const hemo::Array<hemo::Array<signed int,2>,6> & outer_edges_per_neighbour_sign = cellConstants.vertex_outer_edges_per_vertex_sign[vertex];
      
      hemo::Array<T,3> sum_vec_Tij = {0.,0.,0.};
      hemo::Array<hemo::Array<T,3>,3> sum_second_term = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
      hemo::Array<T,3> sum_fourth_term = {0.,0.,0.};

      hemo::Array<hemo::Array<T,3>,3> unit_vectors = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};

      T sum_len_Tij = 0.0;

      for (unsigned int n = 0 ; n < n_neighbours ; n++) {
        const unsigned int n_lower = ((n - 1)+n_neighbours)%n_neighbours;
        const unsigned int n_upper = (n+1)%n_neighbours;
        const hemo::Array<T,3>  n_edge = edges_vector[neighbour_edges[n]]*neighbour_edges_sign[n]*-1.0;
        const hemo::Array<T,3>  n_lower_edge = edges_vector[neighbour_edges[n_lower]]*neighbour_edges_sign[n_lower]*-1.0;
        const hemo::Array<T,3>  n_lower_outer_edge = edges_vector[outer_edges_per_neighbour[n][0]]*outer_edges_per_neighbour_sign[n][0]*-1.0;
        const hemo::Array<T,3>  n_upper_edge = edges_vector[neighbour_edges[n_upper]]*neighbour_edges_sign[n_upper]*-1.0;
        const hemo::Array<T,3>  n_upper_outer_edge = edges_vector[outer_edges_per_neighbour[n][1]]*outer_edges_per_neighbour_sign[n][1]*-1.0;
        const T n_edge_len = norm(n_edge);
        const T n_upper_edge_len = norm(n_upper_edge);
        const T n_lower_edge_len = norm(n_lower_edge);


        T n_lower_adj, n_lower_hyp, n_lower_opp, 
                     n_upper_adj, n_upper_hyp, n_upper_opp;

        computeLengthsPythagoras(n_lower_edge,n_lower_outer_edge,n_lower_adj, n_lower_opp, n_lower_hyp);
        computeLengthsPythagoras(n_upper_edge,n_upper_outer_edge,n_upper_adj, n_upper_opp, n_upper_hyp);


        const T Tij = n_lower_adj/n_lower_opp + 
                           n_upper_adj/n_upper_opp;

        sum_vec_Tij += n_edge*Tij;
        sum_len_Tij += pow(n_edge_len,2)*Tij;

        const T cos_lower = n_lower_adj/n_lower_hyp;
        const T cos_upper = n_upper_adj/n_upper_hyp;

        const hemo::Array<T,3> delta_cos_x_lower = (1.0/(n_edge_len*n_lower_edge_len))*
        (n_lower_edge*-1.0+n_edge*-1.0-(n_lower_edge_len/n_edge_len)*cos_lower*n_edge*-1.0-
        (n_edge_len/n_lower_edge_len)*cos_lower*n_lower_edge*-1.0
        );

        const hemo::Array<T,3> delta_cos_x_upper = (1.0/(n_edge_len*n_upper_edge_len))*
        (n_upper_edge*-1.0+n_edge*-1.0-(n_upper_edge_len/n_edge_len)*cos_upper*n_edge*-1.0-
        (n_edge_len/n_upper_edge_len)*cos_upper*n_upper_edge*-1.0
        );

        const hemo::Array<T,3> delta_Tij = (1.0/pow(1.0-pow(cos_lower,2),1.5))*delta_cos_x_lower +
                                                (1.0/pow(1.0-pow(cos_upper,2),1.5))*delta_cos_x_upper;

        for (unsigned int k = 0 ; k < 3 ; k++) {
          sum_second_term[k] += n_edge*delta_Tij[k] + Tij*unit_vectors[k];
          sum_fourth_term[k] += Tij*2*n_edge[k] + pow(n_edge_len,2)*delta_Tij[k];
        }
        
      }
      hemo::Array<T,3> force_bend;
      for (unsigned int k = 0 ; k < 3 ; k++) {
        force_bend[k]  = (dot(4.0*sum_vec_Tij/sum_len_Tij,sum_second_term[k])-(2.0*dot(sum_vec_Tij,sum_vec_Tij)/pow(sum_len_Tij,2))*sum_fourth_term[k]) * k_bend;
      }
      
      *cell[vertex]->force_bending -= force_bend;

     // const hemo::Array<T,3> force_bend_neg = force_bend/cellConstants.vertex_n_vertexes[vertex];
      
    //  for (unsigned int j = 0 ; j < cellConstants.vertex_n_vertexes[vertex]; j++ ) {
     //  *cell[cellConstants.vertex_vertexes[vertex][j]]->force_bending += force_bend_neg;
    //  }              
    }
  }
};

void RbcHighOrderModelnewBending::statistics() {
    pcout << "(Cell-mechanics model) High Order model parameters for " << cellField.name << " cellfield" << std::endl; 
    pcout << "\t k_link:   " << k_link << std::endl; 
    pcout << "\t k_area:   " << k_area << std::endl; 
    pcout << "\t k_bend: : " << k_bend << std::endl; 
    pcout << "\t k_volume: " << k_volume << std::endl;
    pcout << "\t eta_m:    " << eta_m << std::endl;
    pcout << "\t eta_v:    " << eta_v << std::endl;
};


// Provide methods to calculate and scale to coefficients from here

T RbcHighOrderModelnewBending::calculate_etaV(Config & cfg ){
  return cfg["MaterialModel"]["eta_v"].read<T>() * param::dx * param::dt / param::dm; //== dx^2/dN/dt
};

T RbcHighOrderModelnewBending::calculate_etaM(Config & cfg ){
  return cfg["MaterialModel"]["eta_m"].read<T>() * param::dx / param::dt / param::df;
};

T RbcHighOrderModelnewBending::calculate_kBend(Config & cfg, MeshMetrics<T> & meshmetric ){
  return cfg["MaterialModel"]["kBend"].read<T>(); //Not needed because of things // * param::kBT_lbm / eqLength;
};

T RbcHighOrderModelnewBending::calculate_kVolume(Config & cfg, MeshMetrics<T> & meshmetric){
  T kVolume =  cfg["MaterialModel"]["kVolume"].read<T>();
  T eqLength = meshmetric.getMeanLength();
  kVolume *= param::kBT_lbm/(eqLength*eqLength*eqLength);
  //kVolume /= meshmetric.getNumVertices();
  return kVolume;
};

T RbcHighOrderModelnewBending::calculate_kArea(Config & cfg, MeshMetrics<T> & meshmetric){
  T kArea =  cfg["MaterialModel"]["kArea"].read<T>();
  T eqLength = meshmetric.getMeanLength();
  kArea *= param::kBT_lbm/(eqLength*eqLength);
  return kArea;
};

T RbcHighOrderModelnewBending::calculate_kLink(Config & cfg, MeshMetrics<T> & meshmetric){
  T kLink = cfg["MaterialModel"]["kLink"].read<T>();
  T persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  T plc = persistenceLengthFine/param::dx;
  kLink *= param::kBT_lbm/plc;
  return kLink;
};
