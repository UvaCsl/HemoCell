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
#include "rbcFedosovModel.h"
#include "logfile.h"

namespace hemo {

RbcFedosovModel::RbcFedosovModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_, modelCfg_),
                  cellField(cellField_),
                  eta_m( RbcFedosovModel::calculate_etaM(modelCfg_) ),
                  kd( RbcFedosovModel::calculate_kd(modelCfg_,*cellField_.meshmetric) ),
                  ka( RbcFedosovModel::calculate_ka(modelCfg_,*cellField_.meshmetric) ),
                  kv( RbcFedosovModel::calculate_kv(modelCfg_,*cellField_.meshmetric) ),
                  kb( RbcFedosovModel::calculate_kb(modelCfg_,*cellField_.meshmetric) ),
                  x0( RbcFedosovModel::calculate_x0(modelCfg_,*cellField_.meshmetric) ),
                  m( RbcFedosovModel::calculate_m(modelCfg_,*cellField_.meshmetric)) {};

void RbcFedosovModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, size_t ctype) {

  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.

    const int & cid = pair.first;

    vector<HemoCellParticle*> & cell = particles_per_cell[cid];

    if (cell.size() == 0) continue;
    if (cell[0]->sv.celltype != ctype) continue; //only execute on correct particle

    //Calculate Cell Values that need all particles
    T volume = 0.0;
    vector<T> triangle_areas;
    triangle_areas.reserve(cellConstants.triangle_list.size());
    vector<hemo::Array<T,3>> triangle_normals;
    triangle_normals.reserve(cellConstants.triangle_list.size());

    // Volume calculation
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
      volume += (-v210+v120+v201-v021-v102+v012);
    }
    volume *= (1.0/6.0);

    // Calculating the current surface area
    T surface = 0.0;
    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      const hemo::Array<T,3> & v0 = cell[triangle[0]]->sv.position;
      const hemo::Array<T,3> & v1 = cell[triangle[1]]->sv.position;
      const hemo::Array<T,3> & v2 = cell[triangle[2]]->sv.position;
      
      //Area
      T area_tmp = 0.0;
      area_tmp = computeTriangleArea(v0, v1, v2);
      surface += area_tmp;
    }

    // Per-triangle calculations
    int triangle_n = 0;
    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      const hemo::Array<T,3> & v0 = cell[triangle[0]]->sv.position;
      const hemo::Array<T,3> & v1 = cell[triangle[1]]->sv.position;
      const hemo::Array<T,3> & v2 = cell[triangle[2]]->sv.position;
      
      //Area
      T area; 
      hemo::Array<T,3> t_normal;
      computeTriangleAreaAndUnitNormal(v0, v1, v2, area, t_normal);

       // ======================================= AREA =======================================
      T alpha_a = -ka * (surface - cellConstants.surface_eq) / (4 * cellConstants.surface_eq * area);
      T alpha_d = -kd * (area - cellConstants.triangle_area_eq_list[triangle_n]) / (4 * cellConstants.triangle_area_eq_list[triangle_n] * area);

      hemo::Array<T,3> v10 = v1 - v0;
      hemo::Array<T,3> v02 = v0 - v2;
      hemo::Array<T,3> v21 = v2 - v1;

      hemo::Array<T,3> si_k = crossProduct(v10, -v02);

      hemo::Array<T,3> force0localGlobal = (alpha_d + alpha_a) * crossProduct(si_k, v21);
      hemo::Array<T,3> force1localGlobal = (alpha_d + alpha_a) * crossProduct(si_k, v02);
      hemo::Array<T,3> force2localGlobal = (alpha_d + alpha_a) * crossProduct(si_k, v10);

      *cell[triangle[0]]->force_area += force0localGlobal;
      *cell[triangle[1]]->force_area += force1localGlobal;
      *cell[triangle[2]]->force_area += force2localGlobal;

      //Store values necessary later
      triangle_areas.push_back(area);
      triangle_normals.push_back(t_normal);

      // ======================================= VOLUME ======================================= 
      // Calculate centroid of the triangle
      hemo::Array<T,3> tc_k = (v0 + v1 + v2) / 3.0;

      T beta_v = -kv * (volume - cellConstants.volume_eq) / cellConstants.volume_eq;
      hemo::Array<T,3> force0v = beta_v * (si_k / 3.0 + crossProduct(tc_k, v21));
      hemo::Array<T,3> force1v = beta_v * (si_k / 3.0 + crossProduct(tc_k, v02));
      hemo::Array<T,3> force2v = beta_v * (si_k / 3.0 + crossProduct(tc_k, v10));

      *cell[triangle[0]]->force_volume += force0v;
      *cell[triangle[1]]->force_volume += force1v;
      *cell[triangle[2]]->force_volume += force2v;

#ifdef INTERIOR_VISCOSITY
      // Add the normal direction here, always pointing outward
      // const hemo::Array<T, 3> local_normal_dir = (triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      const hemo::Array<T, 3> local_normal_dir = t_normal * (area / cellConstants.area_mean_eq);
      cell[triangle[0]]->normalDirection += local_normal_dir;
      cell[triangle[1]]->normalDirection += local_normal_dir;
      cell[triangle[2]]->normalDirection += local_normal_dir;
#endif

      triangle_n++;
    }

    // Per-edge calculations ===================== LINK ==================================================

    int edge_n=0;
    for (const hemo::Array<plint,2> & edge : cellConstants.edge_list) {

      hemo::Array<T,3> & p0 = cell[edge[0]]->sv.position;
      hemo::Array<T,3> & p1 = cell[edge[1]]->sv.position;

      hemo::Array<T,3> v32 = p1 - p0;
      T v32_mag = norm(v32);
      hemo::Array<T,3> v32_norm = v32 / v32_mag;

      // Fedosov model
      T edge_frac = (v32_mag - cellConstants.edge_length_eq_list[edge_n]) / cellConstants.edge_length_eq_list[edge_n];
      T x_frac = x0 * (edge_frac + 1);

      // Cutting off x_frac to avoid numerical problems
      if (x_frac > 0.999) {
          x_frac = 0.99;
      }

      T edge_force_WLC_mag = cellConstants.edge_kLinkWLC_list[edge_n] * (1.0 / (4.0 * (1.0 - x_frac) * (1.0 - x_frac)) - 1.0 / 4.0 + x_frac);
      T edge_force_POW_mag = cellConstants.edge_kp_list[edge_n] / std::pow(x_frac * cellConstants.edge_length_max_list[edge_n], m);
      T edge_force_mag = edge_force_WLC_mag - edge_force_POW_mag;
      hemo::Array<T,3> edge_force = v32_norm * edge_force_mag;
      *cell[edge[0]]->force_link += edge_force;
      *cell[edge[1]]->force_link -= edge_force;

      // ======================================= Viscosity =======================================
      
      if (eta_m != 0.0) {
        // Membrane viscosity of bilipid layer
        // F = eta * (dv/l) * l. 
        const hemo::Array<T,3> rel_vel = cell[edge[1]]->sv.v - cell[edge[0]]->sv.v;
        const hemo::Array<T,3> rel_vel_projection = dot(rel_vel, v32_norm) * v32_norm;
        hemo::Array<T,3> Fvisc_memb = eta_m * rel_vel_projection;

        // Limit membrane viscosity
        const T Fvisc_memb_mag = norm(Fvisc_memb);
        if (Fvisc_memb_mag > FORCE_LIMIT / 4.0) {
          Fvisc_memb *= (FORCE_LIMIT / 4.0) / Fvisc_memb_mag;
        }

        *cell[edge[0]]->force_visc += Fvisc_memb;
        *cell[edge[1]]->force_visc -= Fvisc_memb; 
      }

      // ======================================= Bending =======================================
      hemo::Array<T,3> v1 = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->sv.position;
      hemo::Array<T,3> v2 = p0;
      hemo::Array<T,3> v3 = p1;
      hemo::Array<T,3> v4 = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->sv.position;

      // hemo::Array<T,3> v32 = v3 - v2;  // Already calculated above
      hemo::Array<T,3> v13 = v1 - v3;
      hemo::Array<T,3> v34 = v3 - v4;
      hemo::Array<T,3> v42 = v4 - v2;
      hemo::Array<T,3> v21 = v2 - v1;

      // to be used for theta sign determination
      hemo::Array<T,3> v14 = v1 - v4;

      hemo::Array<T,3> si = crossProduct(v21, -v13);
      hemo::Array<T,3> zeta = crossProduct(v34, -v42);

      T si_mag = norm(si);
      T zeta_mag = norm(zeta);

      hemo::Array<T,3> si_norm = si / si_mag;
      hemo::Array<T,3> zeta_norm = zeta / zeta_mag;

      // USERMESO implementation for stability
      T cosine_theta = dot(si_norm, zeta_norm);
      if (cosine_theta > 1.0) cosine_theta = 1.0;
      if (cosine_theta < -1.0) cosine_theta = -1.0;
      T sine_theta = std::sqrt(1 - cosine_theta * cosine_theta);
      if (sine_theta < T(1.0e-6)) sine_theta = 1.0e-6;
      T mx = dot((si_norm - zeta_norm), v14);
      if (mx < 0.0) sine_theta = -sine_theta;
      T theta = std::atan2(sine_theta, cosine_theta);

      // Fedosov method
      T beta = kb * (sine_theta * cellConstants.edge_cos_angle_eq_list[edge_n] - cosine_theta * cellConstants.edge_sin_angle_eq_list[edge_n]) / 
                      sine_theta;
      T dTheta = theta - cellConstants.edge_angle_eq_list[edge_n];

      T b11 = -beta * cosine_theta / (si_mag * si_mag);
      T b12 = beta / (si_mag * zeta_mag);
      T b22 = -beta * cosine_theta / (zeta_mag * zeta_mag);

      hemo::Array<T,3> force1 = b11 * crossProduct(si, v32) + b12 * crossProduct(zeta, v32);
      hemo::Array<T,3> force2 = b11 * crossProduct(si, v13) + b12 * (crossProduct(si, v34) + crossProduct(zeta, v13)) + b22 * crossProduct(zeta, v34);
      hemo::Array<T,3> force3 = b11 * crossProduct(si, v21) + b12 * (crossProduct(si, v42) + crossProduct(zeta, v21)) + b22 * crossProduct(zeta, v42);
      hemo::Array<T,3> force4 = b12 * crossProduct(si, -v32) + b22 * crossProduct(zeta, -v32);

      // Apply forces to the vertices
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_bending += force1;
      *cell[edge[0]]->force_bending += force2;
      *cell[edge[1]]->force_bending += force3;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_bending += force4;

      edge_n++;
    }
  } 
};

void RbcFedosovModel::statistics() {
    hlog << "(Cell-mechanics model) Fedosov model parameters for " << cellField.name << " cellfield" << std::endl; 
    hlog << "\t eta_m:    " << eta_m << std::endl;
    hlog << "\t kd:   " << kd << std::endl; 
    hlog << "\t ka:   " << ka << std::endl; 
    hlog << "\t kv:   " << kv << std::endl; 
    hlog << "\t kb:   " << kb << std::endl; 
    hlog << "\t mean_edge:" << cellConstants.edge_mean_eq << std::endl;
    hlog << "\t N faces:  " << cellConstants.triangle_list.size() << std::endl;
};

}
