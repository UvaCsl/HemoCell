#include "pltSimpleModel.h"

PltSimpleModel::PltSimpleModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_), 
                  cellField(cellField_),
                  k_volume( PltSimpleModel::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( PltSimpleModel::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( PltSimpleModel::calculate_kBend(modelCfg_) ),
                  eta( PltSimpleModel::calculate_eta(modelCfg_) )
  { };

void PltSimpleModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, pluint ctype) {
  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    if (cell[0]->celltype != ctype) continue; //only execute on correct particle

    //Calculate Cell Values that need all particles (but do it efficiently,
    //tailored to this class)
    double volume = 0.0;
    int triangle_n = 0;
    vector<double> triangle_areas;
    vector<hemo::Array<double,3>> triangle_normals; 

    // Per-triangle calculations
    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      const hemo::Array<double,3> & v0 = cell[triangle[0]]->position;
      const hemo::Array<double,3> & v1 = cell[triangle[1]]->position;
      const hemo::Array<double,3> & v2 = cell[triangle[2]]->position;
      
      //Volume
      const double v210 = v2[0]*v1[1]*v0[2];
      const double v120 = v1[0]*v2[1]*v0[2];
      const double v201 = v2[0]*v0[1]*v1[2];
      const double v021 = v0[0]*v2[1]*v1[2];
      const double v102 = v1[0]*v0[1]*v2[2];
      const double v012 = v0[0]*v1[1]*v2[2];
      volume += (-v210+v120+v201-v021-v102+v012);
      
      //Area
      double area; 
      hemo::Array<double,3> t_normal;
      computeTriangleAreaAndUnitNormal(v0, v1, v2, area, t_normal);
      
      //Store values necessary later
      triangle_areas.push_back(area);
      triangle_normals.push_back(t_normal);

      triangle_n++;
    }

    volume *= (1.0/6.0);

    //Volume
    const double volume_frac = (volume-cellConstants.volume_eq)/cellConstants.volume_eq;
#ifdef FORCE_LIMIT
    const double volume_force = -k_volume * volume_frac/std::fabs(0.01-volume_frac*volume_frac);
#else
    const double volume_force = -k_volume * volume_frac/(0.01-volume_frac*volume_frac);
#endif


    triangle_n = 0;

    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      //Fixed volume force per area
      const hemo::Array<double, 3> local_volume_force = (volume_force*triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      *cell[triangle[0]]->force_volume += local_volume_force;
      *cell[triangle[1]]->force_volume += local_volume_force;
      *cell[triangle[2]]->force_volume += local_volume_force;

      triangle_n++;
    }


//Vertex bending force loop
    for (long unsigned int i = 0 ; i < cell.size() ; i++) {
      Array<double,3> vertexes_sum = {0.,0.,0.};
      Array<double,3> vertices_vel_sum = {0.,0.,0.};

      const Array<plint,6> & edges = cellConstants.vertex_edges[i];
      unsigned int absent = 0;
      for (unsigned int j = 0 ; j < 6 ; j++ ) {
        if (edges[j] == -1) {
          absent++;
          continue;
        }
        vertexes_sum += cell[edges[j]]->position;
        vertices_vel_sum += cell[edges[j]]->v;
      }
      const Array<double,3> vertexes_middle = vertexes_sum/(6.0-absent);
      const Array<double,3> vertices_vavg = vertices_vel_sum/(6.0-absent);

      const Array<double,3> dev_vect = vertexes_middle - cell[i]->position;
      const double dev_n = norm(dev_vect); // absolute distance
      const Array<double, 3> dev_dir = dev_vect / dev_n;

      // Get which side is the vertex on (e.g. inward or outward curve)
      // We dont know how much vertex neighbours exist, but the firts two always has to be present, so get the normal approximation using those
      const Array<double,3> patch_norm_approx = crossProduct( (cell[edges[0]]->position - cell[i]->position),
                                                              (cell[edges[1]]->position - cell[i]->position) );
      const double sign = dot(patch_norm_approx, dev_vect);

      const double dDev = abs(copysign(dev_n, sign) - cellConstants.surface_patch_center_eq_list[i] ) / cellConstants.edge_mean_eq; // Non-dimension
  
      //TODO scale bending force
#ifdef FORCE_LIMIT
      const Array<double,3> bending_force = k_bend * ( dDev + dDev/std::fabs(0.5-dDev*dDev)) * dev_dir; // tau_b comes from the angle limit w. eq.lat.tri. assumptiln
#else
      const Array<double,3> bending_force = k_bend * ( dDev + dDev/std::fabs(0.5-dDev*dDev)) * dev_dir;
#endif      
      // Calculating viscous term
      //const Array<double,3> rel_vel_v = vertices_vavg - cell[i]->v;
      //const Array<double,3> rel_vel_proj = dot(dev_dir, rel_vel_v) * dev_dir;
      //const Array<double,3> Fvisc_vol = eta_v * rel_vel_proj * 0.866 * cellConstants.edge_mean_eq;

      //Apply bending force
      *cell[i]->force_bending += bending_force;
      //*cell[i]->force_visc += Fvisc_vol;
      const Array<double,3> negative_bending_force = -bending_force/(6.0-absent);
      for (unsigned int j = 0 ; j < 6 - absent; j++ ) {
       *cell[edges[j]]->force_bending += negative_bending_force;
      }              
    }

    // Per-edge calculations
    int edge_n=0;
    for (const hemo::Array<plint,2> & edge : cellConstants.edge_list) {
      const hemo::Array<double,3> & v0 = cell[edge[0]]->position;
      const hemo::Array<double,3> & v1 = cell[edge[1]]->position;

      // Link force
      const hemo::Array<double,3> edge_v = v1-v0;
      const double edge_length = sqrt(edge_v[0]*edge_v[0]+edge_v[1]*edge_v[1]+edge_v[2]*edge_v[2]);
      const hemo::Array<double,3> edge_uv = edge_v/edge_length;
      const double edge_frac = (edge_length-cellConstants.edge_length_eq_list[edge_n])/cellConstants.edge_length_eq_list[edge_n];

#ifdef FORCE_LIMIT
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#else
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch     
#endif
      const hemo::Array<double,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_link += force;
      *cell[edge[1]]->force_link -= force;

#if 0
      // calculate triangle normals, this should be in a function

      const plint b0 = cellConstants.edge_bending_triangles_list[edge_n][0];
      const plint b1 = cellConstants.edge_bending_triangles_list[edge_n][1];

      const hemo::Array<double,3> b00 = particles_per_cell[cid][cellField.meshElement.getVertexId(b0,0)]->position;
      const hemo::Array<double,3> b01 = particles_per_cell[cid][cellField.meshElement.getVertexId(b0,1)]->position;
      const hemo::Array<double,3> b02 = particles_per_cell[cid][cellField.meshElement.getVertexId(b0,2)]->position;
      
      const hemo::Array<double,3> b10 = particles_per_cell[cid][cellField.meshElement.getVertexId(b1,0)]->position;
      const hemo::Array<double,3> b11 = particles_per_cell[cid][cellField.meshElement.getVertexId(b1,1)]->position;
      const hemo::Array<double,3> b12 = particles_per_cell[cid][cellField.meshElement.getVertexId(b1,2)]->position;

      const hemo::Array<double,3> V1 = computeTriangleNormal(b00,b01,b02, false);
      const hemo::Array<double,3> V2 = computeTriangleNormal(b10,b11,b12, false);

     
      const hemo::Array<double,3> x2 = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->position;

      //calculate angle
      // double angle = angleBetweenVectors(V1, V2);
      // const plint sign = dot(x2-v0, V2) >= 0 ? 1 : -1;
      // if (sign <= 0) {
      //   angle = 2 * PI - angle;
      // }
      // if (angle > PI) {
      //   angle = angle - 2*PI; 
      // }

      double angle = getAngleBetweenFaces(V1, V2, edge_uv);
      
      //calculate resulting bending force
      const double angle_frac = cellConstants.edge_angle_eq_list[edge_n] - angle;

#ifdef FORCE_LIMIT
      const double force_magnitude = - k_bend * (angle_frac + angle_frac / std::fabs(2.467 - angle_frac * angle_frac) ); // tau_b = pi/2
#else
      const double force_magnitude = - k_bend * (angle_frac + angle_frac / (2.467 - angle_frac * angle_frac) ); // tau_b = pi/2
#endif
      //TODO Make bending force differ with area!
      const hemo::Array<double,3> bending_force = force_magnitude*(V1 + V2)*0.5;
      *cell[edge[0]]->force_bending += bending_force;
      *cell[edge[1]]->force_bending += bending_force;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_bending -= bending_force;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_bending -= bending_force;


      // Viscosity based on relative vertex velocity
      const hemo::Array<double,3> outer_end_rel_vel = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->v 
                                              - cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->v;
      const hemo::Array<double,3> section = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->position 
                                    - cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->position;
      const hemo::Array<double,3> section_dir = section / norm(section);
      const double norm_vel = dot(outer_end_rel_vel, section_dir); // relative velocity magn. between the points
      const double visc_mag = eta * norm_vel * 0.5;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_visc +=  visc_mag * section_dir;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_visc -=  visc_mag * section_dir;   
#endif
      edge_n++;
    }

  } 
}

void PltSimpleModel::statistics() {
    pcout << "(Cell-mechanics model) Reduced-model parameters for " << cellField.name << " cellfield" << std::endl;
    pcout << "\t k_link:   " << k_link << std::endl; 
    pcout << "\t k_bend: : " << k_bend << std::endl; 
    pcout << "\t k_volume: " << k_volume << std::endl; 
    pcout << "\t eta:      " << eta << std::endl;
};


double PltSimpleModel::calculate_eta(Config & cfg ){
  return cfg["MaterialModel"]["eta"].read<double>() * param::dx * param::dt / param::dm;
};

double PltSimpleModel::calculate_kBend(Config & cfg ){
  return cfg["MaterialModel"]["kBend"].read<double>() * param::kBT_lbm;
};

double PltSimpleModel::calculate_kVolume(Config & cfg, MeshMetrics<double> & meshmetric){
  double kVolume =  cfg["MaterialModel"]["kVolume"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kVolume *= param::kBT_lbm/(eqLength*eqLength*eqLength);
  return kVolume;
};

double PltSimpleModel::calculate_kLink(Config & cfg, MeshMetrics<double> & meshmetric){
  double kLink = cfg["MaterialModel"]["kLink"].read<double>();
  double persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  double plc = persistenceLengthFine/param::dx; //* sqrt((meshmetric.getNumVertices()-2.0) / (23867-2.0)); //Kaniadakis magic
  kLink *= param::kBT_lbm/plc;
  return kLink;
};