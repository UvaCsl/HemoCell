#include "rbcHighOrderModel.h"
//TODO Make all inner array variables constant as well


RbcHighOrderModel::RbcHighOrderModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_),
                  cellField(cellField_),
                  k_volume( RbcHighOrderModel::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ),
                  k_area( RbcHighOrderModel::calculate_kArea(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( RbcHighOrderModel::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( RbcHighOrderModel::calculate_kBend(modelCfg_) ),
                  eta( RbcHighOrderModel::calculate_eta(modelCfg_) )
    {};

void RbcHighOrderModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> particles_per_cell, map<int,bool> lpc, pluint ctype) {

  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    if (cell[0]->celltype != ctype) continue; //only execute on correct particle

    //Calculate Cell Values that need all particles (but do it most efficient
    //tailored to this class)
    double volume = 0.0;
    int triangle_n = 0;
    vector<double> triangle_areas;
    vector<Array<double,3>> triangle_normals; 


    for (const Array<plint,3> & triangle : cellConstants.triangle_list) {
      const Array<double,3> & v0 = cell[triangle[0]]->position;
      const Array<double,3> & v1 = cell[triangle[1]]->position;
      const Array<double,3> & v2 = cell[triangle[2]]->position;
      
      //Volume
      const double v210 = v2[0]*v1[1]*v0[2];
      const double v120 = v1[0]*v2[1]*v0[2];
      const double v201 = v2[0]*v0[1]*v1[2];
      const double v021 = v0[0]*v2[1]*v1[2];
      const double v102 = v1[0]*v0[1]*v2[2];
      const double v012 = v0[0]*v1[1]*v2[2];
      volume += (1.0/6.0)*(-v210+v120+v201-v021-v102+v012);
      
      //Area
      double area; 
      Array<double,3> t_normal;
      computeTriangleAreaAndUnitNormal(v0, v1, v2, area, t_normal);

      const double areaRatio = (area - /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n])
                               / /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n];      
       
      //area force magnitude
      const double afm = k_area * (areaRatio+areaRatio/(0.04-areaRatio*areaRatio));

        Array<double,3> centroid;
        centroid[0] = (v0[0]+v1[0]+v2[0])/3.0;
        centroid[1] = (v0[1]+v1[1]+v2[1])/3.0;
        centroid[2] = (v0[2]+v1[2]+v2[2])/3.0;
        Array<double,3> av0 = centroid - v0;
        Array<double,3> av1 = centroid - v1;
        Array<double,3> av2 = centroid - v2;

      // if(afm > 0) { // We need to contract the face
        *cell[triangle[0]]->force_area += afm*av0;
        *cell[triangle[1]]->force_area += afm*av1;
        *cell[triangle[2]]->force_area += afm*av2;
      // }
      // else { // We need to expand it, so we reverse their length, the farthest vertex is forced the least
      //   // One possible expansion
      //   double l0 = norm(av0);
      //   double l1 = norm(av1);
      //   double l2 = norm(av2);
      //   double sum_l = l0+l1+l2;  
      //   *cell[triangle[0]]->force_area += afm / l0 * (sum_l-(l1+l2)) * av0;
      //   *cell[triangle[1]]->force_area += afm / l1 * (sum_l-(l0+l2)) * av1;
      //   *cell[triangle[2]]->force_area += afm / l2 * (sum_l-(l0+l1)) * av2;
        
        // Another version
        // double l0 = norm(v1-v2);
        // double l1 = norm(v2-v0);
        // double l2 = norm(v0-v1);
        // double sum_l = l0+l1+l2;
        // *cell[triangle[0]]->force_area += afm / sum_l * l0 * av0;
        // *cell[triangle[1]]->force_area += afm / sum_l * l1 * av1;
        // *cell[triangle[2]]->force_area += afm / sum_l * l2 * av2;

        // Yet another version
        // Array<double,3> av0;
        // crossProduct(v1-v2,t_normal,av0);
        // Array<double,3> av1;
        // crossProduct(v2-v0,t_normal,av1);
        // Array<double,3> av2;
        // crossProduct(v0-v1,t_normal,av2);
        // *cell[triangle[0]]->force_area -= afm*av1*0.5; //av1 edge force is divided over two neighbouring edges
        // *cell[triangle[0]]->force_area -= afm*av2*0.5;
        // *cell[triangle[1]]->force_area -= afm*av0*0.5;
        // *cell[triangle[1]]->force_area -= afm*av2*0.5;
        // *cell[triangle[2]]->force_area -= afm*av0*0.5;
        // *cell[triangle[2]]->force_area -= afm*av1*0.5; //av1 edge force is divided over two neighbouring edges
      //}


      //Store values necessary later
      triangle_areas.push_back(area);
      triangle_normals.push_back(t_normal);

      triangle_n++;
    }
    
    //Volume
    const double volume_frac = (volume-cellConstants.volume_eq)/cellConstants.volume_eq;
    const double volume_force = -k_volume * volume_frac/(0.01-volume_frac*volume_frac);

    triangle_n = 0;

    for (const Array<plint,3> & triangle : cellConstants.triangle_list) {
      //Fixed volume force per area
      const Array<double, 3> local_volume_force = (volume_force*triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      *cell[triangle[0]]->force_volume += local_volume_force;
      *cell[triangle[1]]->force_volume += local_volume_force;
      *cell[triangle[2]]->force_volume += local_volume_force;

      triangle_n++;
    }


    // Edges
    int edge_n=0;
    for (const Array<plint,2> & edge : cellConstants.edge_list) {
      const Array<double,3> & v0 = cell[edge[0]]->position;
      const Array<double,3> & v1 = cell[edge[1]]->position;

      // Link force
      const Array<double,3> edge_v = v1-v0;
      const double edge_length = norm(edge_v);
      const Array<double,3> edge_uv = edge_v/edge_length;
      const double edge_frac = (edge_length - /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n])
                               / /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n];
      
      //if (edge_frac > 0) {
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/(1.0-edge_frac*edge_frac));   // allows at max. 120% stretch
      const Array<double,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_link += force;
      *cell[edge[1]]->force_link -= force;
      //TODO: for proper contraction a different tau_link might be needed (e.g. > 1.0)
      //      instead of making it softer
      // } else{
      //    // less stiff compression resistance -> let compression be dominated
      //    // by area conservation force
      //   const double edge_force_scalar = k_link * edge_frac * edge_frac * edge_frac;
      //   const Array<double,3> force = edge_uv*edge_force_scalar;
      //   *cell[edge[0]]->force_link += force;
      //   *cell[edge[1]]->force_link -= force;
      // }

      // calculate triangle normals, this should be in a function

      const plint b0 = cellConstants.edge_bending_triangles_list[edge_n][0];
      const plint b1 = cellConstants.edge_bending_triangles_list[edge_n][1];

      const Array<double,3> b00 = particles_per_cell[cid][cellField.meshElement.getVertexId(b0,0)]->position;
      const Array<double,3> b01 = particles_per_cell[cid][cellField.meshElement.getVertexId(b0,1)]->position;
      const Array<double,3> b02 = particles_per_cell[cid][cellField.meshElement.getVertexId(b0,2)]->position;
      
      const Array<double,3> b10 = particles_per_cell[cid][cellField.meshElement.getVertexId(b1,0)]->position;
      const Array<double,3> b11 = particles_per_cell[cid][cellField.meshElement.getVertexId(b1,1)]->position;
      const Array<double,3> b12 = particles_per_cell[cid][cellField.meshElement.getVertexId(b1,2)]->position;

      const Array<double,3> V1 = plb::computeTriangleNormal(b00,b01,b02, false);
      const Array<double,3> V2 = plb::computeTriangleNormal(b10,b11,b12, false);

     
      const Array<double,3> x2 = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->position;

      //calculate angle
      double angle = angleBetweenVectors(V1, V2);
      const plint sign = dot(x2-v0, V2) >= 0 ? 1 : -1;
      if (sign <= 0) {
        angle = 2 * PI - angle;
      }
      if (angle > PI) {
        angle = angle - 2*PI; 
      }

      //calculate resulting bending force
      const double angle_frac = cellConstants.edge_angle_eq_list[edge_n]/*cellConstants.angle_mean_eq*/ - angle;

      const double force_magnitude = - k_bend * (angle_frac + angle_frac / ( 0.62 - (angle_frac * angle_frac)));

      //TODO bending force differs with area - That is intentional, and necessary!
      const Array<double,3> v1v2 = (V1 + V2)*0.5; 
      const Array<double,3> bending_force = force_magnitude*v1v2;
      *cell[edge[0]]->force_bending += bending_force;
      *cell[edge[1]]->force_bending += bending_force;
      //TODO Negate the force with 4 point bending
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_bending -= force_magnitude * V1;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_bending -= force_magnitude * V2;


      // Viscosity based on relative vertex velocity
      const Array<double,3> outer_end_rel_vel = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->v 
                                              - cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->v;
      const Array<double,3> section = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->position 
                                    - cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->position;
      const Array<double,3> section_dir = section / norm(section);
      const double norm_vel = dot(outer_end_rel_vel, section_dir); // relative velocity magn. between the points
      const double visc_mag = eta * norm_vel * 0.5;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_visc +=  visc_mag * section_dir;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_visc -=  visc_mag * section_dir;                            

      edge_n++;
    }

  } 
};

void RbcHighOrderModel::statistics() {
    pcout << "(Cell-mechanics model) High Order model parameters for " << cellField.name << " cellfield" << std::endl; 
    pcout << "\t k_link:   " << k_link << std::endl; 
    pcout << "\t k_area:   " << k_area << std::endl; 
    pcout << "\t k_bend: : " << k_bend << std::endl; 
    pcout << "\t k_volume: " << k_volume << std::endl;
    pcout << "\t eta:      " << eta << std::endl;
};


// Provide methods to calculate and scale to coefficients from here

double RbcHighOrderModel::calculate_eta(Config & cfg ){
  return cfg["MaterialModel"]["eta"].read<double>() * param::dx * param::dt / param::dm;
};

double RbcHighOrderModel::calculate_kBend(Config & cfg ){
  return cfg["MaterialModel"]["kBend"].read<double>() * param::kBT_lbm;
};

double RbcHighOrderModel::calculate_kVolume(Config & cfg, MeshMetrics<double> & meshmetric){
  double kVolume =  cfg["MaterialModel"]["kVolume"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kVolume *= param::kBT_lbm/(eqLength*eqLength*eqLength);
  //kVolume /= meshmetric.getNumVertices();
  return kVolume;
};

double RbcHighOrderModel::calculate_kArea(Config & cfg, MeshMetrics<double> & meshmetric){
  double kArea =  cfg["MaterialModel"]["kArea"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kArea *= param::kBT_lbm/(eqLength*eqLength);
  return kArea;
}

double RbcHighOrderModel::calculate_kLink(Config & cfg, MeshMetrics<double> & meshmetric){
  double kLink = cfg["MaterialModel"]["kLink"].read<double>();
  double persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  double plc = persistenceLengthFine/param::dx; //* sqrt((meshmetric.getNumVertices()-2.0) / (23867-2.0)); //Kaniadakis magic
  kLink *= param::kBT_lbm/plc;
  return kLink;
}
