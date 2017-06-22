#include "rbcHighOrderModel.h"
//TODO Make all inner array variables constant as well


RbcHighOrderModel::RbcHighOrderModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_),
                  cellField(cellField_),
                  k_volume( RbcHighOrderModel::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ),
                  k_area( RbcHighOrderModel::calculate_kArea(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( RbcHighOrderModel::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( RbcHighOrderModel::calculate_kBend(modelCfg_) ),
                  eta_m( RbcHighOrderModel::calculate_etaM(modelCfg_) ),
                  eta_v( RbcHighOrderModel::calculate_etaV(modelCfg_) )
    {};

void RbcHighOrderModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, pluint ctype) {

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

    // Per-triangle calculations
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
      volume += (-v210+v120+v201-v021-v102+v012); // the factor of 1/6 moved to after the summation -> saves a few flops
      
      //Area
      double area; 
      Array<double,3> t_normal;
      computeTriangleAreaAndUnitNormal(v0, v1, v2, area, t_normal);

      const double areaRatio = (area - /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n])
                               / /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n];      
       
      //area force magnitude
#ifdef FORCE_LIMIT
      const double afm = k_area * (areaRatio+areaRatio/std::fabs(0.09-areaRatio*areaRatio));
#else
      const double afm = k_area * (areaRatio+areaRatio/(0.09-areaRatio*areaRatio));
#endif

      Array<double,3> centroid;
      centroid[0] = (v0[0]+v1[0]+v2[0])/3.0;
      centroid[1] = (v0[1]+v1[1]+v2[1])/3.0;
      centroid[2] = (v0[2]+v1[2]+v2[2])/3.0;
      Array<double,3> av0 = centroid - v0;
      Array<double,3> av1 = centroid - v1;
      Array<double,3> av2 = centroid - v2;

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
    const double volume_frac = (volume-cellConstants.volume_eq)/cellConstants.volume_eq;
#ifdef FORCE_LIMIT
    const double volume_force = -k_volume * volume_frac/std::fabs(0.01-volume_frac*volume_frac);
#else
    const double volume_force = -k_volume * volume_frac/(0.01-volume_frac*volume_frac);    
#endif
    triangle_n = 0;

    for (const Array<plint,3> & triangle : cellConstants.triangle_list) {
      // Scale volume force with local face area
      const Array<double, 3> local_volume_force = (volume_force*triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      *cell[triangle[0]]->force_volume += local_volume_force;
      *cell[triangle[1]]->force_volume += local_volume_force;
      *cell[triangle[2]]->force_volume += local_volume_force;

      triangle_n++;
    }


    // Per-edge calculations
    int edge_n=0;
    for (const Array<plint,2> & edge : cellConstants.edge_list) {
      const Array<double,3> & p0 = cell[edge[0]]->position;
      const Array<double,3> & p1 = cell[edge[1]]->position;

      // Link force
      const Array<double,3> edge_vec = p1-p0;
      const double edge_length = norm(edge_vec);
      const Array<double,3> edge_uv = edge_vec/edge_length;
      const double edge_frac = (edge_length - /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n])
                               / /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n];

#ifdef FORCE_LIMIT
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#else
      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
#endif
      const Array<double,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_link += force;
      *cell[edge[1]]->force_link -= force;

      // Membrane viscosity of bilipid layer
      // F = eta * (dv/l) * l. 
      const Array<double,3> rel_vel = cell[edge[1]]->v - cell[edge[0]]->v;
      const Array<double,3> rel_vel_projection = dot(rel_vel, edge_uv) * edge_uv;
      const Array<double,3> Fvisc_memb = eta_m * rel_vel_projection;
      *cell[edge[0]]->force_visc += Fvisc_memb;
      *cell[edge[1]]->force_visc -= Fvisc_memb; 


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


      double angle = getAngleBetweenFaces(V1, V2, edge_uv);

      //calculate resulting bending force
      const double angle_frac = cellConstants.edge_angle_eq_list[edge_n]/*cellConstants.angle_mean_eq*/ - angle;
#ifdef FORCE_LIMIT
      const double force_magnitude = - k_bend * (angle_frac + angle_frac / std::fabs(2.467 - angle_frac * angle_frac)); // tau_b = pi/2
#else
      const double force_magnitude = - k_bend * (angle_frac + angle_frac / (2.467 - angle_frac * angle_frac)); // tau_b = pi/2      
#endif

      //TODO make bending force differ with area, V1 and V2 are unit vectors right now!
      const Array<double,3> v1v2 = (V1 + V2)*0.5; 
      const Array<double,3> bending_force = force_magnitude*v1v2;
      *cell[edge[0]]->force_bending += bending_force;
      *cell[edge[1]]->force_bending += bending_force;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_bending -= force_magnitude * V1;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_bending -= force_magnitude * V2;


      // Volume viscosity of cytoplasm based on relative outer vertex velocity
      // F = eta * (dv/2l) * area. | area = sqrt(3)*l^2/4 => F = eta * dv * sqrt(3)/8 * l
      const Array<double,3> outer_end_rel_vel = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->v 
                                              - cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->v;
      const Array<double,3> section = cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->position 
                                    - cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->position;
      const Array<double,3> section_dir = section / norm(section);
      const Array<double,3> vel_projection = dot(outer_end_rel_vel, section_dir) * section_dir; // relative velocity magn. between the points projected on the line between them
      const Array<double,3> vel_rejection = outer_end_rel_vel - vel_projection;
      const Array<double,3> Fvisc_vol = eta_v * vel_rejection * 0.2165 * cellConstants.edge_mean_eq; // assuming similar triangle areas
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_visc +=  Fvisc_vol * 0.5 ;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_visc -=  Fvisc_vol * 0.5;                            

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
    pcout << "\t eta_m:    " << eta_m << std::endl;
    pcout << "\t eta_v:    " << eta_v << std::endl;
};


// Provide methods to calculate and scale to coefficients from here

double RbcHighOrderModel::calculate_etaV(Config & cfg ){
  return cfg["MaterialModel"]["eta_v"].read<double>() * param::dx * param::dt / param::dm; //== dx^2/dN/dt
};

double RbcHighOrderModel::calculate_etaM(Config & cfg ){
  return cfg["MaterialModel"]["eta_m"].read<double>() * param::dx / param::dt / param::df;
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
};

double RbcHighOrderModel::calculate_kLink(Config & cfg, MeshMetrics<double> & meshmetric){
  double kLink = cfg["MaterialModel"]["kLink"].read<double>();
  double persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  double plc = persistenceLengthFine/param::dx;
  kLink *= param::kBT_lbm/plc;
  return kLink;
};