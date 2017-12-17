#include "wbcHighOrderModel.h"
//TODO Make all inner hemo::Array variables constant as well


WbcHighOrderModel::WbcHighOrderModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_, modelCfg_),
                  cellField(cellField_),
                  k_volume( WbcHighOrderModel::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ),
                  k_area( WbcHighOrderModel::calculate_kArea(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( WbcHighOrderModel::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( WbcHighOrderModel::calculate_kBend(modelCfg_,*cellField_.meshmetric) ),
                  eta_m( WbcHighOrderModel::calculate_etaM(modelCfg_) ),
                  eta_v( WbcHighOrderModel::calculate_etaV(modelCfg_) ),
                  k_inner_rigid( WbcHighOrderModel::calculate_kInnerRigid(modelCfg_) ),
                  k_cytoskeleton( WbcHighOrderModel::calculate_kCytoskeleton(modelCfg_) ),
                  core_radius(WbcHighOrderModel::calculate_coreRadius(modelCfg_)),
                  radius(WbcHighOrderModel::calculate_radius(modelCfg_))
    {};

void WbcHighOrderModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, size_t ctype) {

  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    if (cell[0]->celltype != ctype) continue; //only execute on correct particle

    //Calculate Cell Values that need all particles (but do it most efficient
    //tailored to this class)
    double volume = 0.0;
    int triangle_n = 0;
    vector<double> triangle_areas;
    triangle_areas.reserve(cellConstants.triangle_list.size());
    vector<hemo::Array<double,3>> triangle_normals;
    triangle_normals.reserve(cellConstants.triangle_list.size());

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
      volume += (-v210+v120+v201-v021-v102+v012); // the factor of 1/6 moved to after the summation -> saves a few flops
      
      //Area
      double area; 
      hemo::Array<double,3> t_normal;
      computeTriangleAreaAndUnitNormal(v0, v1, v2, area, t_normal);

      const double areaRatio = (area - /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n])
                               / /*cellConstants.area_mean_eq*/ cellConstants.triangle_area_eq_list[triangle_n];      
       
      //area force magnitude
      const double afm = k_area * (areaRatio+areaRatio/std::fabs(0.09-areaRatio*areaRatio));

      hemo::Array<double,3> centroid;
      centroid[0] = (v0[0]+v1[0]+v2[0])/3.0;
      centroid[1] = (v0[1]+v1[1]+v2[1])/3.0;
      centroid[2] = (v0[2]+v1[2]+v2[2])/3.0;
      hemo::Array<double,3> av0 = centroid - v0;
      hemo::Array<double,3> av1 = centroid - v1;
      hemo::Array<double,3> av2 = centroid - v2;

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
    const double volume_force = -k_volume * volume_frac/std::fabs(0.01-volume_frac*volume_frac);
    triangle_n = 0;

//Volume force loop
    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      // Scale volume force with local face area
      const hemo::Array<double, 3> local_volume_force = (volume_force*triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      *cell[triangle[0]]->force_volume += local_volume_force;
      *cell[triangle[1]]->force_volume += local_volume_force;
      *cell[triangle[2]]->force_volume += local_volume_force;

      triangle_n++;
    }

//Per-vertex bending force loop
    for (long unsigned int i = 0 ; i < cell.size() ; i++) {
      hemo::Array<double,3> vertexes_sum = {0.,0.,0.};
      hemo::Array<double,3> vertices_vel_sum = {0.,0.,0.};

      for(unsigned int j = 0; j < cellConstants.vertex_n_vertexes[i]; j++) {
        vertexes_sum += cell[cellConstants.vertex_vertexes[i][j]]->position;
        vertices_vel_sum += cell[cellConstants.vertex_vertexes[i][j]]->v;
      }
      const hemo::Array<double,3> vertexes_middle = vertexes_sum/cellConstants.vertex_n_vertexes[i];
      const hemo::Array<double,3> vertices_vavg = vertices_vel_sum/cellConstants.vertex_n_vertexes[i];

      const hemo::Array<double,3> dev_vect = vertexes_middle - cell[i]->position;
      
      
      // Get the local surface normal
      hemo::Array<double,3> patch_normal = {0.,0.,0.};
      for(unsigned int j = 0; j < cellConstants.vertex_n_vertexes[i]-1; j++) {
        hemo::Array<double,3> triangle_normal = crossProduct(cell[cellConstants.vertex_vertexes[i][j]]->position - cell[i]->position, 
                                                             cell[cellConstants.vertex_vertexes[i][j+1]]->position - cell[i]->position);
        triangle_normal /= norm(triangle_normal);  
        patch_normal += triangle_normal;                                                   
      }
      hemo::Array<double,3> triangle_normal = crossProduct(cell[cellConstants.vertex_vertexes[i][cellConstants.vertex_n_vertexes[i]-1]]->position - cell[i]->position, 
                                                           cell[cellConstants.vertex_vertexes[i][0]]->position - cell[i]->position);
      triangle_normal /= norm(triangle_normal);
      patch_normal += triangle_normal;
 
      patch_normal /= norm(patch_normal);
              
      const double ndev = dot(patch_normal, dev_vect); // distance along patch normal

#ifdef PRECALCULATED_ANGLES
      const double dDev = (ndev - cellConstants.surface_patch_center_dist_eq_list[i] ) / cellConstants.edge_mean_eq; // Non-dimensional
#else 
      const double dDev = ndev / cellConstants.edge_mean_eq; // Non-dimensional
#endif

      //TODO scale bending force
      const hemo::Array<double,3> bending_force = k_bend * ( dDev + dDev/std::fabs(0.055-dDev*dDev)) * patch_normal; // tau_b comes from the angle limit w. eq.lat.tri. assumptiln
      
      // Calculating viscous term
      const hemo::Array<double,3> rel_vel_v = vertices_vavg - cell[i]->v;
      const hemo::Array<double,3> rel_vel_proj = dot(patch_normal, rel_vel_v) * patch_normal;
      hemo::Array<double,3> Fvisc_vol = eta_v * rel_vel_proj * 0.866 * cellConstants.edge_mean_eq; // last term is triangle area from equilateral approx. x ratio of #triangle/#vertices -> surface area belonging to a vertex

      //Apply bending force
      *cell[i]->force_bending += bending_force;

      // Limit volume viscosity
      const double Fvisc_vol_mag = norm(Fvisc_vol);
      if (Fvisc_vol_mag > FORCE_LIMIT / 4.0) {
        Fvisc_vol *= (FORCE_LIMIT / 4.0) / Fvisc_vol_mag;
      }

      // Apply volume viscosity
      *cell[i]->force_visc += Fvisc_vol;
      const hemo::Array<double,3> negative_bending_force = -bending_force/cellConstants.vertex_n_vertexes[i];
      const hemo::Array<double,3> negative_bending_viscous_force = -Fvisc_vol/cellConstants.vertex_n_vertexes[i];
          
      for (unsigned int j = 0 ; j < cellConstants.vertex_n_vertexes[i]; j++ ) {
       *cell[cellConstants.vertex_vertexes[i][j]]->force_bending += negative_bending_force;
       *cell[cellConstants.vertex_vertexes[i][j]]->force_visc += negative_bending_viscous_force;
      }              
    }

    // Per-edge calculations
    int edge_n=0;
    for (const hemo::Array<plint,2> & edge : cellConstants.edge_list) {
      const hemo::Array<double,3> & p0 = cell[edge[0]]->position;
      const hemo::Array<double,3> & p1 = cell[edge[1]]->position;

      // Link force
      const hemo::Array<double,3> edge_vec = p1-p0;
      const double edge_length = norm(edge_vec);
      const hemo::Array<double,3> edge_uv = edge_vec/edge_length;
      const double edge_frac = (edge_length - /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n])
                               / /*cellConstants.edge_mean_eq*/ cellConstants.edge_length_eq_list[edge_n];

      const double edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
      const hemo::Array<double,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_link += force;
      *cell[edge[1]]->force_link -= force;

      // Membrane viscosity of bilipid layer
      // F = eta * (dv/l) * l. 
      const hemo::Array<double,3> rel_vel = cell[edge[1]]->v - cell[edge[0]]->v;
      const hemo::Array<double,3> rel_vel_projection = dot(rel_vel, edge_uv) * edge_uv;
      hemo::Array<double,3> Fvisc_memb = eta_m * rel_vel_projection;

      // Limit membrane viscosity
      const double Fvisc_memb_mag = norm(Fvisc_memb);
      if (Fvisc_memb_mag > FORCE_LIMIT / 4.0) {
        Fvisc_memb *= (FORCE_LIMIT / 4.0) / Fvisc_memb_mag;
      }

      *cell[edge[0]]->force_visc += Fvisc_memb;
      *cell[edge[1]]->force_visc -= Fvisc_memb; 

      edge_n++;
    }
    
    // Enforce rigid inner core size
    for (const hemo::Array<plint,2> & edge : cellConstants.inner_edge_list) {
      const hemo::Array<double,3> & p0 = cell[edge[0]]->position;
      const hemo::Array<double,3> & p1 = cell[edge[1]]->position;

      // Inner link forces
      const hemo::Array<double,3> edge_vec = p1-p0;
      const double edge_length = norm(edge_vec);

      if (edge_length < 2*radius && edge_length >= 2*core_radius){
        const hemo::Array<double,3> edge_uv = edge_vec/edge_length;
        const hemo::Array<double,3> force = edge_uv*(2*radius-edge_length)*k_cytoskeleton;
        *cell[edge[0]]->force_inner_link -= force;
        *cell[edge[1]]->force_inner_link += force;
      }

      if (edge_length < 2*core_radius){
        const hemo::Array<double,3> edge_uv = edge_vec/edge_length;
        const hemo::Array<double,3> force = edge_uv*(2*core_radius-edge_length)*k_inner_rigid;
        *cell[edge[0]]->force_inner_link -= force;
        *cell[edge[1]]->force_inner_link += force;
      }
    }
  } 
};

void WbcHighOrderModel::statistics() {
    pcout << "(Cell-mechanics model) High Order model parameters for " << cellField.name << " cellfield" << std::endl; 
    pcout << "\t k_link:   " << k_link << std::endl; 
    pcout << "\t k_area:   " << k_area << std::endl; 
    pcout << "\t k_bend: : " << k_bend << std::endl; 
    pcout << "\t k_volume: " << k_volume << std::endl;
    pcout << "\t k_cytoskeleton: " << k_cytoskeleton<< std::endl;
    pcout << "\t k_inner_rigid: " << k_inner_rigid << std::endl;
    pcout << "\t eta_m:    " << eta_m << std::endl;
    pcout << "\t eta_v:    " << eta_v << std::endl;
    pcout << "\t wbc_radius:    " << radius << std::endl;
    pcout << "\t core_radius:    " << core_radius << std::endl;
};


// Provide methods to calculate and scale to coefficients from here

double WbcHighOrderModel::calculate_etaV(Config & cfg ){
  return cfg["MaterialModel"]["eta_v"].read<double>() * param::dx * param::dt / param::dm; //== dx^2/dN/dt
};

double WbcHighOrderModel::calculate_etaM(Config & cfg ){
  return cfg["MaterialModel"]["eta_m"].read<double>() * param::dx / param::dt / param::df;
};

double WbcHighOrderModel::calculate_kBend(Config & cfg, MeshMetrics<double> & meshmetric ){
  double eqLength = meshmetric.getMeanLength();
  return cfg["MaterialModel"]["kBend"].read<double>() * param::kBT_lbm / eqLength;
};

double WbcHighOrderModel::calculate_kVolume(Config & cfg, MeshMetrics<double> & meshmetric){
  double kVolume =  cfg["MaterialModel"]["kVolume"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kVolume *= param::kBT_lbm/(eqLength*eqLength*eqLength);
  //kVolume /= meshmetric.getNumVertices();
  return kVolume;
};

double WbcHighOrderModel::calculate_kArea(Config & cfg, MeshMetrics<double> & meshmetric){
  double kArea =  cfg["MaterialModel"]["kArea"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kArea *= param::kBT_lbm/(eqLength*eqLength);
  return kArea;
};

double WbcHighOrderModel::calculate_kLink(Config & cfg, MeshMetrics<double> & meshmetric){
  double kLink = cfg["MaterialModel"]["kLink"].read<double>();
  double persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  //TODO: It should scale with the number of surface points!
  double plc = persistenceLengthFine/param::dx;
  kLink *= param::kBT_lbm/plc;
  return kLink;
};

double WbcHighOrderModel::calculate_kInnerRigid(Config & cfg){
  double kInnerRigid = cfg["MaterialModel"]["kInnerRigid"].read<double>();
  //TODO: convert to proper dimension
  return kInnerRigid;
};

double WbcHighOrderModel::calculate_kCytoskeleton(Config & cfg){
  double kCytoskeleton = cfg["MaterialModel"]["kCytoskeleton"].read<double>();
  return kCytoskeleton;
};

double WbcHighOrderModel::calculate_coreRadius(Config & cfg){
  double coreRadius = cfg["MaterialModel"]["coreRadius"].read<double>();
  return coreRadius/param::dx;
};

double WbcHighOrderModel::calculate_radius(Config & cfg){
  double radius = cfg["MaterialModel"]["radius"].read<double>();
  return radius/param::dx;
};
