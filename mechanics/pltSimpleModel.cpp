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
#include "pltSimpleModel.h"
#include "logfile.h"
#include "octree.h"
#include "mollerTrumbore.h"

#include "palabos3D.h"
#include "palabos3D.hh"


namespace hemo {
PltSimpleModel::PltSimpleModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(cellField_, modelCfg_), 
                  cellField(cellField_),
                  k_volume( PltSimpleModel::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ),
                  k_area( PltSimpleModel::calculate_kArea(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( PltSimpleModel::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( PltSimpleModel::calculate_kBend(modelCfg_,*cellField_.meshmetric) ),
                  eta_m( PltSimpleModel::calculate_etaM(modelCfg_))
  { };

void PltSimpleModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> & particles_per_cell, const map<int,bool> & lpc, pluint ctype) {
  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    if (cell.size() == 0) continue;
    if (cell[0]->sv.celltype != ctype) continue; //only execute on correct particle

    //Calculate Cell Values that need all particles (but do it efficiently,
    //tailored to this class)
    T volume = 0.0;
    int triangle_n = 0;
    vector<T> triangle_areas;
    vector<hemo::Array<T,3>> triangle_normals; 

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
      volume += (-v210+v120+v201-v021-v102+v012);
      
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

    for (const hemo::Array<plint,3> & triangle : cellConstants.triangle_list) {
      //Fixed volume force per area
      const hemo::Array<T, 3> local_volume_force = (volume_force*triangle_normals[triangle_n])*(triangle_areas[triangle_n]/cellConstants.area_mean_eq);
      *cell[triangle[0]]->force_volume += local_volume_force;
      *cell[triangle[1]]->force_volume += local_volume_force;
      *cell[triangle[2]]->force_volume += local_volume_force;

      triangle_n++;
    }


    // Per-edge calculations
    int edge_n=0;
    for (const hemo::Array<plint,2> & edge : cellConstants.edge_list) {
      const hemo::Array<T,3> & v0 = cell[edge[0]]->sv.position;
      const hemo::Array<T,3> & v1 = cell[edge[1]]->sv.position;

      // Link force
      const hemo::Array<T,3> edge_v = v1-v0;
      const T edge_length = sqrt(edge_v[0]*edge_v[0]+edge_v[1]*edge_v[1]+edge_v[2]*edge_v[2]);
      const hemo::Array<T,3> edge_uv = edge_v/edge_length;
      const T edge_frac = (edge_length-cellConstants.edge_length_eq_list[edge_n])/cellConstants.edge_length_eq_list[edge_n];

      const T edge_force_scalar = k_link * ( edge_frac + edge_frac/std::fabs(9.0-edge_frac*edge_frac));   // allows at max. 300% stretch
      
      const hemo::Array<T,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_link += force;
      *cell[edge[1]]->force_link -= force;

      // Membrane viscosity of bilipid layer
      // F = eta * (dv/l) * l. 
      const hemo::Array<T,3> rel_vel = cell[edge[1]]->sv.v - cell[edge[0]]->sv.v;
      const hemo::Array<T,3> rel_vel_projection = dot(rel_vel, edge_uv) * edge_uv;
      hemo::Array<T,3> Fvisc_memb = eta_m * rel_vel_projection;

      // Limit membrane viscosity
      const T Fvisc_memb_mag = norm(Fvisc_memb);
      if (Fvisc_memb_mag > FORCE_LIMIT / 4.0) {
        Fvisc_memb *= (FORCE_LIMIT / 4.0) / Fvisc_memb_mag;
      }

      *cell[edge[0]]->force_visc += Fvisc_memb;
      *cell[edge[1]]->force_visc -= Fvisc_memb; 


      const plint b0 = cellConstants.edge_bending_triangles_list[edge_n][0];
      const plint b1 = cellConstants.edge_bending_triangles_list[edge_n][1];

      const hemo::Array<T,3> b00 = particles_per_cell[cid][cellField.triangle_list[b0][0]]->sv.position;
      const hemo::Array<T,3> b01 = particles_per_cell[cid][cellField.triangle_list[b0][1]]->sv.position;
      const hemo::Array<T,3> b02 = particles_per_cell[cid][cellField.triangle_list[b0][2]]->sv.position;
      
      const hemo::Array<T,3> b10 = particles_per_cell[cid][cellField.triangle_list[b1][0]]->sv.position;
      const hemo::Array<T,3> b11 = particles_per_cell[cid][cellField.triangle_list[b1][1]]->sv.position;
      const hemo::Array<T,3> b12 = particles_per_cell[cid][cellField.triangle_list[b1][2]]->sv.position;

      const hemo::Array<T,3> V1 = computeTriangleNormal(b00,b01,b02, false);
      const hemo::Array<T,3> V2 = computeTriangleNormal(b10,b11,b12, false);

      T angle = getAngleBetweenFaces(V1, V2, edge_uv);
      
      //calculate resulting bending force
      const T angle_frac = angle - cellConstants.edge_angle_eq_list[edge_n];

      const T force_magnitude = k_bend * (angle_frac + angle_frac / std::fabs(2.467 - angle_frac * angle_frac) ); // tau_b = pi/2
      
      //TODO Make bending force differ with area!
      const hemo::Array<T,3> bending_force = force_magnitude*(V1 + V2)*0.5;
      *cell[edge[0]]->force_bending += bending_force;
      *cell[edge[1]]->force_bending += bending_force;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][0]]->force_bending -= bending_force;
      *cell[cellConstants.edge_bending_triangles_outer_points[edge_n][1]]->force_bending -= bending_force;

      edge_n++;
    }

    // Per-inner-edge caluclations
    int inner_edge_n=0;
    for (const hemo::Array<plint,2> & edge : cellConstants.inner_edge_list) {
      const hemo::Array<T,3> & v0 = cell[edge[0]]->sv.position;
      const hemo::Array<T,3> & v1 = cell[edge[1]]->sv.position;

      // Link force
      const hemo::Array<T,3> edge_v = v1-v0;
      const T edge_length = sqrt(edge_v[0]*edge_v[0]+edge_v[1]*edge_v[1]+edge_v[2]*edge_v[2]);
      const hemo::Array<T,3> edge_uv = edge_v/edge_length;
      const T edge_frac = (edge_length-cellConstants.inner_edge_length_eq_list[inner_edge_n])/cellConstants.inner_edge_length_eq_list[inner_edge_n];

      const T edge_force_scalar = k_link * 5.0 * edge_frac; // Keep the linear part only for stability  
      
      const hemo::Array<T,3> force = edge_uv*edge_force_scalar;
      *cell[edge[0]]->force_inner_link += force;
      *cell[edge[1]]->force_inner_link -= force;
      inner_edge_n++;
    }

  } 
}

#ifdef SOLIDIFY_MECHANICS
void PltSimpleModel::solidifyMechanics(const std::map<int,std::vector<int>>& ppc,std::vector<HemoCellParticle>& particles,plb::BlockLattice3D<T,DESCRIPTOR> * fluid,plb::BlockLattice3D<T,CEPAC_DESCRIPTOR> * CEPAC, pluint ctype) {
  hemo::Array<T,3> * pos;
  Dot3D const& location_CEPAC = CEPAC->getLocation();
  Dot3D const& location_fluid = fluid->getLocation();
  double threshold = cfg["MaterialModel"]["solidifyThreshold"].read<double>();
  double shear_threshold = cfg["MaterialModel"]["shearThreshold"].read<double>();
  double wall_distance = cfg["MaterialModel"]["distanceThreshold"].read<double>();

  //For all cells
  for (auto & pair : ppc) {
    bool broken = false;
    const std::vector<int> & cell = pair.second;
    //For all particles of cell
    for (const int & particle : cell ) {
      //Skip non-complete and non-platelets
      if (particle == -1) { broken = true; break; }
      if (particles[particle].sv.celltype != ctype) { broken = true; break; }
    }
    if (broken) { continue; }
    
    bool solidify = false;
    // Complete and Correct Type, do solidify mechanics:
    for (const int & particle : cell) {
      //Firstly check if any particle should be solidified
      if (particles[particle].sv.solidify) {
        solidify = true;
          break;
      }
    }

    // If it was tagged last round, solidify it now
    if (solidify) {
        //pcout << particles[cell[0]].sv.cellId << endl;
      hemo::OctreeStructCell octCell(3, 1, 30, cellConstants.triangle_list, particles, cell);
 
      vector<Array<plint,3>> innerNodes;
      octCell.findInnerNodes(fluid,particles,cell,innerNodes);
      //pcout << innerNodes[0][0] << endl;      
      for (Array<plint,3> & node : innerNodes) {
            //pcout << node[0] << endl;
          if (!fluid->get(node[0],node[1],node[2]).getDynamics().isBoundary()) {
          defineDynamics(*fluid,node[0],node[1],node[2],new BounceBack<T,DESCRIPTOR>(1.));
          
        }
      }

      vector<Array<plint,3>> innerNodesCEPAC;
      octCell.findInnerNodes(CEPAC,particles,cell,innerNodesCEPAC);
      for (Array<plint,3> & node : innerNodesCEPAC) {
        if (!CEPAC->get(node[0],node[1],node[2]).getDynamics().isBoundary()) {
          defineDynamics(*CEPAC,node[0],node[1],node[2],new BounceBack<T,CEPAC_DESCRIPTOR>(1.));
        }
      }


      for (const int & particle : cell) {
        particles[particle].tag = 1; //tag for removal
      }
    } //else {
      //Otherwise, see if we have to tag for solidify
//      for (const int & particle : cell) {
//        pos = &particles[particle].sv.position;
//      
//        int x = pos->operator[](0)-location_CEPAC.x+0.5;
//        int y = pos->operator[](1)-location_CEPAC.y+0.5;
//        int z = pos->operator[](2)-location_CEPAC.z+0.5;
//        if ((x >= 0) && (x < CEPAC->getNx()) &&
//            (y >= 0) && (y < CEPAC->getNy()) &&
//            (z >= 0) && (z < CEPAC->getNz()) ) {
 //         if (CEPAC->get(x,y,z).computeDensity() > threshold) {
//            particles[particle].sv.solidify = true;
//            //pcout << particles[particle].sv.cellId << endl;
//          }
//        }
//      }
//      for (const int & particle : cell) {
//        pos = &particles[particle].sv.position;
//        
//        int x_f = pos->operator[](0)-location_fluid.x+0.5;
//        int y_f = pos->operator[](1)-location_fluid.y+0.5;
//        int z_f = pos->operator[](2)-location_fluid.z+0.5;
//        
//         if ( (x_f > fluid->getBoundingBox().x0-0.5) && (x_f <= fluid->getBoundingBox().x1+0.5) &&
//	   (y_f > fluid->getBoundingBox().y0-0.5) && (y_f <= fluid->getBoundingBox().y1+0.5) &&
//	   (z_f > fluid->getBoundingBox().z0-0.5) && (z_f <= fluid->getBoundingBox().z1+0.5) );
//        { continue;
//         }
//           
//        T x_f_r = pos->operator[](0)-location_fluid.x;
//        T y_f_r = pos->operator[](1)-location_fluid.y;
//        T z_f_r = pos->operator[](2)-location_fluid.z;
//        
//         if ( (x_f_r > fluid->getBoundingBox().x0-0.5) && (x_f_r <= fluid->getBoundingBox().x1+0.5) &&
//	   (y_f_r > fluid->getBoundingBox().y0-0.5) && (y_f_r <= fluid->getBoundingBox().y1+0.5) &&
//	   (z_f_r > fluid->getBoundingBox().z0-0.5) && (z_f_r <= fluid->getBoundingBox().z1+0.5) );
//        { continue;
//         }
//        
//        //plb::Array<T,6> shearrate;
//        double avg_shearrate;
//        //fluid->get(x_f,y_f,z_f).computeShearRate(*fluid,shearrate,x_f,y_f,z_f);
//        //avg_shearrate=sqrt(shearrate[1]*shearrate[1]+shearrate[2]*shearrate[2]+shearrate[4]*shearrate[4])/param::dt;
//        
//        std::auto_ptr<TensorField3D<T,6> > strainrate (computeStrainRateFromStress(*fluid));
//        std::auto_ptr<ScalarField3D<T> > shearrate_scalar (computeSymmetricTensorNorm(*strainrate));
//        avg_shearrate = computeAverage(*shearrate_scalar);
//        //plb::Array<T,6> strainrate;
//        //computeStrainRateFromStress(*fluid,strainrate,fluid->getBoundingBox());
        
//        if ((avg_shearrate > shear_threshold )) {
//            for (int xx = x_f-6; xx<=x_f+6; xx++) {
//                for (int yy = y_f-6; yy<= y_f+6; yy++) {
//                    for (int zz = z_f-6; zz<=z_f+6; zz++) {
//                        if (fluid->get(xx,yy,zz).getDynamics().isBoundary()) {
//                            if(sqrt((xx-x_f_r)*(xx-x_f_r) + (yy-y_f_r)*(yy-y_f_r)+(zz-z_f_r)*(zz-z_f_r)) > wall_distance ) {
//                                particles[particle].sv.solidify = true;
//                            }
//                        }
//                    }
//                }
//                
//            }
//        }
//            
//        }
//    }
//  }
}
}
#endif

void PltSimpleModel::statistics() {
    hlog << "(Cell-mechanics model) Reduced-model parameters for " << cellField.name << " cellfield" << std::endl;
    hlog << "\t k_link:   " << k_link << std::endl; 
    hlog << "\t k_bend: : " << k_bend << std::endl; 
    hlog << "\t k_volume: " << k_volume << std::endl; 
    hlog << "\t eta_m:      " << eta_m << std::endl;
};
}
