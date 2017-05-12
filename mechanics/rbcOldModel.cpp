#include "rbcOldModel.h"

vector<HemoCellParticle*> * glob_cell;
HemoCellField * glob_cf;

RbcOldModel::RbcOldModel(Config & modelCfg_, HemoCellField & cellField_) : CellMechanics(),
                  cellConstants(CommonCellConstants::CommonCellConstantsConstructor(cellField_)), cellField(cellField_),
                  k_volume( RbcOldModel::calculate_kVolume(modelCfg_,*cellField_.meshmetric) ),
                  k_area( RbcOldModel::calculate_kArea(modelCfg_,*cellField_.meshmetric) ), 
                  k_link( RbcOldModel::calculate_kLink(modelCfg_,*cellField_.meshmetric) ), 
                  k_bend( RbcOldModel::calculate_kBend(modelCfg_) )
    {};


inline Array<double,3> computeVolumeConservationForce(Array<double,3> const& x1, Array<double,3> const& x2, 
                               Array<double,3> const& x3, double cVolume) {
  Array<double,3> tmp;
  crossProduct(x2,x3,tmp);
  return -cVolume * 1.0/60 * tmp;
}

inline Array<double,3> computeHighOrderLocalAreaConservationForce(
        Array<double,3> const& x1, Array<double,3> const& x2, Array<double,3> const& x3,
        Array<double,3> const& triangleNormal,  double triangleArea, double eqArea,
        double cArea) {

    Array<double,3> tmp;
    crossProduct(triangleNormal,
                 x3 - x2,
                 tmp);
    Array<double,3> dAdx = 0.5 *tmp;

    double dAreaRatio = (triangleArea - eqArea) / eqArea;

    // Only allow a stretch/compression of a maximum of 20% -> this usally results in a 11% stretch at max
    return -cArea * (dAreaRatio + dAreaRatio / (0.04-dAreaRatio*dAreaRatio) ) * dAdx;

}

inline Array<double,3> computeInPlaneHighOrderForce(Array<double,3> const& x1, Array<double,3> const& x2, double eqLength, double k_inPlane) {

    Array<double,3> vL = (x1 - x2);
    double L = norm(vL);
    Array<double,3> eij = vL/L;

    double dL = (L-eqLength)/eqLength;

    /* Spectrin links that are somewhat compressible */
    double force1D;
    if(dL > 0)
        force1D = - k_inPlane * ( dL + dL/(0.64-dL*dL) );   // allows at max. 80% stretch
    else
        force1D = - k_inPlane * dL * dL * dL;   // less stiff compression resistance -> let compression be dominated by area conservation force

    Array<double,3> force = eij * force1D;

    return force;
}

inline Array<double,3> getVertex(plint id) {
  return (*glob_cell)[id]->position;
}

vector<plint> getAdjacentTriangleIds(plint iV, plint jV) {
  return (*glob_cf).meshElement.getAdjacentTriangleIds(iV, jV);
}

inline Array<double,3> computeTriangleNormal(Array<double,3> v0, Array<double,3> v1, Array<double,3> v2) {
  Array<double,3> e01 = v1 - v0;
  Array<double,3> e02 = v2 - v0;
  Array<double,3> n;
  crossProduct(e01,e02,n);
  n /= norm(n);
  return n;
}

inline double computeSignedAngle(plint iVertex, plint jVertex, plint & kVertex, plint & lVertex, bool& found, plint edgeId,const vector<Array<plint,3>> & triangles) {
  found = true;
    Array<double,3> x1 = getVertex(iVertex), x2(0.,0.,0.), x3(0.,0.,0.), x4(0.,0.,0.);

    std::vector<plint> adjacentTriangles = getAdjacentTriangleIds(iVertex, jVertex);
    plint iTriangle=adjacentTriangles[0], jTriangle=adjacentTriangles[1];
    x3 = getVertex(jVertex);
    double foundVertices=0;
    for (pluint id = 0; id < 3; ++id) {
        kVertex = triangles[iTriangle][id];
        if ( (kVertex != iVertex) && (kVertex != jVertex)) {
            x2 = getVertex(kVertex);
            foundVertices += 1;
            break;
        }
    }
    for (pluint id = 0; id < 3; ++id) {
        lVertex = triangles[jTriangle][id];
        if ( (lVertex != iVertex) && (lVertex != jVertex)) {
            x4 = getVertex(lVertex);
            foundVertices += 1;
            break;
        }
    }
  found = (foundVertices == 2); //Assert if some particles are outside of the domain
//What does this signedAngles vector do? circumvent the double edge problem? please... not like
//this.. TODO fixit
  if (not found) { /*signedAngles.erase(edgeId);*/  return 0.0; }
  //if (signedAngles.count(edgeId) == 0) {
      Array<double,3> V1 = computeTriangleNormal(getVertex(triangles[iTriangle][0]),getVertex(triangles[iTriangle][1]),getVertex(triangles[iTriangle][2]));
      Array<double,3> V2 = computeTriangleNormal(getVertex(triangles[jTriangle][0]),getVertex(triangles[jTriangle][1]),getVertex(triangles[jTriangle][2]));
      double angle = angleBetweenVectors(V1, V2);
    plint sign = dot(x2-x1, V2) >= 0?1:-1;

    if (sign <= 0) {
      angle = 2*PI-angle;
    }
    angle = (angle > PI)?angle-2*PI:angle;
   // signedAngles[edgeId] = angle;
  //}
  return angle;
}

inline Array<double,3> computeBendingForce4p (Array<double,3> const& xi, Array<double,3> const& xj, 
                                Array<double,3> const& xk, Array<double,3> const& xl, 
                                Array<double,3> const& nTk, Array<double,3> const& nTl, 
                                //T Ai, T Aj, 
                                double eqArea, double eqLength, double eqAngle, double k, 
                                Array<double,3> & iFx, Array<double,3> & jFx, Array<double,3> & kFx, Array<double,3> & lFx) { 
/* The most messy force! 
 * Triangles are: 
 *      (i, j, k) and (l, k, j). These triangles share 
 *      (x2, x1, x3) and (x4, x3, x1). These triangles share 
 *      the common edge j-k or 1-3. 
 * 
 *      crossProduct(jPosition - iPosition, kPosition - jPosition, nijk); 
 *      crossProduct(kPosition - lPosition, jPosition - kPosition, nlkj); 
 *      crossProduct(x1 - x2, x3 - x1, nijk); 
 *      crossProduct(x3 - x4, x1 - x3, nlkj); 
*/ 
    double dAngle; 
    double edgeAngle = angleBetweenVectors(nTk, nTl); 
         
    plint sign = dot(xk-xi, nTl) > 0?1:-1; 
    if (sign <= 0) { 
        edgeAngle = 2*PI-edgeAngle; 
    } 
    edgeAngle = (edgeAngle > PI)?edgeAngle-2*PI:edgeAngle; 
    eqAngle = (eqAngle > 2*PI)?eqAngle-2*PI:eqAngle; 
    eqAngle = (eqAngle > PI)?eqAngle-2*PI:eqAngle; 
 
    // WARNING: sign is mixed up! (Clockwise vs. anti-clockwise problem) 
    //dAngle = (edgeAngle-eqAngle); 
    dAngle = eqAngle-edgeAngle; 
 
    // Linear force 
    //T force = -k*(dAngle) * (eqLength*0.5/eqArea); 
    double force = -k * dAngle; 
 
    // Force based on sphere-curvature model 
    // T force = -k * sin( dAngle ); 
 
    #if HEMOCELL_MEMBRANE_BENDING == 2 
        kFx = -force*nTk*0.5; // Reduced to half to avoid self-cancellation. 
        lFx = -force*nTl*0.5; 
        iFx = -(kFx+lFx);    
        jFx = iFx; 
    #elif HEMOCELL_MEMBRANE_BENDING == 3 
        kFx = -force*nTk; // These on a closed and regular surface will mostly cancel eachother out 
        lFx = -force*nTl; 
        iFx = -(kFx+lFx)*0.5;    
        jFx = iFx; 
    #endif 
 
    return iFx; 
} 


 


void RbcOldModel::ParticleMechanics(map<int,vector<HemoCellParticle *>> particles_per_cell, map<int,bool> lpc, pluint ctype) {

  for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
    const int & cid = pair.first;
    vector<HemoCellParticle*> & cell = particles_per_cell[cid];
    glob_cell = &cell;
    glob_cf = &cellField;
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
      
      //Triangle Normal and area
      //TODO this is so wrong
      const Array<double,3> e01 = v1 - v0;
      const Array<double,3> e02 = v2 - v0;
      Array<double,3> triangle_normal;
      crossProduct(e01, e02, triangle_normal);
      const double normN = norm(triangle_normal);
      triangle_normal /= normN;
      const double triangle_area = normN * 0.5;
        
      //Area
      //push back area force
      *cell[triangle[0]]->force_area += computeHighOrderLocalAreaConservationForce(v0,v1,v2,triangle_normal,triangle_area,cellConstants.triangle_area_eq_list[triangle_n],k_area);
      *cell[triangle[1]]->force_area += computeHighOrderLocalAreaConservationForce(v1,v2,v0,triangle_normal,triangle_area,cellConstants.triangle_area_eq_list[triangle_n],k_area);
      *cell[triangle[2]]->force_area += computeHighOrderLocalAreaConservationForce(v2,v0,v1,triangle_normal,triangle_area,cellConstants.triangle_area_eq_list[triangle_n],k_area);

      triangle_n++;
    }

    //Volume
    const double volume_frac = (volume-cellConstants.volume_eq)/cellConstants.volume_eq;
    const double volume_force = k_volume * volume_frac/(0.01-volume_frac*volume_frac);

    for (const Array<plint,3> & triangle : cellConstants.triangle_list) {
      //TODO volume force per area
      const Array<double,3> & v0 = cell[triangle[0]]->position;
      const Array<double,3> & v1 = cell[triangle[1]]->position;
      const Array<double,3> & v2 = cell[triangle[2]]->position;
      *cell[triangle[0]]->force_volume += computeVolumeConservationForce(v0,v1,v2,volume_force);
      *cell[triangle[1]]->force_volume += computeVolumeConservationForce(v1,v2,v0,volume_force);
      *cell[triangle[2]]->force_volume += computeVolumeConservationForce(v2,v0,v1,volume_force);

    }


    // Edges
    int edge_n=0;
    for (const Array<plint,2> & edge : cellConstants.edge_list) {
      const Array<double,3> & v0 = cell[edge[0]]->position;
      const Array<double,3> & v1 = cell[edge[1]]->position;

      // Link force
      const Array<double,3> l_force = computeInPlaneHighOrderForce(v0,v1,cellConstants.edge_length_eq_list[edge_n],k_link);
      *cell[edge[0]]->force_link += l_force;
      *cell[edge[1]]->force_link -= l_force;
      
			//Bending
			bool angleFound;
			plint kVertex, lVertex;
			const double edgeAngle = computeSignedAngle(edge[0],edge[1],kVertex, lVertex, angleFound,edge_n,cellConstants.triangle_list);

			if (angleFound) {
				Array<double,3> iNormal = computeTriangleNormal(getVertex(edge[0]),getVertex(edge[1]),getVertex(kVertex));
				Array<double,3> jNormal = computeTriangleNormal(getVertex(edge[0]),getVertex(edge[1]),getVertex(lVertex));
        double Ai = computeTriangleArea(getVertex(edge[0]),getVertex(edge[1]),getVertex(kVertex));
        double Aj = computeTriangleArea(getVertex(edge[0]),getVertex(edge[1]),getVertex(lVertex));

        Array<double,3> fi,fj,fk,fl;
        //OKAY... eqArea is actually the mean area of all triangles, this is
        //like MEGA-WRONG. lol
        fi = computeBendingForce4p (v0,v1,cell[kVertex]->position, cell[lVertex]->position, iNormal, jNormal, 
                                    cellConstants.area_mean_eq, cellConstants.edge_length_eq_list[edge_n], 
                                    cellConstants.edge_angle_eq_list[edge_n], k_bend, fi, fj, fk, fl);

                                    
        *cell[edge[0]]->force_bending += fi;
        *cell[edge[1]]->force_bending += fj;
        *cell[kVertex]->force_bending += fk;
        *cell[lVertex]->force_bending += fl;
      }
      edge_n++;
    }

  } 
};

void RbcOldModel::statistics() {
    pcout << "High Order (old model) forces for " << cellField.name << " cellfield" << std::endl;
    pcout << "k_volume: " << k_volume << std::endl; 
    pcout << "k_area:   " << k_area << std::endl; 
    pcout << "k_link:   " << k_link << std::endl; 
    pcout << "k_bend: : " << k_bend << std::endl; 
};


// Provide methods to calculate and scale to coefficients from here

double RbcOldModel::calculate_kBend(Config & cfg ){
  return cfg["MaterialModel"]["kBend"].read<double>() * param::kBT_lbm;
};

double RbcOldModel::calculate_kVolume(Config & cfg, MeshMetrics<double> & meshmetric){
  double kVolume =  cfg["MaterialModel"]["kVolume"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kVolume *= param::kBT_lbm/(eqLength*eqLength*eqLength);
  return kVolume;
};

double RbcOldModel::calculate_kArea(Config & cfg, MeshMetrics<double> & meshmetric){
  double kArea =  cfg["MaterialModel"]["kArea"].read<double>();
  double eqLength = meshmetric.getMeanLength();
  kArea *= param::kBT_lbm/(eqLength*eqLength);
  return kArea;
}

double RbcOldModel::calculate_kLink(Config & cfg, MeshMetrics<double> & meshmetric){
  double kLink = cfg["MaterialModel"]["kLink"].read<double>();
  double persistenceLengthFine = 7.5e-9; // In meters -> this is a biological value
  // TODO: this is a fixed number, no need to calculate it like this
  double plc = persistenceLengthFine/param::dx * sqrt((meshmetric.getNumVertices()-2.0) / (23867-2.0)); //Kaniadakis magic
  kLink *= param::kBT_lbm/(4.0*plc);
  return kLink;
}
