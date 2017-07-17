#ifndef Hemo_CommonCellConstants_CPP
#define Hemo_CommonCellConstants_CPP
#include "commonCellConstants.h"
#include "hemoCellParticleType.h"

CommonCellConstants::CommonCellConstants(HemoCellField & cellField_,
                      vector<hemo::Array<plint,3>> triangle_list_,
                      vector<hemo::Array<plint,2>> edge_list_,
                      vector<double> edge_length_eq_list_,
                      vector<double> edge_angle_eq_list_,
                      vector<hemo::Array<plint,2>> edge_bending_triangles_list_,
                      vector<hemo::Array<plint,2>> edge_bending_triangles_outer_points_,
                      vector<double> triangle_area_eq_list_,
                      double volume_eq_, double area_mean_eq_, 
                      double edge_mean_eq_, double angle_mean_eq_) :
    cellField(cellField_),
    triangle_list(triangle_list_),
    edge_list(edge_list_),
    edge_length_eq_list(edge_length_eq_list_),
    edge_angle_eq_list(edge_angle_eq_list_),
    edge_bending_triangles_list(edge_bending_triangles_list_),
    edge_bending_triangles_outer_points(edge_bending_triangles_outer_points_),
    triangle_area_eq_list(triangle_area_eq_list_),
    volume_eq(volume_eq_),
    area_mean_eq(area_mean_eq_),
    edge_mean_eq(edge_mean_eq_),
    angle_mean_eq(angle_mean_eq_)
  {};

CommonCellConstants CommonCellConstants::CommonCellConstantsConstructor(HemoCellField & cellField_) {
    HemoCellField & cellField = cellField_;
    //Calculate triangles
    vector<hemo::Array<plint,3>> triangle_list_ = cellField.triangle_list;

    //Calculate edges //TSM uses "directed edges", We dont care about memory, we
    //just store them, so calculate all vertices here
    //TODO, every edge is in here twice, clean that?
    vector<hemo::Array<plint,2>> edge_list_;
    for (const hemo::Array<plint,3> & triangle : triangle_list_) {
 
    //FIX by allowing only incremental edges
      if(triangle[0] < triangle[1]) {
        edge_list_.push_back({triangle[0],triangle[1]});
      }
      if(triangle[1] < triangle[2]) {
        edge_list_.push_back({triangle[1],triangle[2]});
      }
      if (triangle[2] < triangle[0]) {
        edge_list_.push_back({triangle[2],triangle[0]});
      }
    }

    //Calculate eq edges
    vector<double> edge_length_eq_list_;
    for (const hemo::Array<plint,2> & edge : edge_list_) {
      edge_length_eq_list_.push_back(cellField.meshElement.computeEdgeLength(edge[0],edge[1]));
    }

    //Calculate eq edges angles
    vector<double> edge_angle_eq_list_;
    hemo::Array<double,3> x2 = {0.0,0.0,0.0};

    for (const hemo::Array<plint,2> & edge : edge_list_) {
      const vector<plint> adjacentTriangles = cellField.meshElement.getAdjacentTriangleIds(edge[0], edge[1]);
      //const hemo::Array<double,3> x1 = cellField.meshElement.getVertex(edge[0]);

      for (pluint id = 0 ; id < 3 ; id ++ ) {
        const plint kVertex = cellField.meshElement.getVertexId(adjacentTriangles[0],id);
        if (kVertex != edge[0] && kVertex != edge[1]) {
          x2 = cellField.meshElement.getVertex(kVertex);
          break;
        }
      }

      const hemo::Array<double,3> V1 = cellField.meshElement.computeTriangleNormal(adjacentTriangles[0]);
      const hemo::Array<double,3> V2 = cellField.meshElement.computeTriangleNormal(adjacentTriangles[1]);

      // double angle = angleBetweenVectors(V1, V2);

      // ///TODO does this gives us a correct "sign"?
      // const plint sign = dot(x2-x1, V2) >= 0 ? 1 : -1;
      // if (sign <= 0) {
      //   angle = 2 * PI - angle;
      // }
      
      // if (angle > PI) {
      //   angle = angle - 2*PI;
      // }

      const hemo::Array<double,3> p0 = cellField.meshElement.getVertex(edge[0]);
      const hemo::Array<double,3> p1 = cellField.meshElement.getVertex(edge[1]);

      const hemo::Array<double,3> edge_vec = p1-p0;
      const double edge_length = norm(edge_vec);
      const hemo::Array<double,3> edge_uv = edge_vec/edge_length;
      double angle = getAngleBetweenFaces(V1, V2, edge_uv);
      
      edge_angle_eq_list_.push_back(angle);
    }

    //Calculate triangle eq
    vector<double> triangle_area_eq_list_;
    for (pluint iTriangle = 0 ; iTriangle < cellField.triangle_list.size(); iTriangle++) {
      triangle_area_eq_list_.push_back(cellField.meshElement.computeTriangleArea(iTriangle));
    }

    
    double volume_eq_ = MeshMetrics<double>(cellField.meshElement).getVolume();


    //store important points for bending calculation
    vector<hemo::Array<plint,2>> edge_bending_triangles_;
    vector<hemo::Array<plint,2>> edge_bending_triangles_outer_points_;
    for (const hemo::Array<plint,2> & edge : edge_list_) {
      const vector<plint> adjacentTriangles = cellField.meshElement.getAdjacentTriangleIds(edge[0], edge[1]);
      edge_bending_triangles_.push_back({adjacentTriangles[0],adjacentTriangles[1]});
      hemo::Array<plint,2> op;
      for (int i = 0; i < 3; i++) {
        if ((triangle_list_[adjacentTriangles[0]][i] != edge[0]) && (triangle_list_[adjacentTriangles[0]][i] != edge[1])) {
          op[0] = triangle_list_[adjacentTriangles[0]][i];
        }
        if ((triangle_list_[adjacentTriangles[1]][i] != edge[0]) && (triangle_list_[adjacentTriangles[1]][i] != edge[1])) {
          op[1] = triangle_list_[adjacentTriangles[1]][i];
        }
      }
      edge_bending_triangles_outer_points_.push_back(op);
    }

    //Calculate mean face areas
    double mean_area_eq_ = 0.;
    for (const double & area : triangle_area_eq_list_) {
      mean_area_eq_ += area;
    }
    mean_area_eq_ /= triangle_area_eq_list_.size();
    
    // Calculate mean edge lengths
    double mean_edge_eq_ = 0;
    for(const double & edge : edge_length_eq_list_) {
      mean_edge_eq_ += edge;
    }
    mean_edge_eq_ /= edge_length_eq_list_.size();

    // Calculate mean edge angles
    double mean_angle_eq_ = 0;
    for(const double & angle : edge_angle_eq_list_) {
      mean_angle_eq_ += angle;
    }
    mean_angle_eq_ /= edge_angle_eq_list_.size();

    CommonCellConstants CCC(cellField_,triangle_list_,edge_list_,edge_length_eq_list_,edge_angle_eq_list_,edge_bending_triangles_,edge_bending_triangles_outer_points_,triangle_area_eq_list_,volume_eq_,mean_area_eq_, mean_edge_eq_, mean_angle_eq_);
    return CCC;
};

#endif
