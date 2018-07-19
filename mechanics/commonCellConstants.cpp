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
#include "commonCellConstants.h"

namespace hemo {
 using namespace std; 
CommonCellConstants::CommonCellConstants(HemoCellField & cellField_,
                      vector<hemo::Array<plint,3>> triangle_list_,
                      vector<hemo::Array<plint,2>> edge_list_,
                      vector<T> edge_length_eq_list_,
                      vector<T> edge_angle_eq_list_,
                      vector<T> surface_patch_center_dist_eq_list_,
                      vector<hemo::Array<plint,2>> edge_bending_triangles_list_,
                      vector<hemo::Array<plint,2>> edge_bending_triangles_outer_points_,
                      vector<T> triangle_area_eq_list_,
                      vector<hemo::Array<plint,6>> vertex_vertexes_,
                      vector<hemo::Array<plint,6>> vertex_edges_,
                      vector<hemo::Array<signed int,6>> vertex_edges_sign_,
                      vector<unsigned int> vertex_n_vertexes_,
                      vector<hemo::Array<hemo::Array<plint,2>,6>> vertex_outer_edges_per_vertex_,
                      vector<hemo::Array<hemo::Array<signed int,2>,6>> vertex_outer_edges_per_vertex_sign_,
                      T volume_eq_, T area_mean_eq_, 
                      T edge_mean_eq_, T angle_mean_eq_,
                      vector<hemo::Array<plint,2>> inner_edge_list_,
                      vector<T> inner_edge_length_eq_list_) :
    cellField(cellField_),
    triangle_list(triangle_list_),
    edge_list(edge_list_),
    edge_length_eq_list(edge_length_eq_list_),
    edge_angle_eq_list(edge_angle_eq_list_),
    surface_patch_center_dist_eq_list(surface_patch_center_dist_eq_list_),
    edge_bending_triangles_list(edge_bending_triangles_list_),
    edge_bending_triangles_outer_points(edge_bending_triangles_outer_points_),
    triangle_area_eq_list(triangle_area_eq_list_),
    vertex_vertexes(vertex_vertexes_),
    vertex_edges(vertex_edges_),
    vertex_edges_sign(vertex_edges_sign_),
    vertex_n_vertexes(vertex_n_vertexes_),
    vertex_outer_edges_per_vertex(vertex_outer_edges_per_vertex_),
    vertex_outer_edges_per_vertex_sign(vertex_outer_edges_per_vertex_sign_),
    volume_eq(volume_eq_),
    area_mean_eq(area_mean_eq_),
    edge_mean_eq(edge_mean_eq_),
    angle_mean_eq(angle_mean_eq_),
    inner_edge_list(inner_edge_list_),
    inner_edge_length_eq_list(inner_edge_length_eq_list_)
  {};

CommonCellConstants CommonCellConstants::CommonCellConstantsConstructor(HemoCellField & cellField_, Config & modelCfg_) {
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
    vector<T> edge_length_eq_list_;
    for (const hemo::Array<plint,2> & edge : edge_list_) {
      edge_length_eq_list_.push_back(cellField.meshElement.computeEdgeLength(edge[0],edge[1]));
    }

    //Calculate eq edges angles
    vector<T> edge_angle_eq_list_;
    hemo::Array<T,3> x2 = {0.0,0.0,0.0};

    for (const hemo::Array<plint,2> & edge : edge_list_) {
      const vector<plint> adjacentTriangles = cellField.meshElement.getAdjacentTriangleIds(edge[0], edge[1]);
      //const hemo::Array<T,3> x1 = cellField.meshElement.getVertex(edge[0]);

      for (pluint id = 0 ; id < 3 ; id ++ ) {
        const plint kVertex = cellField.meshElement.getVertexId(adjacentTriangles[0],id);
        if (kVertex != edge[0] && kVertex != edge[1]) {
          x2 = cellField.meshElement.getVertex(kVertex);
          break;
        }
      }

      const hemo::Array<T,3> V1 = cellField.meshElement.computeTriangleNormal(adjacentTriangles[0]);
      const hemo::Array<T,3> V2 = cellField.meshElement.computeTriangleNormal(adjacentTriangles[1]);

      // T angle = angleBetweenVectors(V1, V2);

      // ///TODO does this gives us a correct "sign"?
      // const plint sign = dot(x2-x1, V2) >= 0 ? 1 : -1;
      // if (sign <= 0) {
      //   angle = 2 * PI - angle;
      // }
      
      // if (angle > PI) {
      //   angle = angle - 2*PI;
      // }
      const hemo::Array<T,3> p0 = cellField.meshElement.getVertex(edge[0]);
      const hemo::Array<T,3> p1 = cellField.meshElement.getVertex(edge[1]);

      const hemo::Array<T,3> edge_vec = p1-p0;
      const T edge_length = norm(edge_vec);
      const hemo::Array<T,3> edge_uv = edge_vec/edge_length;
      T angle = getAngleBetweenFaces(V1, V2, edge_uv);
      
      edge_angle_eq_list_.push_back(angle);
    }

    // Define opposite points TODO: make it automatic / or add the 33 links version ;) 
    vector<hemo::Array<plint,2>> inner_edge_list_; 
    try {
      int v1 = 0, v2 = 0;
      tinyxml2::XMLElement * ie = modelCfg_["MaterialModel"]["InnerEdges"].getOrig();
      for (tinyxml2::XMLElement* edge = ie->FirstChildElement(); edge != NULL; edge = edge->NextSiblingElement())
      {
        if (sscanf(edge->GetText(), "%d %d", &v1, &v2) != 2 ) {
         pcout << "Inner Edges not read, somethings wrong" << endl; 
        }
        // read out integers
        inner_edge_list_.push_back({v1,v2});
      }
    } catch (std::invalid_argument & exeption) {}

    //Calculate eq edges lengths
    vector<T> inner_edge_length_eq_list_;
    for (const hemo::Array<plint,2> & edge : inner_edge_list_) {
      inner_edge_length_eq_list_.push_back(cellField.meshElement.computeEdgeLength(edge[0],edge[1]));
    }

    //Calculate triangle eq
    vector<T> triangle_area_eq_list_;
    for (pluint iTriangle = 0 ; iTriangle < cellField.triangle_list.size(); iTriangle++) {
      triangle_area_eq_list_.push_back(cellField.meshElement.computeTriangleArea(iTriangle));
    }

    
    T volume_eq_ = MeshMetrics<T>(cellField.meshElement).getVolume();


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
    T mean_area_eq_ = 0.;
    for (const T & area : triangle_area_eq_list_) {
      mean_area_eq_ += area;
    }
    mean_area_eq_ /= triangle_area_eq_list_.size();
    
    // Calculate mean edge lengths
    T mean_edge_eq_ = 0;
    for(const T & edge : edge_length_eq_list_) {
      mean_edge_eq_ += edge;
    }
    mean_edge_eq_ /= edge_length_eq_list_.size();

    // Calculate mean edge angles
    T mean_angle_eq_ = 0;
    for(const T & angle : edge_angle_eq_list_) {
      mean_angle_eq_ += angle;
    }
    mean_angle_eq_ /= edge_angle_eq_list_.size();

    //Pre-init, making things easy
    vector<hemo::Array<plint,6>> vertex_vertexes_(cellField.meshElement.getNumVertices(),{-1,-1,-1,-1,-1,-1});
    //Calculate triangles neighbouring vertices
    for (const hemo::Array<plint,2> & edge:  edge_list_) {
      for (int k = 0 ; k < 6 ; k++) {
        if (vertex_vertexes_[edge[0]][k] == -1) {
          vertex_vertexes_[edge[0]][k] = edge[1];
          break;
        }
      }
     for (int k = 0 ; k < 6 ; k++) {
        if (vertex_vertexes_[edge[1]][k] == -1) {
          vertex_vertexes_[edge[1]][k] = edge[0];
          break;
        }
      }
    }
    
    //Number of neighbouring vertices per vertex
    vector<unsigned int> vertex_n_vertexes_(cellField.meshElement.getNumVertices(),0);
    for (unsigned int i = 0 ; i < vertex_vertexes_.size() ; i++) {
      for (unsigned int j = 0 ; j < vertex_vertexes_[i].size() ; j++ ) {
        if (vertex_vertexes_[i][j] != -1) {
          vertex_n_vertexes_[i]++;
        }
      }
    }
    
    
    //Sort vertexes in order around vertex
    for (unsigned int vertex = 0 ; vertex < vertex_vertexes_.size() ; vertex++) {
      plint n_vertex = vertex_vertexes_[vertex][0];
      plint nextVertex = -1;
      for (unsigned int n = 1 ; n < vertex_n_vertexes_[vertex] ; n++) {
        vector<plint> triangleIDS = cellField.meshElement.getAdjacentTriangleIds(vertex, n_vertex);
        plint v0 = cellField.meshElement.getVertexId(triangleIDS[0],0);
        plint v1 = cellField.meshElement.getVertexId(triangleIDS[0],1);
        plint v2 = cellField.meshElement.getVertexId(triangleIDS[0],2);
        plint w0 = cellField.meshElement.getVertexId(triangleIDS[1],0);
        plint w1 = cellField.meshElement.getVertexId(triangleIDS[1],1);
        plint w2 = cellField.meshElement.getVertexId(triangleIDS[1],2);

        if ((v0 == vertex && v1 == n_vertex) || (v1 == vertex && v2 == n_vertex) || (v2 == vertex && v0 == n_vertex)) {
          if (v0 != vertex && v0 != n_vertex) {
           nextVertex = v0;
          }
          if (v1 != vertex && v1 != n_vertex) {
           nextVertex = v1;
          }
          if (v2 != vertex && v2 != n_vertex) {
           nextVertex = v2;
          }
        } else {
          if ((w0 == vertex && w1 == n_vertex) || (w1 == vertex && w2 == n_vertex) || (w2 == vertex && w0 == n_vertex)) {
            if (w0 != vertex && w0 != n_vertex) {
             nextVertex = w0;
            }
            if (w1 != vertex && w1 != n_vertex) {
             nextVertex = w1;
            }
            if (w2 != vertex && w2 != n_vertex) {
             nextVertex = w2;
            }
          }
        }
        n_vertex = nextVertex;
        vertex_vertexes_[vertex][n] = n_vertex;
      }
    }
    
    
    // Calculate center point deviation for surface patches around each vertex
    vector<T> surface_patch_center_dist_eq_list_;
    for (long signed int i = 0 ; i < cellField.meshElement.getNumVertices() ; i++) {
      
      hemo::Array<T,3> vertexes_sum = {0.,0.,0.};
      for(unsigned int j = 0; j < vertex_n_vertexes_[i]; j++) 
        vertexes_sum += cellField.meshElement.getVertex(vertex_vertexes_[i][j]);  
      
      const hemo::Array<T,3> vertexes_middle = vertexes_sum/vertex_n_vertexes_[i];
      const hemo::Array<T,3> localVertex = cellField.meshElement.getVertex(i);

      const hemo::Array<T,3> dev_vect = vertexes_middle - localVertex;
      
      // Get the local surface normal
      hemo::Array<T,3> patch_normal = {0.,0.,0.};
      for(unsigned int j = 0; j < vertex_n_vertexes_[i]-1; j++) {
        hemo::Array<T,3> triangle_normal = crossProduct(hemo::Array<T,3> (cellField.meshElement.getVertex(vertex_vertexes_[i][j])) - localVertex, 
                                                             hemo::Array<T,3> (cellField.meshElement.getVertex(vertex_vertexes_[i][j+1])) - localVertex);
        triangle_normal /= norm(triangle_normal);  
        patch_normal += triangle_normal;                                                   
      }
      hemo::Array<T,3> triangle_normal = crossProduct(hemo::Array<T,3> (cellField.meshElement.getVertex(vertex_vertexes_[i][vertex_n_vertexes_[i]-1])) -localVertex, 
                                                           hemo::Array<T,3> (cellField.meshElement.getVertex(vertex_vertexes_[i][0])) -localVertex);
      triangle_normal /= norm(triangle_normal);
      patch_normal += triangle_normal;

      patch_normal /= norm(patch_normal);
              
      const T ndev = dot(patch_normal, dev_vect); // distance along patch normal

      surface_patch_center_dist_eq_list_.push_back(ndev);
    }

    
    //Get edges around vertex
    vector<hemo::Array<plint,6>> vertex_edges_(cellField.meshElement.getNumVertices(),{-1,-1,-1,-1,-1,-1});
    vector<hemo::Array<signed int,6>> vertex_edges_sign_(cellField.meshElement.getNumVertices(), {0,0,0,0,0,0});
    plint v0,v1;
    for (unsigned int i = 0 ; i < vertex_vertexes_.size() ; i++) {
      for (unsigned int j = 0 ; j < vertex_n_vertexes_[i] ; j++) {
        v0 = i;
        v1 = vertex_vertexes_[v0][j];
        for (unsigned int e = 0 ; e < edge_list_.size() ; e++) {
          if ( edge_list_[e][0] == v0 && edge_list_[e][1] == v1) {
            vertex_edges_[i][j] = e;
            vertex_edges_sign_[i][j] = 1;
            break;
          }
          if ( edge_list_[e][1] == v0 && edge_list_[e][0] == v1) {
            vertex_edges_[i][j] = e;
            vertex_edges_sign_[i][j] = -1;
            break;
          }
        }  
      }
    }
    
    vector<hemo::Array<hemo::Array<plint,2>,6>> vertex_outer_edges_per_vertex(cellField.meshElement.getNumVertices(),{{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}});
    vector<hemo::Array<hemo::Array<signed int, 2>,6>> vertex_outer_edges_per_vertex_sign(cellField.meshElement.getNumVertices());
    plint lower, upper;
    for (unsigned int i = 0 ; i < vertex_vertexes_.size() ; i++){
      for (unsigned int j = 0 ; j < vertex_n_vertexes_[i] ; j++){
        
        //first edge
        lower = vertex_vertexes_[i][j];
        upper = vertex_vertexes_[i][((j-1)+vertex_n_vertexes_[i])%vertex_n_vertexes_[i]];
        if (lower > upper) { 
          swap(lower,upper); 
          vertex_outer_edges_per_vertex_sign[i][j][0] = -1;
        } else {
          vertex_outer_edges_per_vertex_sign[i][j][0] = 1;
        }

        //Find corresponding edge
        for (unsigned int k = 0 ; k < edge_list_.size() ; k++) {
          if (edge_list_[k][0] == lower and edge_list_[k][1] == upper) {
          vertex_outer_edges_per_vertex[i][j][0] = k;
          break;
          }
        }
        
        //Second edge
        lower = vertex_vertexes_[i][j];
        upper = vertex_vertexes_[i][(j+1)%vertex_n_vertexes_[i]];
        if (lower > upper) { 
          swap(lower,upper); 
          vertex_outer_edges_per_vertex_sign[i][j][1] = -1;
        } else {
          vertex_outer_edges_per_vertex_sign[i][j][1] = 1;
        }

        //Find corresponding edge
        for (unsigned int k = 0 ; k < edge_list_.size() ; k++) {
          if (edge_list_[k][0] != lower) continue;
          if (edge_list_[k][1] != upper) continue;
          
          vertex_outer_edges_per_vertex[i][j][1] = k;
          break;
        }
        
      }
    }
    
    
    CommonCellConstants CCC(cellField_,
            triangle_list_,
            edge_list_,
            edge_length_eq_list_,
            edge_angle_eq_list_,
            surface_patch_center_dist_eq_list_,
            edge_bending_triangles_,
            edge_bending_triangles_outer_points_,
            triangle_area_eq_list_,
            vertex_vertexes_,
            vertex_edges_,
            vertex_edges_sign_,
            vertex_n_vertexes_,
            vertex_outer_edges_per_vertex,
            vertex_outer_edges_per_vertex_sign,
            volume_eq_,
            mean_area_eq_, 
            mean_edge_eq_, 
            mean_angle_eq_,
            inner_edge_list_,
            inner_edge_length_eq_list_);
    return CCC;
};

}