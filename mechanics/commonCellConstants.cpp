class CommonCellConstants {
private:
  CommonCellConstants(HemoCellField & cellField_,
                      vector<Array<plint,3>> triangle_list_,
                      vector<Array<plint,2>> edge_list_,
                      vector<double> edge_length_eq_list_,
                      vector<double> edge_angle_eq_list_,
                      vector<double> triangle_area_eq_list_,
                      double volume_eq_) :
    cellField(cellField_),
    triangle_list(triangle_list_),
    edge_list(edge_list_),
    edge_length_eq_list(edge_length_eq_list_),
    edge_angle_eq_list(edge_angle_eq_list_),
    triangle_area_eq_list(triangle_area_eq_list_),
    volume_eq(volume_eq_)
  {};

public:
  static CommonCellConstants CommonCellConstantsConstructor(HemoCellField & cellField_) {
    HemoCellField & cellField = cellField_;
    //Calculate triangles
    vector<Array<plint,3>> triangle_list_ = cellField.triangle_list;

    //Calculate edges
    vector<Array<plint,2>> edge_list_;
    const TriangularSurfaceMesh<double> constMesh = cellField.meshElement;
    for (const Edge & edge : constMesh.edges()) {
      edge_list_.push_back({edge.pv,edge.ne});
    }

    //Calculate eq edges
    vector<double> edge_length_eq_list_;
    for (const Array<plint,2> & edge : edge_list_) {
      edge_length_eq_list_.push_back(cellField.meshElement.computeEdgeLength(edge[0],edge[1]));
    }

    //Calculate eq edges angles
    vector<double> edge_angle_eq_list_;
    Array<double,3> x2;
    double angle;
    for (const Array<plint,2> & edge : edge_list_) {
      const vector<plint> adjacentTriangles = cellField.meshElement.getAdjacentTriangleIds(edge[0], edge[1]);
      const Array<double,3> x1 = cellField.meshElement.getVertex(edge[0]);

      for (pluint id = 0 ; id < 3 ; id ++ ) {
        const plint kVertex = cellField.meshElement.getVertexId(adjacentTriangles[0],id);
        if (kVertex != edge[0] && kVertex != edge[1]) {
          x2 = cellField.meshElement.getVertex(kVertex);
          break;
        }
      }

      const Array<double,3> V1 = cellField.meshElement.computeTriangleNormal(adjacentTriangles[0]);
      const Array<double,3> V2 = cellField.meshElement.computeTriangleNormal(adjacentTriangles[1]);

      angle = angleBetweenVectors(V1, V2);
      const plint sign = dot(x2-x1, V2) >= 0 ? 1 : -1;
      if (sign <= 0) {
        angle = 2 * pi - angle;
      }
      edge_angle_eq_list_.push_back((angle > pi) ? angle - 2 * pi : angle);
    }

    //Calculate triangle eq
    vector<double> triangle_area_eq_list_;
    for (pluint iTriangle ; iTriangle < cellField.triangle_list.size(); iTriangle++) {
      triangle_area_eq_list_.push_back(cellField.meshElement.computeTriangleArea(iTriangle));
    }

    
    double volume_eq_ = MeshMetrics<double>(cellField.meshElement).getVolume();

    CommonCellConstants CCC(cellField_,triangle_list_,edge_list_,edge_length_eq_list_,edge_angle_eq_list_,triangle_area_eq_list_,volume_eq_);
    return CCC;

    
  };

public:
  HemoCellField & cellField;
  const vector<Array<plint,3>> triangle_list;
  const vector<Array<plint,2>> edge_list;

  const vector<double> edge_length_eq_list;
  const vector<double> edge_angle_eq_list;
  const vector<double> triangle_area_eq_list;

  const double volume_eq;

};
