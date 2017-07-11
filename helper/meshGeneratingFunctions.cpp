#include "meshGeneratingFunctions.h"

TriangularSurfaceMesh<double> * constructStringMeshFromConfig(Config & cfg) {
  vector<plb::Array<double,3>> * vertexList = new vector<plb::Array<double,3>>();
  vector<plint> * emanatingEdgesList = new vector<plint>();
  vector<Edge> * edgeList = new vector<Edge>();
  TriangularSurfaceMesh<double> * mesh = new TriangularSurfaceMesh<double>(*vertexList,*emanatingEdgesList,*edgeList);
  return mesh;
}
