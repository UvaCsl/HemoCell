#ifndef MESH_METRICS_H
#define MESH_METRICS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <map>
#include "computeCellForces3D.h"
#include <vector>

using namespace std;
using namespace plb;

template<typename T>
class ElementsOfTriangularSurfaceMesh {
public:
    std::vector<Array<T,3> > vertexList;
    std::vector<plint> emanatingEdgeList;
    std::vector<Edge> edgeList;
};

template<typename T>
TriangularSurfaceMesh<T> * copyTriangularSurfaceMesh(TriangularSurfaceMesh<T> const& mesh, ElementsOfTriangularSurfaceMesh<T> & emptyEoTSM) {
    emptyEoTSM.vertexList = mesh.vertices();
    emptyEoTSM.emanatingEdgeList = mesh.emanatingEdges();
    emptyEoTSM.edgeList= mesh.edges();
    TriangularSurfaceMesh<T> * newMesh = new TriangularSurfaceMesh<T>(emptyEoTSM.vertexList, emptyEoTSM.emanatingEdgeList, emptyEoTSM.edgeList);
    return newMesh;
}


template<typename T>
class MeshMetrics
{
public:
    MeshMetrics(MeshMetrics<T> const& rhs);
    MeshMetrics(TriangleBoundary3D<T> const& Cells);
    MeshMetrics(TriangularSurfaceMesh<T>  const& mesh_);
    ~MeshMetrics();
    void write(plb_ofstream & meshFile);
    void init();
    void write() { plb_ofstream meshQualityFile((global::directories().getLogOutDir() + "plbMeshQuality.log").c_str());  this->write(meshQualityFile); } ;
    void set_dx(T dx_) { dx = dx_; } ;
    void set_dt(T dt_) { dt = dt_; } ;
    void set_dm(T dm_) { dm = dm_; } ;
    void set_dxdtdm(T dx_, T dt_, T dm_) { dx = dx_; dt = dt_; dm = dm_;} ;
    TriangularSurfaceMesh<T>  const& getMesh() { return mesh; }
    T getMeanLength() { return length; }
    T getMaxLength() { return maxLength; }
    T getMinLength() { return minLength; }

    T getMeanAngle() { return angle; }
    T getMaxAngle() { return maxAngle; }
    T getMinAngle() { return minAngle; }

    T getMeanArea() { return area; }
    T getMaxArea() { return maxArea; }
    T getMinArea() { return minArea; }
    // Computed as the maximum dimension from the BoundingBox
    T getRadius() { return cellRadius; }
    T getSurface() { return Nt*area; }
    T getVolume() { return volume; }
    Array<T,3> getMeanVertexPosition() { return meanVertexPosition; }

    T getNumVertices() { return Nv; }
    T getNumTriangles() { return Nt; }

private:
    TriangularSurfaceMesh<T>  const& mesh;
    T Nv, Nt, Nn, Nn6, Nn5, Nn7;
    T area, length, angle, volume;
    Array<T,3> meanVertexPosition;
    T sigmaArea, sigmaLength, sigmaAngle, sigmaNn;
    T minArea, minLength, minAngle, minNn;
    T maxArea, maxLength, maxAngle, maxNn;
    T dx,dt,dm;
    T cellRadius;
};



template<typename T>
void writeSurfaceMeshAsciiSTL(TriangularSurfaceMesh<T> const& mesh, std::string fname);



#include "meshMetrics.hh"
#endif  // MESH_METRICS_H
