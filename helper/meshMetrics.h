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
class MeshMetrics
{
public:
    MeshMetrics(TriangleBoundary3D<T> const& Cells);
    MeshMetrics(TriangularSurfaceMesh<T>  const& mesh);
    ~MeshMetrics();
    void write(plb_ofstream & meshFile);
    void init(TriangularSurfaceMesh<T>  const& mesh);
    void write() { this->write(meshQualityFile); } ;
    void set_dx(T dx_) { dx = dx_; } ;
    void set_dt(T dt_) { dt = dt_; } ;
    void set_dm(T dm_) { dm = dm_; } ;
    void set_dxdtdm(T dx_, T dt_, T dm_) { dx = dx_; dt = dt_; dm = dm_;} ;
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
    T getVolume() { return Nt*area; }

    T getNumVertices() { return Nv; }
    T getNumTriangles() { return Nt; }

private:
    T Nv, Nt, Nn, Nn6, Nn5, Nn7;
    T area, length, angle, volume;
    T sigmaArea, sigmaLength, sigmaAngle, sigmaNn;
    T minArea, minLength, minAngle, minNn;
    T maxArea, maxLength, maxAngle, maxNn;
    T dx,dt,dm;
    T cellRadius;
    plb_ofstream meshQualityFile;
};





#include "meshMetrics.hh"
#endif  // MESH_METRICS_H
