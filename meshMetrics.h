#ifndef MESH_METRICS_H
#define MESH_METRICS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <map>
#include "computeCellForces3D.h"


template<typename T>
class MeshMetrics
{
public:
    MeshMetrics(TriangleBoundary3D<T> const& Cells);
    ~MeshMetrics();
    void write(plb_ofstream & meshFile);
    void write() { this->write(meshQualityFile); } ;
    void set_dx(T dx_) { dx = dx_; } ;
    void set_dt(T dt_) { dt = dt_; } ;
    void set_dm(T dm_) { dm = dm_; } ;
    void set_dxdtdm(T dx_, T dt_, T dm_) { dx = dx_; dt = dt_; dm = dm_;} ;
private:
    T Nv, Nt, Nn, Nn6, Nn5, Nn7;
    T area, length, angle;
    T sigmaArea, sigmaLength, sigmaAngle, sigmaNn;
    T minArea, minLength, minAngle, minNn;
    T maxArea, maxLength, maxAngle, maxNn;
    T dx,dt,dm;
    plb_ofstream meshQualityFile;
};





#include "meshMetrics.hh"
#endif  // MESH_METRICS_H
