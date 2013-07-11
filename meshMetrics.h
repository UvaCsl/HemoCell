#ifndef MESH_METRICS_H
#define MESH_METRICS_H

#include "palabos3D.h"
#include "palabos3D.hh"


template<typename T>
class MeshMetrics
{
public:
    MeshMetrics(TriangleBoundary3D<T> const& Cells);
    ~MeshMetrics() {};
    void write(plb_ofstream & meshFile) {
        meshFile << "# Deviation in %, defined as 100*sigma(l)/mean(l), sl =  0" << std::endl;

        meshFile << "Number of vertices, Nv =  " << Nv << std::endl;
        meshFile << "Number of triangles, Nt =  " << Nt << std::endl;
        meshFile << std::endl;
        meshFile << "Mean Area per face, A =  " << area << std::endl;
        meshFile << "Deviation of Area %, sA =  " << 100*sigmaArea/area << std::endl;
        meshFile << "Max Area of face, maxA =  " << maxArea << std::endl;
        meshFile << "Min Area of face, minA =  " << minArea << std::endl;
        meshFile << std::endl;
        meshFile << "Mean Length per vertex, L =  " << length << std::endl;
        meshFile << "Deviation of Length %, sL =  " << 100*sigmaLength/length << std::endl;
        meshFile << "Max Length of edge, maxA =  " << maxLength << std::endl;
        meshFile << "Min Length of edge, minA =  " << minLength << std::endl;
        meshFile << std::endl;
        meshFile << "Mean Angle per vertex, theta =  " << angle << std::endl;
        meshFile << "Deviation of Angle %, sTheta =  " << fabs(100*sigmaAngle/angle) << std::endl;
        meshFile << "Max Angle of edge, maxA =  " << maxAngle << std::endl;
        meshFile << "Min Angle of edge, minA =  " << minAngle << std::endl;
        meshFile << std::endl;
        meshFile << "Mean number of neighbours, Nn =  " << Nn << std::endl;
        meshFile << "Deviation for number of Neighbours %, sNn =  " << 100*sigmaNn/Nn << std::endl;
        meshFile << "Number of 5 neighbour in %, Nn5 =  " << 100*Nn5/Nv << std::endl;
        meshFile << "Number of 6 neighbour in %, Nn6 =  " << 100*Nn6/Nv << std::endl;
        meshFile << "Number of 7 neighbour in %, Nn7 =  " << 100*Nn7/Nv << std::endl;
    };
    void write() { this->write(*meshQualityFile); } ;
    void setResultFile(plb_ofstream & meshFile_) { meshQualityFile = &meshFile_; } ;
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
    plb_ofstream* meshQualityFile;
};





#include "meshMetrics.hh"

#endif  // MESH_METRICS_H
