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
#ifndef MESH_METRICS_H
#define MESH_METRICS_H

namespace plb{
template<typename T>
class MeshMetrics;
}

#include <map>
#include <vector>

#include "core/array.h"
#include "offLattice/triangularSurfaceMesh.hh"
#include "offLattice/triangleBoundary3D.hh"
#include "array.h"

namespace plb {
using namespace std;

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
    hemo::Array<T,3> getMeanVertexPosition() { return meanVertexPosition; }

    T getNumVertices() { return Nv; }
    T getNumTriangles() { return Nt; }

private:
    TriangularSurfaceMesh<T>  const& mesh;
    T Nv, Nt, Nn, Nn6, Nn5, Nn7;
    T area, length, angle, volume;
    hemo::Array<T,3> meanVertexPosition;
    T sigmaArea, sigmaLength, sigmaAngle, sigmaNn;
    T minArea, minLength, minAngle, minNn;
    T maxArea, maxLength, maxAngle, maxNn;
    T dx,dt,dm;
    T cellRadius;
};

template<typename T>
MeshMetrics<T>::MeshMetrics(MeshMetrics<T> const& rhs) : mesh(rhs.mesh)  {
    init();
}


template<typename T>
MeshMetrics<T>::MeshMetrics(TriangleBoundary3D<T> const& Cells) : mesh(Cells.getMesh())  {
    init();
}


template<typename T>
MeshMetrics<T>::MeshMetrics(TriangularSurfaceMesh<T> const& mesh_) : mesh(mesh_)  {
    init();
}

template<typename T>
void MeshMetrics<T>::init()    {
    minArea=std::numeric_limits<T>::max(); minLength=std::numeric_limits<T>::max();
    minAngle=std::numeric_limits<T>::max(); minNn=std::numeric_limits<T>::max();
    maxArea=0; maxLength=0; maxAngle=0; maxNn=0;
    meanVertexPosition.resetToZero();

    Nv = mesh.getNumVertices();
    Nt = mesh.getNumTriangles();
    plb::Array<T,2> xRange;
    plb::Array<T,2> yRange;
    plb::Array<T,2> zRange;
    mesh.computeBoundingBox (xRange, yRange, zRange);
    cellRadius = max(xRange[1] - xRange[0], yRange[1] - yRange[0]);
    cellRadius = max(cellRadius , zRange[1] - zRange[0]) * 0.5;

    Nn=0; Nn6=0; Nn5=0; Nn7=0;
    area=0; length=0; angle=0;
    T varArea=0, varLength=0, varAngle=0, varNn=0;
    T tmp;
    // Vertices
    pluint NEdges = 0;
    volume = 0.0;
    for (int iV = 0; iV < Nv; ++iV) {
        meanVertexPosition += mesh.getVertex(iV);
        std::vector<plint> nvid = mesh.getNeighborVertexIds(iV);
        T NumNeighbors = nvid.size();
        Nn += nvid.size();
        minNn = minNn>NumNeighbors?NumNeighbors:minNn;
        maxNn = maxNn<NumNeighbors?NumNeighbors:maxNn;
        for (int ijV = 0; ijV < NumNeighbors; ++ijV) {
            int jV = nvid[ijV];
            tmp = mesh.computeEdgeLength(iV, jV);
            minLength = minLength>tmp?tmp:minLength;
            maxLength = maxLength<tmp?tmp:maxLength;
            length += tmp;
            tmp = calculateSignedAngle(mesh, iV, jV);
            minAngle = minAngle>tmp?tmp:minAngle;
            maxAngle = maxAngle<tmp?tmp:maxAngle;
            angle += tmp;
            NEdges++;
        }
        if (NumNeighbors == 5) {
            Nn5 += 1.0;
        } else if (NumNeighbors == 6) {
            Nn6 += 1.0;
        } else if (NumNeighbors == 7) {
            Nn7 += 1.0;
        }
        std::vector<plint> neighborTriangleIds = mesh.getNeighborTriangleIds(iV);
        for (pluint iB = 0; iB < neighborTriangleIds.size(); ++iB) {
            plint iTriangle = neighborTriangleIds[iB];
            hemo::Array<T,3> v0 = mesh.getVertex(iTriangle, 0);
            hemo::Array<T,3> v1 = mesh.getVertex(iTriangle, 1);
            hemo::Array<T,3> v2 = mesh.getVertex(iTriangle, 2);
            hemo::Array<T,3> tmp;
            crossProduct(v1, v2, tmp);
            T triangleVolumeT6 =  dot(v0,tmp);
            volume += triangleVolumeT6/6.0/3.0; // every volume is evaluated 3 times
        }
    }
    Nn /= Nv; length /= NEdges; angle /= NEdges;
    meanVertexPosition /= Nv;

// Compute vars
    for (int iV = 0; iV < Nv; ++iV) {
        std::vector<plint> nvid = mesh.getNeighborVertexIds(iV);
        T NumNeighbors = nvid.size();
        tmp = nvid.size()-Nn;
        varNn += tmp*tmp;
        for (int ijV = 0; ijV < NumNeighbors; ++ijV) {
            int jV = nvid[ijV];
            tmp=(mesh.computeEdgeLength(iV, jV)-length);
            varLength += tmp*tmp;
            tmp = (calculateSignedAngle(mesh, iV, jV) - angle);
            varAngle += tmp*tmp;
        }
    }
// Triangles
    area=0;
    for (int iT = 0; iT < Nt; ++iT) {
        tmp = mesh.computeTriangleArea(iT);
        minArea = minArea>tmp?tmp:minArea;
        maxArea = maxArea<tmp?tmp:maxArea;
        area += tmp;
    }
    area /= Nt; 
    for (int iT = 0; iT < Nt; ++iT) {
        tmp = mesh.computeTriangleArea(iT) - area;
        varArea += tmp*tmp;
    }
    sigmaNn = sqrt(varNn/Nv);
    sigmaArea = sqrt(varArea/Nt);
    sigmaLength = sqrt(varLength/NEdges);
    sigmaAngle = sqrt(varAngle/NEdges);
}

template<typename T>
MeshMetrics<T>::~MeshMetrics() {

};

template<typename T>
void MeshMetrics<T>::write(plb_ofstream & meshFile) {
    meshFile << "# Deviation in %, defined as 100*sigma(l)/mean(l), sl =  0" << std::endl;

    meshFile << "Number of vertices, Nv =  " << Nv << std::endl;
    meshFile << "Number of triangles, Nt =  " << Nt << std::endl;
    meshFile << "Surface, S =  " << getSurface() << std::endl;
    meshFile << "Volume, V =  " << getVolume() << std::endl;
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

template<typename T>
void writeSurfaceMeshAsciiSTL(TriangularSurfaceMesh<T> const& mesh, std::string fname);


}
#include "meshMetrics.hh"

#endif  // MESH_METRICS_H
