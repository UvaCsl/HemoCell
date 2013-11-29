#ifndef MESH_METRICS_HH
#define MESH_METRICS_HH

#include "meshMetrics.h"


template<typename T>
MeshMetrics<T>::MeshMetrics(TriangleBoundary3D<T> const& Cells) {
    minArea=10000000; minLength=10000000; minAngle=10000000; minNn=10000000;
    maxArea=0; maxLength=0; maxAngle=0; maxNn=0;

    TriangularSurfaceMesh<T>  mesh = Cells.getMesh();
    Nv = mesh.getNumVertices();
    Nt = mesh.getNumTriangles();

    Nn=0; Nn6=0; Nn5=0; Nn7=0;
    area=0; length=0; angle=0;
    sigmaArea=0; sigmaLength=0; sigmaAngle=0; sigmaNn=0;
    T tmp;
    // Compute Mean Values
    pluint NEdges = 0;
    for (int iV = 0; iV < Nv; ++iV) {
        tmp = mesh.computeVertexArea(iV);
        minArea = minArea>tmp?tmp:minArea;
        maxArea = maxArea<tmp?tmp:maxArea;
        area += tmp;
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
            tmp = calculateSignedAngle(mesh, iV, jV) * 180/3.14159;
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
    }
    for (int iT = 0; iT < Nt; ++iT) {
        area += mesh.computeTriangleArea(iT);
    }
    Nn /= Nv; area /= Nt; length /= NEdges; angle /= NEdges;
    // Compute Sigmas
    for (int iV = 0; iV < Nv; ++iV) {
        tmp = mesh.computeVertexArea(iV) - area;
        sigmaArea += tmp*tmp;
        std::vector<plint> nvid = mesh.getNeighborVertexIds(iV);
        T NumNeighbors = nvid.size();
        tmp = nvid.size()-Nn;
        sigmaNn += tmp*tmp;
        for (int ijV = 0; ijV < NumNeighbors; ++ijV) {
            int jV = nvid[ijV];
            tmp=(mesh.computeEdgeLength(iV, jV)-length);
            sigmaLength += tmp*tmp;
            tmp = (calculateSignedAngle(mesh, iV, jV) * 180/3.14159 - angle);
            sigmaAngle += tmp*tmp;
        }
    }
    for (int iT = 0; iT < Nt; ++iT) {
        tmp = mesh.computeTriangleArea(iT) - area;
        sigmaArea += tmp*tmp;
    }
    sigmaNn = sqrt(sigmaNn/Nv);
    sigmaArea = sqrt(sigmaArea/Nt);
    sigmaLength = sqrt(sigmaLength/NEdges);
    sigmaAngle = sqrt(sigmaAngle/NEdges);
}





#endif  // MESH_METRICS_HH
