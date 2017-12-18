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
#ifndef MESH_METRICS_HH
#define MESH_METRICS_HH

#include "meshMetrics.h"
/*
 * Helper function, calculates the angle between -pi and pi.
 * The edge is iVertex-jVertex.
 */
template<typename T>
T calculateSignedAngle(TriangularSurfaceMesh<T> const& mesh, plint iVertex, plint jVertex, plint & kVertex, plint & lVertex) {
    hemo::Array<T,3> x1 = mesh.getVertex(iVertex), x2({0.,0.,0.}), x3({0.,0.,0.}), x4({0.,0.,0.});

    std::vector<plint> adjacentTriangles = mesh.getAdjacentTriangleIds(iVertex, jVertex);
	plint iTriangle=adjacentTriangles[0], jTriangle=adjacentTriangles[1];
    x3 = mesh.getVertex(jVertex);
    T foundVertices=0;
    for (pluint id = 0; id < 3; ++id) {
        kVertex = mesh.getVertexId(iTriangle,id);
        if ( (kVertex != iVertex) && (kVertex != jVertex) ) {
            x2 = mesh.getVertex(kVertex);
            foundVertices += 1;
            break;
        }
    }
    for (pluint id = 0; id < 3; ++id) {
        lVertex = mesh.getVertexId(jTriangle,id);
        if ( (lVertex != iVertex) && (lVertex != jVertex) ) {
            x4 = mesh.getVertex(lVertex);
            foundVertices += 1;
            break;
        }
    }
    PLB_ASSERT(foundVertices == 2); //Assert if some particles are outside of the domain

    hemo::Array<T,3> V1 = mesh.computeTriangleNormal(iTriangle);
    hemo::Array<T,3> V2 = mesh.computeTriangleNormal(jTriangle);
    T angle = angleBetweenVectors(V1, V2);
	plint sign = dot(x2-x1, V2) >= 0?1:-1;
	if (sign <= 0) {
		angle = 2*PI-angle;
	}
	angle = (angle > PI)?angle-2*PI:angle;
	return angle;
}


/*
 *  * Helper function, calculates the angle between -pi and pi
 *   */
template<typename T>
T calculateSignedAngle(TriangularSurfaceMesh<T> const& mesh, plint iVertex, plint jVertex) {
    plint kVertex, lVertex;
      return calculateSignedAngle(mesh, iVertex, jVertex, kVertex, lVertex);
}

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
void writeSurfaceMeshAsciiSTL(TriangularSurfaceMesh<T> const& mesh, std::string fname)
{
	T dx = 1;
    // Output only from one MPI process.

    FILE *fp = fopen(fname.c_str(), "w");
    PLB_ASSERT(fp != NULL);

    char fmt1[64] = "  facet normal ";
    char fmt2[64] = "      vertex ";
    if (sizeof(T) == sizeof(long double)) {
        strcat(fmt1, "% Le % Le % Le\n");
        strcat(fmt2, "% Le % Le % Le\n");
    }
    else if (sizeof(T) == sizeof(float) ||
             sizeof(T) == sizeof(double)) {
        strcat(fmt1, "% e % e % e\n");
        strcat(fmt2, "% e % e % e\n");
    }
    else {
        PLB_ASSERT(false);
    }

    fprintf(fp, "solid surface\n");
    for (plint i = 0; i < mesh.getNumTriangles(); i++) {
        hemo::Array<T,3> n = mesh.computeTriangleNormal(i);
        hemo::Array<T,3> v;
        fprintf(fp, fmt1, n[0], n[1], n[2]);
        fprintf(fp, "    outer loop\n");
        v = dx * mesh.getVertex(i, 0);
        fprintf(fp, fmt2, v[0], v[1], v[2]);
        v = dx * mesh.getVertex(i, 1);
        fprintf(fp, fmt2, v[0], v[1], v[2]);
        v = dx * mesh.getVertex(i, 2);
        fprintf(fp, fmt2, v[0], v[1], v[2]);
        fprintf(fp, "    endloop\n");
        fprintf(fp, "  endfacet\n");
    }
    fprintf(fp, "endsolid surface\n");

    fclose(fp);
}


#endif  // MESH_METRICS_HH
