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

namespace plb {
  
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

}
#endif  // MESH_METRICS_HH
