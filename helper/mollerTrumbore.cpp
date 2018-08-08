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
#include "mollerTrumbore.h"
namespace hemo {
int MollerTrumbore(const hemo::Array<double,3> v0, const hemo::Array<double,3> v1,
	  const hemo::Array<double,3> v2, hemo::Array<double, 3> rayVector, 
	  hemo::Array<plint, 3> latticeSite, const double EPSILON) {
  
    float det,invDet,u,v; // Some floats
    hemo::Array<double, 3> edge1, edge2, pvec, svec, qvec;

    // Define edges
    edge1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};  
    edge2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]}; 

    pvec = hemo::crossProduct(rayVector, edge2);
    det = hemo::dot(edge1, pvec);

    if (det > -EPSILON && det < EPSILON) {
      return 0;
    }

    invDet = 1/det;

    // Construct s
    svec = {latticeSite[0] - v0[0], latticeSite[1] - v0[1], latticeSite[2] - v0[2]};
    u = invDet*hemo::dot(svec, pvec);

    if (u < 0.0 || u > 1.0) {
      return 0;
    }

    qvec = hemo::crossProduct(svec, edge1);
    v = invDet * hemo::dot(rayVector, qvec);

    if (v < 0.0 || u + v > 1.0) {
      return 0;
    }

    // Check where the intersection point is
    float t = invDet*hemo::dot(edge2, qvec);

    if (t > EPSILON) {
      return 1;
    }

    return 0;
  }
}