/* This file is part of the Palabos library.
 * Copyright (C) 2009, 2010 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CELL_MODEL_3D_HH
#define CELL_MODEL_3D_HH

#include <cmath>
#include "cellModel3D.h"

namespace plb {

template<typename T>
CellModel3D<T>::CellModel3D (
        T density_, T k_stretch_, T k_shear_, T k_bend_,
        T eqArea_, T eqLength_, T eqAngle_)
    : ShellModel3D<T>(density_),
      k_stretch(k_stretch_),
      k_shear(k_shear_),
      k_bend(k_bend_),
      eqLength(eqLength_),
      eqArea(eqArea_),
      eqAngle(eqAngle_)
{ }


template<typename T>
Array<T,3> CellModel3D<T>::computeElasticForce (
        TriangleBoundary3D<T> const& boundary,
        plint iVertex )
{
    // Select dynamic, open mesh.
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();

    Array<T,3> force; force.resetToZero();
    static T eps = (sizeof(T) == sizeof(float) ?
            100.0 * std::numeric_limits<T>::epsilon() :
            std::numeric_limits<float>::epsilon());

    Array<T,3> vertex = dynMesh.getVertex(iVertex);

    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, dynMesh, k_stretch,
                k_shear, k_bend,
                eqArea, eqLength, eqAngle);
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, dynMesh, k_stretch,
                k_shear, k_bend,
                eqArea, eqLength, eqAngle);
        force[i] = -(up-um) / (2.0*eps);
    }
    return force;
}

template<typename T>
CellModel3D<T>* CellModel3D<T>::clone() const {
    return new CellModel3D<T>(*this);
}


namespace shellModelHelper3D {

namespace cellModelHelper3D {

template<typename T>
T computePotential(plint iVertex, Array<T,3> const& iPosition,
                   TriangularSurfaceMesh<T> const& dynMesh, 
                   T k_stretch, T k_shear, T k_bend,
                   T eqArea, T eqLength, T eqAngle)
{
    T u = 0.0;

    // Membrane streching mode

    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);

    pluint sz = neighborVertexIds.size();

    for (pluint i = 0; i < sz; i++) {
        plint jVertex = neighborVertexIds[i];
        u += computeStretchPotential(iPosition, dynMesh.getVertex(jVertex),
                eqLength, k_stretch);
    }

    // Membrane shearing mode

    pluint iMax = (dynMesh.isBoundaryVertex(iVertex) ? sz-1 : sz);

    for (pluint i = 0; i < iMax; i++) {
        plint jVertex = neighborVertexIds[i];
        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
        u += computeShearPotential(iPosition,
                                   dynMesh.getVertex(jVertex),
                                   dynMesh.getVertex(kVertex),
                                   eqArea, k_shear);
    }

    // Shell bending mode

    for (pluint i = 0; i < iMax; i++) {
        plint jVertex = neighborVertexIds[i];
        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
        if (dynMesh.isInteriorEdge(iVertex, jVertex)) {
            std::vector<plint> v = dynMesh.getNeighborVertexIds(iVertex, jVertex);
            PLB_ASSERT(v.size() == 2);
            PLB_ASSERT(kVertex == v[0] || kVertex == v[1]);
            plint lVertex = (kVertex==v[0] ? v[1] : v[0]);
            u += computeBendPotential(dynMesh.getVertex(kVertex),
                                      iPosition,
                                      dynMesh.getVertex(jVertex),
                                      dynMesh.getVertex(lVertex),
                                      eqAngle, eqLength, eqArea,
                                      k_bend);
        }
        if (dynMesh.isInteriorEdge(jVertex, kVertex)) {
            std::vector<plint> v = dynMesh.getNeighborVertexIds(jVertex, kVertex);
            PLB_ASSERT(v.size() == 2);
            PLB_ASSERT(iVertex == v[0] || iVertex == v[1]);
            plint lVertex = (iVertex==v[0] ? v[1] : v[0]);
            u += computeBendPotential(iPosition,
                                      dynMesh.getVertex(jVertex),
                                      dynMesh.getVertex(kVertex),
                                      dynMesh.getVertex(lVertex),
                                      eqAngle, eqLength, eqArea,
                                      k_bend);
        }
    }

    return u;
}

template<typename T>
T computeStretchPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                          T eqLength,
                          T k)
{
    T length   = norm(iPosition   - jPosition);
    return 0.5*k*(length - eqLength)*(length - eqLength);
}

template<typename T>
T computeShearPotential(Array<T,3> const& iPosition,
                        Array<T,3> const& jPosition,
                        Array<T,3> const& kPosition,
                        T eqArea,
                        T k)
{
    T area   = computeTriangleArea(iPosition,   jPosition,   kPosition);
    return 0.5*k*(area - eqArea)*(area - eqArea);
}

template<typename T>
T computeBendPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                       Array<T,3> const& kPosition, Array<T,3> const& lPosition,
                       T eqAngle, T eqLength, T eqArea,
                       T k)
{
    // The position arrays specify two oriented triangles, namely 
    // (i, j, k) and (l, k, j). These triangles share
    // the common edge j-k.

    // Compute quantities related to the equilibrium configuration

    Array<T,3> ijEdge ;
    Array<T,3> jkEdge ;
    Array<T,3> kjEdge ;
    Array<T,3> lkEdge ;

    Array<T,3> nijk;
    Array<T,3> nlkj;

    // Compute quantities related to the current configuration

    ijEdge = jPosition - iPosition;
    jkEdge = kPosition - jPosition;
    kjEdge = -jkEdge;
    lkEdge = kPosition - lPosition;

    crossProduct(ijEdge, jkEdge, nijk);

    crossProduct(lkEdge, kjEdge, nlkj);

    T angle = angleBetweenVectors(nijk, nlkj);
    T Dangle = angle - eqAngle;
    T Dangle2 = Dangle * Dangle ;
    return k*Dangle2*( 0.5 + Dangle2/24.0 + Dangle2*Dangle2/720.0) ;
    T he = (eqArea*2.0)/(3.0*eqLength);
    return 0.5*k*(angle - eqAngle)*(angle - eqAngle)*eqLength / he;
}

}  // namespace cellModelHelper3D


}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_HH
