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


#ifndef KBT__
#define KBT__
const double kB_p = 1.3806503e-23; // In SI, m2 kg s-2 K-1 (or J/K)
const double kBT_p = 4.100531391e-21; // In SI, m2 kg s-2 (or J)
double kB, kBT;
#endif  // KBT__

namespace plb {

template<typename T>
CellModel3D<T>::CellModel3D (
        T density_, T k_shear_, T k_bend_, T k_WLC_, T k_elastic_,
        T k_volume_, T k_surface_,
        T eqArea_, T eqLength_, T eqAngle_,
        T eqVolume_, T eqSurface_,
        T maxLength_, T persistenceLength_)
    : ShellModel3D<T>(density_),
      eqLength(eqLength_),
      eqArea(eqArea_),
      eqAngle(eqAngle_),
      eqVolume(eqVolume_),
      eqSurface(eqSurface_),
      maxLength(maxLength_),
      persistenceLength(persistenceLength_)
{
    T x0 = eqLength*1.0/maxLength;

    k_bend = k_bend_  * kBT ;
    k_shear = k_shear_ * kBT / (eqLength*eqLength);
    k_surface = k_surface_* kBT / (eqLength*eqLength);
    k_volume = k_volume_* kBT / (eqLength*eqLength*eqLength);

    k_WLC = k_WLC_ * kBT;
    C_WLC = k_WLC * maxLength/(4.0*persistenceLength);
    k_elastic = k_elastic_ * kBT ;
    C_elastic = k_elastic * 3.0 * sqrt(3.0)*
        (maxLength*maxLength*maxLength)*
        (x0*x0*x0*x0)/(64.0*persistenceLength)*
        (4*x0*x0 - 9*x0 + 6) /
        (1-x0)*(1-x0);
}


template<typename T>
Array<T,3> CellModel3D<T>::computeCellForce (
        TriangleBoundary3D<T> const& boundary,
        T cellVolume, T cellSurface,
        plint iVertex )
{
    // Force calculation according to KrugerThesis, Appendix C
    Array<T,3> volumeForce; volumeForce.resetToZero();
    Array<T,3> surfaceForce; surfaceForce.resetToZero();
    T volumeCoefficient = k_volume * (cellVolume - eqVolume)*1.0/eqVolume;
    T surfaceCoefficient= k_surface * (cellSurface - eqSurface)*1.0/eqSurface;
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    Array<T,3> xi = dynMesh.getVertex(iVertex);
    std::vector<plint> neighbors = dynMesh.getNeighborTriangleIds(iVertex);

    Array<T,3> x[3], triangleNormal, tmp, dAdx, dVdx;
    plint iTriangle, localVertex, localVertexPlus1, localVertexPlus2;
    for (pluint iB = 0; iB < neighbors.size(); ++iB) {
        iTriangle = neighbors[iB];
        triangleNormal = dynMesh.computeTriangleNormal(iTriangle);
        x[0] = dynMesh.getVertex(iTriangle,0);
        x[1] = dynMesh.getVertex(iTriangle,1);
        x[2] = dynMesh.getVertex(iTriangle,2);
        localVertex =  plint( (xi==x[1])*1.0 + (xi==x[2])*2.0 );
        localVertexPlus1 = (localVertex + 1)%3;
        localVertexPlus2 = (localVertex + 2)%3;
        crossProduct(triangleNormal,
                        x[localVertexPlus2] - x[localVertexPlus1],
                        tmp);
        dAdx = 0.5 *tmp;
        crossProduct(x[localVertexPlus1], x[localVertexPlus2], tmp);
        dVdx = 1.0/6.0 * tmp;
        volumeForce  += -volumeCoefficient  * dVdx;
        surfaceForce += -surfaceCoefficient * dAdx;
    }
    return volumeForce + surfaceForce;
}
//    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
//    pluint sz = neighborVertexIds.size();
//    Array<T,3> dAdx, dVdx, tmp;
//    for (pluint i = 0; i < sz; i++) {
//        plint jVertex = neighborVertexIds[i];
//        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
//        Array<T,3> xj = dynMesh.getVertex(jVertex);
//        Array<T,3> xk = dynMesh.getVertex(kVertex);
//        Array<T,3> triangleNormal = dynMesh.computeTriangleNormal(iVertex, jVertex, kVertex);
//        crossProduct(triangleNormal,
//                        xk - xj,
//                        tmp);
//        dAdx = 0.5 *tmp;
//        crossProduct(xk, xj, tmp);
//        dVdx = 1.0/6.0 * tmp;
//        volumeForce  += -volumeCoefficient  * dVdx;
//        surfaceForce += -surfaceCoefficient * dAdx;
//    }
//    return volumeForce + surfaceForce;
//}


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
                iVertex, iPosition, dynMesh,
                k_shear, k_bend, C_WLC,
                maxLength, eqArea, eqLength, eqAngle, C_elastic);
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, dynMesh,
                k_shear, k_bend, C_WLC,
                maxLength, eqArea, eqLength, eqAngle, C_elastic);
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
                   T k_shear, T k_bend, T C_WLC,
                   T maxLength, T eqArea, T eqLength, T eqAngle, T C_elastic)
{
    T u = 0.0;
    T area = 0;

    // Membrane streching mode
    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
    plint sz = neighborVertexIds.size();
    for (plint i = 0; i < sz; i++) {
        plint jVertex = neighborVertexIds[i];
        plint kVertex = ((i+1)==sz ? neighborVertexIds[0] : neighborVertexIds[i+1]);
        plint lVertex = ((i-1)==-1 ? neighborVertexIds[sz-1] : neighborVertexIds[i-1]);

        Array<T,3> jPosition = dynMesh.getVertex(jVertex);
        Array<T,3> kPosition = dynMesh.getVertex(kVertex);
        Array<T,3> lPosition = dynMesh.getVertex(lVertex);
        // Membrane shearing mode
        u += computeShearPotential(iPosition, jPosition, kPosition,
                                   eqArea, area, k_shear) * 1.0/3.0 ; //  (x 1.0 / 3.0),it refers to one vertex of the triangle
        // In-Plane potential
        u +=  computeInPlanePotential(iPosition, jPosition, maxLength, C_WLC) * 0.5  ; // ( x 0.5),it refers to one vertex of the edge
        u += (C_elastic*1.0)/area * 1.0/3.0 ; //  (x 1.0 / 3.0),it refers to one vertex of the triangle
        u += computeBendPotential(lPosition,
                                  iPosition,
                                  jPosition,
                                  kPosition,
                                  eqAngle, eqLength, eqArea,
                                  k_bend) * 0.5;  // ( x 0.5),it refers to one vertex of the edge
    }
    return u;
}
//        if (dynMesh.isInteriorEdge(iVertex, jVertex)) {
//            std::vector<plint> v = dynMesh.getNeighborVertexIds(iVertex, jVertex);
//            PLB_ASSERT(v.size() == 2);
//            PLB_ASSERT(kVertex == v[0] || kVertex == v[1]);
//            plint lVertex = (kVertex==v[0] ? v[1] : v[0]);
//            u += computeBendPotential(dynMesh.getVertex(kVertex),
//                                      iPosition,
//                                      dynMesh.getVertex(jVertex),
//                                      dynMesh.getVertex(lVertex),
//                                      eqAngle, eqLength, eqArea,
//                                      k_bend);
//        }

template<typename T>
T computeInPlanePotential(Array<T,3> const& iPosition,
                          Array<T,3> const& jPosition,
                          T maxLength, T k)
{
    T x = norm(iPosition - jPosition)*1.0/maxLength;
    return k*(3*x*x - 2*x*x*x)/(1-x);
}


template<typename T>
T computeStretchPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                          T eqLength,
                          T k)
{
    T length   = norm(iPosition   - jPosition);
    return 0.5*k*(length - eqLength)*(length - eqLength)/eqLength;
//    return 0.5*k*(length - eqLength)*(length - eqLength);
}

template<typename T>
T computeShearPotential(Array<T,3> const& iPosition,
                        Array<T,3> const& jPosition,
                        Array<T,3> const& kPosition,
                        T eqArea,
                        T k)
{
    T area   = computeTriangleArea(iPosition,   jPosition,   kPosition);
    return 0.5*k*(area - eqArea)*(area - eqArea)/eqArea;
//    return 0.5*k*(area - eqArea)*(area - eqArea);
}

template<typename T>
T computeShearPotential(Array<T,3> const& iPosition,
                        Array<T,3> const& jPosition,
                        Array<T,3> const& kPosition,
                        T eqArea, T & area,
                        T k)
{
    area   = computeTriangleArea(iPosition,   jPosition,   kPosition);
    return 0.5*k*(area - eqArea)*(area - eqArea)/eqArea;
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
    //pcout << "angle: " << angle << std::endl;
    T Dangle = angle - eqAngle;
    T Dangle2 = Dangle * Dangle ;
//    return k*Dangle2*( 0.5 + Dangle2/24.0 + Dangle2*Dangle2/720.0) ;
    return k*( 1 - cos(Dangle) ) ;
//    T he = (eqArea*2.0)/(3.0*eqLength);
//    return 0.5*k*(angle - eqAngle)*(angle - eqAngle)*eqLength / he;
}

}  // namespace cellModelHelper3D


}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_HH
