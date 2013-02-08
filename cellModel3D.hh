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
#include <map>
#include "cellModel3D.h"


#ifndef KBT__
#define KBT__
const double kB_p = 1.3806503e-23; // In SI, m2 kg s-2 K-1 (or J/K)
const double kBT_p = 4.100531391e-21; // In SI, m2 kg s-2 (or J) for T=300
double kB=0, kBT=0;
double dNewton=0;
#endif  // KBT__

#ifndef PI__
#define PI__
const double pi = 4.*atan(1.);
#endif  // PI__


namespace plb {

template<typename T>
CellModel3D<T>::CellModel3D (
        T density_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
        T k_volume_, T k_surface_,
        T eqArea_, T eqLength_, T eqAngle_,
        T eqVolume_, T eqSurface_, T eqTileSpan_,
        T eta_m_,
        T maxLength_, T persistenceLength_)
    : ShellModel3D<T>(density_),
      k_shear(k_shear_),
      k_bend(k_bend_),
      k_stretch(k_stretch_),
      k_elastic(k_elastic_),
      k_surface(k_surface_),
      k_volume(k_volume_),
      eqLength(eqLength_),
      eqArea(eqArea_),
      eqAngle(eqAngle_),
      eqVolume(eqVolume_),
      eqSurface(eqSurface_),
      eqTileSpan(eqTileSpan_),
      eta_m(eta_m_),
      maxLength(maxLength_),
      persistenceLength(persistenceLength_)
{
    T x0 = eqLength*1.0/maxLength;

//    pcout << "K_surf " << k_surface* kBT *1.0/ (eqLength*eqLength) << std::endl;
//    k_surface = k_surface* kBT / (eqLength*eqLength);
//    k_volume = k_volume* kBT / (eqLength*eqLength*eqLength);

    gamma_T = eta_m * 4 * sqrt(3.0) / 13.0;
    gamma_C = gamma_T/3.0;

    k_WLC = k_WLC_ * kBT * maxLength/(4.0*persistenceLength);
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
 /* Force calculation according to KrugerThesis, Appendix C */
{
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    Array<T,3> x1 = dynMesh.getVertex(iVertex), x2, x3, x4;

    T volumeCoefficient = k_volume * (cellVolume - eqVolume)*1.0/eqVolume;
    T surfaceCoefficient= k_surface * (cellSurface - eqSurface)*1.0/eqSurface;
    /* Force initializations */
    Array<T,3> inPlaneForce; inPlaneForce.resetToZero();
    Array<T,3> elasticForce; elasticForce.resetToZero();
    Array<T,3> bendingForce; bendingForce.resetToZero();
    Array<T,3> surfaceForce; surfaceForce.resetToZero();
    Array<T,3> shearForce; shearForce.resetToZero();
    Array<T,3> volumeForce; volumeForce.resetToZero();


    Array<T,3> dAdx, dVdx;
    Array<T,3> dl, eij; dl.resetToZero(); eij.resetToZero();
    Array<T,3> triangleNormal;
    std::map<plint, T> trianglesArea;
    std::map<plint, Array<T,3> > trianglesNormal;

    plint iTriangle, iX1, iX2, iX3;
    T r;
    Array<T,3> xi[3], tmp;
    plint jTriangle;
    /* Run through all the neighboring faces of iVertex and calculate:
     *
     x Volume conservation force
     x Surface conservation force
     x Elastic  WLC force
     *
     */
    std::vector<plint> neighborTriangles = dynMesh.getNeighborTriangleIds(iVertex);
    for (pluint iB = 0; iB < neighborTriangles.size(); ++iB) {
        iTriangle = neighborTriangles[iB];
        for (pluint id = 0; id < 3; ++id) { xi[id] = dynMesh.getVertex(iTriangle,id); }
        if (xi[0] ==x1) { iX1 = 0;}
        else if (xi[1] ==x1) { iX1 = 1;}
        else { iX1 = 2;}
        iX2 = (iX1 + 1)%3;
        iX3 = (iX1 + 2)%3;
        trianglesArea[iTriangle] = dynMesh.computeTriangleArea(iTriangle);
        triangleNormal = trianglesNormal[iTriangle] = dynMesh.computeTriangleNormal(iTriangle);
        /* Surface conservation force */
        crossProduct(triangleNormal,
                    xi[iX3] - xi[iX2],
                    tmp);
        dAdx = 0.5 *tmp;
        surfaceForce += -surfaceCoefficient * dAdx;
        shearForce += - k_shear * (trianglesArea[iTriangle] - eqArea)*1.0/eqArea * dAdx;
        // /* Elastice Force */
        // elasticForce += - (C_elastic*1.0)/(trianglesArea[iTriangle]*trianglesArea[iTriangle]) * dAdx;
        /* Volume conservation force */
        crossProduct(xi[iX3], xi[iX2], tmp);
        dVdx = 1.0/6.0 * tmp;
        volumeForce  += -volumeCoefficient  * dVdx;
    }

    /* Run through all the neighboring vertices of iVertex and calculate:
     *
     x In plane force
     x Bending force
     *
     */
    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
    for (pluint jV = 0; jV < neighborVertexIds.size(); jV++) {
        plint jVertex = neighborVertexIds[jV];
        x3 = dynMesh.getVertex(jVertex);
        /* In Plane Force (WLC) */
        dl = (x3 - x1)*1.0;
        r = norm(dl)/maxLength;
        eij = dl/(r*maxLength);

        inPlaneForce += k_WLC / maxLength * eij * (-1.0 +  4.0*r + 1.0/((1.0-r)*(1.0-r)) );
        inPlaneForce += 0.0;

        /*  Bending Forces Calculations */
        std::vector<plint> triangles = dynMesh.getAdjacentTriangleIds(iVertex, jVertex);
        iTriangle = triangles[0];
        jTriangle = triangles[1];
        plint kVertex, lVertex;
        for (pluint id = 0; id < 3; ++id) {
        	kVertex = dynMesh.getVertexId(iTriangle,id);
        	if ( (kVertex != iVertex) && (kVertex != jVertex) ) {
            	x2 = dynMesh.getVertex(kVertex);
        		break;
        	}
        }
        for (pluint id = 0; id < 3; ++id) {
        	lVertex = dynMesh.getVertexId(jTriangle,id);
        	if ( (lVertex != iVertex) && (lVertex != jVertex) ) {
            	x4 = dynMesh.getVertex(lVertex);
        		break;
        	}
        }
        // Bending Force
        bendingForce += computeBendingForce (x1, x2, x3, x4,
                            trianglesNormal[iTriangle], trianglesNormal[jTriangle],
                            trianglesArea[iTriangle], trianglesArea[jTriangle],
                            eqAngle, k_bend);
    }
    return volumeForce + surfaceForce + shearForce + inPlaneForce + elasticForce + bendingForce;
}


template<typename T>
Array<T,3> CellModel3D<T>::computeElasticForce (
        TriangleBoundary3D<T> const& boundary,
        plint iVertex )
{
    TriangularSurfaceMesh<T> const& mesh = boundary.getMesh();

    Array<T,3> force; force.resetToZero();
    static T eps = (sizeof(T) == sizeof(float) ?
            100.0 * std::numeric_limits<T>::epsilon() :
            std::numeric_limits<float>::epsilon());

    Array<T,3> vertex = mesh.getVertex(iVertex);

    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, mesh,
                k_WLC, k_elastic, k_shear, k_bend,
                eqLength, maxLength, eqArea, eqAngle, eqTileSpan);
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, mesh,
                k_WLC, k_elastic, k_shear, k_bend,
                eqLength, maxLength, eqArea, eqAngle, eqTileSpan);
        force[i] = -(up-um) / (2.0*eps);
    }
    return force;
}


template<typename T>
CellModel3D<T>* CellModel3D<T>::clone() const {
    return new CellModel3D<T>(*this);
}


template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj, T eqAngle, T k)
{
    Array<T,3> dAngledx; dAngledx.resetToZero();
    Array<T,3> dninjdx; dninjdx.resetToZero();
    Array<T,3> tmp; tmp.resetToZero();
    T angle = angleBetweenVectors(ni,nj);
    T ninj = dot(ni, nj);
    crossProduct(
            x2-x3,
            nj - ninj*ni,
            tmp);
    dninjdx = 0.5/Ai*tmp;
    crossProduct(
            x3-x4,
            ni - ninj*nj,
            tmp);
    dninjdx += 0.5/Aj*tmp;
    dAngledx = -1/sqrt(1 - ninj*ninj) * dninjdx;
    return -k*(angle - eqAngle)*dAngledx;
}


namespace shellModelHelper3D {

namespace cellModelHelper3D {

template<typename T>
T computePotential(plint iVertex, Array<T,3> const& iPosition,
                   TriangularSurfaceMesh<T> const& mesh,
                   T k_WLC, T k_elastic, T k_shear, T k_bend,
                   T eqLength, T maxLength, T eqArea, T eqAngle, T eqTileSpan
				   )
{
    T u = 0.0;

    // Membrane In Plane mode

    std::vector<plint> neighborVertexIds = mesh.getNeighborVertexIds(iVertex);

    pluint sz = neighborVertexIds.size();

    for (pluint i = 0; i < sz; i++) {
        plint jVertex = neighborVertexIds[i];
        u += computeInPlanePotential(iPosition, mesh.getVertex(jVertex),
        		k_WLC, maxLength);
    }

    pluint iMax = (mesh.isBoundaryVertex(iVertex) ? sz-1 : sz);

    // Membrane shearing mode

    for (pluint i = 0; i < iMax; i++) {
        plint jVertex = neighborVertexIds[i];
        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
        u += computeShearPotential(iPosition,
                                   mesh.getVertex(jVertex),
                                   mesh.getVertex(kVertex),
								   eqArea, k_shear);
    }



    // Cell bending mode

    for (pluint i = 0; i < iMax; i++) {
        plint jVertex = neighborVertexIds[i];
        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
        if (mesh.isInteriorEdge(iVertex, jVertex)) {
            std::vector<plint> v = mesh.getNeighborVertexIds(iVertex, jVertex);
            PLB_ASSERT(v.size() == 2);
            PLB_ASSERT(kVertex == v[0] || kVertex == v[1]);
            plint lVertex = (kVertex==v[0] ? v[1] : v[0]);
            u += computeBendPotential(mesh.getVertex(kVertex),
                                      iPosition,
                                      mesh.getVertex(jVertex),
                                      mesh.getVertex(lVertex),
                                      k_bend, eqAngle, eqLength, eqTileSpan);
        }
    }

    return u;
}

template<typename T>
T computeInPlanePotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                          T k, T maxLength)
{
    T x   = norm(iPosition   - jPosition) / maxLength;
    return k * (x*x) * (3.0 - 2 *x) / (1 - x);
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
}

template<typename T>
T computeBendPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                       Array<T,3> const& kPosition, Array<T,3> const& lPosition,
                       T k, T eqAngle, T eqLength, T eqTileSpan)
{
    // The position arrays specify two oriented triangles, namely
    // (i, j, k) and (l, k, j). These triangles share
    // the common edge j-k.

    // Compute quantities related to the current configuration
    Array<T,3> ijEdge = jPosition - iPosition;
    Array<T,3> jkEdge = kPosition - jPosition;
    Array<T,3> kjEdge = -jkEdge;
    Array<T,3> lkEdge = kPosition - lPosition;
	Array<T,3> nijk, nlkj;

    crossProduct(ijEdge, jkEdge, nijk);

    crossProduct(lkEdge, kjEdge, nlkj);

    T angle = angleBetweenVectors(nijk, nlkj);

    return 0.5*k*(angle - eqAngle)*(angle - eqAngle)* eqLength / eqTileSpan;
}

}  // namespace cellModelHelper3D


}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_HH
