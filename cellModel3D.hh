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
        T density_, T k_rest_, T k_stretch_, T k_shear_, T k_bend_, T k_vol_, T k_surface_)
    : ShellModel3D<T>(density_),
      k_rest(k_rest_),
      k_stretch(k_stretch_),
      k_shear(k_shear_),
      k_bend(k_bend_),
      k_vol(k_vol_),
      k_surface(k_surface_)
{ }


template<typename T>
Array<T,3> CellModel3D<T>::computeElasticForce (
        TriangleBoundary3D<T> const& boundary,
        plint iVertex )
{
    // Select dynamic, open mesh.
    boundary.pushSelect(0,1);
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    boundary.popSelect();
    // Select static, open mesh.
    boundary.pushSelect(0,0);
    TriangularSurfaceMesh<T> const& eqMesh = boundary.getMesh();
    boundary.popSelect();

    Array<T,3> force; force.resetToZero();
    static T eps = (sizeof(T) == sizeof(float) ?
            100.0 * std::numeric_limits<T>::epsilon() :
            std::numeric_limits<float>::epsilon());

    Array<T,3> vertex = dynMesh.getVertex(iVertex);

    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, dynMesh, eqMesh, k_rest, k_stretch,
                k_shear, k_bend, k_vol, k_surface );
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, dynMesh, eqMesh, k_rest, k_stretch,
                k_shear, k_bend, k_vol, k_surface );
        force[i] = -(up-um) / (2.0*eps);
    }
    return force;
}

template<typename T>
CellModel3D<T>* CellModel3D<T>::clone() const {
    return new CellModel3D<T>(*this);
}

//template<typename T>
//MembraneBendingModel3D<T>::MembraneBendingModel3D (
//        T density_, T youngModulus_, T poissonRatio_, T thickness_, T k_rest_ )
//    : CellModel3D<T>(density_),
//      youngModulus(youngModulus_),
//      poissonRatio(poissonRatio_),
//      thickness(thickness_),
//      k_rest(k_rest_)
//{ }
//
//template<typename T>
//Array<T,3> MembraneBendingModel3D<T>::computeElasticForce (
//        TriangleBoundary3D<T> const& boundary,
//        plint iVertex )
//{
//    // Select dynamic, open mesh.
//    boundary.pushSelect(0,1);
//    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
//    boundary.popSelect();
//    // Select static, open mesh.
//    boundary.pushSelect(0,0);
//    TriangularSurfaceMesh<T> const& eqMesh = boundary.getMesh();
//    boundary.popSelect();
//
//    Array<T,3> force; force.resetToZero();
//    static T eps = (sizeof(T) == sizeof(float) ?
//            100.0 * std::numeric_limits<T>::epsilon() :
//            std::numeric_limits<float>::epsilon());
//
//    Array<T,3> vertex = dynMesh.getVertex(iVertex);
//
//    for (int i = 0; i < 3; i++) {
//        Array<T,3> iPosition = vertex;
//        iPosition[i] += eps;
//        T up = shellModelHelper3D::membraneBendingModelHelper3D::computePotential (
//                iVertex, iPosition, dynMesh, eqMesh, youngModulus,
//                poissonRatio, thickness, k_rest );
//        iPosition = vertex;
//        iPosition[i] -= eps;
//        T um = shellModelHelper3D::membraneBendingModelHelper3D::computePotential (
//                iVertex, iPosition, dynMesh, eqMesh, youngModulus,
//                poissonRatio, thickness, k_rest );
//        force[i] = -(up-um) / (2.0*eps);
//    }
//    return force;
//}
//
//template<typename T>
//MembraneBendingModel3D<T>* MembraneBendingModel3D<T>::clone() const {
//    return new MembraneBendingModel3D<T>(*this);
//}

namespace shellModelHelper3D {

namespace cellModelHelper3D {

template<typename T>
T computePotential(plint iVertex, Array<T,3> const& iPosition,
                   TriangularSurfaceMesh<T> const& dynMesh, 
                   TriangularSurfaceMesh<T> const& eqMesh, 
                   T k_rest, T k_stretch, T k_shear, T k_bend,
                   T k_vol, T k_surface)
{
    T u = 0.0;

    // Resting mode

    u += computeRestPotential(iPosition, eqMesh.getVertex(iVertex), k_rest);

    // Membrane streching mode

    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);

    pluint sz = neighborVertexIds.size();

    for (pluint i = 0; i < sz; i++) {
        plint jVertex = neighborVertexIds[i];
        u += computeStretchPotential(iPosition, dynMesh.getVertex(jVertex),
                                     eqMesh.getVertex(iVertex), eqMesh.getVertex(jVertex),
                                     k_stretch);
    }

    // Membrane shearing mode

    pluint iMax = (dynMesh.isBoundaryVertex(iVertex) ? sz-1 : sz);

    for (pluint i = 0; i < iMax; i++) {
        plint jVertex = neighborVertexIds[i];
        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
        u += computeShearPotential(iPosition,
                                   dynMesh.getVertex(jVertex),
                                   dynMesh.getVertex(kVertex),
                                   eqMesh.getVertex(iVertex),
                                   eqMesh.getVertex(jVertex),
                                   eqMesh.getVertex(kVertex), k_shear);
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
            T he = eqMesh.computeEdgeTileSpan(iVertex, jVertex);
            u += computeBendPotential(dynMesh.getVertex(kVertex),
                                      iPosition,
                                      dynMesh.getVertex(jVertex),
                                      dynMesh.getVertex(lVertex),
                                      eqMesh.getVertex(kVertex),
                                      eqMesh.getVertex(iVertex),
                                      eqMesh.getVertex(jVertex),
                                      eqMesh.getVertex(lVertex),
                                      k_bend, he);
        }
        if (dynMesh.isInteriorEdge(jVertex, kVertex)) {
            std::vector<plint> v = dynMesh.getNeighborVertexIds(jVertex, kVertex);
            PLB_ASSERT(v.size() == 2);
            PLB_ASSERT(iVertex == v[0] || iVertex == v[1]);
            plint lVertex = (iVertex==v[0] ? v[1] : v[0]);
            T he = eqMesh.computeEdgeTileSpan(jVertex, kVertex);
            u += computeBendPotential(iPosition,
                                      dynMesh.getVertex(jVertex),
                                      dynMesh.getVertex(kVertex),
                                      dynMesh.getVertex(lVertex),
                                      eqMesh.getVertex(iVertex),
                                      eqMesh.getVertex(jVertex),
                                      eqMesh.getVertex(kVertex),
                                      eqMesh.getVertex(lVertex),
                                      k_bend, he);
        }
    }

    return u;
}

template<typename T>
T computeStretchPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                          Array<T,3> const& iEqPosition, Array<T,3> const& jEqPosition,
                          T k)
{
    T length   = norm(iPosition   - jPosition);
    T eqLength = norm(iEqPosition - jEqPosition);
    return 0.5*k*(length - eqLength)*(length - eqLength);
}

template<typename T>
T computeShearPotential(Array<T,3> const& iPosition,
                        Array<T,3> const& jPosition,
                        Array<T,3> const& kPosition,
                        Array<T,3> const& iEqPosition,
                        Array<T,3> const& jEqPosition,
                        Array<T,3> const& kEqPosition,
                        T k)
{
    T area   = computeTriangleArea(iPosition,   jPosition,   kPosition);
    T eqArea = computeTriangleArea(iEqPosition, jEqPosition, kEqPosition);
    return 0.5*k*(area - eqArea)*(area - eqArea);
}

template<typename T>
T computeBendPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                       Array<T,3> const& kPosition, Array<T,3> const& lPosition,
                       Array<T,3> const& iEqPosition, Array<T,3> const& jEqPosition,
                       Array<T,3> const& kEqPosition, Array<T,3> const& lEqPosition,
                       T k, T he)
{
    // The position arrays specify two oriented triangles, namely 
    // (i, j, k) and (l, k, j). These triangles share
    // the common edge j-k.

    // Compute quantities related to the equilibrium configuration

    Array<T,3> ijEdge = jEqPosition - iEqPosition;
    Array<T,3> jkEdge = kEqPosition - jEqPosition;
    Array<T,3> kjEdge = -jkEdge;
    Array<T,3> lkEdge = kEqPosition - lEqPosition;

    Array<T,3> nijk;
    crossProduct(ijEdge, jkEdge, nijk);

    Array<T,3> nlkj;
    crossProduct(lkEdge, kjEdge, nlkj);

    T eqAngle = angleBetweenVectors(nijk, nlkj);
    T eqLength = norm(jkEdge);

    // Compute quantities related to the current configuration

    ijEdge = jPosition - iPosition;
    jkEdge = kPosition - jPosition;
    kjEdge = -jkEdge;
    lkEdge = kPosition - lPosition;

    crossProduct(ijEdge, jkEdge, nijk);

    crossProduct(lkEdge, kjEdge, nlkj);

    T angle = angleBetweenVectors(nijk, nlkj);

    return 0.5*k*(angle - eqAngle)*(angle - eqAngle)*eqLength / he;
}

}  // namespace cellModelHelper3D

//namespace membraneBendingModelHelper3D {
//
//template<typename T>
//T computePotential(plint iVertex, Array<T,3> const& iPosition,
//                   TriangularSurfaceMesh<T> const& dynMesh,
//                   TriangularSurfaceMesh<T> const& eqMesh,
//                   T Y, T nu, T h, T k)
//{
//    T u = 0.0;
//
//    // Resting mode
//
//    u += computeRestPotential(iPosition, eqMesh.getVertex(iVertex), k);
//
//    // Compute local vertex area
//
//    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
//    pluint sz = neighborVertexIds.size();
//    pluint iMax = (dynMesh.isBoundaryVertex(iVertex) ? sz-1 : sz);
//
//    T area = 0.0;
//    for (pluint i = 0; i < iMax; i++) {
//        plint jVertex = neighborVertexIds[i];
//        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
//        area += computeTriangleArea(iPosition, dynMesh.getVertex(jVertex),
//                                    dynMesh.getVertex(kVertex));
//    }
//    area /= 3.0;
//
//    // Membrane mode
//
//    T ud = 0.0;
//    Array<Array<T,3>,3> E;
//
//    for (pluint i = 0; i < iMax; i++) {
//        plint jVertex = neighborVertexIds[i];
//        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
//        computeMembraneStrain(iPosition, dynMesh.getVertex(jVertex), dynMesh.getVertex(kVertex),
//                              eqMesh.getVertex(iVertex),
//                              eqMesh.getVertex(jVertex),
//                              eqMesh.getVertex(kVertex), E);
//        ud += computeMembranePotentialDensity(E, Y, nu, h);
//    }
//
//    // Bending mode
//
//    for (pluint i = 0; i < iMax; i++) {
//        plint jVertex = neighborVertexIds[i];
//        plint kVertex = (i==sz-1 ? neighborVertexIds[0] : neighborVertexIds[i+1]);
//
//        Array<T,3> jPosition = dynMesh.getVertex(jVertex);
//        Array<T,3> kPosition = dynMesh.getVertex(kVertex);
//        Array<T,3> lPosition, mPosition, nPosition;
//
//        Array<T,3> const* iPositionP = &iPosition;
//        Array<T,3> const* jPositionP = &jPosition;
//        Array<T,3> const* kPositionP = &kPosition;
//        Array<T,3> *lPositionP, *mPositionP, *nPositionP;
//
//        Array<T,3> iEqPosition = eqMesh.getVertex(iVertex);
//        Array<T,3> jEqPosition = eqMesh.getVertex(jVertex);
//        Array<T,3> kEqPosition = eqMesh.getVertex(kVertex);
//        Array<T,3> lEqPosition, mEqPosition, nEqPosition;
//
//        Array<T,3> *iEqPositionP = &iEqPosition;
//        Array<T,3> *jEqPositionP = &jEqPosition;
//        Array<T,3> *kEqPositionP = &kEqPosition;
//        Array<T,3> *lEqPositionP, *mEqPositionP, *nEqPositionP;
//
//        if (dynMesh.isInteriorEdge(iVertex, jVertex)) {
//            std::vector<plint> v = dynMesh.getNeighborVertexIds(iVertex, jVertex);
//            PLB_ASSERT(v.size() == 2);
//            PLB_ASSERT(kVertex == v[0] || kVertex == v[1]);
//            plint nVertex = (kVertex==v[0] ? v[1] : v[0]);
//            nPosition = dynMesh.getVertex(nVertex);
//            nPositionP = &nPosition;
//            nEqPosition = eqMesh.getVertex(nVertex);
//            nEqPositionP = &nEqPosition;
//        }
//        else {
//            nPositionP = nEqPositionP = NULL;
//        }
//
//        if (dynMesh.isInteriorEdge(jVertex, kVertex)) {
//            std::vector<plint> v = dynMesh.getNeighborVertexIds(jVertex, kVertex);
//            PLB_ASSERT(v.size() == 2);
//            PLB_ASSERT(iVertex == v[0] || iVertex == v[1]);
//            plint lVertex = (iVertex==v[0] ? v[1] : v[0]);
//            lPosition = dynMesh.getVertex(lVertex);
//            lPositionP = &lPosition;
//            lEqPosition = eqMesh.getVertex(lVertex);
//            lEqPositionP = &lEqPosition;
//        }
//        else {
//            lPositionP = lEqPositionP = NULL;
//        }
//
//        if (dynMesh.isInteriorEdge(kVertex, iVertex)) {
//            std::vector<plint> v = dynMesh.getNeighborVertexIds(kVertex, iVertex);
//            PLB_ASSERT(v.size() == 2);
//            PLB_ASSERT(jVertex == v[0] || jVertex == v[1]);
//            plint mVertex = (jVertex==v[0] ? v[1] : v[0]);
//            mPosition = dynMesh.getVertex(mVertex);
//            mPositionP = &mPosition;
//            mEqPosition = eqMesh.getVertex(mVertex);
//            mEqPositionP = &mEqPosition;
//        }
//        else {
//            mPositionP = mEqPositionP = NULL;
//        }
//
//        computeBendingStrain(iPositionP, jPositionP, kPositionP,
//                             lPositionP, mPositionP, nPositionP,
//                             iEqPositionP, jEqPositionP, kEqPositionP,
//                             lEqPositionP, mEqPositionP, nEqPositionP,
//                             E);
//
//        ud += computeBendingPotentialDensity(E, Y, nu, h);
//
//        if (dynMesh.isInteriorEdge(jVertex, kVertex)) {
//            std::vector<plint> v = dynMesh.getNeighborVertexIds(jVertex, kVertex); // The next 8 lines are
//            PLB_ASSERT(v.size() == 2);                                             // not strictly needed!
//            PLB_ASSERT(iVertex == v[0] || iVertex == v[1]);
//            plint lVertex = (iVertex==v[0] ? v[1] : v[0]);
//            lPosition = dynMesh.getVertex(lVertex);
//            lPositionP = &lPosition;
//            lEqPosition = eqMesh.getVertex(lVertex);
//            lEqPositionP = &lEqPosition;
//
//            if (dynMesh.isInteriorEdge(jVertex, lVertex)) {
//                std::vector<plint> v = dynMesh.getNeighborVertexIds(jVertex, lVertex);
//                PLB_ASSERT(v.size() == 2);
//                PLB_ASSERT(kVertex == v[0] || kVertex == v[1]);
//                plint mVertex = (kVertex==v[0] ? v[1] : v[0]);
//                mPosition = dynMesh.getVertex(mVertex);
//                mPositionP = &mPosition;
//                mEqPosition = eqMesh.getVertex(mVertex);
//                mEqPositionP = &mEqPosition;
//            }
//            else {
//                mPositionP = mEqPositionP = NULL;
//            }
//
//            if (dynMesh.isInteriorEdge(lVertex, kVertex)) {
//                std::vector<plint> v = dynMesh.getNeighborVertexIds(lVertex, kVertex);
//                PLB_ASSERT(v.size() == 2);
//                PLB_ASSERT(jVertex == v[0] || jVertex == v[1]);
//                plint nVertex = (jVertex==v[0] ? v[1] : v[0]);
//                nPosition = dynMesh.getVertex(nVertex);
//                nPositionP = &nPosition;
//                nEqPosition = eqMesh.getVertex(nVertex);
//                nEqPositionP = &nEqPosition;
//            }
//            else {
//                nPositionP = nEqPositionP = NULL;
//            }
//
//            computeBendingStrain(kPositionP, jPositionP, lPositionP,
//                                 mPositionP, nPositionP, iPositionP,
//                                 kEqPositionP, jEqPositionP, lEqPositionP,
//                                 mEqPositionP, nEqPositionP, iEqPositionP,
//                                 E);
//
//            ud += computeBendingPotentialDensity(E, Y, nu, h);
//        }
//    }
//
//    return (u+ud*area);
//}
//
//template<typename T>
//void computeMembraneStrain(Array<T,3> const& iPosition,
//                           Array<T,3> const& jPosition,
//                           Array<T,3> const& kPosition,
//                           Array<T,3> const& iEqPosition,
//                           Array<T,3> const& jEqPosition,
//                           Array<T,3> const& kEqPosition,
//                           Array<Array<T,3>,3>& E)
//{
//    Array<Array<T,3>,3> v;
//    v[0] = jPosition - iPosition;
//    v[1] = kPosition - jPosition;
//    v[2] = iPosition - kPosition;
//
//    Array<Array<T,3>,3> vEq;
//    vEq[0] = jEqPosition - iEqPosition;
//    vEq[1] = kEqPosition - jEqPosition;
//    vEq[2] = iEqPosition - kEqPosition;
//
//    Array<T,3> nEq;
//    crossProduct(vEq[0], vEq[1], nEq);
//    T norm_nEq = norm(nEq);
//    nEq /= norm_nEq;
//    T aEq = 0.5 * norm_nEq;
//
//    Array<Array<T,3>,3> t;
//    crossProduct(vEq[0], nEq, t[0]);
//    crossProduct(vEq[1], nEq, t[1]);
//    crossProduct(vEq[2], nEq, t[2]);
//
//    E[0].resetToZero();
//    E[1].resetToZero();
//    E[2].resetToZero();
//
//    for (int i = 0; i < 3; i++) { // Loop over local edges
//        int j = (i == 0 ? 1 : (i == 1 ? 2 : 0));
//        int k = (j == 0 ? 1 : (j == 1 ? 2 : 0));
//
//        T l   = norm(v[i]);
//        T lEq = norm(vEq[i]);
//
//        T coef = (l*l - lEq*lEq) / (8.0*aEq*aEq);
//
//        for (int p = 0; p < 3; p++)
//            for (int q = 0; q < 3; q++)
//                E[p][q] += coef*(t[j][p]*t[k][q] + t[k][p]*t[j][q]);
//    }
//}
//
//template<typename T>
//void computeBendingStrain(Array<T,3> const* iPositionP, Array<T,3> const* jPositionP,
//                          Array<T,3> const* kPositionP,
//                          Array<T,3> const* lPositionP, Array<T,3> const* mPositionP,
//                          Array<T,3> const* nPositionP,
//                          Array<T,3> const* iEqPositionP, Array<T,3> const* jEqPositionP,
//                          Array<T,3> const* kEqPositionP,
//                          Array<T,3> const* lEqPositionP, Array<T,3> const* mEqPositionP,
//                          Array<T,3> const* nEqPositionP,
//                          Array<Array<T,3>,3>& E)
//{
//    // The position arrays specify at most four oriented triangles.
//    // (i, j, k) is the basic triangle under consideration.
//    // The rest of the triangles are the neighboring triangles of the basic one
//    // so some of them might not exist. If one of them does not exist,
//    // then the corresponding pointer argument of its vertex (which is
//    // not a vertex of the basic triangle) is NULL.
//    // The oriented neighboring triangles are:
//    // (l, k, j) which shares the edge k-j with the basic triangle,
//    // (m, i, k) which shares the edge i-k with the basic triangle,
//    // (n, j, i) which shares the edge j-i with the basic triangle.
//
//    E[0].resetToZero();
//    E[1].resetToZero();
//    E[2].resetToZero();
//
//    Array<T,3> tmp0 = *jEqPositionP - *iEqPositionP;
//    Array<T,3> tmp1 = *kEqPositionP - *iEqPositionP;
//    Array<T,3> nEq;
//    crossProduct(tmp0, tmp1, nEq);
//    T norm_nEq = norm(nEq);
//    nEq /= norm_nEq;
//    T aEq = 0.5 * norm_nEq;
//
//    for (int i = 0; i < 3; i++) { // Loop over local edges
//        Array<T,3> const *pos3 = (i == 0 ? nPositionP : (i == 1 ? lPositionP : mPositionP));
//        if (pos3 != NULL) {
//            Array<T,3> const *pos0 = (i == 0 ? kPositionP : (i == 1 ? iPositionP : jPositionP));
//            Array<T,3> const *pos1 = (i == 0 ? iPositionP : (i == 1 ? jPositionP : kPositionP));
//            Array<T,3> const *pos2 = (i == 0 ? jPositionP : (i == 1 ? kPositionP : iPositionP));
//
//            Array<T,3> const *eqPos0 = (i == 0 ? kEqPositionP : (i == 1 ? iEqPositionP : jEqPositionP));
//            Array<T,3> const *eqPos1 = (i == 0 ? iEqPositionP : (i == 1 ? jEqPositionP : kEqPositionP));
//            Array<T,3> const *eqPos2 = (i == 0 ? jEqPositionP : (i == 1 ? kEqPositionP : iEqPositionP));
//            Array<T,3> const *eqPos3 = (i == 0 ? nEqPositionP : (i == 1 ? lEqPositionP : mEqPositionP));
//
//            Array<T,3> edge01, edge12, edge21, edge32, n012, n321;
//
//            // Compute quantities related to the equilibrium configuration
//
//            edge12 = *eqPos2 - *eqPos1;
//            edge21 = -edge12;
//            edge32 = *eqPos2 - *eqPos3;
//
//            crossProduct(edge32, edge21, n321);
//
//            T eqAngle = angleBetweenVectors(nEq, n321);
//
//            T lEq = norm(edge12);
//            Array<T,3> t;
//            crossProduct(edge12, nEq, t);
//
//            // Compute quantities related to the current configuration
//
//            edge01 = *pos1 - *pos0;
//            edge12 = *pos2 - *pos1;
//            edge21 = -edge12;
//            edge32 = *pos2 - *pos3;
//
//            crossProduct(edge01, edge12, n012);
//            crossProduct(edge32, edge21, n321);
//
//            T angle = angleBetweenVectors(n012, n321);
//
//            T coef = 2.0 * (tan(0.5*angle) - tan(0.5*eqAngle)) / (2.0*aEq*lEq);
//
//            // Compute the strain tensor
//
//            for (int p = 0; p < 3; p++)
//                for (int q = 0; q < 3; q++)
//                    E[p][q] += coef*t[p]*t[q];
//        }
//    }
//}
//
//template<typename T>
//T computeMembranePotentialDensity(Array<Array<T,3>,3> const& E, T Y, T nu, T h)
//{
//    T trE = E[0][0] + E[1][1] + E[2][2];
//
//    T trE2 = E[0][0]*E[0][0] + E[0][1]*E[1][0] + E[0][2]*E[2][0] +
//             E[1][0]*E[0][1] + E[1][1]*E[1][1] + E[1][2]*E[2][1] +
//             E[2][0]*E[0][2] + E[2][1]*E[1][2] + E[2][2]*E[2][2];
//
//    return (Y*h* ((1.0-nu)*trE2 + nu*trE*trE) / (2.0*(1.0 - nu*nu)));
//}
//
//template<typename T>
//T computeBendingPotentialDensity(Array<Array<T,3>,3> const& E, T Y, T nu, T h)
//{
//    T trE = E[0][0] + E[1][1] + E[2][2];
//
//    T trE2 = E[0][0]*E[0][0] + E[0][1]*E[1][0] + E[0][2]*E[2][0] +
//             E[1][0]*E[0][1] + E[1][1]*E[1][1] + E[1][2]*E[2][1] +
//             E[2][0]*E[0][2] + E[2][1]*E[1][2] + E[2][2]*E[2][2];
//
//    return (Y*h*h*h* ((1.0-nu)*trE2 + nu*trE*trE) / (24.0*(1.0 - nu*nu)));
//}
//
//}  // namespace membraneBendingModelHelper3D
//
//template<typename T>
//T computeRestPotential(Array<T,3> const& iPosition, Array<T,3> const& iEqPosition, T k)
//{
//    //T length = norm(iPosition - iEqPosition);
//    // ATTENTION: This code is in debugging state. Must be undone.
//    T length = fabs(iPosition[0] - iEqPosition[0]);
//    return 0.5*k*length*length;
//}

// template<typename T>
// T computeTriangleArea(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
//                       Array<T,3> const& kPosition)
// {
//     Array<T,3> ijEdge = jPosition - iPosition;
//     Array<T,3> ikEdge = kPosition - iPosition;
//     Array<T,3> n;
//     crossProduct(ijEdge, ikEdge, n);
//     return 0.5*norm(n);
// }

}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_HH
