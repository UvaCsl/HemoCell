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



namespace plb {

template<typename T>
CellModel3D<T>::CellModel3D (
        T density_, T k_rest_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_rep_, T k_elastic_,
        T k_volume_, T k_surface_, T eta_m_,
        T eqArea_, T eqLength_, T eqAngle_,
        T eqVolume_, T eqSurface_, T eqTileSpan_,
        T maxLength_, T persistenceLengthFine, pluint Nv)
    : ShellModel3D<T>(density_),
      k_rest(k_rest_),
      k_shear(k_shear_),
      k_bend(k_bend_),
      k_stretch(k_stretch_),
      k_elastic(k_elastic_),
      k_surface(k_surface_),
      k_volume(k_volume_),
      eta_m(eta_m_),
      eqLength(eqLength_),
      eqArea(eqArea_),
      eqAngle(eqAngle_),
      eqVolume(eqVolume_),
      eqSurface(eqSurface_),
      eqTileSpan(eqTileSpan_),
      maxLength(maxLength_)
{
    T x0 = eqLength*1.0/maxLength;
    persistenceLengthCoarse = persistenceLengthFine ;//* sqrt( (Nv-2.0) / (23867-2.0)) ;
    /* Use dimensionless coefficients */
    k_volume *= kBT/pow(eqLength,3);
    k_surface *= kBT/pow(eqLength,2);
    k_shear *= kBT/pow(eqLength,2);

    /* Dissipative term coefficients from FedosovCaswellKarniadakis2010 */
    gamma_T = (eta_m * 4 * sqrt(3.0) / 13.0);
    gamma_C = (gamma_T/3.0);
    /* Multiplying with eqLength because the units on the paper are wrong */
    gamma_T *= eqLength;
    gamma_C *= eqLength;

    k_WLC = k_WLC_ * kBT * maxLength/(4.0*persistenceLengthCoarse);
    /* Solving f_WLC + f_rep =0 for x=eqLength, f_rep = k_rep/L^m, m=2. */
    k_rep = k_rep_*(k_WLC*maxLength*pow(x0,3)*(6 - 9*x0 + 4*pow(x0,2)))/pow(-1 + x0,2);


    T forceSum =
            (k_WLC*x0*(-6 + (9 - 4*x0)*x0))/(maxLength*pow(-1 + x0,2)) // WLC part
           + k_rep/(eqLength*eqLength); // repulsive part
    C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
             * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
             / (64.0*persistenceLengthCoarse)
             * (4*x0*x0 - 9*x0 + 6)
             / (1-x0)*(1-x0);
    pcout << std::endl;
    pcout << " ============================================= " << std::endl;
    pcout << "k_WLC: " << k_WLC << ",\t eqLength: " << eqLength << std::endl;
    pcout << "k_rep: " << k_rep << ",\t forceSum: " << forceSum << std::endl;
    pcout << "k_bend: " << k_bend << ",\t eqAngle (degrees): " << eqAngle*180.0/pi << std::endl;
    pcout << "k_volume: " << k_volume << ",\t eqVolume: " << eqVolume << std::endl;
    pcout << "k_surface: " << k_surface << ",\t eqSurface: " << eqSurface << std::endl;
    pcout << "k_shear: " << k_shear << ",\t eqArea: " << eqArea << std::endl;
    pcout << "eta_m: " << eta_m << ",\t x0: " << x0 << std::endl;
    pcout << "gamma_T: " << gamma_T << ",\t persistenceLengthCoarse: " << persistenceLengthCoarse << std::endl;
    pcout << "gamma_C: " << gamma_C << ",\t maxLength: " << maxLength << std::endl;
    pcout << "k_rest: " << k_rest << ",\t 0 : " << 0 << std::endl;
    pcout << "k_stretch: " << k_stretch << ",\t eqTileSpan: " << eqTileSpan << std::endl;
    pcout << "k_elastic: " << k_elastic << ",\t eqLength: " << eqLength << std::endl;
    pcout << "* k_bend: " << k_bend/kBT << std::endl;
    pcout << "* k_volume: " <<  k_volume/(kBT/pow(eqLength,3)) <<  std::endl;
    pcout << "* k_surface: " << k_surface/(kBT/pow(eqLength,2)) <<  std::endl;
    pcout << "* k_shear: " << k_shear/(kBT/pow(eqLength,2)) <<  std::endl;
    pcout << "* eqLength from eqArea: " << sqrt(4*eqArea/sqrt(3.0)) << ",\t eqLength: " << eqLength << std::endl;
    pcout << " ============================================= " << std::endl;
}


template<typename T>
Array<T,3> CellModel3D<T>::computeCellForce (
        TriangleBoundary3D<T> const& boundary,
        T cellVolume, T cellSurface, T & iSurface,
        std::map< plint, Array<T,3> > & particleVelocity,
        std::map< plint, Array<T,3> > & particleForces,
        plint iVertex,
        Array<T,3> & f_wlc, Array<T,3> & f_bending, Array<T,3> & f_volume,
        Array<T,3> & f_surface, Array<T,3> & f_shear, Array<T,3> & f_viscosity)
 /* Some force calculations are according to KrugerThesis, Appendix C */
{
    /* Initializations from input */
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    Array<T,3> x1 = dynMesh.getVertex(iVertex), x2, x3, x4;
    Array<T,3> iVelocity = particleVelocity[iVertex];
    /* ===================== In case of quasi-rigid object =====================
     *
     * If this is a boundary element (k_rest != 0), get the reference locations
     * of iVertex and calculate and return the force for quasi-rigid objects.
     *          (FengMichaelides2004, J.Comp.Phys. 195(2))
     *
     * */
    if (k_rest != 0.0) {
        Array<T,3> x1ref;
        boundary.pushSelect(0,1);
        x1ref = boundary.getMesh().getVertex(iVertex);
        boundary.popSelect();
        Array<T,3> dx = x1-x1ref;
        return (-k_rest*dx) + (-gamma_C*dot(iVelocity,dx)/eqLength * dx); // Dissipative term from Dupin2007
    }
    /* Force initializations */
    Array<T,3> inPlaneForce; inPlaneForce.resetToZero();
    Array<T,3> elasticForce; elasticForce.resetToZero();
    Array<T,3> repulsiveForce; repulsiveForce.resetToZero();
    Array<T,3> bendingForce; bendingForce.resetToZero();
    Array<T,3> surfaceForce; surfaceForce.resetToZero();
    Array<T,3> shearForce; shearForce.resetToZero();
    Array<T,3> volumeForce; volumeForce.resetToZero();
    Array<T,3> stretchForce; stretchForce.resetToZero();
    Array<T,3> dissipativeForce; dissipativeForce.resetToZero();
    /* Calculate cell coefficients */
    T volumeCoefficient = k_volume * (cellVolume - eqVolume)*1.0/eqVolume;
    T surfaceCoefficient = k_surface * (cellSurface - eqSurface)*1.0/eqSurface;
    T areaCoefficient = k_shear/eqArea ;
    iSurface = 0.0;

    /* Run through all the neighbouring faces of iVertex and calculate:
         x Volume conservation force
         x Surface conservation force
         x Shear force
     */
    Array<T,3> dAdx, dVdx;
    std::map<plint, T> trianglesArea;
    T triangleArea;
    std::map<plint, Array<T,3> > trianglesNormal;
    Array<T,3> triangleNormal;

    plint iTriangle, jTriangle;
    plint iX1, iX2, iX3;
    Array<T,3> xi[3], tmp;

    std::vector<plint> neighborTriangles = dynMesh.getNeighborTriangleIds(iVertex);
    for (pluint iB = 0; iB < neighborTriangles.size(); ++iB) {
        iTriangle = neighborTriangles[iB];
        for (pluint id = 0; id < 3; ++id) { xi[id] = dynMesh.getVertex(iTriangle,id); }
        if (xi[0] ==x1) { iX1 = 0;}
        else if (xi[1] ==x1) { iX1 = 1;}
        else { iX1 = 2;}
        iX2 = (iX1 + 1)%3;
        iX3 = (iX1 + 2)%3;
        triangleArea = trianglesArea[iTriangle] = dynMesh.computeTriangleArea(iTriangle);
        iSurface += triangleArea/3.0;
        trianglesNormal[iTriangle] = triangleNormal = dynMesh.computeTriangleNormal(iTriangle);
        /* Surface conservation force */
        crossProduct(triangleNormal,
                    xi[iX3] - xi[iX2],
                    tmp);
        dAdx = 0.5 *tmp;
        surfaceForce += -surfaceCoefficient * dAdx;
        /* Shear force */
        shearForce += - areaCoefficient * (triangleArea - eqArea)*1.0 * dAdx;
        /* Volume conservation force */
        crossProduct(xi[iX3], xi[iX2], tmp);
        dVdx = - 1.0/6.0 * tmp;
        volumeForce  += -volumeCoefficient  * dVdx;
        /* Elastice Force */
        elasticForce += (C_elastic*1.0)/(triangleArea*triangleArea) * dAdx;
    }

    /* Run through all the neighbouring vertices of iVertex and calculate:
         x In plane (WLC) force
         x Repulsive force
         o Stretch force
         x Dissipative force
         x Bending force
     */
    T L, r;
    Array<T,3> dL, eij; dL.resetToZero(); eij.resetToZero();

    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
    for (pluint jV = 0; jV < neighborVertexIds.size(); jV++) {
        plint jVertex = neighborVertexIds[jV];
        x3 = dynMesh.getVertex(jVertex);
        dL = (x3 - x1)*1.0;
        L = norm(dL);
        r = L/maxLength;
        eij = dL/L;
        /* In Plane Force (WLC) */
        inPlaneForce += eij * (k_WLC*r*(-6 + (9 - 4*r)*r))/(maxLength*pow(-1 + r,2)); // inPlaneForce += k_WLC / maxLength * eij * (-1.0 +  4.0*r + 1.0/((1.0-r)*(1.0-r)) );
        /* Repulsive Force */
        repulsiveForce +=  eij * k_rep/(L*L);
        /* Stretch force */
//        stretchForce += - eij * k_stretch * (L/eqLength - 1.0);
        /*  Dissipative Forces Calculations */
        Array<T,3> vij = iVelocity - particleVelocity[jVertex];
        dissipativeForce += -gamma_T*vij -gamma_C*dot(vij,eij)*eij;
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
        Array<T,3> iVertexBendingForce = computeBendingForce (x1, x2, x3, x4,
                            trianglesNormal[iTriangle], trianglesNormal[jTriangle],
                            trianglesArea[iTriangle], trianglesArea[jTriangle],
                            eqTileSpan, eqLength, eqAngle, k_bend);
        bendingForce += iVertexBendingForce;
        T iTriangleBendingCoefficient = trianglesArea[iTriangle] / (trianglesArea[iTriangle] + trianglesArea[jTriangle]);
        T jTriangleBendingCoefficient = trianglesArea[jTriangle] / (trianglesArea[iTriangle] + trianglesArea[jTriangle]);
        particleForces[kVertex] += -iTriangleBendingCoefficient * iVertexBendingForce/2.0;
        particleForces[lVertex] += -jTriangleBendingCoefficient * iVertexBendingForce/2.0;
    }
    f_wlc = inPlaneForce + repulsiveForce;
    f_bending = bendingForce;
    f_volume = volumeForce;
    f_surface = surfaceForce;
    f_shear = shearForce;
    f_viscosity = dissipativeForce;
    return (inPlaneForce + elasticForce + repulsiveForce) + bendingForce + (volumeForce + surfaceForce + shearForce) + dissipativeForce;// + stretchForce;
}


template<typename T>
Array<T,3> CellModel3D<T>::computeElasticForce (
        TriangleBoundary3D<T> const& boundary,
        plint iVertex )
{
    TriangularSurfaceMesh<T> const& mesh = boundary.getMesh();
    return Array<T,3>(0,0,0);

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
                eqLength, maxLength, eqArea, eqAngle, eqTileSpan,
                k_WLC, k_elastic, k_shear, k_bend);
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = shellModelHelper3D::cellModelHelper3D::computePotential (
                iVertex, iPosition, mesh,
                eqLength, maxLength, eqArea, eqAngle, eqTileSpan,
                k_WLC, k_elastic, k_shear, k_bend);
        force[i] = -(up-um) / (2.0*eps);
    }
    return force;
}

template<typename T>
CellModel3D<T>* CellModel3D<T>::clone() const {
    return new CellModel3D<T>(*this);
}

//template<typename T>
//Array<T,3> computeStrainForce (Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
//                               T Ai, T eqArea, T eqLength, T k_s, T k_a)
//{
//    T I1, dI1dx;
//    T I2, dI2dx;
//
//    x21 = x1 - x2;
//    x31 = x1 - x3;
//
//    T sin_phi = 2*Ai/(norm(x21)*norm(x31));
//    T cos_phi = sqrt(1 - sin_phi*sin_phi);
//
//    T sin_phi0 = 2*eqArea/(eqLength*eqLength);
//    T cos_phi0 = sqrt(1 - sin_phi0*sin_phi0);
//
//    Dxx = x21/eqLength;
//    dDxxdx = 1.0/eqLength;
//
//
//
//    I1 = Dxx*Dxx + Dyy*Dyy + Dxy*Dxy - 2.0;
//    I2 = Dxx*Dxx*Dyy*Dyy - 1.0;
//
//
//    T desdx = k_s/12.0*(2*I1*dI1dx + 2*dI1dx - 2*dI2dx);
//    desdx += k_a/12.0*(2*I2*dI2dx);
//}


template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k)
{
    // (i, j, k) and (l, k, j). These triangles share
    // (x2, x1, x3) and (x4, x3, x1). These triangles share
    // the common edge j-k.
    // crossProduct(jPosition - iPosition, kPosition - jPosition, nijk);
    // crossProduct(kPosition - lPosition, jPosition - kPosition, nlkj);

    // crossProduct(x1 - x2, x3 - x1, nijk);
    // crossProduct(x3 - x4, x1 - x3, nlkj);

    Array<T,3> v1, v2, D1, D2, dAngledx,tmp;
    Array<T,3> x32=x2-x3, x43=x3-x4;
//    =============================================================
//     Probably the following computations are the correct ones,
//     but they do not agree with the potential calculation.
//    =============================================================
//    v1 = 2 * Ai * ni;
//    v2 = 2 * Aj * nj;
//    pcout << "aN1 = (" << v1[0] << ", " << v1[1] << ", " << v1[2] << ") " << std::endl;
//    pcout << "aN2 = (" << v2[0] << ", " << v2[1] << ", " << v2[2] << ") " << std::endl;
//    =============================================================
//     Probably the following computations are incorrect,
//     but they do agree with the potential calculation.
//    =============================================================
    crossProduct(x1 - x2, x32, v1);
    crossProduct(x1 - x3, x43, v2);
//    pcout << "cP1 = (" << v1[0] << ", " << v1[1] << ", " << v1[2] << ") " << std::endl;
//    pcout << "cP2 = (" << v2[0] << ", " << v2[1] << ", " << v2[2] << ") " << std::endl;
    T angle = angleBetweenVectors(v1, v2);
    T cosAngle=cos(angle),overSinAngle=1.0/sin(angle);
    T nv1 = 2 * Ai, nv2 = 2 * Aj;
    dAngledx.resetToZero();
    for (int var = 0; var < 3; ++var) {
        D1.resetToZero(); D2.resetToZero(); tmp.resetToZero();
        tmp[var] = 1.0;
        crossProduct(tmp, x32, D1);
        crossProduct(tmp, x43, D2);
        dAngledx[var] = (dot(D1,v2) + dot(D2,v1))/nv1/nv2;
        dAngledx[var] -= cosAngle*(dot(D1,v1)/nv1/nv1 + dot(D2,v2)/nv2/nv2);
        dAngledx[var] *= -overSinAngle;
    }
    return -k * sin(angle - eqAngle) * dAngledx ; //* eqLength / eqTileSpan;
//    pcout << "fangle: " << angle << " " << eqAngle << std::endl;
//    T dAngle = (angle - eqAngle);
//    return -k* (dAngle*(1 - dAngle*dAngle/6.0)) * dAngledx ; //* eqLength / eqTileSpan;
}


namespace shellModelHelper3D {

namespace cellModelHelper3D {

template<typename T>
T computePotential(plint iVertex, Array<T,3> const& iPosition,
                   TriangularSurfaceMesh<T> const& mesh,
                   T eqLength, T maxLength, T eqArea, T eqAngle, T eqTileSpan,
                   T k_WLC, T k_elastic, T k_shear, T k_bend)
{
    T u = 0.0;

    // Membrane In Plane mode

    std::vector<plint> neighborVertexIds = mesh.getNeighborVertexIds(iVertex);

    pluint sz = neighborVertexIds.size();

    for (pluint i = 0; i < sz; i++) {
        plint jVertex = neighborVertexIds[i];
        u += computeInPlanePotential(iPosition, mesh.getVertex(jVertex),
                maxLength, k_WLC);
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
                                      eqTileSpan, eqLength, eqAngle, k_bend);
        }
    }

    return u;
}

template<typename T>
T computeInPlanePotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                             T maxLength, T k)
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
                       T eqTileSpan, T eqLength, T eqAngle, T k)
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

    return 0.5*k*(angle - eqAngle)*(angle - eqAngle);// * eqLength / eqTileSpan;
}

}  // namespace cellModelHelper3D


}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_HH
