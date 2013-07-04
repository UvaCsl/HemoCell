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

#ifndef COMPUTE_CELL_FORCES3D_HH
#define COMPUTE_CELL_FORCES3D_HH

#include <cmath>
#include <map>
#include "cellModel3D.h"


namespace plb {


template<typename T>
Array<T,3> computeInPlaneForce(Array<T,3> const& x1, Array<T,3> const& x2, T maxLength, T k_WLC, T k_rep) {
/*
 *  Computes In-Plane forces based on Worm-like chain forces and a repulsive potential.
 *
 *    k_WLC = kBT * maxLength/(4.0*persistenceLengthCoarse);
 *    A repulsive potential is added to equilibrate the forces at point x0.
 *    // Solving f_WLC + f_rep =0 for x=eqLength, f_rep = k_rep/L^m, m=2. //
 *    k_rep = k_rep_*(k_WLC*maxLength*pow(x0,3)*(6 - 9*x0 + 4*pow(x0,2)))/pow(-1 + x0,2);
 *
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b, Pivkin2008]

*/
    Array<T,3> dL = (x1 - x2)*1.0;
    T L = norm(dL);
    Array<T,3> eij = dL/L;
    T r = L*1.0/maxLength;
    /* In Plane Force (WLC) */
    Array<T,3> tmpForce = eij * (k_WLC*r*(-6 + (9 - 4*r)*r))/(maxLength*pow(-1 + r,2));
    /* Repulsive Force */
    tmpForce += eij * k_rep/(L*L);
    return tmpForce;
}


template<typename T>
Array<T,3> computeInPlaneExplicitForce(Array<T,3> const& x1, Array<T,3> const& x2, T eqLengthRatio, T eqLength, T k_inPlane) {
/*
 *  Computes In-Plane forces based on Worm-like chain forces and a repulsive potential.
 *  Has the same function as computeInPlaneForce but with different arguments.
 *
 *    eqLengthRatio = eqLength/maxLength ; // Used to be x0
 *    k_inPlane = kBT /(4.0*persistenceLengthCoarse);
 *
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b, Pivkin2008]

*/
    Array<T,3> dL = (x1 - x2)*1.0;
    T L = norm(dL);
    T r = L/(eqLength*eqLengthRatio), r0=1.0/eqLengthRatio; // T Lmax = eqLength*eqLengthRatio;
    Array<T,3> eij = dL/L;
    /* In Plane Force (WLC) and Repulsive Force */
    Array<T,3> tmpForce =  eij * k_inPlane * r * (-6 + (9 - 4*r)*r)/( (r-1)*(r-1) );
    T k_rep = -(eqLength*eqLength)* k_inPlane * r0 * (-6 + (9 - 4*r0)*r0)/( (r0-1)*(r0-1) );
    tmpForce +=  eij * k_rep / (L*L) ;
//    Array<T,3> tmpForce =  eij * k_inPlane *
//                    (1 - (4*L*eqLengthRatio)/eqLength +
//                            pow(eqLength,2)*((-1 + pow(-1 + eqLengthRatio,-2) + 4*eqLengthRatio)/pow(L,2) - pow(eqLength - L*eqLengthRatio,-2)));
    return tmpForce;
}


template<typename T>
Array<T,3> computeDissipativeForce(Array<T,3> const& x1, Array<T,3> const& x2,
                                   Array<T,3> const& v1, Array<T,3> const& v2,
                                   T gamma_T, T gamma_C) {
/*
 *
//    Dissipative term coefficients from FedosovCaswellKarniadakis2010 //
        gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
        gamma_C = (gamma_T/3.0);
 *  Where eta_m is the 2D membrane viscosity.
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b]
*/
    Array<T,3> dL = (x1 - x2)*1.0;
    Array<T,3> eij = dL/norm(dL);
    Array<T,3> vij = v1 - v2;
    return -gamma_T*vij -gamma_C*dot(vij,eij)*eij;
}


template<typename T>
Array<T,3> computeVolumeConservationForce(
                       Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3, T cVolume) {
/*
 * Global volume conservation force, acting on a vertex.
 *
        cVolume = k_volume * kBT/pow(eqLength,3) * (cellVolume - eqVolume)*1.0/eqVolume;

 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b]
*/
    Array<T,3> tmp;
    crossProduct(x3, x2, tmp);
    return -cVolume  * (- 1.0/6.0 * tmp);
}

template<typename T>
Array<T,3> computeSurfaceConservationForce(
                Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
                Array<T,3> triangleNormal, T cSurface) {
/*
 * Global surface conservation force, acting on a vertex.
        cSurface = k_surface * kBT/pow(eqLength,2) * (cellSurface - eqSurface)*1.0/eqSurface;
 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b]
*/
    Array<T,3> tmp;
    return computeSurfaceConservationForce(x1, x2, x3,
                triangleNormal, cSurface, tmp);
}

template<typename T>
Array<T,3> computeSurfaceConservationForce(
                Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
                Array<T,3> const& triangleNormal, T cSurface, Array<T,3> & dAdx) {
/*
 * Global surface conservation force, acting on a vertex.
        cSurface = k_surface * kBT/pow(eqLength,2) * (cellSurface - eqSurface)*1.0/eqSurface;
 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b]
*/
    Array<T,3> tmp;
    crossProduct(triangleNormal,
                x3 - x2,
                tmp);
    dAdx = 0.5 *tmp; tmp.resetToZero();
    return -cSurface * dAdx;
}

template<typename T>
Array<T,3> computeLocalAreaConservationForce(
                Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
                Array<T,3> const& triangleNormal,  T triangleArea, T eqArea,
                T cArea) {
/*
 * Local surface conservation force.
        cArea = k_shear*kBT/pow(eqLength,2)/eqArea;
 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b]
*/
    Array<T,3> tmp;
    crossProduct(triangleNormal,
                x3 - x2,
                tmp);
    Array<T,3> dAdx = 0.5 *tmp;
    return computeLocalAreaConservationForce(dAdx, triangleArea, eqArea, cArea);
}

template<typename T>
Array<T,3> computeLocalAreaConservationForce(
                Array<T,3> const& dAdx, T triangleArea, T eqArea,
                T cArea) {
/*
 * Local surface conservation force.
        cArea = k_shear*kBT/pow(eqLength,2)/eqArea;
 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b]
*/
    return -cArea * (triangleArea - eqArea)*1.0 * dAdx;
}


/* Elastic force to couple with the WLC force.
    C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
             * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
             / (64.0*persistenceLengthCoarse)
             * (4*x0*x0 - 9*x0 + 6)
             / (1-x0)*(1-x0);
 *  Related publications: [Pivkin2008, FedosovCaswellKarniadakis2010] */
template<typename T>
Array<T,3> computeElasticRepulsiveForce(
        Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
        Array<T,3> const& triangleNormal,  T triangleArea, T cElastic) {
    Array<T,3> tmp;
    crossProduct(triangleNormal,
                x3 - x2,
                tmp);
    Array<T,3> dAdx = 0.5 *tmp;
    return computeElasticRepulsiveForce(dAdx, triangleArea, cElastic);
}


/* Elastic force to couple with the WLC force.
    C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
             * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
             / (64.0*persistenceLengthCoarse)
             * (4*x0*x0 - 9*x0 + 6)
             / (1-x0)*(1-x0);
 *  Related publications: [Pivkin2008, FedosovCaswellKarniadakis2010] */
template<typename T>
Array<T,3> computeElasticRepulsiveForce(Array<T,3> const& dAdx, T triangleArea, T cElastic) {
    return (cElastic*1.0)/(triangleArea*triangleArea) * dAdx;
}


template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k)
{
/*
 * The most messy force!
 *
 *
 * Triangles are:
 *      (i, j, k) and (l, k, j). These triangles share
 *      (x2, x1, x3) and (x4, x3, x1). These triangles share
 *      the common edge j-k.
 *
 *      crossProduct(jPosition - iPosition, kPosition - jPosition, nijk);
 *      crossProduct(kPosition - lPosition, jPosition - kPosition, nlkj);
 *      crossProduct(x1 - x2, x3 - x1, nijk);
 *      crossProduct(x3 - x4, x1 - x3, nlkj);
*/

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


template<typename T>
Array<T,3> computeBendingForce_Krueger (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k,
                                Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4)
{
    // (i, j, k) and (l, k, j). These triangles share
    // (x2, x1, x3) and (x4, x3, x1). These triangles share
    // the common edge j-k.

    fx2.resetToZero(); fx3.resetToZero(); fx4.resetToZero();

    T ninj = dot(ni,nj);
    T edgeAngle = angleBetweenVectors(ni, nj);
    T factor = k * (edgeAngle - eqAngle) * (-1.0/(1.0 - sqrt(ninj*ninj)) );
    Array<T,3> fx1(0,0,0), dninjdx(0,0,0);
    Array<T,3> dx23(0,0,0), dx34(0,0,0);
    Array<T,3> dx31(0,0,0), dx13(0,0,0);

    // Calculation of force for x1
    crossProduct(x2-x3, nj - ninj*ni, dx23);
    crossProduct(x3-x4, ni - ninj*nj, dx34);
    dninjdx = 0.5/Ai*dx23 + 0.5/Aj*dx34;
    fx1 = -factor*dninjdx;

    // Calculation of force for x2
    crossProduct(x3-x1, nj - ninj*ni, dx31);
    dninjdx = 0.5/Ai*dx31;
    fx2 = -factor*dninjdx;

    // Calculation of force for x4
    crossProduct(x1-x3, ni - ninj*nj, dx13);
    dninjdx = 0.5/Aj*dx13;
    fx4 = -factor*dninjdx;

    // Calculation of force for x3
    fx3 = - fx1 - fx2 - fx4;
    return fx1;
}


template<typename T>
T computeBendingPotential (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                T eqTileSpan, T eqLength, T eqAngle, T k)
{
    Array<T,3> ni(0.,0.,0.),  nj(0.,0.,0.);
    crossProduct(x3-x1,  x4-x1, nj);
    crossProduct(x2-x1,  x3-x1, ni);
    T edgeAngle = angleBetweenVectors(ni, nj);
    return k * (0.5*(edgeAngle - eqAngle)*(edgeAngle - eqAngle));
//    return k * (1-cos(edgeAngle - eqAngle));
}


template<typename T>
Array<T,3> computeBendingForceFromPotential (
        Array<T,3> const& x1, Array<T,3> const& x2,
        Array<T,3> const& x3, Array<T,3> const& x4,
        T eqTileSpan, T eqLength, T eqAngle, T k,
        Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4)
{
    Array<T,3> fx1, vertex;
    static T eps = (sizeof(T) == sizeof(float) ?
            100.0 * std::numeric_limits<T>::epsilon() :
            std::numeric_limits<float>::epsilon());

    vertex = x1;
    fx1.resetToZero();
    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = computeBendingPotential (
                iPosition, x2, x3, x4,
                eqTileSpan, eqLength, eqAngle, k
                );
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = computeBendingPotential (
                iPosition, x2, x3, x4,
                eqTileSpan, eqLength, eqAngle, k);
        fx1[i] = -(up-um) / (2.0*eps);
    }

    vertex = x2;
    fx2.resetToZero();
    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = computeBendingPotential (
                x1, iPosition, x3, x4,
                eqTileSpan, eqLength, eqAngle, k
                );
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = computeBendingPotential (
                x1, iPosition, x3, x4,
                eqTileSpan, eqLength, eqAngle, k);
        fx2[i] = -(up-um) / (2.0*eps);
    }

    vertex = x3;
    fx3.resetToZero();
    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = computeBendingPotential (
                x1, x2, iPosition, x4,
                eqTileSpan, eqLength, eqAngle, k
                );
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = computeBendingPotential (
                x1, x2, iPosition, x4,
                eqTileSpan, eqLength, eqAngle, k);
        fx3[i] = -(up-um) / (2.0*eps);
    }

    vertex = x4;
    fx4.resetToZero();
    for (int i = 0; i < 3; i++) {
        Array<T,3> iPosition = vertex;
        iPosition[i] += eps;
        T up = computeBendingPotential (
                x1, x2, x3, iPosition,
                eqTileSpan, eqLength, eqAngle, k
                );
        iPosition = vertex;
        iPosition[i] -= eps;
        T um = computeBendingPotential (
                x1, x2, x3, iPosition,
                eqTileSpan, eqLength, eqAngle, k);
        fx4[i] = -(up-um) / (2.0*eps);
    }
    Array<T,3> df = fx1 + fx2 + fx3 + fx4;
    if (norm(df) > 1e-10) {
        pcout << "!!! Something's wrong with the angles !!!" << norm(df) << std::endl;
        fx1 = fx1 - df/4.;
        fx2 = fx2 - df/4.;
        fx3 = fx3 - df/4.;
        fx4 = fx4 - df/4.;
    }
    return fx1;
}



}  // namespace plb

#endif  // COMPUTE_CELL_FORCES3D_HH
