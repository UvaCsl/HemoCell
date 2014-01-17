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

#ifndef COMPUTE_CELL_FORCES3D_H
#define COMPUTE_CELL_FORCES3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
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


/*  Computes In-Plane forces based on Worm-like chain forces and a repulsive potential.
 *
 *    k_WLC = kBT * maxLength/(4.0*persistenceLengthCoarse);
 *    A repulsive potential is added to equilibrate the forces at point x0.
 *    // Solving f_WLC + f_rep =0 for x=eqLength, f_rep = k_rep/L^m, m=2. //
 *    k_rep = k_rep_*(k_WLC*maxLength*pow(x0,3)*(6 - 9*x0 + 4*pow(x0,2)))/pow(-1 + x0,2);
 *
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b, Pivkin2008] */
template<typename T>
Array<T,3> computeInPlaneForce(Array<T,3> const& x1, Array<T,3> const& x2, T maxLength, T k_WLC, T k_rep);


/*  Computes In-Plane forces based on Worm-like chain forces and a repulsive potential.
 *  Has the same function as computeInPlaneForce but with different arguments.
 *
 *    eqLengthRatio = eqLength/maxLength ; // Used to be x0
 *    k_inPlane = kBT /(4.0*persistenceLengthCoarse);
 *
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b, Pivkin2008] */
template<typename T>
Array<T,3> computeInPlaneExplicitForce(Array<T,3> const& x1, Array<T,3> const& x2, T eqLengthRatio, T eqLength, T k_inPlane);


/* Dissipative term coefficients from FedosovCaswellKarniadakis2010
        gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
        gamma_C = (gamma_T/3.0);
 *  Where eta_m is the 2D membrane viscosity.
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b] */
template<typename T>
Array<T,3> computeDissipativeForce(Array<T,3> const& x1, Array<T,3> const& x2,
                                   Array<T,3> const& v1, Array<T,3> const& v2,
                                   T gamma_T, T gamma_C);


/* Global volume conservation force, acting on a vertex.
    cVolume = k_volume * kBT/pow(eqLength,3) * (cellVolume - eqVolume)*1.0/eqVolume;
*  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b] */
template<typename T>
Array<T,3> computeVolumeConservationForce(
                       Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3, T cVolume);


/* Global surface conservation force, acting on a vertex.
       cSurface = k_surface * kBT/pow(eqLength,2) * (cellSurface - eqSurface)*1.0/eqSurface;
*  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b] */
template<typename T>
Array<T,3> computeSurfaceConservationForce(
                Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
                Array<T,3> triangleNormal, T cSurface);

/* Global surface conservation force, acting on a vertex.
       cSurface = k_surface * kBT/pow(eqLength,2) * (cellSurface - eqSurface)*1.0/eqSurface;
*  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b] */
template<typename T>
Array<T,3> computeSurfaceConservationForce(
                Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
                Array<T,3> const& triangleNormal, T cSurface, Array<T,3> & dAdx);


/* Local surface conservation force.
        cArea = k_shear*kBT/pow(eqLength,2)/eqArea;
 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b] */
template<typename T>
Array<T,3> computeLocalAreaConservationForce(
                Array<T,3> const& x1, Array<T,3> const& x2, Array<T,3> const& x3,
                Array<T,3> const& triangleNormal,  T triangleArea, T eqArea,
                T cArea);


/* Local surface conservation force.
        cArea = k_shear*kBT/pow(eqLength,2)/eqArea;
 *  Related publications: [Pivkin2008] and secondary [FedosovCaswellKarniadakis2010, FedosovCaswell2010b] */
template<typename T>
Array<T,3> computeLocalAreaConservationForce(
                Array<T,3> const& dAdx, T triangleArea, T eqArea,
                T cArea);


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
        Array<T,3> const& triangleNormal,  T triangleArea, T cElastic);


/* Elastic force to couple with the WLC force.
    C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
             * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
             / (64.0*persistenceLengthCoarse)
             * (4*x0*x0 - 9*x0 + 6)
             / (1-x0)*(1-x0);
 *  Related publications: [Pivkin2008, FedosovCaswellKarniadakis2010] */
template<typename T>
Array<T,3> computeElasticRepulsiveForce(Array<T,3> const& dAdx, T triangleArea, T cElastic);

/* Calculates the bending force for the triangles formed by the vertices:
 * (x1, x2, x3) and (x1,x3,x4) with the common edge (x1,x3).
 * eqAngle is expected to be between [-pi,pi].  */
/* The most messy force! */
template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k,
								Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4);

/*
 * Calculates the bending forces.
 * Angles are expected to be between [-pi,pi].
 */
/* The most messy force! */
template<typename T>
Array<T,3> computeBendingForces (T edgeAngle, T eqAngle, T k,
								Array<T,3> const& ni, Array<T,3> const& nj,
								Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4) ;

template<typename T>
Array<T,3> computeBendingForceEdge (T edgeAngle, T eqAngle, T k,
								Array<T,3> const& ni, Array<T,3> const& nj) ;

template<typename T>
Array<T,3> computeBendingForceLateral (T edgeAngle, T eqAngle, T k,
								Array<T,3> const& ni) ;

/*
 * Calculates the bending potential.
 * Angles are expected to be between [-pi,pi].
 */
template<typename T>
T computeBendingPotential (T edgeAngle, T eqAngle, T k);

/*
 * Helper function, calculates the angle between -pi and pi
 */
template<typename T>
T calculateSignedAngle(TriangularSurfaceMesh<T> const& mesh, plint iVertex, plint jVertex, plint & kVertex, plint & lVertex) ;

/*
 * Helper function, calculates the angle between -pi and pi
 */
template<typename T>
T calculateSignedAngle(TriangularSurfaceMesh<T> const& mesh, plint iVertex, plint jVertex) ;

}  // namespace plb

#include "computeCellForces3D.hh"

#endif  // COMPUTE_CELL_FORCES3D_H







/*
 * Legacy Code, working only for convex surfaces.
 *
 *
 *


template<typename T>
T computeBendingPotential (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                T eqTileSpan, T eqLength, T eqAngle, T k)
{
    Array<T,3> ni(0.,0.,0.),  nj(0.,0.,0.);
    crossProduct(x2-x1,  x3-x2, nj);
    crossProduct(x3-x1,  x4-x3, ni);
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
        pcout << "!!! X-Axis (" << x1[0] << ",\t" << x2[0] << ",\t" << x3[0] << ",\t" << x4[0] << ")" << std::endl;
        pcout << "!!! Y-Axis (" << x1[1] << ",\t" << x2[1] << ",\t" << x3[1] << ",\t" << x4[1] << ")" << std::endl;
        pcout << "!!! Z-Axis (" << x1[2] << ",\t" << x2[2] << ",\t" << x3[2] << ",\t" << x4[2] << ")" << std::endl;
        pcout << "!!! Fx-Axis (" << fx1[0] << ",\t" << fx2[0] << ",\t" << fx3[0] << ",\t" << fx4[0] << ")" << std::endl;
        pcout << "!!! Fy-Axis (" << fx1[1] << ",\t" << fx2[1] << ",\t" << fx3[1] << ",\t" << fx4[1] << ")" << std::endl;
        pcout << "!!! Fz-Axis (" << fx1[2] << ",\t" << fx2[2] << ",\t" << fx3[2] << ",\t" << fx4[2] << ")" << std::endl;
        Array<T,3> ni(0.,0.,0.),  nj(0.,0.,0.);
        crossProduct(x1-x3,  x2-x3, ni);
        crossProduct(x3-x1,  x4-x1, nj);
        T edgeAngle = angleBetweenVectors(ni, nj);
        pcout << "!!! Angle " << edgeAngle*180/pi << ", eqAngle " << eqAngle*180/pi <<std::endl;
    }
//    fx1 = fx1 - df/4.;
//    fx2 = fx2 - df/4.;
//    fx3 = fx3 - df/4.;
//    fx4 = fx4 - df/4.;
    return fx1;
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
    // T factor = k * sin(edgeAngle - eqAngle) * (-1.0/(1.0 - sqrt(ninj*ninj)) );
    Array<T,3> fx1(0,0,0), dninjdx(0,0,0);
    Array<T,3> dx23(0,0,0), dx34(0,0,0);
    Array<T,3> dx31(0,0,0), dx13(0,0,0);
    Array<T,3> dx12(0,0,0), dx41(0,0,0);

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
//    fx3 = - fx1 - fx2 - fx4;
    crossProduct(x1-x2, nj - ninj*ni, dx12);
    crossProduct(x4-x1, ni - ninj*nj, dx41);
    dninjdx = 0.5/Ai*dx23 + 0.5/Aj*dx34;
    fx3 = -factor*dninjdx;

    return fx1;
}




*/
