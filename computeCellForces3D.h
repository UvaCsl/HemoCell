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

#include "core/globalDefs.h"
#include "offLattice/triangleBoundary3D.h"
#include "shellModel3D.h"


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


template<typename T>
Array<T,3> computeBendingForce_Krueger (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k,
                                Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4);


/* The most messy force! */
template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k);

}  // namespace plb

#endif  // COMPUTE_CELL_FORCES3D_H
