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

#include "computeCellForces3D.h"


namespace plb {


template<typename T>
Array<T,3> computeInPlaneForce(Array<T,3> const& x1, Array<T,3> const& x2, T maxLength, T k_WLC, T k_rep, T & potential) {
/*
 * ** DOES NOT COMPUTE POTENTIAL **
 *
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

    potential = 0.0;
//    potential = k_WLC * ( 3.0*pow(r, 2) - 2.0*pow(r, 3) ) * 1.0 / (1.0 - r);
//    potential += - k_rep / pow(L,1);


    return tmpForce;
}


template<typename T>
Array<T,3> computeInPlaneExplicitForce(Array<T,3> const& x1, Array<T,3> const& x2, T eqLengthRatio, T eqLength, T k_inPlane, T & potential) {
/*
 *  Computes In-Plane forces based on Worm-like chain forces and a repulsive potential.
 *  Has the same function as computeInPlaneForce but with different arguments.
 *
 *    eqLengthRatio = eqLength/maxLength ; // Used to be x0
 *    k_inPlane = kBT /(4.0*persistenceLengthCoarse);
 *
 *  Related publications: [FedosovCaswellKarniadakis2010, FedosovCaswell2010b, Pivkin2008]

*/
    T Lmax = eqLength*eqLengthRatio;
    Array<T,3> dL = (x1 - x2)*1.0;
    T L = norm(dL);
    T r = L/Lmax, r0=1.0/eqLengthRatio; // T Lmax = eqLength*eqLengthRatio;
    Array<T,3> eij = dL/L;
    /* In Plane Force (WLC) and Repulsive Force */
    Array<T,3> tmpForce =  eij * k_inPlane * r * (-6 + (9 - 4*r)*r)/( (r-1)*(r-1) );
    T k_rep = (eqLength*eqLength)* k_inPlane * r0 * (-6 + (9 - 4*r0)*r0)/( (r0-1)*(r0-1) );
    tmpForce +=  - eij * k_rep / (L*L) ;
//    Array<T,3> tmpForce =  eij * k_inPlane *
//                    (1 - (4*L*eqLengthRatio)/eqLength +
//                            pow(eqLength,2)*((-1 + pow(-1 + eqLengthRatio,-2) + 4*eqLengthRatio)/pow(L,2) - pow(eqLength - L*eqLengthRatio,-2)));
    potential = k_inPlane * Lmax * ( 3.0*pow(r, 2) - 2.0*pow(r, 3) ) * 1.0 / (1.0 - r);
    potential +=  k_rep * 1.0 / L;
    potential -= k_inPlane * Lmax * ( 3.0*pow(r0, 2) - 2.0*pow(r0, 3) ) * 1.0 / (1.0 - r0) + k_rep * 1.0 / eqLength; // Energy residual from potential calculation (potential at rest state)
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
        T Lmax = eqLength*eqLengthRatio;
        Array<T,3> dL = (x1 - x2)*1.0;
        T L = norm(dL);
        T r = L/Lmax, r0=1.0/eqLengthRatio; // T Lmax = eqLength*eqLengthRatio;
        Array<T,3> eij = dL/L;
        /* In Plane Force (WLC) and Repulsive Force */
        Array<T,3> tmpForce =  eij * k_inPlane * r * (-6 + (9 - 4*r)*r)/( (r-1)*(r-1) );
        T k_rep = (eqLength*eqLength)* k_inPlane * r0 * (-6 + (9 - 4*r0)*r0)/( (r0-1)*(r0-1) );
        tmpForce +=  - eij * k_rep / (L*L) ;
    return tmpForce;
}

template<typename T>
Array<T,3> computeInPlaneForce(Array<T,3> const& x1, Array<T,3> const& x2, T maxLength, T k_WLC, T k_rep) {
/*
 *  Wrapper when calling without potential

*/
    T potential = 0.0;
    return computeInPlaneExplicitForce(x1,x2, maxLength, k_WLC, k_rep, potential);
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


/* The most messy force! */
/* Calculates the bending force for the triangles formed by the vertices:
 * (x1, x2, x3) and (x1,x3,x4) with the common edge (x1,x3).
 * eqAngle is expected to be between [-pi,pi].  */
template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqArea, T eqLength, T eqAngle, T k,
								Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4) {
/* The most messy force!
 * Triangles are:
 *      (i, j, k) and (l, k, j). These triangles share
 *      (x2, x1, x3) and (x4, x3, x1). These triangles share
 *      the common edge j-k or 1-3.
 *
 *      crossProduct(jPosition - iPosition, kPosition - jPosition, nijk);
 *      crossProduct(kPosition - lPosition, jPosition - kPosition, nlkj);
 *      crossProduct(x1 - x2, x3 - x1, nijk);
 *      crossProduct(x3 - x4, x1 - x3, nlkj);
*/
	Array<T,3> fx1, tmp;
	T dAngle;
	T edgeAngle = angleBetweenVectors(ni, nj);
	plint sign = dot(x2-x1, nj) > 0?1:-1;
	if (sign <= 0) {
		edgeAngle = 2*pi-edgeAngle;
	}
	edgeAngle = (edgeAngle > pi)?edgeAngle-2*pi:edgeAngle;
	eqAngle = (eqAngle > 2*pi)?eqAngle-2*pi:eqAngle;
	eqAngle = (eqAngle > pi)?eqAngle-2*pi:eqAngle;

	dAngle = (edgeAngle-eqAngle);
//	dAngle = (edgeAngle);
//	dAngle = pow(dAngle, 3);
//	dAngle = dAngle + pow(dAngle, 3) + tan(edgeAngle/2)/50;
	fx2 = -k*dAngle*ni * (eqLength*0.5/eqArea);
    fx4 = -k*dAngle*nj * (eqLength*0.5/eqArea);
    fx1 = -(fx2+fx4)*0.5;
    fx3 = fx1;
    return fx1;
}

/* The most messy force! */
/* Calculates the bending force for the triangles formed by the vertices:
 * (x1, x2, x3) and (x1,x3,x4) with the common edge (x1,x3).
 * eqAngle is expected to be between [-pi,pi].  */
template<typename T>
Array<T,3> computeBendingForce_Krueger (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqArea, T eqLength, T eqAngle, T k,
                                Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4) {
/* The most messy force!
 * Triangles are:
 *      (i, j, k) and (l, k, j). These triangles share
 *      the common edge j-k.
 *      (x2, x1, x3) and (x4, x3, x1). These triangles share
 *      the common edge x1-x3.
 *
 *      crossProduct(jPosition - iPosition, kPosition - jPosition, nijk);
 *      crossProduct(kPosition - lPosition, jPosition - kPosition, nlkj);
 *      crossProduct(x1 - x2, x3 - x1, nijk);
 *      crossProduct(x3 - x4, x1 - x3, nlkj);
*/
    Array<T,3> fx1, tmp;
    T dAngle;
    T edgeAngle = angleBetweenVectors(ni, nj);
    plint sign = dot(x2-x1, nj) > 0?1:-1;
    if (sign <= 0) {
        edgeAngle = 2*pi-edgeAngle;
    }
    edgeAngle = (edgeAngle > pi)?edgeAngle-2*pi:edgeAngle;
    eqAngle = (eqAngle > 2*pi)?eqAngle-2*pi:eqAngle;
    eqAngle = (eqAngle > pi)?eqAngle-2*pi:eqAngle;
    dAngle = (edgeAngle-eqAngle);

    T niDotnj = dot(ni,nj);
    T dAngledxkFactor = -1.0/(sqrt(1.0 - niDotnj*niDotnj ) );

    // f1
    crossProduct(x2 - x3, nj - ni*niDotnj, tmp);
    fx1 = -k * dAngle * dAngledxkFactor * (0.5/Ai) * tmp;
    crossProduct(x3 - x4, ni - nj*niDotnj, tmp);
    fx1 = fx1 + (-k * dAngle * dAngledxkFactor * (0.5/Aj) * tmp);
    // f2
    crossProduct(x3 - x1, nj - ni*niDotnj, tmp);
    fx2 = -k * dAngle * dAngledxkFactor * (0.5/Ai) * tmp;
    // f4
    crossProduct(x1 - x3, ni - nj*niDotnj, tmp);
    fx4 = -k * dAngle * dAngledxkFactor * (0.5/Aj) * tmp;
    // f3
    fx3 = -(fx1+fx2+fx3);
    return fx1;
}



/* Calculates the bending forces.
 * Angles are expected to be between [-pi,pi]. */
/* The most messy force! */
template<typename T>
Array<T,3> computeBendingForces (T edgeAngle, T eqAngle, T k,
								Array<T,3> const& ni, Array<T,3> const& nj,
								Array<T,3> & fx2, Array<T,3> & fx3, Array<T,3> & fx4) {
	Array<T,3> fx1;
	edgeAngle = (edgeAngle > pi)?edgeAngle-2*pi:edgeAngle;
	eqAngle = (eqAngle > 2*pi)?eqAngle-2*pi:eqAngle;
	eqAngle = (eqAngle > pi)?eqAngle-2*pi:eqAngle;
	T dAngle = (edgeAngle-eqAngle);
	fx2 = -k*dAngle*ni;
    fx4 = -k*dAngle*nj;
    fx1 = -(fx2+fx4)*0.5;
    fx3 = fx1;
    return fx1;
}

template<typename T>
Array<T,3> computeBendingForceEdge (T edgeAngle, T eqAngle, T k,
								Array<T,3> const& ni, Array<T,3> const& nj) {
//	Array<T,3> fx1, fx2, fx3, fx4;
	edgeAngle = (edgeAngle > pi)?edgeAngle-2*pi:edgeAngle;
	eqAngle = (eqAngle > 2*pi)?eqAngle-2*pi:eqAngle;
	eqAngle = (eqAngle > pi)?eqAngle-2*pi:eqAngle;
	T dAngle = (edgeAngle-eqAngle);
	return k*dAngle*(ni+nj)*0.5;
}


template<typename T>
Array<T,3> computeBendingForceLateral (T edgeAngle, T eqAngle, T k,
								Array<T,3> const& ni) {
	Array<T,3> fx1;
	edgeAngle = (edgeAngle > pi)?edgeAngle-2*pi:edgeAngle;
	eqAngle = (eqAngle > 2*pi)?eqAngle-2*pi:eqAngle;
	eqAngle = (eqAngle > pi)?eqAngle-2*pi:eqAngle;
	T dAngle = (edgeAngle-eqAngle);
	return -k*dAngle*ni;
}


/*
 * Calculates the bending potential.
 * Angles are expected to be between [-pi,pi].
 */
template<typename T>
T computeBendingPotential (T edgeAngle, T eqAngle, T k) {
	edgeAngle = (edgeAngle > pi)?edgeAngle-2*pi:edgeAngle;
	eqAngle = (eqAngle > 2*pi)?eqAngle-2*pi:eqAngle;
	eqAngle = (eqAngle > pi)?eqAngle-2*pi:eqAngle;
	return 0.5*k*(edgeAngle-eqAngle)*(edgeAngle-eqAngle);
}


/*
 * Helper function, calculates the angle between -pi and pi.
 * The edge is iVertex-jVertex.
 */
template<typename T>
T calculateSignedAngle(TriangularSurfaceMesh<T> const& mesh, plint iVertex, plint jVertex, plint & kVertex, plint & lVertex) {
    Array<T,3> x1 = mesh.getVertex(iVertex), x2(0.,0.,0.), x3(0.,0.,0.), x4(0.,0.,0.);

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

    Array<T,3> V1 = mesh.computeTriangleNormal(iTriangle);
    Array<T,3> V2 = mesh.computeTriangleNormal(jTriangle);
    T angle = angleBetweenVectors(V1, V2);
	plint sign = dot(x2-x1, V2) >= 0?1:-1;
	if (sign <= 0) {
		angle = 2*pi-angle;
	}
	angle = (angle > pi)?angle-2*pi:angle;
	return angle;
}


/*
 * Helper function, calculates the angle between -pi and pi
 */
template<typename T>
T calculateSignedAngle(TriangularSurfaceMesh<T> const& mesh, plint iVertex, plint jVertex) {
	plint kVertex, lVertex;
	return calculateSignedAngle(mesh, iVertex, jVertex, kVertex, lVertex);
}


}  // namespace plb

#endif  // COMPUTE_CELL_FORCES3D_HH
