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
#include "computeCellForces3D.hh"
#include "cellModel3D.h"


//ShapeMemoryModel3D<T>::ShapeMemoryModel3D (
//        T density_, T k_rest_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
//        T k_volume_, T k_surface_, T eta_m_,
//        vector<T> eqArea_, map<plint,T> eqLength_, map<plint,T> eqAngle_,
//        T eqVolume_, T eqSurface_, T eqTileSpan_,
//        T persistenceLengthFine, T eqLengthRatio_, pluint cellNumTriangles_, pluint cellNumVertices_)

namespace plb {

template<typename T>
CellModel3D<T>::CellModel3D (
        T density_, T k_rest_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
        T k_volume_, T k_surface_, T eta_m_,
        T eqArea_, T eqLength_, T eqAngle_,
        T eqVolume_, T eqSurface_, T eqTileSpan_,
        T persistenceLengthFine, T eqLengthRatio_, pluint cellNumTriangles_, pluint cellNumVertices_)
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
//      maxLength(maxLength_),
      eqLengthRatio(eqLengthRatio_),
      cellNumTriangles(cellNumTriangles_),
      cellNumVertices(cellNumVertices_)
{
    T x0 = eqLengthRatio;
    T maxLength = eqLength*eqLengthRatio;
    persistenceLengthCoarse = persistenceLengthFine * sqrt( (cellNumVertices-2.0) / (23867-2.0)) ;
    /* Use dimensionless coefficients */
    k_volume *= kBT/pow(eqLength,3);
    k_surface *= kBT/pow(eqLength,2);
    k_shear *= kBT/pow(eqLength,2);
    k_bend *= kBT;
    k_inPlane = k_WLC_ * kBT /(4.0*persistenceLengthCoarse);
    // Calculating eqAngle and eqLength according to FedosovCaswellKarniadakis2010
    eqAngle = acos( (sqrt(3.)*(cellNumVertices-2.0) - 5*pi)/(sqrt(3.)*(cellNumVertices-2.0) - 3*pi) );
//    eqLength = sqrt( 2.0*eqArea*1.0/sqrt(3) );
    /* Dissipative term coefficients from FedosovCaswellKarniadakis2010 */
    gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
    gamma_C = (gamma_T/3.0);
    /* The units on the paper are wrong, should have been fixed on config.xml */
    // gamma_T *= eqLength;
    // gamma_C *= eqLength;

    T k_WLC = k_WLC_ * kBT * maxLength/(4.0*persistenceLengthCoarse);
    /* Solving f_WLC + f_rep =0 for x=eqLength, f_rep = k_rep/L^m, m=2. */
    T k_rep = (k_WLC*maxLength*pow(x0,3)*(6 - 9*x0 + 4*pow(x0,2)))/pow(-1 + x0,2);

    Array<T,3> x1(0.,0.,0.), x3(0.,0.,0.);
    x3[0] = eqLength;
    T forceSum = norm(computeInPlaneExplicitForce(x1, x3, eqLengthRatio, eqLength, k_inPlane));
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
        std::map< plint, Array<T,3>* > & particleForces,
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
        return (-k_rest*dx) + (-gamma_T*iVelocity); // Dissipative term from Dupin2007
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
    Array<T,3> tmpForce; tmpForce.resetToZero();
    /* Calculate cell coefficients */
    T volumeCoefficient = k_volume * (cellVolume - eqVolume)*1.0/eqVolume;
    T surfaceCoefficient = k_surface * (cellSurface - eqSurface)*1.0/eqSurface;
    T eqMeanArea = eqSurface/cellNumTriangles;
    T areaCoefficient = k_shear/eqMeanArea ;
    iSurface = 0.0;

    /* Run through all the neighbouring faces of iVertex and calculate:
         x Volume conservation force
         x Surface conservation force
         x Shear force
     */
    Array<T,3> dAdx, dVdx, tmp(0,0,0);
    std::map<plint, T> trianglesArea;
    std::map<plint, Array<T,3> > trianglesNormal;
    T triangleArea;
    Array<T,3> triangleNormal;

    plint iTriangle;
    plint jVertex, kVertex, lVertex;
    plint vertexIds[3];

    std::vector<plint> neighborTriangles = dynMesh.getNeighborTriangleIds(iVertex);
    for (pluint iB = 0; iB < neighborTriangles.size(); ++iB) {
        iTriangle = neighborTriangles[iB];
        for (pluint id = 0; id < 3; ++id) {
            vertexIds[id] = dynMesh.getVertexId(iTriangle,id);
        }
        plint localVertexId = (iVertex == vertexIds[1])*1 + (iVertex == vertexIds[2])*2;
        jVertex = vertexIds[(localVertexId + 1)%3];
        kVertex = vertexIds[(localVertexId + 2)%3];
        x2 = dynMesh.getVertex(jVertex);
        x3 = dynMesh.getVertex(kVertex);
        triangleNormal = trianglesNormal[iTriangle] = dynMesh.computeTriangleNormal(iTriangle);
        triangleArea = trianglesArea[iTriangle] = dynMesh.computeTriangleArea(iTriangle);
        iSurface += triangleArea/3.0;
        /* Surface conservation force */
        surfaceForce += computeSurfaceConservationForce(x1, x2, x3, triangleNormal, surfaceCoefficient, dAdx);
        /* Shear force */
        shearForce += computeLocalAreaConservationForce(dAdx, triangleArea, eqArea, areaCoefficient);
        /* Elastice Force */
        elasticForce += computeElasticRepulsiveForce(dAdx, triangleArea, C_elastic);
        /* Volume conservation force */
        volumeForce  += computeVolumeConservationForce(x1, x2, x3, volumeCoefficient);
        }

    /* Run through all the neighbouring vertices of iVertex and calculate:
         x In plane (WLC) force
         x Repulsive force
         o Stretch force
         x Dissipative force
         x Bending force
     */
    std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
    for (pluint jV = 0; jV < neighborVertexIds.size(); jV++) {
        jVertex = neighborVertexIds[jV];
        x3 = dynMesh.getVertex(jVertex);
        /* In Plane (WLC) and repulsive forces*/
        inPlaneForce += computeInPlaneExplicitForce(x1, x3, eqLengthRatio, eqLength, k_inPlane);
        /*  Dissipative Forces Calculations */
        dissipativeForce += computeDissipativeForce(x1, x3, iVelocity, particleVelocity[jVertex], gamma_T, gamma_C);
        /*  Bending Forces Calculations */
        T edgeAngle = calculateSignedAngle(dynMesh, iVertex, jVertex, kVertex, lVertex); //edge is iVertex, jVertex
        Array<T,3> ftmp2, ftmp3, ftmp4;
        Array<T,3> iNormal = dynMesh.computeTriangleNormal(iVertex, jVertex, kVertex);
        Array<T,3> jNormal = dynMesh.computeTriangleNormal(iVertex, jVertex, lVertex);
		tmpForce = computeBendingForce (edgeAngle, eqAngle, k_bend,
				iNormal, jNormal,
                ftmp2, ftmp3, ftmp4);
//        tmpForce = computeBendingForceFromPotential (x1, x2, x3, x4,
//                            eqTileSpan, eqLength, eqAngle, k_bend,
//                            tmp2, tmp3, tmp4);
        bendingForce += tmpForce * 0.5; // Multiplied by 0.5, because this force is calculated twice (x2 also calculates this force)
        particleForces[kVertex][1] += ftmp2 * 0.5;
        particleForces[jVertex][1] += ftmp3 * 0.5;
        particleForces[lVertex][1] += ftmp4 * 0.5;

        T t1 = computeBendingPotential (edgeAngle, eqAngle, k_bend);
        particleForces[iVertex][6][0] += t1;
        particleForces[kVertex][6][0] += t1;
        particleForces[jVertex][6][0] += t1;
        particleForces[lVertex][6][0] += t1;
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
    return Array<T,3>(0,0,0);
}

template<typename T>
CellModel3D<T>* CellModel3D<T>::clone() const {
    return new CellModel3D<T>(*this);
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
