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

#ifndef SHAPE_MEMORY_MODEL_3D_HH
#define SHAPE_MEMORY_MODEL_3D_HH

#include <cmath>
#include <map>
#include "shapeMemoryModel3D.h"
#include "computeCellForces3D.hh"



namespace plb {

template<typename T>
ShapeMemoryModel3D<T>::ShapeMemoryModel3D (
        T density_, T k_rest_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
        T k_volume_, T k_surface_, T eta_m_,
        vector<T> eqArea_, map<plint,T> eqLength_, map<plint,T> eqAngle_,
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
      eqAreaPerTriangle(eqArea_),
      eqLengthPerEdge(eqLength_),
      eqAnglePerEdge(eqAngle_),
      eqVolume(eqVolume_),
      eqSurface(eqSurface_),
      eqTileSpan(eqTileSpan_),
      eqLengthRatio(eqLengthRatio_),
      cellNumTriangles(cellNumTriangles_),
      cellNumVertices(cellNumVertices_)
{
    persistenceLengthCoarse = persistenceLengthFine * sqrt( (cellNumVertices-2.0) / (23867-2.0)) ;

    T eqMeanArea = eqSurface/cellNumTriangles;
    T eqLength = sqrt(4*eqMeanArea/sqrt(3.0));
    T eqAngle=0.0;

    typename map<plint,T>::reverse_iterator iter = eqAnglePerEdge.rbegin();
    for (iter = eqAnglePerEdge.rbegin(); iter != eqAnglePerEdge.rend(); ++iter) {
       eqAngle = iter->second;
    }
    eqAngle /= eqAnglePerEdge.size();

    /* Use dimensionless coefficients */

    /* Use Mean eqLength for stiffness coefficients of:
     *  Volume, Surface and Shear.
     *  Shear also uses mean eqArea in the calculation.
     *  */
    k_volume *= kBT/pow(eqLength,3);
    k_surface *= kBT/pow(eqLength,2);
    k_shear *= kBT/pow(eqLength,2);
    k_bend *= kBT;
    /* In plane coefficient initialization */
    k_inPlane = k_WLC_ * kBT /(4.0*persistenceLengthCoarse);
    C_elastic =  kBT*k_elastic * (3.0*sqrt(3.0)*pow(eqLength,3) * eqLengthRatio *
                    (6 - 9*eqLengthRatio + 4*pow(eqLengthRatio,2))) / (64.*persistenceLengthCoarse);
    /* Dissipative term coefficients from FedosovCaswellKarniadakis2010 */
    gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
    gamma_C = (gamma_T/3.0);

    pcout << std::endl;
    pcout << " ============================================= " << std::endl;
    pcout << "k_bend: " << k_bend << ",\t eqAngle (degrees): " << eqAngle*180.0/pi << std::endl;
    pcout << "k_volume: " << k_volume << ",\t eqVolume: " << eqVolume << std::endl;
    pcout << "k_surface: " << k_surface << ",\t eqSurface: " << eqSurface << std::endl;
    pcout << "k_shear: " << k_shear << ",\t eqMeanArea: " << eqMeanArea << std::endl;
    pcout << "eta_m: " << eta_m << ",\t eqLengthRatio: " << eqLengthRatio << std::endl;
    pcout << "gamma_T: " << gamma_T << ",\t persistenceLengthCoarse: " << persistenceLengthCoarse << std::endl;
    pcout << "gamma_C: " << gamma_C << ",\t 0: " << 0 << std::endl;
    pcout << "k_rest: " << k_rest << ",\t 0 : " << 0 << std::endl;
    pcout << "k_stretch: " << k_stretch << ",\t eqTileSpan: " << eqTileSpan << std::endl;
    pcout << "k_elastic: " << k_elastic << ",\t eqLength: " << eqLength << std::endl;
    pcout << "* k_bend: " << k_bend/kBT << std::endl;
    pcout << "* k_volume: " <<  k_volume/(kBT/pow(eqLength,3)) <<  std::endl;
    pcout << "* k_surface: " << k_surface/(kBT/pow(eqLength,2)) <<  std::endl;
    pcout << "* k_shear: " << k_shear/(kBT/pow(eqLength,2)) <<  std::endl;
    pcout << "* eqLength from eqArea: " << sqrt(4*eqMeanArea/sqrt(3.0)) << ",\t eqLength: " << eqLength << std::endl;
    pcout << " ============================================= " << std::endl;
}


template<typename T>
plint ShapeMemoryModel3D<T>::getTriangleId(plint iTriangle)
{
    return iTriangle % cellNumTriangles;
};


template<typename T>
plint ShapeMemoryModel3D<T>::getEdgeId(plint iVertex, plint jVertex) {
    iVertex = iVertex % cellNumVertices;
    jVertex = jVertex % cellNumVertices;
    if (iVertex > jVertex){
        return (iVertex*(iVertex - 1))/2 + jVertex;
    } else if (iVertex < jVertex) {
        return (jVertex*(jVertex - 1))/2 + iVertex;
    }
    return -1;
};


template<typename T>
Array<T,3> ShapeMemoryModel3D<T>::computeCellForce (
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

    plint iTriangle, jTriangle;
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
        T eqArea = eqAreaPerTriangle[getTriangleId(iTriangle)];
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
        T eqLength = eqLengthPerEdge[getEdgeId(iVertex, jVertex)];
        T eqAngle = eqAnglePerEdge[getEdgeId(iVertex, jVertex)];
        /* In Plane (WLC) and repulsive forces*/
        inPlaneForce += computeInPlaneExplicitForce(x1, x3, eqLengthRatio, eqLength, k_inPlane);

        /*  Dissipative Forces Calculations */
        dissipativeForce += computeDissipativeForce(x1, x3, iVelocity, particleVelocity[jVertex], gamma_T, gamma_C);
        /*  Bending Forces Calculations */
        std::vector<plint> triangles = dynMesh.getAdjacentTriangleIds(iVertex, jVertex);
        iTriangle = triangles[0]; jTriangle = triangles[1];
        plint foundVertices = 0;
        for (pluint id = 0; id < 3; ++id) {
            kVertex = dynMesh.getVertexId(iTriangle,id);
            if ( (kVertex != iVertex) && (kVertex != jVertex) ) {
                x2 = dynMesh.getVertex(kVertex);
                foundVertices += 1;
                break;
            }
        }
        for (pluint id = 0; id < 3; ++id) {
            lVertex = dynMesh.getVertexId(jTriangle,id);
            if ( (lVertex != iVertex) && (lVertex != jVertex) ) {
                x4 = dynMesh.getVertex(lVertex);
                foundVertices += 1;
                break;
            }
        }
        PLB_ASSERT(foundVertices == 2); //Assert if some particles are outside of the domain

//        tmpForce = computeBendingForce (x1, x2, x3, x4,
//                            trianglesNormal[iTriangle], trianglesNormal[jTriangle],
//                            trianglesArea[iTriangle], trianglesArea[jTriangle],
//                            eqTileSpan, eqLength, eqAngle, k_bend);
//        T iTriangleBendingCoefficient = trianglesArea[iTriangle] / (trianglesArea[iTriangle] + trianglesArea[jTriangle]);
//        T jTriangleBendingCoefficient = trianglesArea[jTriangle] / (trianglesArea[iTriangle] + trianglesArea[jTriangle]);
//        bendingForce += tmpForce;
//        particleForces[kVertex][1] += -iTriangleBendingCoefficient * tmpForce;
//        particleForces[lVertex][1] += -jTriangleBendingCoefficient * tmpForce;
        Array<T,3> tmp2, tmp3, tmp4;
        tmpForce = computeBendingForce_Krueger (x1, x2, x3, x4,
                            trianglesNormal[iTriangle], trianglesNormal[jTriangle],
                            trianglesArea[iTriangle], trianglesArea[jTriangle],
                            eqTileSpan, eqLength, eqAngle, k_bend,
                            tmp2, tmp3, tmp4);
        bendingForce += tmpForce * 0.5; // Multiplied by 0.5, because this force is calculated twice (x2 also calculates this force)
        particleForces[kVertex][1] += tmp2 * 0.5;
        particleForces[jVertex][1] += tmp3 * 0.5;
        particleForces[lVertex][1] += tmp4 * 0.5;
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
Array<T,3> ShapeMemoryModel3D<T>::computeElasticForce (
        TriangleBoundary3D<T> const& boundary,
        plint iVertex )
{
    return Array<T,3>(0,0,0);
}

template<typename T>
ShapeMemoryModel3D<T>* ShapeMemoryModel3D<T>::clone() const {
    return new ShapeMemoryModel3D<T>(*this);
}



}  // namespace plb

#endif  // SHAPE_MEMORY_MODEL_3D_HH
