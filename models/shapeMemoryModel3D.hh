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

#include "shapeMemoryModel3D.h"


namespace plb {

template<typename T, template<typename U> class Descriptor>
ShapeMemoryModel3D<T,Descriptor>::ShapeMemoryModel3D (
        T density_, T k_rest_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
        T k_volume_, T k_surface_, T eta_m_,
        vector<T> eqArea_, map<plint,T> eqLength_, map<plint,T> eqAngle_,
        T eqVolume_, T eqSurface_, T eqTileSpan_,
        T persistenceLengthFine, T eqLengthRatio_,
        std::map<plint, Particle3D<T,Descriptor>*> & tagToParticle3D_,
        pluint cellNumTriangles_, pluint cellNumVertices_)
    : RBCModel3D<T>(density_),
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
      tagToParticle3D(tagToParticle3D_),
      cellNumTriangles(cellNumTriangles_),
      cellNumVertices(cellNumVertices_)
{
    persistenceLengthCoarse = persistenceLengthFine * sqrt( (cellNumVertices-2.0) / (23867-2.0)) ;

    T eqMeanArea = eqArea = eqSurface/cellNumTriangles;
    eqAngle=0.0;

    typename map<plint,T>::reverse_iterator iter = eqAnglePerEdge.rbegin();
    for (iter = eqAnglePerEdge.rbegin(); iter != eqAnglePerEdge.rend(); ++iter) {
       eqAngle += iter->second;
    }
    eqAngle /= eqAnglePerEdge.size();

    iter = eqLengthPerEdge.rbegin();
    for (iter = eqLengthPerEdge.rbegin(); iter != eqLengthPerEdge.rend(); ++iter) {
       eqLength += iter->second;
    }
    eqLength /= eqLengthPerEdge.size();

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
//    C_elastic =  kBT*k_elastic * (3.0*sqrt(3.0)*pow(eqLength,3) * eqLengthRatio *
//                    (6 - 9*eqLengthRatio + 4*pow(eqLengthRatio,2))) / (64.*persistenceLengthCoarse);
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


template<typename T, template<typename U> class Descriptor>
plint ShapeMemoryModel3D<T,Descriptor>::getTriangleId(plint iTriangle)
{
    return iTriangle % cellNumTriangles;
};


template<typename T, template<typename U> class Descriptor>
plint ShapeMemoryModel3D<T,Descriptor>::getEdgeId(plint iVertex, plint jVertex) {
    iVertex = iVertex % cellNumVertices;
    jVertex = jVertex % cellNumVertices;
    if (iVertex > jVertex){
        return (iVertex*(iVertex - 1))/2 + jVertex;
    } else if (iVertex < jVertex) {
        return (jVertex*(jVertex - 1))/2 + iVertex;
    }
    return -1;
};


template<typename T, template<typename U> class Descriptor>
Array<T,3> ShapeMemoryModel3D<T,Descriptor>::computeCellForce (
        TriangleBoundary3D<T> const& boundary,
        T cellVolume, T cellSurface, T & iSurface,
        plint iVertex)
 /* Some force calculations are according to KrugerThesis, Appendix C */
{
    if (not ((cellVolume > 0) and (cellSurface > 0))) {
        cout << "iVertex: " << iVertex
             << ", processor: " << dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (tagToParticle3D[iVertex])->get_processor()
             << ", cellId: " << dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (tagToParticle3D[iVertex])->get_cellId()
             << ", volume: " << cellVolume
             << ", surface: " << cellSurface
             << endl;
        PLB_PRECONDITION( (cellVolume > 0) and (cellSurface > 0) );
    }
    /* Initializations from input */
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    Array<T,3> x1 = dynMesh.getVertex(iVertex), x2, x3, x4;
    ImmersedCellParticle3D<T,Descriptor>* iParticle = dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (tagToParticle3D[iVertex]);
    Array<T,3> iVelocity = iParticle->get_v();
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

    iParticle->get_E_volume() = 0.5*volumeCoefficient*(cellVolume - eqVolume)*1.0/cellNumVertices;
    iParticle->get_E_area() = 0.5*surfaceCoefficient*(cellSurface - eqSurface)*1.0/cellNumVertices;

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
        T eqArea = eqAreaPerTriangle[getTriangleId(iTriangle)];
        shearForce += computeLocalAreaConservationForce(dAdx, triangleArea, eqArea, areaCoefficient);
        iParticle->get_E_area() += 0.5*areaCoefficient*pow(eqArea - triangleArea, 2)/3.0; // 3 vertices on each triangle

        /* Elastice Force */
//        elasticForce += computeElasticRepulsiveForce(dAdx, triangleArea, C_elastic);
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
        plint edgeId = getEdgeId(iVertex, jVertex);
        T eqLength = eqLengthPerEdge[edgeId];
        /* In Plane (WLC) and repulsive forces*/
        T potential;
        inPlaneForce += computeInPlaneExplicitForce(x1, x3, eqLengthRatio, eqLength, k_inPlane, potential);
        iParticle->get_E_inPlane() += potential / 2.0; // Each spring consists of 2 particles.

        /*  Dissipative Forces Calculations */
//        dissipativeForce += computeDissipativeForce(x1, x3, iVelocity, particleVelocity[jVertex], gamma_T, gamma_C);
        /*  Bending Forces Calculations */
        /*  Bending Forces Calculations */
        T eqAngle = eqAnglePerEdge[edgeId];
        T edgeAngle = calculateSignedAngle(dynMesh, iVertex, jVertex, kVertex, lVertex); //edge is iVertex, jVertex
        Array<T,3> iNormal = dynMesh.computeTriangleNormal(iVertex, jVertex, kVertex);
        Array<T,3> jNormal = dynMesh.computeTriangleNormal(iVertex, jVertex, lVertex);

        // Compute bending force for the vertex as part of the main edge
        bendingForce += computeBendingForceEdge (edgeAngle, eqAngle, k_bend,
				iNormal, jNormal);
        iParticle->get_E_bending() += computeBendingPotential (edgeAngle, eqAngle, k_bend);
        // Compute bending force for the vertex as a lateral vertex
        edgeId = getEdgeId(jVertex, kVertex);
        eqAngle = eqAnglePerEdge[edgeId];
        edgeAngle = calculateSignedAngle(dynMesh, jVertex, kVertex);
        bendingForce += computeBendingForceLateral (edgeAngle, eqAngle, k_bend, iNormal)*0.5; // Is calculated twice, for the other triangle as well
        iParticle->get_E_bending() += computeBendingPotential (edgeAngle, eqAngle, k_bend)*0.5; // Is calculated twice, for the other triangle as well
        // Compute bending force for the vertex as a lateral vertex
        edgeId = getEdgeId(jVertex, lVertex);
        eqAngle = eqAnglePerEdge[edgeId];
        edgeAngle = calculateSignedAngle(dynMesh, jVertex, lVertex);
        bendingForce += computeBendingForceLateral (edgeAngle, eqAngle, k_bend, jNormal)*0.5; // Is calculated twice, for the other triangle as well
        iParticle->get_E_bending() += computeBendingPotential (edgeAngle, eqAngle, k_bend)*0.5; // Is calculated twice, for the other triangle as well
    }
    iParticle->get_f_wlc() = inPlaneForce + repulsiveForce;
    iParticle->get_f_bending() = bendingForce;
    iParticle->get_f_volume() = volumeForce;
    iParticle->get_f_surface() = surfaceForce;
    iParticle->get_f_shear() = shearForce;
    iParticle->get_f_viscosity() = dissipativeForce;
    return (inPlaneForce + elasticForce + repulsiveForce) + bendingForce + (volumeForce + surfaceForce + shearForce) + dissipativeForce;// + stretchForce;
}


template<typename T, template<typename U> class Descriptor>
Array<T,3> ShapeMemoryModel3D<T,Descriptor>::computeElasticForce (
        TriangleBoundary3D<T> const& boundary,
        plint iVertex )
{
    return Array<T,3>(0,0,0);
}

template<typename T, template<typename U> class Descriptor>
ShapeMemoryModel3D<T, Descriptor>* ShapeMemoryModel3D<T,Descriptor>::clone() const {
    return new ShapeMemoryModel3D<T,Descriptor>(*this);
}


template<typename T>
void getCellShapeQuantitiesFromMesh(TriangleBoundary3D<T>& boundary,
                            vector<T> & eqAreaPerTriangle, map<plint,T> & eqLengthPerEdge, map<plint,T> & eqAnglePerEdge,
                            plint cellNumTriangles, plint cellNumPartsPerCell) {
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    for (plint iVertex = 0; iVertex < cellNumPartsPerCell; ++iVertex) {
        std::vector<plint> neighborVertexIds = dynMesh.getNeighborVertexIds(iVertex);
        for (pluint jV = 0; jV < neighborVertexIds.size(); jV++) {
            plint jVertex = neighborVertexIds[jV];
            if (iVertex > jVertex){
                plint edgeId = (iVertex*(iVertex - 1))/2 + jVertex;
                eqLengthPerEdge[edgeId] = dynMesh.computeEdgeLength(iVertex, jVertex);
                eqAnglePerEdge[edgeId] = calculateSignedAngle(dynMesh, iVertex, jVertex);
            }
        }
    }
    for (plint iTriangle = 0; iTriangle < cellNumTriangles; ++iTriangle) {
        eqAreaPerTriangle[iTriangle] = dynMesh.computeTriangleArea(iTriangle);
    }
}

}  // namespace plb

#endif  // SHAPE_MEMORY_MODEL_3D_HH
