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
ShapeMemoryModel3D<T,Descriptor>::ShapeMemoryModel3D (T density_, T k_rest_,
        T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
        T k_volume_, T k_surface_, T eta_m_,
        T persistenceLengthFine, T eqLengthRatio_,
        T dx_, T dt_, T dm_,
        TriangularSurfaceMesh<T> const& meshElement)
    : ConstitutiveModel<T,Descriptor>(density_),
      k_rest(k_rest_),
      k_shear(k_shear_),
      k_bend(k_bend_),
      k_stretch(k_stretch_),
      k_elastic(k_elastic_),
      k_surface(k_surface_),
      k_volume(k_volume_),
      eta_m(eta_m_),
      eqLengthRatio(eqLengthRatio_),
      dx(dx_), dt(dt_), dm(dm_)
{
    T dNewton = (dm*dx/(dt*dt)) ;
    k_WLC_ *= 1.0;     k_elastic *= 1.0;     k_bend *= 1.0;
    k_volume *= 1.0;     k_surface *= 1.0;     k_shear *= 1.0;
    eta_m /= dNewton*dt/dx;     k_stretch /= dNewton;    k_rest /= dNewton/dx;

    T x0 = eqLengthRatio;
    syncRequirements.insert(volumeAndSurfaceReductions);

    cellNumTriangles = meshElement.getNumTriangles();
    T cellNumPartsPerCell = meshElement.getNumVertices();
    getCellShapeQuantitiesFromMesh(meshElement, eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, cellNumTriangles, cellNumPartsPerCell);

    persistenceLengthCoarse = persistenceLengthFine * sqrt( (cellNumVertices-2.0) / (23867-2.0)) ;

    T eqMeanArea = eqArea = eqSurface/cellNumTriangles;
    eqAngle=0.0;

    typename map<plint,T>::reverse_iterator iter = eqAnglePerEdge.rbegin();
    for (iter = eqAnglePerEdge.rbegin(); iter != eqAnglePerEdge.rend(); ++iter) {
       eqAngle += iter->second;
    }
    eqAngle /= eqAnglePerEdge.size();

    maxLength = 0;
    iter = eqLengthPerEdge.rbegin();
    for (iter = eqLengthPerEdge.rbegin(); iter != eqLengthPerEdge.rend(); ++iter) {
       eqLength += iter->second;
       maxLength = max(iter->second, maxLength);
    }
    eqLength /= eqLengthPerEdge.size();
    /* Calculate cell Radius */
    Array<T,2> xRange;     Array<T,2> yRange;     Array<T,2> zRange;
    meshElement.computeBoundingBox (xRange, yRange, zRange);
    cellRadiusLU = max(xRange[1] - xRange[0], yRange[1] - yRange[0]);
    cellRadiusLU = max(cellRadiusLU , zRange[1] - zRange[0]);

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
void ShapeMemoryModel3D<T, Descriptor>::computeCellForce (Cell3D<T,Descriptor> & cell, T ratio) {
     /* Some force calculations are according to KrugerThesis, Appendix C */
    T cellVolume = cell.getVolume();
    T cellSurface = cell.getSurface();
    if (not ((cellVolume > 0) and (cellSurface > 0))) {
        cout << ", processor: " << cell.getMpiProcessor()
             << ", cellId: " << cell.get_cellId()
             << ", volume: " << cellVolume
             << ", surface: " << cellSurface
             << endl;
        PLB_PRECONDITION( (cellVolume > 0) and (cellSurface > 0) );
    }
    std::vector<plint> const& triangles = cell.getTriangles();
    std::vector<Array<plint,2> > const& edges = cell.getEdges();
    std::vector<plint > const& vertices = cell.getVertices();
    for (int iV = 0; iV < vertices.size(); ++iV) {
        castParticleToICP3D(cell.getParticle3D(vertices[iV]))->resetForces();
    }
    plint iTriangle;
    plint iVertex, jVertex, kVertex, lVertex;

    /* Run through all the edges and calculate:
         x In plane (WLC and repulsive) force
         x Dissipative force
         x Bending force
         o Stretch force
     */
    Array<T,3> force1, force2;
    T potential;
    for (int iE = 0; iE < edges.size(); ++iE) {
        iVertex = edges[iE][0];  jVertex = edges[iE][1];
        plint edgeId = getEdgeId(iVertex, jVertex);
        Array<T,3> const& iX = cell.getVertex(iVertex);
        Array<T,3> const& jX = cell.getVertex(jVertex);
        ImmersedCellParticle3D<T,Descriptor>* iParticle = castParticleToICP3D(cell.getParticle3D(iVertex));
        ImmersedCellParticle3D<T,Descriptor>* jParticle = castParticleToICP3D(cell.getParticle3D(jVertex));
          /* ------------------------------------*/
         /* In Plane forces (WLC and repulsive) */
        /* ------------------------------------*/
        force1 = computeInPlaneExplicitForce(iX, jX, eqLengthRatio, eqLengthPerEdge[edgeId], k_inPlane, potential);
        iParticle->get_force() += force1;
        jParticle->get_force() -= force1;

        iParticle->get_E_inPlane() += potential;
        jParticle->get_E_inPlane() += potential;
        iParticle->get_f_wlc() += force1;
        jParticle->get_f_wlc() -= force1;
          /* ------------------------------------*/
         /*    Dissipative Forces Calculations  */
        /* ------------------------------------*/
        if (gamma_T>0.0) {
            force1 = computeDissipativeForce(iX, jX, iParticle->get_v(), jParticle->get_v(), gamma_T, gamma_C);
            iParticle->get_force() += force1;
            jParticle->get_force() -= force1;

            iParticle->get_f_viscosity() += force1;
            jParticle->get_f_viscosity() -= force1;
        }
        /* -------------------------------------------*/
        /*    Stretch (Hookean) Forces Calculations  */
        /* -----------------------------------------*/
        if (k_stretch>0.0) {

        }
          /* ------------------------------------*/
         /*    Bending Forces Calculations      */
        /* ------------------------------------*/
        bool angleFound;
        plint kVertex, lVertex;
        T edgeAngle = cell.computeSignedAngle(iVertex, jVertex, kVertex, lVertex, angleFound); //edge is iVertex, jVertex
        if (angleFound) {
            Array<T,3> iNormal = cell.computeTriangleNormal(iVertex, jVertex, kVertex);
            Array<T,3> jNormal = cell.computeTriangleNormal(iVertex, jVertex, lVertex);
            ImmersedCellParticle3D<T,Descriptor>* kParticle = castParticleToICP3D(cell.getParticle3D(kVertex));
            ImmersedCellParticle3D<T,Descriptor>* lParticle = castParticleToICP3D(cell.getParticle3D(lVertex));
            Array<T,3> const& kX = cell.getVertex(kVertex);
            Array<T,3> const& lX = cell.getVertex(lVertex);

            /*== Compute bending force for the vertex as part of the main edge ==*/
            force1 = computeBendingForceEdge (edgeAngle, eqAnglePerEdge[edgeId], k_bend, iNormal, jNormal);
            Array<T,3> fi, fk, fj, fl;
            T Ai=0, Aj=0; // Not necessary for this calculation
            computeBendingForce (iX, kX, jX, lX, iNormal, jNormal, Ai, Aj, eqTileSpan, eqLengthPerEdge[edgeId], eqAnglePerEdge[edgeId], k_bend, fk, fj, fl);

            iParticle->get_force() += fi;
            jParticle->get_force() += fj;
            kParticle->get_force() += fk;
            lParticle->get_force() += fl;

            potential = computeBendingPotential (edgeAngle, eqAnglePerEdge[edgeId], k_bend);
            iParticle->get_E_bending() += potential;
            jParticle->get_E_bending() += potential;
            kParticle->get_E_bending() += potential;
            lParticle->get_E_bending() += potential;
        }
    }

    /* ===================== In case of quasi-rigid object =====================
     *
     * If this is a boundary element (k_rest != 0), get the reference locations
     * of iVertex and calculate and return the force for quasi-rigid objects.
     *          (FengMichaelides2004, J.Comp.Phys. 195(2))
     * // CURRENTLY UNAVAILABLE
     * */

    /* Calculate cell coefficients */
    T volumeCoefficient = k_volume * (cellVolume - eqVolume)*1.0/eqVolume;
    T surfaceCoefficient = k_surface * (cellSurface - eqSurface)*1.0/eqSurface;
    T eqMeanArea = eqSurface/cellNumTriangles;
    T areaCoefficient = k_shear/eqMeanArea ;

//    iParticle->get_E_volume() = 0.5*volumeCoefficient*(cellVolume - eqVolume)*1.0/cellNumVertices;
//    iParticle->get_E_area() = 0.5*surfaceCoefficient*(cellSurface - eqSurface)*1.0/cellNumVertices;


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
    plint vertexIds[3];
    for (int iT = 0; iT < triangles.size(); ++iT) {
        iTriangle = triangles[iT];
        triangleNormal = cell.computeTriangleNormal(iTriangle);
        triangleArea = cell.computeTriangleArea(iTriangle);
        iVertex = cell.getVertexId(iTriangle,0);
        jVertex = cell.getVertexId(iTriangle,1);
        kVertex = cell.getVertexId(iTriangle,2);
        Array<T,3> const& x1 = cell.getVertex(iVertex);
        Array<T,3> const& x2 = cell.getVertex(jVertex);
        Array<T,3> const& x3 = cell.getVertex(kVertex);
        ImmersedCellParticle3D<T,Descriptor>* iParticle = castParticleToICP3D(cell.getParticle3D(iVertex));
        ImmersedCellParticle3D<T,Descriptor>* jParticle = castParticleToICP3D(cell.getParticle3D(jVertex));
        ImmersedCellParticle3D<T,Descriptor>* kParticle = castParticleToICP3D(cell.getParticle3D(kVertex));

        /* Surface and local area conservation forces */
        force1  = computeSurfaceConservationForce(x1, x2, x3, triangleNormal, surfaceCoefficient, dAdx);
        force1 += computeLocalAreaConservationForce(dAdx, triangleArea, eqAreaPerTriangle[iTriangle], areaCoefficient);
        force2  = computeSurfaceConservationForce(x2, x3, x1, triangleNormal, surfaceCoefficient, dAdx);
        force2 += computeLocalAreaConservationForce(dAdx, triangleArea, eqAreaPerTriangle[iTriangle], areaCoefficient);
        iParticle->get_force() += force1;
        jParticle->get_force() += force2;
        kParticle->get_force() -= (force1+force2);
        /* Volume conservation forces */
        force1  = computeVolumeConservationForce(x1, x2, x3, volumeCoefficient);
        force2  = computeVolumeConservationForce(x2, x3, x1, volumeCoefficient);
        iParticle->get_force() += force1;
        jParticle->get_force() += force2;
        kParticle->get_force() -= (force1+force2);
    }

}


template<typename T, template<typename U> class Descriptor>
ShapeMemoryModel3D<T, Descriptor>* ShapeMemoryModel3D<T,Descriptor>::clone() const {
    return new ShapeMemoryModel3D<T,Descriptor>(*this);
}


template<typename T>
void getCellShapeQuantitiesFromMesh(TriangularSurfaceMesh<T> const& dynMesh,
                            vector<T> & eqAreaPerTriangle, map<plint,T> & eqLengthPerEdge, map<plint,T> & eqAnglePerEdge,
                            plint cellNumTriangles, plint cellNumPartsPerCell) {
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
