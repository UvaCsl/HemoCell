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

#include "cellModel3D.h"


//ShapeMemoryModel3D<T>::ShapeMemoryModel3D (
//        T density_, T k_rest_, T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
//        T k_volume_, T k_surface_, T eta_m_,
//        vector<T> eqArea_, map<plint,T> eqLength_, map<plint,T> eqAngle_,
//        T eqVolume_, T eqSurface_, T eqTileSpan_,
//        T persistenceLengthFine, T eqLengthRatio_, pluint cellNumTriangles_, pluint cellNumVertices_)

namespace plb {

template<typename T, template<typename U> class Descriptor>
CellModel3D<T, Descriptor>::CellModel3D(T density_, T k_rest_,
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
    MeshMetrics<T> meshmetric(meshElement);

    cellNumVertices = meshmetric.getNumVertices();
    cellNumTriangles = meshmetric.getNumTriangles();
    cellRadiusLU = meshmetric.getRadius();
    eqLength = meshmetric.getMeanLength();
    maxLength = eqLength*eqLengthRatio;
    eqArea = meshmetric.getMeanArea();
//    eqAngle = meshmetric.getMeanAngle();
    eqVolume = meshmetric.getVolume();
    eqSurface = meshmetric.getSurface();
    eqTileSpan = 0.0;

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
//    C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
//             * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
//             / (64.0*persistenceLengthCoarse)
//             * (4*x0*x0 - 9*x0 + 6)
//             / (1-x0)*(1-x0);
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


template<typename T, template<typename U> class Descriptor>
void CellModel3D<T, Descriptor>::computeCellForce (Cell3D<T,Descriptor> & cell, T ratio) {
    T cellVolume = cell.getVolume();
    T cellSurface = cell.getSurface();
    plint iVertex=1;
     /* Some force calculations are according to KrugerThesis, Appendix C */
    if (not ((cellVolume > 0) and (cellSurface > 0))) {
        cout << "iVertex: " << iVertex
//             << ", processor: " << dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (tagToParticle3D[iVertex])->get_processor()
//             << ", cellId: " << dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (tagToParticle3D[iVertex])->get_cellId()
             << ", volume: " << cellVolume
             << ", surface: " << cellSurface
             << endl;
        PLB_PRECONDITION( (cellVolume > 0) and (cellSurface > 0) );
    }
    /* Initializations from input */
//    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
    Array<T,3> x1 = cell.getVertex(iVertex), x2, x3, x4;
//    ImmersedCellParticle3D<T,Descriptor>* iParticle = dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (tagToParticle3D[iVertex]);
    Array<T,3> iVelocity;// = iParticle->get_v();
    /* ===================== In case of quasi-rigid object =====================
     *
     * If this is a boundary element (k_rest != 0), get the reference locations
     * of iVertex and calculate and return the force for quasi-rigid objects.
     *          (FengMichaelides2004, J.Comp.Phys. 195(2))
     * // CURRENTLY UNAVAILABLE
     * */
    Array<T,3> restForce; restForce.resetToZero();
//    if (k_rest != 0.0) {
//        Array<T,3> x1ref;
//        boundary.pushSelect(0,1);
//        x1ref = boundary.getMesh().getVertex(iVertex);
//        boundary.popSelect();
//        Array<T,3> dx = x1-x1ref;
//        restForce = (-k_rest*dx) + (-gamma_T*iVelocity); // Dissipative term from Dupin2007
//    }
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

    plint iTriangle;
    plint jVertex, kVertex, lVertex;
    plint vertexIds[3];

    std::vector<plint> neighborTriangles = cell.getNeighborTriangleIds(iVertex);
    for (pluint iB = 0; iB < neighborTriangles.size(); ++iB) {
        iTriangle = neighborTriangles[iB];
        for (pluint id = 0; id < 3; ++id) {
            vertexIds[id] = cell.getVertexId(iTriangle,id);
        }
        plint localVertexId = (iVertex == vertexIds[1])*1 + (iVertex == vertexIds[2])*2;
        jVertex = vertexIds[(localVertexId + 1)%3];
        kVertex = vertexIds[(localVertexId + 2)%3];
        Array<T,3> x2 = cell.getVertex(jVertex);
        Array<T,3> x3 = cell.getVertex(kVertex);
        triangleNormal = cell.computeTriangleNormal(iTriangle);
        triangleArea = cell.computeTriangleArea(iTriangle);

        /* Surface conservation force */
        surfaceForce += computeSurfaceConservationForce(x1, x2, x3, triangleNormal, surfaceCoefficient, dAdx);
        /* Shear force */
        shearForce += computeLocalAreaConservationForce(dAdx, triangleArea, eqArea, areaCoefficient);
//        iParticle->get_E_area() += 0.5*areaCoefficient*pow(eqArea - triangleArea, 2)/3.0; // 3 vertices on each triangle

        /* Elastice Force, disabled */
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
    std::vector<plint> neighborVertexIds = cell.getNeighborVertexIds(iVertex);
    for (pluint jV = 0; jV < neighborVertexIds.size(); jV++) {
        jVertex = neighborVertexIds[jV];
        Array<T,3> x3 = cell.getVertex(jVertex);
        /* In Plane (WLC) and repulsive forces*/
        T potential;
        inPlaneForce += computeInPlaneExplicitForce(x1, x3, eqLengthRatio, eqLength, k_inPlane, potential);
//        iParticle->get_E_inPlane() += potential / 2.0; // Each spring consists of 2 particles.

        /*  Dissipative Forces Calculations */
//        dissipativeForce += computeDissipativeForce(x1, x3, iVelocity, particleVelocity[jVertex], gamma_T, gamma_C);
        /*  Bending Forces Calculations */
        T edgeAngle = cell.computeSignedAngle(iVertex, jVertex, kVertex, lVertex); //edge is iVertex, jVertex
        Array<T,3> iNormal = cell.computeTriangleNormal(iVertex, jVertex, kVertex);
        Array<T,3> jNormal = cell.computeTriangleNormal(iVertex, jVertex, lVertex);

        /*== Compute bending force for the vertex as part of the main edge ==*/
        bendingForce += computeBendingForceEdge (edgeAngle, eqAngle, k_bend, iNormal, jNormal);
//        iParticle->get_E_bending() += computeBendingPotential (edgeAngle, eqAngle, k_bend);
        // Compute bending force for the vertex as a lateral vertex
        edgeAngle = cell.calculateSignedAngle(kVertex, jVertex);
        bendingForce += computeBendingForceLateral (edgeAngle, eqAngle, k_bend, iNormal)*0.5; // Is calculated twice, for the other triangle as well
//        iParticle->get_E_bending() += computeBendingPotential (edgeAngle, eqAngle, k_bend);
        // Compute bending force for the vertex as a lateral vertex
        edgeAngle = cell.calculateSignedAngle(jVertex, lVertex);
        bendingForce += computeBendingForceLateral (edgeAngle, eqAngle, k_bend, jNormal)*0.5; // Is calculated twice, for the other triangle as well
//        iParticle->get_E_bending() += computeBendingPotential (edgeAngle, eqAngle, k_bend);
    }
//    iParticle->get_f_wlc() = inPlaneForce + repulsiveForce;
//    iParticle->get_f_bending() = bendingForce;
//    iParticle->get_f_volume() = volumeForce;
//    iParticle->get_f_surface() = surfaceForce;
//    iParticle->get_f_shear() = shearForce;
//    iParticle->get_f_viscosity() = dissipativeForce;
    return (inPlaneForce + elasticForce + repulsiveForce) + bendingForce + (volumeForce + surfaceForce + shearForce) + dissipativeForce + restForce;// + stretchForce;
}



template<typename T, template<typename U> class Descriptor>
CellModel3D<T, Descriptor>* CellModel3D<T, Descriptor>::clone() const {
    return new CellModel3D<T, Descriptor>(*this);
}




}  // namespace plb

#endif  // CELL_MODEL_3D_HH
