#ifndef REST_MODEL_3D_HH
#define REST_MODEL_3D_HH

#include "restModel3D.h"
#include "shapeMemoryModel3D.h"

template<typename T, template<typename U> class Descriptor>
RestModel3D<T, Descriptor>::RestModel3D(T density_, T k_rest_,
            T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
            T k_volume_, T k_surface_, T eta_m_,
            T persistenceLengthFine_, T eqLengthRatio_,
            T dx_, T dt_, T dm_, 
            TriangularSurfaceMesh<T> const& meshElement,
            std::map<plint, TriangularSurfaceMesh<T>* > & meshes_)
    : ShellModel3D<T>(density_),
      meshmetric(meshElement),
      k_rest(k_rest_),
      k_shear(k_shear_),
      k_bend(k_bend_),
      k_stretch(k_stretch_),
      k_elastic(k_elastic_),
      k_surface(k_surface_),
      k_volume(k_volume_),
      eta_m(eta_m_),
      eqLengthRatio(eqLengthRatio_),
      dx(dx_), dt(dt_), dm(dm_),
      persistenceLengthFine(persistenceLengthFine_),
      syncRequirements(),
      meshes(meshes_)
{
    T dNewton = (dm*dx/(dt*dt)) ;
    T kBT = kBT_p / ( dm * dx*dx/(dt*dt) );

    k_WLC_ *= 1.0;     k_elastic *= 1.0;     k_bend *= 1.0;
    k_volume *= 1.0;     k_surface *= 1.0;     k_shear *= 1.0;
    eta_m /= dNewton*dt/dx;     k_stretch /= dNewton;    k_rest /= dNewton/dx;


    syncRequirements.insert(volumeAndSurfaceReductions);

    cellNumVertices = meshmetric.getNumVertices();
    cellNumTriangles = meshmetric.getNumTriangles();
    eqLength = meshmetric.getMeanLength();
    maxLength = meshmetric.getMaxLength()*eqLengthRatio;
    T eqMeanArea = eqArea = meshmetric.getMeanArea();
//    eqAngle = meshmetric.getMeanAngle();
    eqVolume = meshmetric.getVolume();
    eqSurface = meshmetric.getSurface();
    eqTileSpan = 0.0;

    getCellShapeQuantitiesFromMesh(meshElement, eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, cellNumTriangles, cellNumVertices);

    persistenceLengthCoarse = persistenceLengthFine/dx * sqrt( (cellNumVertices-2.0) / (23867-2.0)) ;

    eqAngle=0.0;

    typename map<plint,T>::reverse_iterator iter = eqAnglePerEdge.rbegin();
    for (iter = eqAnglePerEdge.rbegin(); iter != eqAnglePerEdge.rend(); ++iter) {
       eqAngle += iter->second;
    }
    eqAngle /= eqAnglePerEdge.size();


    /* Calculate cell Radius */
    cellRadiusLU = meshmetric.getRadius();

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

    T forceSum = 0; // norm(computeInPlaneExplicitForce(x1, x3, eqLengthRatio, eqLength, k_inPlane));
//    C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
//             * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
//             / (64.0*persistenceLengthCoarse)
//             * (4*x0*x0 - 9*x0 + 6)
//             / (1-x0)*(1-x0);
#ifdef PLB_DEBUG // Less Calculations
    pcout << std::endl;
    pcout << " ============================================= " << std::endl;
    pcout << "eqLength: " << eqLength << std::endl;
    pcout << "forceSum: " << forceSum << std::endl;
    pcout << "k_bend: " << k_bend << ",\t eqAngle (degrees): " << eqAngle*180.0/pi << std::endl;
    pcout << "k_volume: " << k_volume << ",\t eqVolume: " << eqVolume << std::endl;
    pcout << "k_surface: " << k_surface << ",\t eqSurface: " << eqSurface << std::endl;
    pcout << "k_shear: " << k_shear << ",\t eqArea: " << eqArea << std::endl;
    pcout << "eta_m: " << eta_m << " ";
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
    pcout << "# mu_0 = " << getMembraneShearModulus()*dNewton/dx << std::endl;
    pcout << "# K = " << getMembraneElasticAreaCompressionModulus()*dNewton/dx << std::endl;
    pcout << "# YoungsModulus = " << getYoungsModulus()*dNewton/dx << std::endl;
    pcout << "# Poisson ratio = " << getPoissonRatio() << std::endl;
    pcout << " ============================================= " << std::endl;
#endif // PLB_DEBUG // Less Calculations

}


template<typename T, template<typename U> class Descriptor>
void RestModel3D<T, Descriptor>::computeCellForce (Cell3D<T,Descriptor> * cell) {
     /* Some force calculations are according to KrugerThesis, Appendix C */
//	pcout << "RestModel3D<T, Descriptor>::computeCellForce" << std::endl;
    T cellVolume = cell->getVolume();
    T cellSurface = cell->getSurface();

    if (not ((cellVolume > 0) and (cellSurface > 0))) {
        cout << "processor: " << cell->getMpiProcessor()
             << ", cellId: " << cell->get_cellId()
             << ", volume: " << cellVolume
             << ", surface: " << cellSurface
             << ", cellNumVertices: " << cellNumVertices
             << endl;
        PLB_PRECONDITION( (cellVolume > 0) and (cellSurface > 0) );
    }
    std::vector<plint> const& triangles = cell->getTriangles();
    std::vector<Array<plint,2> > const& edges = cell->getEdges();
    std::vector<plint > const& vertices = cell->getVertices();
    for (pluint iV = 0; iV < vertices.size(); ++iV) {
        castParticleToICP3D(cell->getParticle3D(vertices[iV]))->resetForces();
    }
    plint iTriangle;
    plint iVertex, jVertex, kVertex=-1, lVertex=-1;
//    T maxDr = 0;
    for (pluint iV = 0; iV < vertices.size(); ++iV) {
    	iVertex = vertices[iV];
        SurfaceParticle3D<T,Descriptor>* iParticle = castParticleToICP3D(cell->getParticle3D(iVertex));
    	plint cellId = iParticle->get_cellId();
		Array<T,3> dr = cell->getVertex(iVertex) - meshes[cellId]->getVertex(iVertex);
		iParticle->get_force() += (-k_rest*dr); // Dissipative term from Dupin2007
//		maxDr = norm(dr)>maxDr?norm(dr):maxDr;
    }
//    pcout << "maxDr " << maxDr << std::endl;
//    return;

    /* Run through all the edges and calculate:
         x In plane (WLC and repulsive) force
         x Dissipative force
         x Bending force
         o Stretch force
     */
    Array<T,3> force1, force2, force3;
    T potential;
    for (pluint iE = 0; iE < edges.size(); ++iE) {
        iVertex = edges[iE][0];  jVertex = edges[iE][1];
        plint edgeId = getEdgeId(iVertex, jVertex);
        Array<T,3> const& iX = cell->getVertex(iVertex);
        Array<T,3> const& jX = cell->getVertex(jVertex);
        SurfaceParticle3D<T,Descriptor>* iParticle = castParticleToICP3D(cell->getParticle3D(iVertex));
        SurfaceParticle3D<T,Descriptor>* jParticle = castParticleToICP3D(cell->getParticle3D(jVertex));
          /* ------------------------------------*/
         /* In Plane forces (WLC and repulsive) */
        /* ------------------------------------*/
        force1 = computeInPlaneExplicitForce(iX, jX, eqLengthRatio, eqLengthPerEdge[edgeId], k_inPlane, potential);
        iParticle->get_force() += force1;
        jParticle->get_force() -= force1;
#ifdef PLB_DEBUG // Less Calculations
        iParticle->get_E_inPlane() += potential;
        jParticle->get_E_inPlane() += potential;
        iParticle->get_f_wlc() += force1;
        jParticle->get_f_wlc() -= force1;
#endif
          /* ------------------------------------*/
         /*    Dissipative Forces Calculations  */
        /* ------------------------------------*/
        if (gamma_T>0.0) {
            force1 = computeDissipativeForce(iX, jX, iParticle->get_v(), jParticle->get_v(), gamma_T, gamma_C);
            iParticle->get_force() += force1;
            jParticle->get_force() -= force1;
#ifdef PLB_DEBUG // Less Calculations
            iParticle->get_f_viscosity() += force1;
            jParticle->get_f_viscosity() -= force1;
#endif
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

        T edgeAngle = cell->computeSignedAngle(iVertex, jVertex, kVertex, lVertex, angleFound); //edge is iVertex, jVertex
        if (angleFound) {
            Array<T,3> iNormal = cell->computeTriangleNormal(iVertex, jVertex, kVertex);
            Array<T,3> jNormal = cell->computeTriangleNormal(iVertex, jVertex, lVertex);
            T Ai = cell->computeTriangleArea(iVertex, jVertex, kVertex);
            T Aj = cell->computeTriangleArea(iVertex, jVertex, lVertex);
            SurfaceParticle3D<T,Descriptor>* kParticle = castParticleToICP3D(cell->getParticle3D(kVertex));
            SurfaceParticle3D<T,Descriptor>* lParticle = castParticleToICP3D(cell->getParticle3D(lVertex));
            Array<T,3> const& kX = cell->getVertex(kVertex);
            Array<T,3> const& lX = cell->getVertex(lVertex);

            /*== Compute bending force for the vertex as part of the main edge ==*/
            Array<T,3> fi, fk, fj, fl;
            fi = computeBendingForce (iX, kX, jX, lX, iNormal, jNormal, Ai, Aj, eqArea, eqLengthPerEdge[edgeId], eqAnglePerEdge[edgeId], k_bend, fk, fj, fl);

            iParticle->get_force() += fi;
            jParticle->get_force() += fj;
            kParticle->get_force() += fk;
            lParticle->get_force() += fl;

#ifdef PLB_DEBUG // Less Calculations
            iParticle->get_f_bending() += fi;
            jParticle->get_f_bending() += fj;
            kParticle->get_f_bending() += fk;
            lParticle->get_f_bending() += fl;
            potential = computeBendingPotential (edgeAngle, eqAnglePerEdge[edgeId], k_bend);
            iParticle->get_E_bending() += potential;
            jParticle->get_E_bending() += potential;
            kParticle->get_E_bending() += potential;
            lParticle->get_E_bending() += potential;
#endif
        } else {
            cout << global::mpi().getRank() << " angle not found " << std::endl;
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
    Array<T,3> dAdx1, dAdx2, dAdx3, dVdx, tmp(0,0,0);
    std::map<plint, T> trianglesArea;
    std::map<plint, Array<T,3> > trianglesNormal;
    T triangleArea;
    Array<T,3> triangleNormal;

    for (pluint iT = 0; iT < triangles.size(); ++iT) {
        iTriangle = triangles[iT];
        triangleNormal = cell->computeTriangleNormal(iTriangle);
        triangleArea = cell->computeTriangleArea(iTriangle);
        iVertex = cell->getVertexId(iTriangle,0);
        jVertex = cell->getVertexId(iTriangle,1);
        kVertex = cell->getVertexId(iTriangle,2);
        Array<T,3> const& x1 = cell->getVertex(iVertex);
        Array<T,3> const& x2 = cell->getVertex(jVertex);
        Array<T,3> const& x3 = cell->getVertex(kVertex);
        SurfaceParticle3D<T,Descriptor>* iParticle = castParticleToICP3D(cell->getParticle3D(iVertex));
        SurfaceParticle3D<T,Descriptor>* jParticle = castParticleToICP3D(cell->getParticle3D(jVertex));
        SurfaceParticle3D<T,Descriptor>* kParticle = castParticleToICP3D(cell->getParticle3D(kVertex));

        /* Surface conservation forces */
        force1  = computeSurfaceConservationForce(x1, x2, x3, triangleNormal, surfaceCoefficient, dAdx1);
        force2  = computeSurfaceConservationForce(x2, x3, x1, triangleNormal, surfaceCoefficient, dAdx2);
        force3  = computeSurfaceConservationForce(x3, x1, x2, triangleNormal, surfaceCoefficient, dAdx3);
        iParticle->get_force() += force1;
        jParticle->get_force() += force2;
        kParticle->get_force() += force3;

#ifdef PLB_DEBUG // Less Calculations
        iParticle->get_f_surface() += force1;
        jParticle->get_f_surface() += force2;
        kParticle->get_f_surface() += force3;
        potential = 0.5*surfaceCoefficient*(cellSurface - eqSurface)*1.0/cellNumVertices;
        iParticle->get_E_area() += potential;
        jParticle->get_E_area() += potential;
        kParticle->get_E_area() += potential;
#endif
        /* Local area conservation forces */
        force1 = computeLocalAreaConservationForce(dAdx1, triangleArea, eqAreaPerTriangle[iTriangle], areaCoefficient);
        force2 = computeLocalAreaConservationForce(dAdx2, triangleArea, eqAreaPerTriangle[iTriangle], areaCoefficient);
        force3 = computeLocalAreaConservationForce(dAdx3, triangleArea, eqAreaPerTriangle[iTriangle], areaCoefficient);
        iParticle->get_force() += force1;
        jParticle->get_force() += force2;
        kParticle->get_force() += force3;

#ifdef PLB_DEBUG // Less Calculations
        iParticle->get_f_shear() += force1;
        jParticle->get_f_shear() += force2;
        kParticle->get_f_shear() += force3;
        potential = 0.5*areaCoefficient*pow(eqAreaPerTriangle[iTriangle] - triangleArea, 2)/3.0; // 3 vertices on each triangle
        iParticle->get_E_area() += potential;
        jParticle->get_E_area() += potential;
        kParticle->get_E_area() += potential;
#endif

        /* Volume conservation forces */
        force1  = computeVolumeConservationForce(x1, x2, x3, volumeCoefficient);
        force2  = computeVolumeConservationForce(x2, x3, x1, volumeCoefficient);
        force3  = computeVolumeConservationForce(x3, x1, x2, volumeCoefficient);
        iParticle->get_force() += force1;
        jParticle->get_force() += force2;
        kParticle->get_force() += force3;

#ifdef PLB_DEBUG // Less Calculations
        iParticle->get_f_volume() += force1;
        jParticle->get_f_volume() += force2;
        kParticle->get_f_volume() += force3;
        potential = 0.5*volumeCoefficient*(cellVolume - eqVolume)*1.0/cellNumVertices;
        iParticle->get_E_volume() += potential;
        jParticle->get_E_volume() += potential;
        kParticle->get_E_volume() += potential;
#endif

    }

}


/* ******** GetCellVertexPositions *********************************** */

template<typename T, template<typename U> class Descriptor>
GetCellVertexPositions<T,Descriptor>::GetCellVertexPositions ( std::map<plint, std::map<plint,Array<T,3> > >  & vertexToPosition_)
: vertexToPosition(vertexToPosition_) { }

template<typename T, template<typename U> class Descriptor>
GetCellVertexPositions<T,Descriptor>::GetCellVertexPositions (
            GetCellVertexPositions<T,Descriptor> const& rhs)
    : vertexToPosition(rhs.vertexToPosition)
{ }


template<typename T, template<typename U> class Descriptor>
void GetCellVertexPositions<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{

    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        SurfaceParticle3D<T,Descriptor>* particle =
            dynamic_cast<SurfaceParticle3D<T,Descriptor>*> (particles[iParticle]);
        plint iVertex = particle->getVertexId();
        plint cellId = particle->get_cellId();
        vertexToPosition[cellId][iVertex] = particle->getPosition();
    }

}

template<typename T, template<typename U> class Descriptor>
GetCellVertexPositions<T,Descriptor>*
    GetCellVertexPositions<T,Descriptor>::clone() const
{
    return new GetCellVertexPositions<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT GetCellVertexPositions<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void GetCellVertexPositions<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}

template<typename T, template<typename U> class Descriptor>
void GetCellVertexPositions<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}



#endif  // CELL_MODEL_3D_HH
