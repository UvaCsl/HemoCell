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

#ifndef SHAPE_MEMORY_MODEL_3D_H
#define SHAPE_MEMORY_MODEL_3D_H

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


template<typename T>
void getCellShapeQuantitiesFromMesh(TriangularSurfaceMesh<T> const& dynMesh,
                            vector<T> & eqAreaPerTriangle, map<plint,T> & eqLengthPerEdge, map<plint,T> & eqAnglePerEdge,
                            plint cellNumTriangles, plint cellNumPartsPerCell);



template<typename T, template<typename U> class Descriptor>
class ShapeMemoryModel3D :  public ConstitutiveModel<T,Descriptor>
{
public:
    /* All input should be in dimensionless units */
    ShapeMemoryModel3D(T density_, T k_rest_,
            T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
            T k_volume_, T k_surface_, T eta_m_,
            T persistenceLengthFine, T eqLengthRatio_,
            T dx_, T dt_, T dm_,
            TriangularSurfaceMesh<T> const& meshElement);
    ShapeMemoryModel3D(ShapeMemoryModel3D<T,Descriptor> const& rhs);
    ~ShapeMemoryModel3D() { } ;
    Array<T,3> computeElasticForce (
               TriangleBoundary3D<T> const& boundary,
                plint iVertex ) { return Array<T,3>(0.0, 0.0, 0.0); }
    virtual void computeCellForce (Cell3D<T,Descriptor> * cell);
    virtual plint getMaximumEdgeExtensionLengthLU() { return maximum(ceil(2*maxLength + 0.5),4); };
    virtual plint getMaxCellDiameterLU() { return maximum(ceil(4*cellRadiusLU),4); };
    virtual T getDx() { return dx; };
    virtual T getDt() { return dt; };
    virtual T getDm() { return dm; };

    virtual SyncRequirements & getSyncRequirements() {return syncRequirements;} ;
    virtual SyncRequirements const& getSyncRequirements() const {return syncRequirements;} ;

    virtual ShapeMemoryModel3D<T,Descriptor>* clone() const;
private:
    inline void computeCellForceHighOrder (Cell3D<T,Descriptor> *cell);
    inline void computeCellForceSuresh (Cell3D<T,Descriptor> * cell);

    plint getTriangleId(plint iTriangle);
    plint getEdgeId(plint iVertex, plint jVertex);
private:
    MeshMetrics<T> meshmetric;
    SyncRequirements syncRequirements;
    T maxLength, cellRadiusLU;

    T k_rest, k_shear, k_bend, k_stretch, k_inPlane, k_elastic, k_surface, k_volume;
    T C_elastic;
    T eta_m, gamma_T, gamma_C;
    T eqLength, eqArea, eqAngle;
    vector<T> eqAreaPerTriangle;
    map<plint,T> eqLengthPerEdge, eqAnglePerEdge;
    T eqVolume, eqSurface, eqTileSpan;
    T persistenceLengthCoarse, eqLengthRatio;
    T dx, dt, dm;
    T persistenceLengthFine;
    pluint cellNumTriangles, cellNumVertices;
public:
    /* Computes the equilibrium quantities to correspond to the an inflated cell with
     * 		eqVolume=ratio*eqVolume.
     * Can also be used for deflation. */
    virtual void inflate(T ratio) {
        eqVolume = meshmetric.getVolume() * ratio*ratio*ratio;
        eqSurface = meshmetric.getSurface() * ratio*ratio;

        getCellShapeQuantitiesFromMesh(meshmetric.getMesh(), eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, cellNumTriangles, cellNumVertices);
        persistenceLengthCoarse = persistenceLengthFine/dx * sqrt( (cellNumVertices-2.0) / (23867-2.0)) * ratio;

    	for (pluint i=0; i < eqAreaPerTriangle.size() ; i++) {
    		eqAreaPerTriangle[i] *= ratio*ratio;
    	}
    	typename std::map<plint,T>::iterator iter;
        for (iter = eqLengthPerEdge.begin(); iter != eqLengthPerEdge.end(); ++iter) {
        	eqLengthPerEdge[iter->first] *= ratio;
        }
    	eqLength *= ratio;
    }
public:
    /* Coefficients */
    /* Coefficients */
    virtual T& getRestingStiffness() { return k_rest; }
    virtual T& getBendingStiffness() { return k_bend; }
    virtual T& getStretchingStiffness() { return k_stretch; }
    virtual void setRestingStiffness(T value) { k_rest = value; }
    virtual void setBendingStiffness(T value) { k_bend = value; }
    virtual void setStretchingStiffness(T value) { k_stretch = value; }

    /* TODO: Fix change of coefficients */

    // Units are N/m
    virtual T getMembraneShearModulus() {
        T Lmax = eqLength*eqLengthRatio;
        T x0 = 1.0/eqLengthRatio;
        // k_inPlane = k_WLC_ * kBT /(4.0*persistenceLengthCoarse);
        T kP =  k_inPlane;
        T k_rep = -1.0*((eqLength*eqLength)* kP * x0 * (-6 + (9 - 4*x0)*x0)/( (x0-1)*(x0-1) ));
        T oneMinusX0 = 1-x0;
        T oneMinusX0Square = oneMinusX0*oneMinusX0;
        T oneMinusX0Cube = oneMinusX0*oneMinusX0*oneMinusX0;
        // Fedosov et al.  doi:10.1016/j.bpj.2010.02.002
        // A Multiscale Red Blood Cell Model with Accurate Mechanics, Rheology, and Dynamics  
        T ansFedosov = kP * sqrt(3) / (Lmax*x0) * (x0/(2.0*oneMinusX0Cube) - 1.0/(4.0*oneMinusX0Square) + 1.0/4.0) + 
                3 * sqrt(3) * k_rep/(4*eqLength*eqLength*eqLength);
        // Reasor et al. doi:10.1002/fld.2534
        // Coupling the lattice-Boltzmann and spectrin-link methods for the direct numerical simulation of cellular blood flow
        T ansReasor = sqrt(3) * kP / (Lmax*x0)*( 3.0/( 4*(1-x0)*(1-x0) ) - 3.0/4.0 +4*x0 + x0/(2.0 * (1-x0)*(1-x0)*(1-x0)));
        return ansFedosov;
    }
    // Units are N/m
    virtual T getMembraneElasticAreaCompressionModulus() {
        T Lmax = eqLength*eqLengthRatio;
        T x0 = 1.0/eqLengthRatio;
        // T kP =  kBT /(4.0*persistenceLengthCoarse);
        T kP =  k_inPlane;
        return sqrt(3) * kP / (Lmax*(1-x0)*(1-x0))*(1.5*(6-9*x0+4*x0*x0) + (1+2*(1-x0)*(1-x0)*(1-x0))/(1-x0)  );
    }
    // Units are N/m
    virtual T getYoungsModulus() {
        T mu0 = getMembraneShearModulus();
        T K = getMembraneElasticAreaCompressionModulus();
        return (4*K*mu0)/(K+mu0);
    }
    // Dimensionless number
    virtual T getPoissonRatio() {
        T mu0 = getMembraneShearModulus();
        T K = getMembraneElasticAreaCompressionModulus();
        return (K-mu0)/(K+mu0);
    }
    // Units are N s/m
    virtual T& getMembraneShearViscosity() { return eta_m; }
    virtual void setMembraneShearViscosity(T value) {
        eta_m = value;
        gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
        gamma_C = (gamma_T/3.0);
    }

    virtual T& getDissipativeParameterT() { return gamma_T; }
    virtual void setDissipativeParameterT(T value) {
        gamma_T = value;
    }
    virtual T& getDissipativeParameterC() { return gamma_C; }
    virtual void setDissipativeParameterC(T value) {
        gamma_C = value;
    }
    /* Equilibrium parameters */
    virtual T& getEquilibriumLinkLength() { return eqLength; }
    virtual T& getEquilibriumTriangleArea() { return eqArea; }
    virtual T& getEquilibriumAngle() { return eqAngle; }
    virtual T& getEquilibriumVolume() { return eqVolume; }
    virtual T& getEquilibriumSurface() { return eqSurface; }
    virtual T& getEquilibriumTileSpan() { return eqTileSpan; }
    virtual void setEquilibriumLinkLength(T value) { pcout << "setEquilibriumLinkLength not available for ShapeMemoryModel3D."; }
    virtual void setEquilibriumTriangleArea(T value) { pcout << "setEquilibriumTriangleArea not available for ShapeMemoryModel3D."; }
    virtual void setEquilibriumAngle(T value) { pcout << "setEquilibriumAngle not available for ShapeMemoryModel3D."; }
    virtual void setEquilibriumVolume(T value) { eqVolume = value; }
    virtual void setEquilibriumSurface(T value) { eqSurface = value; }
    virtual void setEquilibriumTileSpan(T value) { eqTileSpan = value; }

    /* State parameters */
    virtual pluint& getNumberOfVertices() { return cellNumVertices; }
};

}

#include "shapeMemoryModel3D.hh"
#endif  // SHAPE_MEMORY_MODEL_3D_H
