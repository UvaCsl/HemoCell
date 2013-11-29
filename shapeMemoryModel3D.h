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

#include "palabos3D.h"
#include "palabos3D.hh"
#include "shellModel3D.h"
#include <cmath>
#include <map>
#include "computeCellForces3D.hh"


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
class ShapeMemoryModel3D : public ShellModel3D<T>
{
public:
    /* All input should be in dimensionless units */
    ShapeMemoryModel3D(T density_, T k_rest_,
                T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
                T k_volume_, T k_surface_, T eta_m_,
                vector<T> eqArea_, map<plint,T> eqLength_, map<plint,T> eqAngle_,
                T eqVolume_, T eqSurface_, T eqTileSpan_,
                T persistenceLengthFine, T eqLengthRatio_, pluint cellNumTriangles_, pluint cellNumVertices_);
    virtual Array<T,3> computeCellForce (
            TriangleBoundary3D<T> const& boundary,
            T cellVolume, T cellSurface, T & iSurface,
            std::map< plint, Array<T,3> > & particleVelocity,
            std::map< plint, Array<T,3>* > & particleForces,
            plint iVertex,
            Array<T,3> & f_wlc, Array<T,3> & f_bending, Array<T,3> & f_volume,
            Array<T,3> & f_surface, Array<T,3> & f_shear, Array<T,3> & f_viscosity);
    virtual Array<T,3> computeElasticForce (
            TriangleBoundary3D<T> const& boundary,
            plint iVertex );
    virtual ShapeMemoryModel3D<T>* clone() const;
private:
    plint getTriangleId(plint iTriangle);
    plint getEdgeId(plint iVertex, plint jVertex);
private:
    T k_rest, k_shear, k_bend, k_stretch, k_inPlane, k_elastic, k_surface, k_volume;
    T C_elastic;
    T eta_m, gamma_T, gamma_C;
    T eqLength, eqAngle;
    vector<T> eqAreaPerTriangle;
    map<plint,T> eqLengthPerEdge, eqAnglePerEdge;
    T eqVolume, eqSurface, eqTileSpan;
    T persistenceLengthCoarse, eqLengthRatio;
    pluint cellNumTriangles, cellNumVertices;
public:
    /* Computes the equilibrium quantities to correspond to the an inflated cell with
     * 		eqVolume=ratio*eqVolume.
     * Can also be used for deflation. */
    void inflate(T ratio, bool scaleCoefficients=true) {
    	eqVolume *= ratio;
    	eqSurface *= pow(ratio,2.0/3.0);
    	eqTileSpan *= pow(ratio,1.0/3.0);
    	for (plint i=0; i < eqAreaPerTriangle.size() ; i++) {
    		eqAreaPerTriangle[i] *= pow(ratio,2.0/3.0);
    	}
    	typename std::map<plint,T>::iterator iter;
        for (iter = eqLengthPerEdge.begin(); iter != eqLengthPerEdge.end(); ++iter) {
        	eqLengthPerEdge[iter->first] *= pow(ratio,1.0/3.0);
        }
    	eqLength *= pow(ratio,1.0/3.0);
    	if (scaleCoefficients) {
    	    k_volume *= 1.0/ratio;
    	    k_surface *= 1.0/pow(ratio,2.0/3.0);
    	    k_shear *= 1.0/pow(ratio,2.0/3.0);
    	}
    }
public:
    /* Coefficients */
    T& getRestingStiffness() { return k_rest; }
    void setRestingStiffness(T value) { k_rest = value; }
    T const& getRestingStiffness() const { return k_rest; }
    T& getBendingStiffness() { return k_bend; }
    void setBendingStiffness(T value) { k_bend = value; }
    T const& getBendingStiffness() const { return k_bend; }
    T& getStretchingStiffness() { return k_stretch; }
    void setStretchingStiffness(T value) { k_stretch = value; }
    T const& getStretchingStiffness() const { return k_stretch; }

    /* TODO: Fix change of coefficients */

    // Units are N/m
    T getMembraneShearModulus() {
        T Lmax = eqLength*eqLengthRatio;
        T x0 = 1.0/eqLengthRatio;
        // T kP =  kBT /(4.0*persistenceLengthCoarse);
        T kP =  k_inPlane;
        return sqrt(3) * kP / (Lmax*x0)*( 3.0/( 4*(1-x0)*(1-x0) ) - 3.0/4.0 +4*x0 + x0/(2.0 * (1-x0)*(1-x0)*(1-x0)));
    }
    // Units are N/m
    T getMembraneElasticAreaCompressionModulus() {
        T Lmax = eqLength*eqLengthRatio;
        T x0 = 1.0/eqLengthRatio;
        // T kP =  kBT /(4.0*persistenceLengthCoarse);
        T kP =  k_inPlane;
        return sqrt(3) * kP / (Lmax*(1-x0)*(1-x0))*(1.5*(6-9*x0+4*x0*x0) + (1+2*(1-x0)*(1-x0)*(1-x0))/(1-x0)  );
    }
    // Units are N/m
    T getYoungsModulus() {
        T mu0 = getMembraneShearModulus();
        T K = getMembraneElasticAreaCompressionModulus();
        return (4*K*mu0)/(K+mu0);
    }
    // Dimensionless number
    T getPoissonRatio() {
        T mu0 = getMembraneShearModulus();
        T K = getMembraneElasticAreaCompressionModulus();
        return (K-mu0)/(K+mu0);
    }
    // Units are N s/m
    T const& getMembraneShearViscosity() const { return eta_m; }
    T& getMembraneShearViscosity() { return eta_m; }
    void setMembraneShearViscosity(T value) {
        eta_m = value;
        gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
        gamma_C = (gamma_T/3.0);
    }

    T& getDissipativeParameterT() { return gamma_T; }
    T const& getDissipativeParameterT() const { return gamma_C; }
    void setDissipativeParameterT(T value) const {
        gamma_T = value;
    }
    T& getDissipativeParameterC() { return gamma_C; }
    T const& getDissipativeParameterC() const { return gamma_C; }
    void setDissipativeParameterC(T value) const {
        gamma_C = value;
    }
    /* Equilibrium parameters */
//    vector<T>& getEquilibriumLinkLength() { return eqLength; }
//    void setEquilibriumLinkLength(vector<T> value) { eqLength = value; }
//    vector<T> const& getEquilibriumLinkLength() const { return eqLength; }
//    vector<T>& getEquilibriumTriangleArea() { return eqArea; }
//    void setEquilibriumTriangleArea(vector<T> value) { eqArea = value; }
//    vector<T> const& getEquilibriumTriangleArea() const { return eqArea; }
//    vector<T>& getEquilibriumAngle() { return eqAngle; }
//    void setEquilibriumAngle(vector<T> value) { eqAngle = value; }
//    vector<T> const& getEquilibriumAngle() const { return eqAngle; }
    T& getEquilibriumVolume() { return eqVolume; }
    void setEquilibriumVolume(T value) { eqVolume = value; }
    T const& getEquilibriumVolume() const { return eqVolume; }
    T& getEquilibriumSurface() { return eqSurface; }
    void setEquilibriumSurface(T value) { eqSurface = value; }
    T const& getEquilibriumSurface() const { return eqSurface; }
    T& getEquilibriumTileSpan() { return eqTileSpan; }
    void setEquilibriumTileSpan(T value) { eqTileSpan = value; }
    T const& getEquilibriumTileSpan() const { return eqTileSpan; }
    /* State parameters */
    pluint& getNumberOfVertices() { return cellNumVertices; }
    void setNumberOfVertices(pluint value) { cellNumVertices = value; }
    pluint const& getNumberOfVertices() const { return cellNumVertices; }
};

}  // namespace plb


#include "shapeMemoryModel3D.hh"
#endif  // SHAPE_MEMORY_MODEL_3D_H
