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

#ifndef CELL_MODEL_3D_H
#define CELL_MODEL_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"

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

template<typename T>
class CellModel3D : public ShellModel3D<T>
{
public:
    /* All input should be in dimensionless units */
    CellModel3D(T density_, T k_rest_,
            T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
            T k_volume_, T k_surface_, T eta_m_,
            T eqArea_, T eqLength_, T eqAngle_,
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
    virtual CellModel3D<T>* clone() const;
private:
    T k_rest, k_shear, k_bend, k_stretch, k_inPlane, k_elastic, k_surface, k_volume;
    T C_elastic;
    T eta_m, gamma_T, gamma_C;
    T eqLength, eqArea, eqAngle;
    T eqVolume, eqSurface, eqTileSpan;
    T persistenceLengthCoarse, eqLengthRatio;
    pluint cellNumTriangles, cellNumVertices;
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
    /*
    T getWormLikeChainCoefficient() { return k_WLC/(kBT * maxLength/(4.0*persistenceLength)); }
    void setWormLikeChainCoefficient(T value) { k_WLC = value*(kBT * maxLength/(4.0*persistenceLength)); }
    T getElasticStiffness() { return k_elastic; }
    void setElasticStiffness(T value) {
        k_elastic = value;
        T x0 = eqLength*1.0/maxLength;
        C_elastic = k_elastic * 3.0 * sqrt(3.0)* kBT
            * (maxLength*maxLength*maxLength) * (x0*x0*x0*x0)
            / (64.0*persistenceLength)
            * (4*x0*x0 - 9*x0 + 6)
            / (1-x0)*(1-x0);
    }
    T const& getElasticStiffness() const { return k_elastic; }
    T getVolumeStiffness() { return k_volume/(kBT/pow(eqLength,3)); }
    void setVolumeStiffness(T value) { k_volume = value*(kBT/pow(eqLength,3)); }
    T getShearingStiffness() { return k_shear/(kBT/pow(eqLength,2)); }
    void setShearingStiffness(T value) { k_shear = value*(kBT/pow(eqLength,2)); }
    T getSurfaceStiffness() { return k_surface/(kBT/pow(eqLength,2)); }
    void setSurfaceStiffness(T value) { k_surface = value*kBT/pow(eqLength,2); }
    T& getPersistenceLength() { return persistenceLength; }
    void setPersistenceLength(T value) { persistenceLength = value; }
    T const& getPersistenceLength() const { return persistenceLength; }

                                                                             */
    T const& getMembraneShearViscosity() const { return eta_m; }
    T& getMembraneShearViscosity() { return eta_m; }
    void setMembraneShearViscosity(T value) {
        eta_m = value;
        gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
        gamma_C = (gamma_T/3.0);
    }
    T& getDissipativeParameterT() { return gamma_T; }
    void setDissipativeParameterT(T value) const { gamma_T = value; }
    T const& getDissipativeParameterT() const { return gamma_C; }
    T& getDissipativeParameterC() { return gamma_C; }
    void setDissipativeParameterC(T value) const { gamma_C = value; }
    T const& getDissipativeParameterC() const { return gamma_C; }
    /* Equilibrium parameters */
    T& getEquilibriumLinkLength() { return eqLength; }
    void setEquilibriumLinkLength(T value) { eqLength = value; }
    T const& getEquilibriumLinkLength() const { return eqLength; }
    T& getEquilibriumTriangleArea() { return eqArea; }
    void setEquilibriumTriangleArea(T value) { eqArea = value; }
    T const& getEquilibriumTriangleArea() const { return eqArea; }
    T& getEquilibriumAngle() { return eqAngle; }
    void setEquilibriumAngle(T value) { eqAngle = value; }
    T const& getEquilibriumAngle() const { return eqAngle; }
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
//    T& getMaximumLinkLength() { return maxLength; }
//    void setMaximumLinkLength(T value) { maxLength = value; }
//    T const& getMaximumLinkLength() const { return maxLength; }
//    pluint& getNumberOfVertices() { return Nv; }
//    void setNumberOfVertices(pluint value) { Nv = value; }
//    pluint const& getNumberOfVertices() const { return Nv; }
};


namespace shellModelHelper3D {

namespace cellModelHelper3D {

template<typename T>
T computePotential(plint iVertex, Array<T,3> const& iPosition,
                   TriangularSurfaceMesh<T> const& mesh,
                   T eqLength, T maxLength, T eqArea, T eqAngle, T eqTileSpan,
                   T k_WLC, T k_elastic, T k_shear, T k_bend );

template<typename T>
T computeInPlanePotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                             T maxLength, T k);

template<typename T>
T computeShearPotential(Array<T,3> const& iPosition,
                        Array<T,3> const& jPosition,
                        Array<T,3> const& kPosition,
                        T eqArea,
                        T k);

template<typename T>
T computeBendPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                       Array<T,3> const& kPosition, Array<T,3> const& lPosition,
                       T eqTileSpan, T eqLength, T eqAngle, T k);

}  // namespace cellModelHelper3D


}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_H
