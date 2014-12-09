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

#ifndef REST_MODEL_3D_H
#define REST_MODEL_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <map>
#include "cellModel3D.h"
#include "cellField3D.h"

using namespace std;


template<typename T, template<typename U> class Descriptor>
class GetCellVertexPositions : public BoxProcessingFunctional3D
{
public:
    GetCellVertexPositions (std::map<plint, std::map<plint,Array<T,3> > > & vertexToPosition_) ;
    ~GetCellVertexPositions() { };
    GetCellVertexPositions(GetCellVertexPositions<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual GetCellVertexPositions<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    std::map<plint, std::map<plint,Array<T,3> > > & vertexToPosition; // cellVertexToPosition[cellId][iVertex] = getPosition();
};


template<typename T, template<typename U> class Descriptor>
class RestModel3D : public ConstitutiveModel<T,Descriptor>
{
public:
    /* All input should be given in dimensional units.
     * Afterwards they are saved in dimensionless LB units.
     * 
     *     meshElement should be ready and in LB units
     * 
     *  */
    RestModel3D(T density_, T k_rest_,
            T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
            T k_volume_, T k_surface_, T eta_m_,
            T persistenceLengthFine, T eqLengthRatio_,
            T dx_, T dt_, T dm_,
            TriangularSurfaceMesh<T> const& meshElement);
    ~RestModel3D() { } ;
    RestModel3D(RestModel3D<T,Descriptor> const& rhs);
    virtual void computeCellForce (Cell3D<T,Descriptor> * cell);
    virtual plint getMaximumEdgeExtensionLengthLU() { return maximum(ceil(2*maxLength + 0.5),4); };
    virtual plint getMaxCellDiameterLU() { return maximum(ceil(2*cellRadiusLU),2); };
    virtual T getDx() { return dx; };
    virtual T getDt() { return dt; };
    virtual T getDm() { return dm; };

    virtual SyncRequirements & getSyncRequirements() {return syncRequirements;} ;
    virtual SyncRequirements const& getSyncRequirements() const {return syncRequirements;} ;

    void freezeVertices(CellField3D<T, Descriptor> & cellField) {
        applyProcessingFunctional ( // compute force applied on the fluid by the particles
                new GetCellVertexPositions<T,Descriptor> (frozenVertices),
                cellField.getBoundingBox(), cellField.getParticleArg() );
    };

    virtual RestModel3D<T, Descriptor>* clone() const;
private:
    MeshMetrics<T> meshmetric;
    T k_rest, k_shear, k_bend, k_stretch, k_inPlane, k_elastic, k_surface, k_volume;
    T C_elastic;
    T eta_m, gamma_T, gamma_C;
    T eqLength, eqArea, eqAngle;
    T eqVolume, eqSurface, eqTileSpan;
    T persistenceLengthCoarse, eqLengthRatio;
    T dx, dt, dm;
    T persistenceLengthFine;
    pluint cellNumTriangles, cellNumVertices;
    plint cellRadiusLU, maxLength;
    SyncRequirements syncRequirements;
    std::map<plint, std::map<plint,Array<T,3> > >  frozenVertices;
public:
    /* Computes the equilibrium quantities to correspond to the an inflated cell with
     * 		eqVolume=ratio*eqVolume.
     * Can also be used for deflation. */
    virtual void inflate(T ratio) {
        eqLength = meshmetric.getMeanLength() * ratio;
        maxLength = meshmetric.getMaxLength()*eqLengthRatio * ratio;
        eqArea = meshmetric.getMeanArea() * ratio * ratio;
        eqVolume = meshmetric.getVolume() * ratio * ratio* ratio;
        eqSurface = meshmetric.getSurface() * ratio * ratio;
        persistenceLengthCoarse = persistenceLengthFine/dx * sqrt( (cellNumVertices-2.0) / (23867-2.0))  * ratio;
    }
public:
    /* Coefficients */
    virtual T& getRestingStiffness() { return k_rest; }
    virtual T& getBendingStiffness() { return k_bend; }
    virtual T& getStretchingStiffness() { return k_stretch; }
    virtual void setRestingStiffness(T value) { k_rest = value; }
    virtual void setBendingStiffness(T value) { k_bend = value; }
    virtual void setStretchingStiffness(T value) { k_stretch = value; }

    /* TODO: Fix change of coefficients */
    // Units are N/m
    // Units are N/m
    virtual T getMembraneShearModulus() {
        T Lmax = eqLength*eqLengthRatio;
        T x0 = 1.0/eqLengthRatio;
        // T kP =  kBT /(4.0*persistenceLengthCoarse);
        T kP =  k_inPlane;
        return sqrt(3) * kP / (Lmax*x0)*( 3.0/( 4*(1-x0)*(1-x0) ) - 3.0/4.0 +4*x0 + x0/(2.0 * (1-x0)*(1-x0)*(1-x0)));
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
    virtual T& getDissipativeParameterC() { return gamma_C; }
    virtual void setDissipativeParameterT(T value) { gamma_T = value; }
    virtual void setDissipativeParameterC(T value) { gamma_C = value; }
    /* Equilibrium parameters */
    virtual T& getEquilibriumLinkLength() { return eqLength; }
    virtual T& getEquilibriumTriangleArea() { return eqArea; }
    virtual T& getEquilibriumAngle() { return eqAngle; }
    virtual T& getEquilibriumVolume() { return eqVolume; }
    virtual T& getEquilibriumSurface() { return eqSurface; }
    virtual T& getEquilibriumTileSpan() { return eqTileSpan; }
    virtual void setEquilibriumLinkLength(T value) { eqLength = value; }
    virtual void setEquilibriumTriangleArea(T value) { eqArea = value; }
    virtual void setEquilibriumAngle(T value) { eqAngle = value; }
    virtual void setEquilibriumVolume(T value) { eqVolume = value; }
    virtual void setEquilibriumSurface(T value) { eqSurface = value; }
    virtual void setEquilibriumTileSpan(T value) { eqTileSpan = value; }
    /* State parameters */
//    T& getMaximumLinkLength() { return maxLength; }
//    void setMaximumLinkLength(T value) { maxLength = value; }
//    T const& getMaximumLinkLength() const { return maxLength; }
//    pluint& getNumberOfVertices() { return Nv; }
//    void setNumberOfVertices(pluint value) { Nv = value; }
//    pluint const& getNumberOfVertices() const { return Nv; }
};




#include "restModel3D.hh"
#endif  // CELL_MODEL_3D_H