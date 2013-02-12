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

#include "core/globalDefs.h"
#include "offLattice/triangleBoundary3D.h"
#include "shellModel3D.h"

namespace plb {

template<typename T>
class CellModel3D : public ShellModel3D<T>
{
public:
    /* All input should be in dimensionless units */
    CellModel3D(T density_,
                T k_shear_, T k_bend_, T k_stretch_, T k_WLC_, T k_elastic_,
                T k_volume_, T k_surface_, T eta_m_,
                T eqArea_, T eqLength_, T eqAngle_,
                T eqVolume_, T eqSurface_, T eqTileSpan_,
                T maxLength_, T persistenceLength_);
    virtual Array<T,3> computeCellForce (
            TriangleBoundary3D<T> const& boundary,
            T cellVolume, T cellSurface,
            std::map< plint, Array<T,3> > particleVelocity,
            plint iVertex );
    virtual Array<T,3> computeElasticForce (
            TriangleBoundary3D<T> const& boundary,
            plint iVertex );
    virtual CellModel3D<T>* clone() const;
    T const& getShearingStiffness() const { return k_shear; }
    T& getShearingStiffness() { return k_shear; }
    T const& getBendingStiffness() const { return k_bend; }
    T& getBendingStiffness() { return k_bend; }
    T const& getSurfaceStiffness() const { return k_surface; }
    T& getSurfaceStiffness() { return k_surface; }
    T const& getVolumeStiffness() const { return k_volume; }
    T& getVolumeStiffness() { return k_volume; }
    T const& getWLCStiffness() const { return k_WLC; }
    T& getWLCStiffness() { return k_WLC; }
    T const& getEquilibriumLength() const { return eqLength; }
    T& getEquilibriumLength() { return eqLength; }
    T const& getEquilibriumAngle() const { return eqAngle; }
    T& getEquilibriumAngle() { return eqAngle; }
    T const& getEquilibriumArea() const { return eqArea; }
    T& getEquilibriumArea() { return eqArea; }
    T const& getEquilibriumVolume() const { return eqVolume; }
    T& getEquilibriumVolume() { return eqVolume; }
    T const& getEquilibriumSurface() const { return eqSurface; }
    T& getEquilibriumSurface() { return eqSurface; }
    T const& getMaximumLength() const { return maxLength; }
    T& getMaximumLength() { return maxLength; }
    T const& getPersistenceLength() const { return persistenceLength; }
    T& getPersistenceLength() { return persistenceLength; }
    void setEquilibriumVolume(T eqVolume_) { eqVolume = eqVolume_; }
    void setEquilibriumSurface(T eqSurface_) { eqSurface = eqSurface_; }
private:
    T k_shear, k_bend, k_stretch, k_WLC, k_elastic, k_surface, k_volume;
    T eta_m, gamma_T, gamma_C;
    T eqLength, eqArea, eqAngle;
    T eqVolume, eqSurface, eqTileSpan;
    T maxLength, persistenceLength, C_elastic;
};


template<typename T>
Array<T,3> computeBendingForce (Array<T,3> const& x1, Array<T,3> const& x2,
                                Array<T,3> const& x3, Array<T,3> const& x4,
                                Array<T,3> const& ni, Array<T,3> const& nj,
                                T Ai, T Aj,
                                T eqTileSpan, T eqLength, T eqAngle, T k);


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
