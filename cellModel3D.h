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
    CellModel3D(T density_,
                T k_stretch_, T k_shear_, T k_bend_,
                T eqArea_, T eqLength_, T eqAngle_
                    );
    virtual Array<T,3> computeElasticForce (
            TriangleBoundary3D<T> const& boundary,
            plint iVertex );
    virtual CellModel3D<T>* clone() const;
    T const& getStretchingStiffness() const { return k_stretch; }
    T& getStretchingStiffness() { return k_stretch; }
    T const& getShearingStiffness() const { return k_shear; }
    T& getShearingStiffness() { return k_shear; }
    T const& getBendingStiffness() const { return k_bend; }
    T& getBendingStiffness() { return k_bend; }
    T const& getEquilibriumLength() const { return eqLength; }
    T& getEquilibriumLength() { return eqLength; }
    T const& getEquilibriumAngle() const { return eqAngle; }
    T& getEquilibriumAngle() { return eqAngle; }
    T const& getEquilibriumArea() const { return eqArea; }
    T& getEquilibriumArea() { return eqArea; }
private:
    T k_stretch, k_shear, k_bend;
    T eqLength, eqArea, eqAngle;
};




namespace shellModelHelper3D {

namespace cellModelHelper3D {

template<typename T>
T computePotential(plint iVertex, Array<T,3> const& iPosition,
                   TriangularSurfaceMesh<T> const& dynMesh, 
                   T k_stretch, T k_shear, T k_bend,
                   T eqArea, T eqLength, T eqAngle);

template<typename T>
T computeStretchPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                          T eqLength,
                          T k);

template<typename T>
T computeShearPotential(Array<T,3> const& iPosition,
                        Array<T,3> const& jPosition,
                        Array<T,3> const& kPosition,
                        T eqArea,
                        T k);

template<typename T>
T computeBendPotential(Array<T,3> const& iPosition, Array<T,3> const& jPosition,
                       Array<T,3> const& kPosition, Array<T,3> const& lPosition,
                       T eqAngle, T eqLength, T eqArea,
                       T k);

}  // namespace cellModelHelper3D


}  // namespace shellModelHelper3D

}  // namespace plb

#endif  // CELL_MODEL_3D_H
