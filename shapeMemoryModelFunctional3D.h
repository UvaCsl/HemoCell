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

#ifndef SHAPE_MEMORY_MODEL_FUNCTIONAL_3D_H
#define SHAPE_MEMORY_MODEL_FUNCTIONAL_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticle3D.h"
#include "immersedBoundaryMethod3D.h"
#include <map>
#include "computeCellForces3D.h"
#include "shapeMemoryModel3D.hh"

namespace plb {


template<typename T, template<typename U> class Descriptor>
class ComputeShapeMemoryModelForce3D : public BoxProcessingFunctional3D
{
public:
    ComputeShapeMemoryModelForce3D (
            TriangleBoundary3D<T> const& triangleBoundary_,
            ShapeMemoryModel3D<T>* cellModel_,
            std::vector<T> const& cellsVolume_, std::vector<T> const& cellsSurface_);
    ~ComputeShapeMemoryModelForce3D();
    ComputeShapeMemoryModelForce3D(ComputeShapeMemoryModelForce3D<T,Descriptor> const& rhs);
    /// Arguments: [0] Particle-field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeShapeMemoryModelForce3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    TriangleBoundary3D<T> const& triangleBoundary;
    ShapeMemoryModel3D<T>* cellModel;
    std::vector<T> const& cellsVolume;
    std::vector<T> const& cellsSurface;
};


template<typename T>
void getCellShapeQuantitiesFromMesh(TriangleBoundary3D<T>& boundary,
                            vector<T> & eqAreaPerTriangle, map<plint,T> & eqLengthPerEdge, map<plint,T> & eqAnglePerEdge,
                            plint cellNumTriangles, plint cellNumPartsPerCell);


}  // namespace plb

#endif  // SHAPE_MEMORY_MODEL_FUNCTIONAL_3D_H

