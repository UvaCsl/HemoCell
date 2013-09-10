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

#ifndef SHAPE_MEMORY_MODEL_FUNCTIONAL_3D_HH
#define SHAPE_MEMORY_MODEL_FUNCTIONAL_3D_HH

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticleFunctional3D.h"
#include "immersedCellParticle3D.h"
#include "immersedBoundaryMethod3D.h"
#include <map>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include "shapeMemoryModelFunctional3D.h"
#include "immersedCellsReductions.hh"

namespace plb {


/* ******** ComputeShapeMemoryModelForce3D *********************************** */

template<typename T, template<typename U> class Descriptor>
ComputeShapeMemoryModelForce3D<T,Descriptor>::ComputeShapeMemoryModelForce3D (
        TriangleBoundary3D<T> const& triangleBoundary_,
        ShapeMemoryModel3D<T>* cellModel_,
        std::vector<T> const& cellsVolume_, std::vector<T> const& cellsSurface_)
    : triangleBoundary(triangleBoundary_),
      cellModel(cellModel_),
      cellsVolume(cellsVolume_),
      cellsSurface(cellsSurface_)
{ }

template<typename T, template<typename U> class Descriptor>
ComputeShapeMemoryModelForce3D<T,Descriptor>::~ComputeShapeMemoryModelForce3D()
{
    delete cellModel;
}

template<typename T, template<typename U> class Descriptor>
ComputeShapeMemoryModelForce3D<T,Descriptor>::ComputeShapeMemoryModelForce3D (
            ComputeShapeMemoryModelForce3D<T,Descriptor> const& rhs)
    : triangleBoundary(rhs.triangleBoundary),
      cellModel(rhs.cellModel->clone()),
      cellsVolume(rhs.cellsVolume),
      cellsSurface(rhs.cellsSurface)

{ }


template<typename T, template<typename U> class Descriptor>
void ComputeShapeMemoryModelForce3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);

    std::map< plint, Array<T,3> > particleVelocity;
//    std::map< plint, Array<T,3> > particleForces;
    std::map< plint, Array<T,3>* > particleForces;
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getTag();
        particleVelocity[vertexId] = particle->get_v();
//        particleForces[vertexId] = Array<T,3>(0., 0., 0.);
        particleForces[vertexId] = new Array<T,3> [7]; // [f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity, E_bending]
        for (pluint var = 0; var < 7; ++var) {
            particleForces[vertexId][var].resetToZero();
        }
    }
    // T eqArea = cellModel->getEquilibriumTriangleArea();
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getTag();
        T iSurface = T(); // one third of the sum of the areas of all triangles that share the given vertex.
        plint cellId = particle->get_cellId();
        if (!isRigid(triangleBoundary.getVertexProperty(vertexId))) {
            Array<T,3> elasticForce = cellModel->computeElasticForce (
                    triangleBoundary, vertexId );
            Array<T,3> f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity;
            f_wlc.resetToZero(); f_bending.resetToZero(); f_volume.resetToZero();
            f_surface.resetToZero(); f_shear.resetToZero(); f_viscosity.resetToZero();
            Array<T,3> cellForce = cellModel->computeCellForce (
                    triangleBoundary, cellsVolume[cellId], cellsSurface[cellId], iSurface, particleVelocity, particleForces, vertexId,
                    f_wlc, f_bending, f_volume, f_surface, f_shear, f_viscosity);
            particle->get_f_wlc() = f_wlc;
            particle->get_f_bending() = f_bending;
            particle->get_f_volume() = f_volume;
            particle->get_f_surface() = f_surface;
            particle->get_f_shear() = f_shear;
            particle->get_f_viscosity() = f_viscosity;
            particle->get_E_bending().resetToZero();

            Array<T,3> force, acc; force.resetToZero(); acc.resetToZero();
            force = elasticForce + cellForce;
            acc = force*1.0 / cellModel->getDensity();
            particle->get_a() = acc;
            particle->get_force() = force;
            particle->get_stress() = force*1.0/iSurface;
        }
    }
    Array<T,3> sforce; sforce.resetToZero();
    Array<T,3> sf_wlc, sf_bending, sf_volume, sf_surface, sf_shear, sf_viscosity;
    sf_wlc.resetToZero(); sf_bending.resetToZero(); sf_volume.resetToZero(); sf_surface.resetToZero();
    sf_shear.resetToZero(); sf_viscosity.resetToZero();
    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypedParticle = found[iParticle];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        plint vertexId = particle->getTag();

        particle->get_f_wlc() += particleForces[vertexId][0];
        particle->get_f_bending() += particleForces[vertexId][1];
        particle->get_f_volume() += particleForces[vertexId][2];
        particle->get_f_surface() += particleForces[vertexId][3];
        particle->get_f_shear() += particleForces[vertexId][4];
        particle->get_f_viscosity() += particleForces[vertexId][5];
        for (pluint var = 0; var < 6; ++var) {
            particle->get_force() += particleForces[vertexId][var];
        }
        particle->get_E_bending() += particleForces[vertexId][6];
        delete [] particleForces[vertexId];

        sforce += particle->get_force();
        sf_wlc += particle->get_f_wlc();
        sf_bending += particle->get_f_bending();
        sf_volume += particle->get_f_volume();
        sf_surface += particle->get_f_surface();
        sf_shear += particle->get_f_shear();
        sf_viscosity += particle->get_f_viscosity();
    }
    bool pcoutForceSum = true;
	#ifdef PLB_MPI_PARALLEL
    	int ntasks;
		MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
		pcoutForceSum = (ntasks == 1);
	#endif
    if (pcoutForceSum and (fabs(sforce[0]) + fabs(sforce[1]) + fabs(sforce[2]) > 1e-10)) {
        pcout << "sforce (" << sforce[0] << ", " << sforce[1] << ", " << sforce[2] << ") " << std::endl;
        pcout << "sf_wlc (" << sf_wlc[0] << ", " << sf_wlc[1] << ", " << sf_wlc[2] << ") " << std::endl;
        pcout << "sf_bending (" << sf_bending[0] << ", " << sf_bending[1] << ", " << sf_bending[2] << ") " << std::endl;
        pcout << "sf_volume (" << sf_volume[0] << ", " << sf_volume[1] << ", " << sf_volume[2] << ") " << std::endl;
        pcout << "sf_surface (" << sf_surface[0] << ", " << sf_surface[1] << ", " << sf_surface[2] << ") " << std::endl;
        pcout << "sf_shear (" << sf_shear[0] << ", " << sf_shear[1] << ", " << sf_shear[2] << ") " << std::endl;
        pcout << "sf_viscosity (" << sf_viscosity[0] << ", " << sf_viscosity[1] << ", " << sf_viscosity[2] << ") " << std::endl;
    }
}

template<typename T, template<typename U> class Descriptor>
ComputeShapeMemoryModelForce3D<T,Descriptor>*
    ComputeShapeMemoryModelForce3D<T,Descriptor>::clone() const
{
    return new ComputeShapeMemoryModelForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeShapeMemoryModelForce3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeShapeMemoryModelForce3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeShapeMemoryModelForce3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}


template<typename T>
void getCellShapeQuantitiesFromMesh(TriangleBoundary3D<T>& boundary,
                            vector<T> & eqAreaPerTriangle, map<plint,T> & eqLengthPerEdge, map<plint,T> & eqAnglePerEdge,
                            plint cellNumTriangles, plint cellNumPartsPerCell) {
    TriangularSurfaceMesh<T> const& dynMesh = boundary.getMesh();
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

#endif  // SHAPE_MEMORY_MODEL_FUNCTIONAL_3D_HH

