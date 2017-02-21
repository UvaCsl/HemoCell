#ifndef CELL_IN_SHEAR_FLOW_3D_HH
#define CELL_IN_SHEAR_FLOW_3D_HH

#include "cellInShearFlow3D.h"

namespace plb {

/* ******** SingleCellInShearFlow *********************************** */

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
SingleCellInShearFlow<T,Descriptor,ParticleFieldT>::SingleCellInShearFlow(TriangleBoundary3D<T> const& Cells_,
        MultiParticleField3D<ParticleFieldT<T,Descriptor> > & particles_, std::vector<plint> cellIds,
        std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume,
        plint numParticlesPerSide_, plint flowType_, T dx_, T dt_, T dNewton_,
        HemoCellField & chq_,
        bool checkpointed_,
        bool store_)
        :
    Cells(Cells_), particles(particles_), numParticlesPerSide(numParticlesPerSide_),
    flowType(flowType_), dx(dx_), dt(dt_), dNewton(dNewton_), dm(dNewton_*dt_*dt_/dx_),
    chq(chq_),
    tagToParticle3D(chq.get_tagToParticle3D()),
    checkpointed(checkpointed_),
    store(store_)
{
    lateralCellParticleTags.push_back(&outerLeftTags);
    lateralCellParticleTags.push_back(&outerRightTags);
    lateralCellParticleTags.push_back(&outerFrontTags);
    lateralCellParticleTags.push_back(&outerBackTags);
    lateralCellParticleTags.push_back(&outerUpTags);
    lateralCellParticleTags.push_back(&outerDownTags);

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    applyProcessingFunctional (
        new FindTagsOfLateralCellParticles3D<T,Descriptor>(numParticlesPerSide, &outerLeftTags, &outerRightTags, 0),
        particles.getBoundingBox(), particleArg );
    applyProcessingFunctional (
        new FindTagsOfLateralCellParticles3D<T,Descriptor>(numParticlesPerSide, &outerUpTags, &outerDownTags, 1),
        particles.getBoundingBox(), particleArg );
    applyProcessingFunctional (
        new FindTagsOfLateralCellParticles3D<T,Descriptor>(numParticlesPerSide, &outerFrontTags, &outerBackTags, 2),
        particles.getBoundingBox(), particleArg );

    std::ostream::openmode mode = std::ostream::out;
    if (not checkpointed) { mode = mode | std::ostream::trunc; }
    else { mode = mode | std::ostream::app; }

    shearResultFile.open((global::directories().getLogOutDir() + "plbShearResults.log").c_str(), mode);
    writeHeader();

    if (flowType==6) {
        updateQuantities(0, cellIds, cellCenters, cellsVolume);
        maxDiameter = max(diameters[0][0], diameters[0][1]);
        maxDiameter = max(maxDiameter,  diameters[0][2]);
    }
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void SingleCellInShearFlow<T,Descriptor,ParticleFieldT>::updateQuantities(plint iteration, std::vector<plint> cellIds,
        std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume)
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    applyProcessingFunctional (
        new MeasureCellStretchDeformation3D<T,Descriptor>(lateralCellParticleTags, &tankTreadingLengths, &tankTreadingAngles, tagToParticle3D),
        particles.getBoundingBox(), particleArg );

    std::vector< std::vector<T> > ellipsoidAngles;
    std::vector< std::vector<T> > ellipsoidSemiAxes;
    std::vector<T> difference;
    std::vector< std::vector<T> > inertia;
    computeEllipsoidFit (Cells, particles, cellIds, cellCenters, cellsVolume,
            ellipsoidAngles, ellipsoidSemiAxes, inertia, difference);
    T currMaxDiameter =  max(max(ellipsoidSemiAxes[0][0], ellipsoidSemiAxes[0][1]),  ellipsoidSemiAxes[0][2]);
    if (maxDiameter <= 0) { maxDiameter = currMaxDiameter; }


    if ( store or (not store and (iterations.size() < 1) ) ) {
        iterations.push_back(iteration);
        deformationIndex.push_back((currMaxDiameter - maxDiameter)*1.0/(currMaxDiameter + maxDiameter));
        diameters.push_back(Array<T,3>(ellipsoidSemiAxes[0][0], ellipsoidSemiAxes[0][1], ellipsoidSemiAxes[0][2]));
        tumblingAngles.push_back(Array<T,3>(ellipsoidAngles[0][0], ellipsoidAngles[0][1], ellipsoidAngles[0][2]));
        symmetryDeviation.push_back(difference[0]);
        inertiaTensor.push_back(inertia[0]);
    } else {
        iterations[0] = iteration;
        deformationIndex[0] = (currMaxDiameter - maxDiameter)*1.0/(currMaxDiameter + maxDiameter);
        diameters[0] = Array<T,3>(ellipsoidSemiAxes[0][0], ellipsoidSemiAxes[0][1], ellipsoidSemiAxes[0][2]);
        tumblingAngles[0] = Array<T,3>(ellipsoidAngles[0][0], ellipsoidAngles[0][1], ellipsoidAngles[0][2]);
        symmetryDeviation[0] = difference[0];
        inertiaTensor[0] =  inertia[0];
    }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void SingleCellInShearFlow<T,Descriptor,ParticleFieldT>::writeHeader(bool writeOutput) {
        if (writeOutput and (not checkpointed)) {
            shearResultFile <<
                    "# time (sec) ; " <<
                    "deformationIndex ; " <<
                    "diameter x (m) ; " <<
                    "diameter y (m) ; " <<
                    "diameter z (m) ; " <<
                    "tumblingAngles x (rad) ; " <<
                    "tumblingAngles y (rad) ; " <<
                    "tumblingAngles z (rad) ; " <<
                    /* tankTreadingAngles */
                    "tankTreadingAngles LR x (rad) ; " <<
                    "tankTreadingAngles LR y (rad) ; " <<
                    "tankTreadingAngles LR z (rad) ; " <<

                    "tankTreadingAngles UD x (rad) ; " <<
                    "tankTreadingAngles UD y (rad) ; " <<
                    "tankTreadingAngles UD z (rad) ; " <<

                    "tankTreadingAngles FB x (rad) ; " <<
                    "tankTreadingAngles FB y (rad) ; " <<
                    "tankTreadingAngles FB z (rad) ; " <<
                    "tankTreadingLengths x (m) ; " <<
                    "tankTreadingLengths y (m) ; " <<
                    "tankTreadingLengths z (m) ; " <<
                    "symmerty deviation ; " <<
                    "inertia tensor Ixx (kg/m2) ; " <<
                    "inertia tensor Ixy (kg/m2) ; " <<
                    "inertia tensor Ixz (kg/m2) ; " <<
                    "inertia tensor Iyx (kg/m2) ; " <<
                    "inertia tensor Iyy (kg/m2) ; " <<
                    "inertia tensor Iyz (kg/m2) ; " <<
                    "inertia tensor Izx (kg/m2) ; " <<
                    "inertia tensor Izy (kg/m2) ; " <<
                    "inertia tensor Izz (kg/m2) " <<
                    std:: endl;
    }
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void SingleCellInShearFlow<T,Descriptor,ParticleFieldT>::write(bool writeOutput) {
        if (writeOutput) {
            std::string delim = " ; ";
            plint N = iterations.size() - 1;
            shearResultFile <<
                iterations[N] * dt<< delim <<
                deformationIndex[N] << delim <<
                diameters[N][0] * dx<< delim <<
                diameters[N][1] * dx<< delim <<
                diameters[N][2] * dx<< delim <<
                tumblingAngles[N][0] << delim <<
                tumblingAngles[N][1] << delim <<
                tumblingAngles[N][2] << delim;
            for (pluint d = 0; d < tankTreadingAngles.size(); ++d) {
                shearResultFile <<
                    tankTreadingAngles[d][0] << delim <<
                    tankTreadingAngles[d][1] << delim <<
                    tankTreadingAngles[d][2] << delim;
            }
            for (pluint  d = 0; d < tankTreadingLengths.size(); ++d) {
                shearResultFile << tankTreadingLengths[d] * dx<< delim ;
            }
            shearResultFile << symmetryDeviation[N] ;
            for (int i = 0; i < 9; ++i) {
                shearResultFile << delim << inertiaTensor[N][i] * dm / (dx*dx);
            }
            shearResultFile << std:: endl;
            std::vector<Array<T,3> > positions;
            std::vector<plint>  tags;
            for (pluint iTag = 0; iTag < lateralCellParticleTags.size(); ++iTag) {
                std::vector<plint>* const pTags = lateralCellParticleTags[iTag];
                for (pluint iVertex = 0; iVertex < pTags->size(); ++iVertex) {
                    SurfaceParticle3D* particle =
                            dynamic_cast<SurfaceParticle3D*> (tagToParticle3D[ (*pTags)[iVertex]]);

                    positions.push_back(particle->get_pbcPosition());
                    tags.push_back(iTag);
                }
            }
            writeImmersedPointsVTK(positions, tags, dx,
                global::directories().getOutputDir()+createFileName("LateralParticles.",iterations[N],10)+".vtk");
        }

}



/* ================================================================================ */
/* ******** InertiaTensorCellReduceFunctional3D *********************************** */
/* ================================================================================ */
template<typename T, template<typename U> class Descriptor>
InertiaTensorCellReduceFunctional3D<T,Descriptor>::InertiaTensorCellReduceFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_,
        std::vector< Array<T,3> >& cellCenters_)
    : triangleBoundary(triangleBoundary_)
{
    numberOfCells = cellIds_.size();
    for (pluint i=0; i< (pluint) (numberOfCells); ++i) {
        quantityIds[ cellIds_[i] ].resetToZero();
        cellCenters[ cellIds_[i] ] = cellCenters_[i];
        Array<plint,9> inertiaTensorId; inertiaTensorId.resetToZero();
        for (pluint j=0; j< (pluint) 9; ++j) {
            inertiaTensorId[j] = this->getStatistics().subscribeSum();
        }
        quantityIds[ cellIds_[i] ] = inertiaTensorId;
    }
}


template<typename T, template<typename U> class Descriptor>
InertiaTensorCellReduceFunctional3D<T,Descriptor>::InertiaTensorCellReduceFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
    : triangleBoundary(triangleBoundary_)
{
    numberOfCells = cellIds_.size();
    for (pluint i=0; i< (pluint) (numberOfCells); ++i) {
        quantityIds[ cellIds_[i] ].resetToZero();
        cellCenters[ cellIds_[i] ] = Array<T,3>(0.,0.,0.);
        Array<plint,9> inertiaTensorId; inertiaTensorId.resetToZero();
        for (pluint j=0; j< (pluint) 9; ++j) {
            inertiaTensorId[j] = this->getStatistics().subscribeSum();
        }
        quantityIds[ cellIds_[i] ] = inertiaTensorId;
    }
}


template<typename T, template<typename U> class Descriptor>
void InertiaTensorCellReduceFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    calculateQuantity(triangleMesh, particles, quantityIds);
}


template<typename T, template<typename U> class Descriptor>
void InertiaTensorCellReduceFunctional3D<T,Descriptor>::calculateQuantity(TriangularSurfaceMesh<T> & triangleMesh,
        std::vector<Particle3D<T,Descriptor>*> & particles, std::map<plint, Array<plint,9> > & quantityIds_)
{
    T rx, ry, rz;
    T Ixx, Ixy, Ixz;
    T Iyx, Iyy, Iyz;
    T Izx, Izy, Izz;
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        SurfaceParticle3D* particle =
                dynamic_cast<SurfaceParticle3D*> (nonTypedParticle);
        plint iVertex = particle->getTag();
        plint cellId = particle->get_cellId();
        Array<T,3> r0 = cellCenters[cellId];
        std::vector<plint> neighbors = triangleMesh.getNeighborTriangleIds(iVertex);
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            T Aj = triangleMesh.computeTriangleArea(neighbors[iB]);
            Array<T,3> nj = triangleMesh.computeTriangleNormal(neighbors[iB]);
            Array<T,3> rj = (triangleMesh.getVertex(neighbors[iB],0) +
                            triangleMesh.getVertex(neighbors[iB],1) + triangleMesh.getVertex(neighbors[iB],2))/3.0;
            rj = rj - r0;
            T dV = 1.0/5.0 * Aj * dot(nj, rj);
            rx = rj[0]; ry = rj[1]; rz = rj[2];
            Ixx = (rz*rz+ry*ry)*dV;
            Iyy = (rz*rz+rx*rx)*dV;
            Izz = (rx*rx+ry*ry)*dV;
            Iyx = Ixy = -rx*ry*dV;
            Izx = Ixz = -rx*rz*dV;
            Izy = Iyz = -ry*rz*dV;


            this->getStatistics().gatherSum(quantityIds_[cellId][0], Ixx/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][1], Ixy/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][2], Ixz/3.0); // every element is evaluated 3 times

            this->getStatistics().gatherSum(quantityIds_[cellId][3], Iyx/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][4], Iyy/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][5], Iyz/3.0); // every element is evaluated 3 times

            this->getStatistics().gatherSum(quantityIds_[cellId][6], Izx/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][7], Izy/3.0); // every element is evaluated 3 times
            this->getStatistics().gatherSum(quantityIds_[cellId][8], Izz/3.0); // every element is evaluated 3 times
        }
    }
}


template<typename T, template<typename U> class Descriptor>
void InertiaTensorCellReduceFunctional3D<T,Descriptor>::getCellQuantityArray(std::vector< std::vector<T> > & cellQuantity, std::vector<plint> cellIds_) const {
    cellQuantity.clear();
    for (pluint i=0; i< (pluint) (cellIds_.size()); ++i) {
        std::vector<T> tmp(9);
        Array<plint,9> inertiaTensorId = this->quantityIds.find(cellIds_[i])->second;
        for (pluint j=0; j < 9; ++j) {
            tmp[j] = this->getStatistics().getSum( inertiaTensorId[j] );
        }
        cellQuantity.push_back(tmp);
    }
}


/* ================================================================================ */
/* ******** computeCellInertia *********************************** */
/* ================================================================================ */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeCellInertia (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,3> >& cellCenters, std::vector< std::vector<T> >& cellInertia) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    InertiaTensorCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds, cellCenters);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellInertia, cellIds);
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeCellInertia (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector< Array<T,9> >& cellInertia) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    InertiaTensorCellReduceFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellQuantityArray(cellInertia, cellIds);
}


/* ================================================================================ */
/* ******** computeEllipsoidFit *********************************** */
/* ================================================================================ */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void computeEllipsoidFit (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, std::vector<plint> cellIds,
                std::vector< Array<T,3> > cellCenters, std::vector<T> cellsVolume,
                std::vector< std::vector<T> > & cellsEllipsoidFitAngles,
                std::vector< std::vector<T> > & cellsEllipsoidFitSemiAxes,
                std::vector< std::vector<T> > & cellInertia,
                std::vector<T> & difference) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    InertiaTensorCellReduceFunctional3D<T,Descriptor> inertiaTensorFunctional(Cells, cellIds, cellCenters);
    applyProcessingFunctional(inertiaTensorFunctional, particles.getBoundingBox(), particleArg);
    inertiaTensorFunctional.getCellQuantityArray(cellInertia, cellIds);

    for (pluint i = 0; i < cellIds.size(); ++i) {
        vector<T> semiAxes, ellipsoidAngles;
        T dif;
        getLambdasAndAngles(cellInertia[i], semiAxes, ellipsoidAngles, dif);

        T f1 = semiAxes[0] * 5 / cellsVolume[i];
        T f2 = semiAxes[1] * 5 / cellsVolume[i];
        T f3 = semiAxes[2] * 5 / cellsVolume[i];
        semiAxes[0] = (sqrt( (f1 + f2 - f3)/2 ));
        semiAxes[1] = (sqrt( (f2 + f3 - f1)/2 ));
        semiAxes[2] = (sqrt( (f3 + f1 - f2)/2 ));
        std::sort(semiAxes.begin(), semiAxes.end());
        cellsEllipsoidFitSemiAxes.push_back(semiAxes);
        cellsEllipsoidFitAngles.push_back(ellipsoidAngles);
        difference.push_back(dif);
    }
}


}
#endif // CELL_IN_SHEAR_FLOW_3D_HH
