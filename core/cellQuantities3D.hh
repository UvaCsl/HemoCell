/* This file is part of the ficsion framework.

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

#ifndef CELL_QUANTITIES_3D_HH
#define CELL_QUANTITIES_3D_HH


#include "cellQuantities3D.h"

// Class storing all the cell measures
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
CellQuantities3D<T,Descriptor,ParticleFieldT>::CellQuantities3D(
        TriangleBoundary3D<T> const& Cells_,
        MultiParticleField3D<ParticleFieldT<T,Descriptor> > &  particles_,
        std::vector<plint> const&  cellIds_, plint numberOfCells_,
        CellField3D<T,Descriptor> & chq_,
        std::string cmhFileName, T dx_, T dt_, plint cellRadiusLU_, bool checkpointed_) :
            Cells(Cells_),
            particles(&particles_),
            cellIds(cellIds_),
            numberOfCells(numberOfCells_),
            chq(chq_),
            tagToParticle3D(chq.get_tagToParticle3D()),
            dx(dx_), dt(dt_), cellRadiusLU(cellRadiusLU_), checkpointed(checkpointed_)
{
        std::ostream::openmode mode = std::ostream::out;
        if (not checkpointed) { mode = mode | std::ostream::trunc; }
        else { mode = mode | std::ostream::app; }
        logFile.open(cmhFileName.c_str(), mode);
        if (not checkpointed_) {
            writeHeader();
        }
        MultiBlockManagement3D const& latticeManagement(particles->getMultiBlockManagement());
        MultiBlockManagement3D reductionParticleManagement (
                latticeManagement.getSparseBlockStructure(),
                latticeManagement.getThreadAttribution().clone(),
                cellRadiusLU*5,
                latticeManagement.getRefinementLevel() );

        reductionParticles = new MultiParticleField3D<DenseParticleField3D<T,Descriptor> >(reductionParticleManagement,
                defaultMultiBlockPolicy3D().getCombinedStatistics() );
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
CellQuantities3D<T,Descriptor,ParticleFieldT>::~CellQuantities3D() { delete reductionParticles; };


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::calculateAll()
{

        cellsMeanTriangleArea.clear(); cellsMeanEdgeDistance.clear();
        cellsMaxEdgeDistance.clear(); cellsMeanAngle.clear(); cellsCenter.clear(); cellsVelocity.clear();
        cellsMeanTileSpan.clear();

        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);
        NumVerticesCellReduceFunctional3D<T,Descriptor> nfunctional(Cells, cellIds, numberOfCells);
        applyProcessingFunctional(nfunctional, particles->getBoundingBox(), particleArg);
        nfunctional.getCellQuantityArray(cellNumVertices, cellIds);

        calculateVolumeAndSurface();
        countCellMeanTriangleArea(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsMeanTriangleArea);
        countCellMeanAngle(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsMeanAngle);
        countCellMeanEdgeDistance(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsMeanEdgeDistance);
        countCellMeanTileSpan(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsMeanTileSpan);
        countCellMaxEdgeDistance(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsMaxEdgeDistance);
        countCellMaxEdgeDistance(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsMaxEdgeDistance);
        countCellCenters(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsCenter, cellNumVertices);
        countCellVelocity(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsVelocity, cellNumVertices);

        plint numVerticesPerCell = cellNumVertices[0];

        chq.clearQuantities();
        std::vector<plint> subscribedQuantities;
        subscribedQuantities.push_back(CCR_POSITION_MEAN);
        subscribedQuantities.push_back(CCR_VOLUME);
        subscribedQuantities.push_back(CCR_SURFACE);
        subscribedQuantities.push_back(CCR_EDGE_DISTANCE_MEAN);
        subscribedQuantities.push_back(CCR_EDGE_DISTANCE_MIN);
        subscribedQuantities.push_back(CCR_EDGE_DISTANCE_MAX);
        subscribedQuantities.push_back(CCR_ANGLE_MEAN);
        subscribedQuantities.push_back(CCR_ANGLE_MIN);
        subscribedQuantities.push_back(CCR_ANGLE_MAX);
        subscribedQuantities.push_back(CCR_ENERGY);

        particleArg.push_back(reductionParticles);
        applyProcessingFunctional(
                new LocalCellReductions3D<T,Descriptor>(Cells, chq, numVerticesPerCell, subscribedQuantities),
                particles->getBoundingBox(), particleArg);
        particleArg.clear();  particleArg.push_back(reductionParticles);
        syncCellQuantities<T,Descriptor>(reductionParticles->getBoundingBox(), particleArg, chq);

        for (plint ci = 0; ci < numberOfCells; ++ci) {
            cout << ci
                    << " vol " << cellsVolume[ci] - chq.getVolume(ci)
                    << " sur " << cellsSurface[ci] - chq.getSurface(ci)
                    << " ang " << cellsMeanAngle[ci] -  chq.getMeanAngle(ci)
                    << " v-v " << cellsVolume[ci] << " " << chq.getVolume(ci)
                                        << std::endl;
        }

}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::calculateVolumeAndSurfaceAndCenters()
{
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);
        cellsVolume.clear(); cellsSurface.clear(); cellsCenter.clear();
        countCellVolume(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsVolume, tagToParticle3D);
        countCellSurface(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsSurface);
        countCellCenters(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsCenter, cellNumVertices);

//        chq.clearQuantities();
//        std::vector<plint> subscribedQuantities;
//        subscribedQuantities.push_back(CCR_POSITION_MEAN);
//        subscribedQuantities.push_back(CCR_VOLUME);
//        subscribedQuantities.push_back(CCR_SURFACE);
//
//        particleArg.push_back(reductionParticles);
//        applyProcessingFunctional(
//                new LocalCellReductions3D<T,Descriptor>(Cells, chq, chq.getNumVerticesPerCell(), subscribedQuantities),
//                particles->getBoundingBox(), particleArg);
//        particleArg.clear();  particleArg.push_back(reductionParticles);
//        syncCellQuantities<T,Descriptor>(reductionParticles->getBoundingBox(), particleArg, chq);
//        for (plint ci = 0; ci < numberOfCells; ++ci) {
//            cellsVolume[ci] = chq.getVolume(ci);
//            cellsSurface[ci] = chq.getSurface(ci);
//            cellsCenter[ci] = chq.getPosition(ci);
//        }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::calculateInertia()
{
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);
        particleArg.push_back(reductionParticles);

        std::vector<plint> subscribedQuantities;
        subscribedQuantities.push_back(CCR_INERTIA);

        applyProcessingFunctional(
                new LocalCellReductions3D<T,Descriptor>(Cells, chq, chq.getNumVerticesPerCell(), subscribedQuantities),
                particles->getBoundingBox(), particleArg);

        particleArg.clear();  particleArg.push_back(reductionParticles);
        syncCellQuantities<T,Descriptor>(reductionParticles->getBoundingBox(), particleArg, chq);
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::calculateVolumeAndSurface()
{
        cellsVolume.clear(); cellsSurface.clear();
        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);
        countCellVolume(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsVolume, tagToParticle3D);
        countCellSurface(Cells, *particles, particles->getBoundingBox(), cellIds, numberOfCells, cellsSurface);
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::writeHeader() {
    std::string delim(" ; ");
    logFile << "# Iteration " << delim
            << " Volume" << delim
            << " Surface" << delim
            << " Mean Triangle Surface" << delim
            << " Mean Edge Distance" << delim
            << " Mean Angle [^o]" << delim
            << " Mean Edge Distance [LU]" << delim
            << " Max Edge Distance [LU]" << delim
            << " x [LU]" << delim
            << " y [LU]" << delim
            << " z [LU]" << delim
            << " vx [LU]" << delim
            << " vy [LU]" << delim
            << " vz [LU]" << delim
            << "00" << std::endl;
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::write(plint iter, T eqVolume, T eqSurface, T eqArea, T eqLength)
{
        std::string delim(" ; ");
        if (iter != 0 or (iter == 0 and not checkpointed)) {
            logFile << iter*1.0 << delim
                    << cellsVolume[0]*100.0/eqVolume - 100 <<  delim
                    << cellsSurface[0]*100.0/eqSurface - 100 << delim
                    << cellsMeanTriangleArea[0]*100.0/eqArea - 100 << delim
                    << cellsMeanEdgeDistance[0]*100.0/eqLength - 100 << delim
                    << cellsMeanAngle[0]*180.0/pi << delim
                    << cellsMeanEdgeDistance[0] << delim
                    << cellsMaxEdgeDistance[0] << delim
                    << cellsCenter[0][0] << delim
                    << cellsCenter[0][1] << delim
                    << cellsCenter[0][2] << delim
                    << cellsVelocity[0][0] << delim
                    << cellsVelocity[0][1] << delim
                    << cellsVelocity[0][2] << delim
                    << "00" << std::endl;
        }
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void CellQuantities3D<T,Descriptor,ParticleFieldT>::print(plint iter, T eqVolume, T eqSurface, T eqArea, T eqLength)
{
        pcout << "=== " << iter << " === " << std::endl;
        pcout << "Volume: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsVolume[iA]*100.0/eqVolume - 100 <<"%, ";
        pcout << std::endl <<"Surface: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsSurface[iA]*100.0/eqSurface - 100 << "%, ";
        pcout << std::endl <<"Mean Edge Distance: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsMeanEdgeDistance[iA]*100.0/eqLength - 100<< "%, ";
        pcout << std::endl <<"Mean Angle [^o]: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsMeanAngle[iA]*180.0/pi << ", ";
        pcout << std::endl <<"Mean Triangle Surface [LU^2]: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsMeanTriangleArea[iA] << ", ";
        pcout << std::endl <<"Mean Edge Distance [LU]: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsMeanEdgeDistance[iA] << ", ";
        pcout << std::endl <<"Max Edge Distance  [LU]: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsMaxEdgeDistance[iA] << ", ";
        pcout << std::endl <<"Volume  [mu{m}^3]: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsVolume[iA]*dx*dx*dx*1e18 <<", ";
        pcout << std::endl <<"Surface [mu{m}^2]: "; for (pluint iA = 0; iA < cellsVolume.size(); ++iA) pcout << cellsSurface[iA]*dx*dx*1e12 << ", ";
        pcout << std::endl;
        for (pluint iA = 0; iA < cellsCenter.size(); ++iA)
            pcout <<"Coordinates: (" << cellsCenter[iA][0] << ", " << cellsCenter[iA][1] << ", " << cellsCenter[iA][2] << ")" << std::endl;
        for (pluint iA = 0; iA < cellsCenter.size(); ++iA)
            pcout <<"Velocity: (" << cellsVelocity[iA][0] << ", " << cellsVelocity[iA][1] << ", " << cellsVelocity[iA][2] << ")" << std::endl;
}


template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsVolume() { return cellsVolume; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsSurface() { return cellsSurface; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsMeanEdgeDistance() { return cellsMeanEdgeDistance; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsMaxEdgeDistance() { return cellsMaxEdgeDistance; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsMeanAngle() { return cellsMeanAngle; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsMeanTriangleArea() { return cellsMeanTriangleArea; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector<T> & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsMeanTileSpan() { return cellsMeanTileSpan; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector< Array<T,3> > & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsCenter() { return cellsCenter; }

template< typename T, template<typename U> class Descriptor,
    template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
std::vector< Array<T,3> > & CellQuantities3D<T,Descriptor,ParticleFieldT>::getCellsVelocity() { return cellsVelocity; }


#endif  // CELL_QUANTITIES_3D_HH


