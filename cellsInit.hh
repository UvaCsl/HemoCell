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

#ifndef CELLS_INIT_HH
#define CELLS_INIT_HH


#include "cellsInit.h"

template<typename T>
void positionCells(plint shape, T radius, plint & npar, IncomprFlowParam<T> const& parameters,
        std::vector<Array<T,3> > & centers, std::vector<T> & radii, plint flowType) {

    std::vector<T> posX, posY, posZ;
    const plint nx = parameters.getNx() - 1 ;
    const plint ny = parameters.getNy()  - 1;
    const plint nz = parameters.getNz()  - 1;
    T dX = 2.1 * radius ;
    T dY = 2.1 * radius * ( (shape==1) ? 0.5 : 1 );
    T dZ = 2.1 * radius;

    plint NdX = (nx-2)*1.0/dX;
    plint NdY = (ny-2)*1.0/dY;
    plint NdZ = (nz-2)*1.0/dZ;
    dX = (nx-2.0) / NdX ;
    dY = (ny-2.0) / NdY ;
    dZ = (nz-2.0) / NdZ ;

    npar = npar<(NdX*NdY*NdZ)?npar:(NdX*NdY*NdZ);
    plint slices = npar/(NdY*NdZ);

    for (plint i = 0; i < slices; ++i) {
        for (plint iy = 0; iy < NdY; ++iy) {
            for (plint iz = 0; iz < NdZ; ++iz) {
                posX.push_back((i+0.5)*dX);
                posY.push_back((iy+0.5)*dY);
                posZ.push_back((iz+0.5)*dZ);
            }
        }
    }
    plint lastSlice = npar%(NdY*NdZ);
    plint rows = lastSlice/NdZ;
    for (plint iy = 1; iy <= rows; ++iy) {
        for (plint iz = 0; iz < NdZ; ++iz) {
            posX.push_back((slices+0.5)*dX);
            posY.push_back(ny * iy*1.0/(rows+1 + 1.0));
            posZ.push_back((iz+0.5)*dZ);
        }
    }

    plint mods = lastSlice%NdZ;
    for (plint iz = 1; iz <= mods; ++iz) {
        posX.push_back((slices+0.5)*dX);
        posY.push_back(ny * (rows+1)*1.0/(rows+1 + 1.0));
        posZ.push_back(nz * iz*1.0/(mods + 1.0));
    }

    T addToX = 0.5;
//    addToX = (NdX - slices) * dX * 0.5 ;


    for (pluint iA = 0; iA < posX.size(); ++iA) {
        centers.push_back(Array<T,3>(posX[iA]+addToX,posY[iA],posZ[iA]));
        radii.push_back(radius);
    }
}




template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void calculateCellMeasures(TriangleBoundary3D<T> Cells, MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles,
                           std::vector<plint> & cellIds, plint numberOfCells,
                           std::vector<T> & cellsVolume, std::vector<T> & cellsSurface, std::vector<T> & cellsMeanTriangleArea,
                           std::vector<T> & cellsMeanEdgeDistance, std::vector<T> & cellsMaxEdgeDistance, std::vector<T> & cellsMeanAngle,
                           std::vector< Array<T,3> > & cellsCenter, std::vector< Array<T,3> > & cellsVelocity,
                           std::vector<T> & cellsMeanTileSpan,
                           std::map <plint, Particle3D<T,Descriptor>*> const & iVertexToParticle3D)
    {
    cellsVolume.clear(); cellsSurface.clear(); cellsMeanTriangleArea.clear(); cellsMeanEdgeDistance.clear();
    cellsMaxEdgeDistance.clear(); cellsMeanAngle.clear(); cellsCenter.clear(); cellsVelocity.clear();
    cellsMeanTileSpan.clear();

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    std::vector<T> cellNumVertices;
    NumVerticesCellReduceFunctional3D<T,Descriptor> nfunctional(Cells, cellIds, numberOfCells);
    applyProcessingFunctional(nfunctional, particles.getBoundingBox(), particleArg);
    nfunctional.getCellQuantityArray(cellNumVertices, cellIds);

    countCellVolume(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsVolume, iVertexToParticle3D);
    countCellSurface(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsSurface);
    countCellMeanTriangleArea(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsMeanTriangleArea);
    countCellMeanAngle(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsMeanAngle);
    countCellMeanEdgeDistance(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsMeanEdgeDistance);
    countCellMeanTileSpan(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsMeanTileSpan);
    countCellMaxEdgeDistance(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsMaxEdgeDistance);
    countCellMaxEdgeDistance(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsMaxEdgeDistance);
    countCellCenters(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsCenter, cellNumVertices);
    countCellVelocity(Cells, particles, particles.getBoundingBox(), cellIds, numberOfCells, cellsVelocity, cellNumVertices);
}


template<typename T>
void printCellMeasures(plint i, TriangleBoundary3D<T> Cells,
                       std::vector<T> & cellsVolume, std::vector<T> & cellsSurface, std::vector<T> & cellsMeanTriangleArea,
                       std::vector<T> & cellsMeanEdgeDistance, std::vector<T> & cellsMaxEdgeDistance, std::vector<T> & cellsMeanAngle,
                       std::vector< Array<T,3> > & cellsCenter, std::vector< Array<T,3> > & cellsVelocity,
                       T eqVolume, T eqSurface, T eqArea, T eqLength, T dx, T dt) {
    pcout << "=== " << i << " === " << std::endl;
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


template<typename T>
void writeCellLog(plint i, plb_ofstream & logFile,
                  std::vector<T> & cellsVolume, std::vector<T> & cellsSurface, std::vector<T> & cellsMeanTriangleArea, std::vector<T> & cellsMeanEdgeDistance,
                  std::vector<T> & cellsMaxEdgeDistance, std::vector<T> & cellsMeanAngle,
                  std::vector< Array<T,3> > & cellsCenter, std::vector< Array<T,3> > & cellsVelocity,
                  T eqVolume, T eqSurface, T eqArea, T eqLength) {
    std::string delim("; ");
    if (i==0) {
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
    logFile << i*1.0 << delim
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

#endif  // CELLS_INIT_HH


