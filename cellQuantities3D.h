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

#ifndef CELL_QUANTITIES_H
#define CELL_QUANTITIES_H

#include "immersedCellsReductions.h"

// Class storing all the cell measures
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
class CellQuantities3D {
public:
        CellQuantities3D(
                TriangleBoundary3D<T> const& Cells_,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> > &  particles_,
                std::vector<plint> const&  cellIds_, plint numberOfCells_,
                std::map <plint, Particle3D<T,Descriptor>*> const & iVertexToParticle3D_,
                std::string cmhFileName, T dx_, T dt_);
        virtual ~CellQuantities3D();
        void calculateAll() ;
        void calculateVolumeAndSurface() ;
        void print(plint iter, T eqVolume, T eqSurface, T eqArea, T eqLength);
        void write(plint iter, T eqVolume, T eqSurface, T eqArea, T eqLength);
public:
        std::vector<T> & getCellsVolume();
        std::vector<T> & getCellsSurface();
        std::vector<T> & getCellsMeanEdgeDistance() ;
        std::vector<T> & getCellsMaxEdgeDistance() ;
        std::vector<T> & getCellsMeanAngle() ;
        std::vector<T> & getCellsMeanTriangleArea() ;
        std::vector<T> & getCellsMeanTileSpan() ;
        std::vector< Array<T,3> > & getCellsCenter() ;
        std::vector< Array<T,3> > & getCellsVelocity() ;
private:
        std::vector<T> cellsVolume, cellsSurface;
        std::vector<T> cellsMeanEdgeDistance, cellsMaxEdgeDistance, cellsMeanAngle, cellsMeanTriangleArea, cellsMeanTileSpan;
        std::vector< Array<T,3> > cellsCenter, cellsVelocity;

        TriangleBoundary3D<T> const & Cells;
        MultiParticleField3D<ParticleFieldT<T,Descriptor> > * particles;
        std::vector<plint> const& cellIds;
        plint numberOfCells;
        std::map <plint, Particle3D<T,Descriptor>*> const & iVertexToParticle3D;
        plb_ofstream logFile;
        T dx, dt;
};




#include "cellQuantities3D.hh"
#endif  // CELLS_INIT_H


