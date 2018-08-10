/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef HEMO_OCTREE_H
#define HEMO_OCTREE_H

#include "hemoCellParticle.h" // Need to make pointers to particle object
#include <vector> 

namespace hemo {
  class OctreeStructCell {     
    private:
      hemo::Array<double, 6> bBox;
      std::vector<hemo::Array<plint,3>> triangle_list; // Keep all the unsorted obs here
      OctreeStructCell * nodes[8];
      std::vector<hemo::Array<T, 6>> octantOfBoundingBox();
      
      plint maxDivisions;
      plint level; // To keep track of how many octants we have already 
      bool stopDividing; // To denote octants which are cool
      bool finalNode = false; // To track if needs more subdividing
      plint limit;
      
    public:
      OctreeStructCell(plint divis, plint l, unsigned int lim, hemo::Array<double, 6> bbox,
                       std::vector<hemo::Array<plint,3>> triangle_list_,
                       std::vector<HemoCellParticle>* part, const std::vector<int>  cell);
      ~OctreeStructCell();
      void constructTree(std::vector<HemoCellParticle>* part,  std::vector<int>  cell, std::vector<hemo::Array<plint,3>> triangle_list_);
      int returnTrianglesAmount();
      void findCrossings(hemo::Array<plint, 3> latticeSite, std::vector<hemo::Array<plint,3>> &);
  };
}

#endif
