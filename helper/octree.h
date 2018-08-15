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
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/blockLattice3D.hh"
#include "mollerTrumbore.h"


namespace hemo {
  class OctreeStructCell {     
    private:
      hemo::Array<T, 6> bBox;
      std::vector<hemo::Array<plint,3>> triangle_list = {}; // Keep all the unsorted obs here
      OctreeStructCell * nodes[8];
      
      plint maxDivisions;
      plint level; // To keep track of how many octants we have already 
      bool finalNode = false; // To track if needs more subdividing
      plint limit;
      
    public:
      OctreeStructCell(plint divis, plint l, unsigned int lim, hemo::Array<double, 6> bbox,
                       std::vector<hemo::Array<plint,3>> triangle_list_,
                       std::vector<HemoCellParticle>& part, const std::vector<int>  cell);
      OctreeStructCell(plint divis, plint l, unsigned int lim,
                       std::vector<hemo::Array<plint,3>> triangle_list_,
                       std::vector<HemoCellParticle>& part, const std::vector<int>  cell);
  private:
      void sharedConstructor(plint divis, plint l, unsigned int lim,
			std::vector<hemo::Array<plint,3>> triangle_list_,
			std::vector<HemoCellParticle> & part, const std::vector<int>  cell);
  public:
      ~OctreeStructCell();
      void constructTree(std::vector<HemoCellParticle>& part,  std::vector<int>  cell, std::vector<hemo::Array<plint,3>> triangle_list_);
      int returnTrianglesAmount();
      void findCrossings(hemo::Array<plint, 3> latticeSite, std::vector<hemo::Array<plint,3>> &);
      
      template<template<typename U> class Descriptor>
      void findInnerNodes(plb::BlockLattice3D<T,Descriptor> * fluid, std::vector<HemoCellParticle> & particles, const std::vector<int> & cell, std::vector<plb::Cell<T,Descriptor>*> & innerNodes) {
        innerNodes.clear();
        hemo::Array<T,6> bbox = bBox;
        //Adjust bbox to fit local atomic block
        bbox[0] = bbox[0] < fluid->getLocation().x ? fluid->getLocation().x : bbox[0];
        bbox[1] = bbox[1] > fluid->getLocation().x + fluid->getNx()-1 ? fluid->getLocation().x + fluid->getNx()-1: bbox[1];
        bbox[2] = bbox[2] < fluid->getLocation().y ? fluid->getLocation().y : bbox[2];
        bbox[3] = bbox[3] > fluid->getLocation().y + fluid->getNy()-1 ? fluid->getLocation().y + fluid->getNy()-1: bbox[3];
        bbox[4] = bbox[4] < fluid->getLocation().z ? fluid->getLocation().z : bbox[4];
        bbox[5] = bbox[5] > fluid->getLocation().z + fluid->getNz()-1 ? fluid->getLocation().z + fluid->getNz()-1: bbox[5];

        // Create a triple for-loop to go over all lattice points in the bounding box of a cell
        for (int x = (int)bbox[0]; x <= (int)bbox[1]+0.5; x++) { 
          for (int y = (int)bbox[2]; y <= (int)bbox[3]+0.5; y++) {
            for (int z = (int)bbox[4]; z <= (int)bbox[5]+0.5; z++) {
              int crossedCounter = 0; // How many triangles are crossed

              hemo::Array<plint, 3> latticeSite = {x, y, z};
              std::vector<hemo::Array<plint,3>> triangles_list;
              findCrossings(latticeSite,triangles_list);

              for (hemo::Array<plint, 3> triangle : triangles_list) {
                // Muller-trumbore intersection algorithm 
                const hemo::Array<double,3> & v0 = particles[cell[triangle[0]]].sv.position;
                const hemo::Array<double,3> & v1 = particles[cell[triangle[1]]].sv.position;
                const hemo::Array<double,3> & v2 = particles[cell[triangle[2]]].sv.position;

                crossedCounter += hemo::MollerTrumbore(v0, v1, v2, latticeSite);
              }

              // Count even-odd crossings
              if (crossedCounter%2) {
                int x_l = x-fluid->getLocation().x;
                int y_l = y-fluid->getLocation().y;
                int z_l = z-fluid->getLocation().z;
                innerNodes.push_back(&fluid->get(x_l,y_l,z_l));
              }
            }
          }
        }
      }
      
      template<template<typename U> class Descriptor>
      void findInnerNodes(plb::BlockLattice3D<T,Descriptor> * fluid, std::vector<HemoCellParticle> & particles, const std::vector<int> & cell, std::vector<Array<plint,3>> & innerNodes) {
        innerNodes.clear();
        hemo::Array<T,6> bbox = bBox;
        //Adjust bbox to fit local atomic block
        bbox[0] = bbox[0] < fluid->getLocation().x ? fluid->getLocation().x : bbox[0];
        bbox[1] = bbox[1] > fluid->getLocation().x + fluid->getNx()-1 ? fluid->getLocation().x + fluid->getNx()-1: bbox[1];
        bbox[2] = bbox[2] < fluid->getLocation().y ? fluid->getLocation().y : bbox[2];
        bbox[3] = bbox[3] > fluid->getLocation().y + fluid->getNy()-1 ? fluid->getLocation().y + fluid->getNy()-1: bbox[3];
        bbox[4] = bbox[4] < fluid->getLocation().z ? fluid->getLocation().z : bbox[4];
        bbox[5] = bbox[5] > fluid->getLocation().z + fluid->getNz()-1 ? fluid->getLocation().z + fluid->getNz()-1: bbox[5];

        // Create a triple for-loop to go over all lattice points in the bounding box of a cell
        for (int x = (int)bbox[0]; x <= (int)bbox[1]+0.5; x++) { 
          for (int y = (int)bbox[2]; y <= (int)bbox[3]+0.5; y++) {
            for (int z = (int)bbox[4]; z <= (int)bbox[5]+0.5; z++) {
              int crossedCounter = 0; // How many triangles are crossed

              hemo::Array<plint, 3> latticeSite = {x, y, z};
              std::vector<hemo::Array<plint,3>> triangles_list;
              findCrossings(latticeSite,triangles_list);

              for (hemo::Array<plint, 3> triangle : triangles_list) {
                // Muller-trumbore intersection algorithm 
                const hemo::Array<double,3> & v0 = particles[cell[triangle[0]]].sv.position;
                const hemo::Array<double,3> & v1 = particles[cell[triangle[1]]].sv.position;
                const hemo::Array<double,3> & v2 = particles[cell[triangle[2]]].sv.position;

                crossedCounter += hemo::MollerTrumbore(v0, v1, v2, latticeSite);
              }

              // Count even-odd crossings
              if (crossedCounter%2) {
                int x_l = x-fluid->getLocation().x;
                int y_l = y-fluid->getLocation().y;
                int z_l = z-fluid->getLocation().z;
                innerNodes.push_back({x_l,y_l,z_l});
              }
            }
          }
        }
      }
  };
}

#endif
