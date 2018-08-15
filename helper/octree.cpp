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
#include "octree.h"
#include "logfile.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "core/geometry3D.h"
#include "mollerTrumbore.h"
#include "atomicBlock/blockLattice3D.hh"


namespace hemo {

using namespace std;


OctreeStructCell::OctreeStructCell(plint divis, plint l, unsigned int lim, hemo::Array<double, 6> bbox,
			vector<hemo::Array<plint,3>> triangle_list_,
			vector<HemoCellParticle> & part, const vector<int>  cell) {
  bBox = bbox;

  sharedConstructor(divis,l,lim,triangle_list_,part,cell);
}

OctreeStructCell::OctreeStructCell(plint divis, plint l, unsigned int lim,
			vector<hemo::Array<plint,3>> triangle_list_,
			vector<HemoCellParticle> & particles, const vector<int>  cell) {
  //The same, but construct bounding box first
  hemo::Array<T,3> * position = &particles[0].sv.position;
  
  bBox[0] = bBox[1] = (*position)[0];
  bBox[2] = bBox[3] = (*position)[1];
  bBox[4] = bBox[5] = (*position)[2];

  for (const int pid : cell ) {

    position = &particles[pid].sv.position;

    bBox[0] = bBox[0] > (*position)[0] ? (*position)[0] : bBox[0];
    bBox[1] = bBox[1] < (*position)[0] ? (*position)[0] : bBox[1];
    bBox[2] = bBox[2] > (*position)[1] ? (*position)[1] : bBox[2];
    bBox[3] = bBox[3] < (*position)[1] ? (*position)[1] : bBox[3];
    bBox[4] = bBox[4] > (*position)[2] ? (*position)[2] : bBox[4];
    bBox[5] = bBox[5] < (*position)[2] ? (*position)[2] : bBox[5];
  }
  
  sharedConstructor(divis,l,lim,triangle_list_,particles,cell);
}

void OctreeStructCell::sharedConstructor(plint divis, plint l, unsigned int lim,
			vector<hemo::Array<plint,3>> triangle_list_,
			vector<HemoCellParticle> & part, const vector<int>  cell) {
  
  maxDivisions = divis;
  level = l;
  limit = lim;

  // Check if subdividing limit is reached?
  if (triangle_list_.size() < lim || level > maxDivisions) {
    finalNode = true;
    triangle_list=triangle_list_;
    return;
  }
  constructTree(part, cell,triangle_list_);
}

OctreeStructCell::~OctreeStructCell() {
  // Check if we need to go over all the children
  if (! finalNode) {
    for (int i = 0; i < 8; i++) {
      delete nodes[i];
    }
  }
}

int OctreeStructCell::returnTrianglesAmount() {
  if (finalNode) {
    return triangle_list.size(); // Will not go deeper
  }     
  int tempSize = triangle_list.size();

  for (int i = 0; i < 8; i++) {     
    int temp = nodes[i]->returnTrianglesAmount();
    tempSize += temp;
  }

  return tempSize;
}

void OctreeStructCell::constructTree(vector<HemoCellParticle> & part, const vector<int> cell,vector<hemo::Array<plint,3>> triangle_list_) {
  // Find the octants of the current bounding box.
  vector<hemo::Array<double, 6>> bBoxes;
  T xHalf = bBox[0] + (bBox[1] - bBox[0])/2;
  T yHalf = bBox[2] + (bBox[3] - bBox[2])/2;
  T zHalf = bBox[4] + (bBox[5] - bBox[4])/2;

  // Not gonna do complicated stuff..
  bBoxes.push_back({bBox[0], xHalf, bBox[2], yHalf, bBox[4], zHalf});
  bBoxes.push_back({xHalf, bBox[1], bBox[2], yHalf, bBox[4], zHalf});
  bBoxes.push_back({xHalf, bBox[1], yHalf, bBox[3], bBox[4], zHalf});
  bBoxes.push_back({bBox[0], xHalf, yHalf, bBox[3], bBox[4], zHalf});
  bBoxes.push_back({bBox[0], xHalf, bBox[2], yHalf, zHalf, bBox[5]});
  bBoxes.push_back({xHalf, bBox[1], bBox[2], yHalf, zHalf, bBox[5]});
  bBoxes.push_back({xHalf, bBox[1], yHalf, bBox[3], zHalf, bBox[5]});
  bBoxes.push_back({bBox[0], xHalf, yHalf, bBox[3], zHalf, bBox[5]}); 
  // Create vectors to store corresponding triangles
  vector<vector<hemo::Array<plint,3>>> list_triangle_list;

  // Fill triangle lists
  for (unsigned int i = 0; i < bBoxes.size(); i++) {
    list_triangle_list.push_back({});
  }
  
  for (hemo::Array<plint,3> & triangle : triangle_list_) {  
    hemo::Array<double,3> & v0 = part[cell[triangle[0]]].sv.position;
    hemo::Array<double,3> & v1 = part[cell[triangle[1]]].sv.position;
    hemo::Array<double,3> & v2 = part[cell[triangle[2]]].sv.position;


    bool broken = false;
    for (unsigned int i = 0; i < bBoxes.size() ; i++) { 
      if ( ( (bBoxes[i][0] <= v0[0] && bBoxes[i][1] >= v0[0]) && 
             (bBoxes[i][2] <= v0[1] && bBoxes[i][3] >= v0[1]) &&
             (bBoxes[i][4] <= v0[2] && bBoxes[i][5] >= v0[2]) ) &&
           ( (bBoxes[i][0] <= v1[0] && bBoxes[i][1] >= v1[0]) && 
             (bBoxes[i][2] <= v1[1] && bBoxes[i][3] >= v1[1]) &&
             (bBoxes[i][4] <= v1[2] && bBoxes[i][5] >= v1[2]) ) &&
           ( (bBoxes[i][0] <= v2[0] && bBoxes[i][1] >= v2[0]) && 
             (bBoxes[i][2] <= v2[1] && bBoxes[i][3] >= v2[1]) &&
             (bBoxes[i][4] <= v2[2] && bBoxes[i][5] >= v2[2]) ) ) {
         list_triangle_list[i].push_back(triangle);
         broken = true;
         break;
      }
    }
    if (!broken) {
      triangle_list.push_back(triangle);
    }  
  }
 
  for (unsigned int i = 0; i < 8 ; i++) {
    nodes[i] = new OctreeStructCell(maxDivisions, level+1, limit, bBoxes[i],
                                    list_triangle_list[i], part, cell);                         
  }
}

void OctreeStructCell::findCrossings(hemo::Array<plint, 3> latticeSite, std::vector<hemo::Array<plint,3>>& output) {

  // The ray passes inside the bboxes, xdir is open, abuse that
  if ( bBox[0] <= latticeSite[1] &&
       bBox[2] <= latticeSite[1] && bBox[3] >= latticeSite[1] &&
       bBox[4] <= latticeSite[2] && bBox[5] >= latticeSite[2]) {
    output.insert(output.end(),triangle_list.begin(),triangle_list.end());
  } else {
return;
  }

  if (finalNode) {
    return;
  }
  
  for (int i = 0; i < 8; i++) {  
    nodes[i]->findCrossings(latticeSite,output);
  }
  return;
}

void OctreeStructCell::findInnerNodes(plb::BlockLattice3D<T,DESCRIPTOR> * fluid, vector<HemoCellParticle> & particles, const vector<int> & cell, vector<plb::Cell<T,DESCRIPTOR>*> & innerNodes) {
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
        vector<hemo::Array<plint,3>> triangles_list;
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

void OctreeStructCell::findInnerNodes(plb::BlockLattice3D<T,DESCRIPTOR> * fluid, vector<HemoCellParticle> & particles, const vector<int> & cell, vector<Array<plint,3>> & innerNodes) {
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
        vector<hemo::Array<plint,3>> triangles_list;
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

}
