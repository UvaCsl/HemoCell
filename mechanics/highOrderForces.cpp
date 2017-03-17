#include "cellMechanics.h"
#include "commonCellConstants.cpp"

class HighOrderForces :public CellMechanics {

  //Variables
  const CommonCellConstants cellConstants;
  HemoCellField & cellField;
  const double k_volume;


  //Constructor
  HighOrderForces(HemoCellField & cellField_, double k_volume_) : cellConstants(CommonCellConstants::CommonCellConstantsConstructor(cellField_)),
                                                                  cellField(cellField_), k_volume(k_volume_)
  { 

  };

  virtual void ParticleMechanics(map<int,vector<SurfaceParticle3D *>> particles_per_cell, map<int,bool> lpc) {
    for (const auto & pair : lpc) { //For all cells with at least one lsp in the local domain.
      const int & cid = pair.first;
      vector<SurfaceParticle3D*> & cell = particles_per_cell[cid];

      //Calculate Cell Values that need all particles (but do it most efficient
      //tailered to this class)
      double volume = 0.0;
      double total_area = 0.0;
      for (const Array<plint,3> & triangle : cellConstants.triangle_list) {
        const Array<double,3> & v0 = cell[triangle[0]]->getPosition();
        const Array<double,3> & v1 = cell[triangle[1]]->getPosition();
        const Array<double,3> & v2 = cell[triangle[2]]->getPosition();
        
        //Volume
        const double v210 = v2[0]*v1[1]*v0[2];
        const double v120 = v1[0]*v2[1]*v0[2];
        const double v201 = v2[0]*v0[1]*v1[2];
        const double v021 = v0[0]*v2[1]*v1[2];
        const double v102 = v1[0]*v0[1]*v2[2];
        const double v012 = v0[0]*v1[1]*v2[2];
        volume += (1.0/6.0)*(-v210+v120+v201-v021-v102+v012);
        
        //Area
        //With herons formula, rewritten for speed
        const double l1 = (v0[0]-v1[0])*(v0[0]-v1[0]) +
                          (v0[1]-v1[1])*(v0[1]-v1[1]) +
                          (v0[2]-v1[2])*(v0[2]-v1[2]);
        const double l2 = (v2[0]-v1[0])*(v2[0]-v1[0]) +
                          (v2[1]-v1[1])*(v2[1]-v1[1]) +
                          (v2[2]-v1[2])*(v2[2]-v1[2]);
        const double l3 = (v0[0]-v2[0])*(v0[0]-v2[0]) +
                          (v0[1]-v2[1])*(v0[1]-v2[1]) +
                          (v0[2]-v2[2])*(v0[2]-v2[2]);
        const double area = (2*l1*l2 + 2*l2*l3 + 2*l1*l3 - l1*l1 - l2*l2 - l3*l3)/16.0


      }



      for (SurfaceParticle3D * particle : particles_per_cell[cid]) {
        //Volume

      }
    } 
  };



};
