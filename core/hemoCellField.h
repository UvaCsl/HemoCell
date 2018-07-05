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
#ifndef HEMOCELLFIELD_H
#define HEMOCELLFIELD_H

class HemoCellField;
#include "hemocell_internal.h"
#include "cellMechanics.h"
#include "meshMetrics.h"
#include "hemoCellFields.h"
#include "readPositionsBloodCells.h"

/*contains information about one particular cellfield, structlike*/
class HemoCellField{
  static vector<int> default_output;
  public:

  HemoCellField(HemoCellFields& cellFields_, TriangularSurfaceMesh<T>& meshElement_, string & name_, unsigned int ctype_);
  T getVolumeFraction();
  //ShellModel3D<T> * model;
  TriangularSurfaceMesh<T> & getMesh();
  std::string name;
  HemoCellFields & cellFields;
  vector<int> desiredOutputVariables;
  TriangularSurfaceMesh<T> & meshElement;
  pluint ctype;
  
  int numVertex = 0;
  T volume = 0;
  T volumeFractionOfLspPerNode = 0;
  T restingCellVolume = 0;
  unsigned int timescale = 1;
  unsigned int minimumDistanceFromSolid = 0;
  bool outputTriangles = false;
  bool deleteIncomplete = true;
  vector<hemo::Array<plint,3>> triangle_list;
  void(*kernelMethod)(BlockLattice3D<T,DESCRIPTOR> &,HemoCellParticle*);
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleField3D();
  MultiBlockLattice3D<T,DESCRIPTOR> * getFluidField3D();
  int getNumberOfCells_Global();
  std::string getIdentifier();
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleArg();
  void setOutputVariables(const vector<int> &);
  CellMechanics * mechanics;
  void statistics();
  /* position is in micrometers, so we still have to convert it*/
  void addSingleCell(hemo::Array<T,3> position, plint cellId);
  hemo::Array<T,6> getOriginalBoundingBox();
  MeshMetrics<T> * meshmetric;
};


#endif /* HEMOCELLFIELD_H */

