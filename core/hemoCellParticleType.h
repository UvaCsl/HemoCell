/* 
 * File:   hemoCellField.h
 * Author: vikko
 *
 * Created on 30 March 2017, 11:24
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

  HemoCellField(HemoCellFields& cellFields_, TriangularSurfaceMesh<double>& meshElement_);
  double getVolumeFraction();
  //ShellModel3D<double> * model;
  TriangularSurfaceMesh<double> & getMesh();
  std::string name;
  pluint ctype;
  int numVertex;
  unsigned int timescale = 1;
  unsigned int minimumDistanceFromSolid = 0;
  bool outputTriangles = false;
  bool outputLines = false;
  bool deleteIncomplete = true;
  HemoCellFields & cellFields;
  vector<int> desiredOutputVariables;
  vector<hemo::Array<plint,3>> triangle_list;
  TriangularSurfaceMesh<double> & meshElement;
  void(*kernelMethod)(BlockLattice3D<double,DESCRIPTOR> const&,HemoCellParticle*);
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleField3D();
  MultiBlockLattice3D<double,DESCRIPTOR> * getFluidField3D();
  int getNumberOfCells_Global();
  std::string getIdentifier();
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleArg();
  void setOutputVariables(const vector<int> &);
  CellMechanics * mechanics;
  void statistics();
  /* position is in micrometers, so we still have to convert it*/
  void addSingleCell(hemo::Array<double,3> position, plint cellId);
  hemo::Array<double,6> getOriginalBoundingBox();
  MeshMetrics<double> * meshmetric;
};


#endif /* HEMOCELLFIELD_H */

