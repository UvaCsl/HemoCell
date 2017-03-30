/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   hemoCellField.h
 * Author: vikko
 *
 * Created on 30 March 2017, 11:24
 */

#ifndef HEMOCELLFIELD_H
#define HEMOCELLFIELD_H
class HemoCellField;
#include "hemocell.h"
#include "immersedBoundaryMethod3D.h"

/*contains information about one particular cellfield, structlike*/
class HemoCellField{
  static vector<int> default_output;
  public:

  HemoCellField(CellFields3D& cellFields_, Cell3D<double,DESCRIPTOR> cell3D_, TriangularSurfaceMesh<double>& meshElement_);
  double getVolumeFraction();
  double hematocrit;
  //ShellModel3D<double> * model;
  TriangularSurfaceMesh<double> & getMesh();
  std::string name;
  int ctype;
  int numVertex;
  bool outputTriangles = false;
  CellFields3D & cellFields;
  vector<int> desiredOutputVariables;
  Cell3D<double,DESCRIPTOR> & cell3D;
  vector<Array<plint,3>> triangle_list;
  TriangularSurfaceMesh<double> & meshElement;
  void(*kernelMethod)(BlockLattice3D<double,DESCRIPTOR> const&,SurfaceParticle3D*);
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleField3D();
  MultiBlockLattice3D<double,DESCRIPTOR> * getFluidField3D();
  int getNumberOfCells_Global();
  std::string getIdentifier();
  MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * getParticleArg();
  std::map<plint, Cell3D<double,DESCRIPTOR>* > getCellIdToCell3D();
  void setOutputVariables(const vector<int> &);
  CellMechanics * mechanics;
  void statistics();
  MeshMetrics<double> * meshmetric;
};


#endif /* HEMOCELLFIELD_H */

