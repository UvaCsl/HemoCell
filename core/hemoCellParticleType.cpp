/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   hemoCellField.cpp
 * Author: vikko
 * 
 * Created on 30 March 2017, 11:24
 */

#include "hemoCellParticleType.h"
#include "immersedBoundaryMethod.h"

//HemoCellField
HemoCellField::HemoCellField(HemoCellFields& cellFields_, TriangularSurfaceMesh<double>& meshElement_)
      :cellFields(cellFields_), desiredOutputVariables(default_output), meshElement(meshElement_) {
         numVertex = meshElement.getNumVertices();
         std::vector<int>::iterator it = std::find(desiredOutputVariables.begin(), desiredOutputVariables.end(),OUTPUT_TRIANGLES);
         if (it != desiredOutputVariables.end()) {
            desiredOutputVariables.erase(it);
            outputTriangles = true;
         }

        for (plint iTriangle = 0; iTriangle < meshElement.getNumTriangles(); iTriangle++) {
          triangle_list.push_back({meshElement.getVertexId(iTriangle,0),
                                   meshElement.getVertexId(iTriangle,1),
                                   meshElement.getVertexId(iTriangle,2) 
                                   });
        }
        meshmetric = new MeshMetrics<double>(meshElement_);

        kernelMethod = interpolationCoefficientsPhi2;

        
       }


void HemoCellField::setOutputVariables(const vector<int> & outputs) { desiredOutputVariables = outputs;
         std::vector<int>::iterator it = std::find(desiredOutputVariables.begin(), desiredOutputVariables.end(),OUTPUT_TRIANGLES);
         if (it != desiredOutputVariables.end()) {
            desiredOutputVariables.erase(it);
            outputTriangles = true;
         } else {
           outputTriangles = false;
         }
  }

/* position is in micrometers, so we still have to convert it*/
void HemoCellField::addSingleCell(Array<double,3> position, plint cellId) {
  pcerr << "(HemoCell) (ParticleType) addSingleCell not implemented, but might definitely be nice to have" << endl;
  exit(0);  
}

void HemoCellField::statistics() {
  pcout << "Cellfield  (+ material model) of " << name << std::endl;

  mechanics->statistics();

}

MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * HemoCellField::getParticleArg() { return cellFields.immersedParticles; }
  int HemoCellField::getNumberOfCells_Global() {return 0;}
  std::string HemoCellField::getIdentifier() {return name;}
vector<int> HemoCellField::default_output ({OUTPUT_POSITION});
 MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * HemoCellField::getParticleField3D() {return cellFields.immersedParticles;}
  MultiBlockLattice3D<double,DESCRIPTOR> * HemoCellField::getFluidField3D() {return cellFields.lattice;}
  TriangularSurfaceMesh<double> & HemoCellField::getMesh() { return meshElement;}
  //double HemoCellField::getVolumeFraction() { return hematocrit;}   // no need to keep track of the hematocrit
