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
#include "hemoCellField.h"
#include "immersedBoundaryMethod.h"
#include "config.h"

//HemoCellField
HemoCellField::HemoCellField(HemoCellFields& cellFields_, TriangularSurfaceMesh<double>& meshElement_, string & name_, unsigned int ctype_)
      :name(name_), cellFields(cellFields_), desiredOutputVariables(default_output), meshElement(meshElement_), ctype(ctype_) {
         numVertex = meshElement.getNumVertices();
         std::vector<int>::iterator it = std::find(desiredOutputVariables.begin(), desiredOutputVariables.end(),OUTPUT_TRIANGLES);
         if (it != desiredOutputVariables.end()) {
            desiredOutputVariables.erase(it);
            outputTriangles = true;
         }
         it = std::find(desiredOutputVariables.begin(), desiredOutputVariables.end(),OUTPUT_LINES);
         if (it != desiredOutputVariables.end()) {
            desiredOutputVariables.erase(it);
            outputLines = true;
         }


        for (plint iTriangle = 0; iTriangle < meshElement.getNumTriangles(); iTriangle++) {
          triangle_list.push_back({meshElement.getVertexId(iTriangle,0),
                                   meshElement.getVertexId(iTriangle,1),
                                   meshElement.getVertexId(iTriangle,2) 
                                   });
        }
        meshmetric = new MeshMetrics<double>(meshElement_);

        kernelMethod = interpolationCoefficientsPhi2;
        
        try {
          string materialXML = name + ".xml";
          hemo::Config materialCfg(materialXML.c_str());
          volume = materialCfg["MaterialModel"]["Volume"].read<double>();
          volumeFractionOfLspPerNode = (volume/numVertex)/pow(param::dx*1e6,3);
        } catch (std::invalid_argument & exeption) {
          pcout << "(HemoCell) (WARNING) (AddCellType) Volume of celltype not present, volume set to zero" << endl;
        }
}

void HemoCellField::setOutputVariables(const vector<int> & outputs) { desiredOutputVariables = outputs;
         std::vector<int>::iterator it = std::find(desiredOutputVariables.begin(), desiredOutputVariables.end(),OUTPUT_TRIANGLES);
         if (it != desiredOutputVariables.end()) {
            desiredOutputVariables.erase(it);
            outputTriangles = true;
         } else {
           outputTriangles = false;
         }
         it = std::find(desiredOutputVariables.begin(), desiredOutputVariables.end(),OUTPUT_LINES);
         if (it != desiredOutputVariables.end()) {
            desiredOutputVariables.erase(it);
            outputLines = true;
         } else {
           outputLines = false;
         }
  }

hemo::Array<double,6> HemoCellField::getOriginalBoundingBox() {
  hemo::Array<double,6> bb;
  bb[0] = bb[1] = meshElement.getVertex(0)[0];
  bb[2] = bb[3] = meshElement.getVertex(0)[1];
  bb[4] = bb[5] = meshElement.getVertex(0)[2];
  
  for (long int i = 0 ; i < meshElement.getNumVertices() ; i++) {
    const hemo::Array<double,3> vertex(meshElement.getVertex(i));
    if (bb[0] > vertex[0]) { bb[0] = vertex[0]; }
    if (bb[1] < vertex[0]) { bb[1] = vertex[0]; }
    if (bb[2] > vertex[1]) { bb[2] = vertex[1]; }
    if (bb[3] < vertex[1]) { bb[3] = vertex[1]; }
    if (bb[4] > vertex[2]) { bb[4] = vertex[2]; }
    if (bb[5] < vertex[2]) { bb[5] = vertex[2]; }    
  }
  return bb;
}

/* position is in micrometers, so we still have to convert it*/
void HemoCellField::addSingleCell(hemo::Array<double,3> position, plint cellId) {
  pcerr << "(HemoCell) (ParticleType) addSingleCell not implemented, but might definitely be nice to have" << endl;
  exit(0);  
}

void HemoCellField::statistics() {
  pcout << "Cellfield  (+ material model) of " << name << std::endl;
  pcout << "  Volume :" << volume << " µm³ VolumeFraction of lsp per fluid node: " << volumeFractionOfLspPerNode << " %" << endl;
  pcout << "  Nvertex: " << numVertex << endl;
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
