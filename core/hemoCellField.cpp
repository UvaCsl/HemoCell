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
#include "logfile.h"
#include "hemoCellParticleField.h"
#include "readPositionsBloodCells.h"

#include <climits>

using namespace hemo;

//HemoCellField
HemoCellField::HemoCellField(HemoCellFields& cellFields_, TriangularSurfaceMesh<T>& meshElement_, string & name_, unsigned int ctype_)
      :name(name_), cellFields(cellFields_), desiredOutputVariables(default_output), meshElement(meshElement_), ctype(ctype_) {
  
  if (ctype_ > UCHAR_MAX) {
    hlog << "(HemoCell) (AddCellType) more celltypes than UCHAR_MAX (255) added, please convert celltype to int or add less celltypes" << endl;
    exit(1);
  }
  
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
        meshmetric = new MeshMetrics<T>(meshElement_);

        kernelMethod = interpolationCoefficientsPhi2;
        
        try {
          string materialXML = name + ".xml";
          hemo::Config materialCfg(materialXML.c_str());
          try{
            volume = materialCfg["MaterialModel"]["Volume"].read<T>();
            volumeFractionOfLspPerNode = (volume/numVertex)/pow(param::dx*1e6,3);
          } catch (std::invalid_argument & exeption) {
              hlog << "(HemoCell) (WARNING) (AddCellType) Volume of celltype " << name << " not present, volume set to zero" << endl;
          }
          try {
            interiorViscosityTau = materialCfg["MaterialModel"]["viscosityRatio"].read<T>()*(param::tau-0.5)+0.5;
          } catch (std::invalid_argument & e) {}
        } catch (std::invalid_argument & e) {}
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

hemo::Array<T,6> HemoCellField::getOriginalBoundingBox() {
  hemo::Array<T,6> bb;
  bb[0] = bb[1] = meshElement.getVertex(0)[0];
  bb[2] = bb[3] = meshElement.getVertex(0)[1];
  bb[4] = bb[5] = meshElement.getVertex(0)[2];
  
  for (long int i = 0 ; i < meshElement.getNumVertices() ; i++) {
    const hemo::Array<T,3> vertex(meshElement.getVertex(i));
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
void HemoCellField::addSingleCell(hemo::Array<T,3> position, plint cellId) {
  pcerr << "(HemoCell) (ParticleType) addSingleCell not implemented, but might definitely be nice to have" << endl;
  exit(1);  
}

void HemoCellField::statistics() {
  hlog << "Cellfield  (+ material model) of " << name << std::endl;
  hlog << "  Volume :" << volume << " µm³ VolumeFraction of lsp per fluid node: " << volumeFractionOfLspPerNode << " %" << endl;
  hlog << "  Nvertex: " << numVertex << endl;
  mechanics->statistics();

}

plb::MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * HemoCellField::getParticleArg() { return cellFields.immersedParticles; }
int HemoCellField::getNumberOfCells_Global() {return 0;}
std::string HemoCellField::getIdentifier() {return name;}
vector<int> HemoCellField::default_output ({OUTPUT_POSITION});
plb::MultiParticleField3D<HEMOCELL_PARTICLE_FIELD> * HemoCellField::getParticleField3D() {return cellFields.immersedParticles;}
plb::MultiBlockLattice3D<T,DESCRIPTOR> * HemoCellField::getFluidField3D() {return cellFields.lattice;}
plb::TriangularSurfaceMesh<T> & HemoCellField::getMesh() { return meshElement;}
