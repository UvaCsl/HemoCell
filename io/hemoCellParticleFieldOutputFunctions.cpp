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
#include "hemoCellParticleField.h"
#include "hemocell.h"

#include <algorithm>

namespace hemo {

void HemoCellParticleField::AddOutputMap() {
  outputFunctionMap[OUTPUT_POSITION] = &HemoCellParticleField::outputPositions;
  outputFunctionMap[OUTPUT_VELOCITY] = &HemoCellParticleField::outputVelocities;
  outputFunctionMap[OUTPUT_FORCE] = &HemoCellParticleField::outputForces;
  outputFunctionMap[OUTPUT_FORCE_VOLUME] = &HemoCellParticleField::outputForceVolume;
  outputFunctionMap[OUTPUT_FORCE_AREA] = &HemoCellParticleField::outputForceArea;
  outputFunctionMap[OUTPUT_FORCE_LINK] = &HemoCellParticleField::outputForceLink;
  outputFunctionMap[OUTPUT_FORCE_INNER_LINK] = &HemoCellParticleField::outputForceInnerLink;
  outputFunctionMap[OUTPUT_FORCE_BENDING] = &HemoCellParticleField::outputForceBending;
  outputFunctionMap[OUTPUT_FORCE_VISC] = &HemoCellParticleField::outputForceVisc;
  outputFunctionMap[OUTPUT_VERTEX_ID] = &HemoCellParticleField::outputVertexId;
  outputFunctionMap[OUTPUT_CELL_ID] = &HemoCellParticleField::outputCellId;
  outputFunctionMap[OUTPUT_FORCE_REPULSION] = &HemoCellParticleField::outputForceRepulsion;
  outputFunctionMap[OUTPUT_RES_TIME] = &HemoCellParticleField::outputResTime;
}

void HemoCellParticleField::passthroughpass(int type, Box3D domain, vector<vector<T>> & output, pluint ctype,
                                            std::string & name) {
  if (outputFunctionMap.find(type) == outputFunctionMap.end()) { return; }
  void (HemoCellParticleField::*badideapointer)(Box3D, vector<vector<T>> &, pluint, std::string &) =
      outputFunctionMap[type];
  (this->*badideapointer)(domain, output, ctype, name);
}

void HemoCellParticleField::getOutputCells(pluint ctype, vector<int> & cellIds) {
  cellIds.clear();
  for (const auto & entry : get_particles_per_cell()) {
    const vector<int> & cell = entry.second;
    if (cell.empty() || cell[0] == -1) { continue; }
    const HemoCellParticle & vertex0 = particles[cell[0]];
    if (vertex0.sv.celltype != ctype) { continue; }
    // Exactly one block holds vertex 0 in its bulk; envelope copies, including
    // the shifted images of cells across a periodic boundary, are skipped.
    if (!isContainedABS(vertex0.sv.position, localDomain)) { continue; }
    if (std::find(cell.begin(), cell.end(), -1) != cell.end()) { continue; }
    cellIds.push_back(entry.first);
  }
}

namespace {
// One output row per vertex of the cells this block writes (getOutputCells()),
// cell by cell in vertex ID order; `row` returns the values of one vertex.
template <class RowFunction>
void outputPerVertex(HemoCellParticleField & field, pluint ctype, vector<vector<T>> & output, RowFunction row) {
  output.clear();
  vector<int> cellIds;
  field.getOutputCells(ctype, cellIds);
  const map<int, vector<int>> & particles_per_cell = field.get_particles_per_cell();
  for (int cellId : cellIds) {
    for (int index : particles_per_cell.at(cellId)) { output.push_back(row(field.particles[index])); }
  }
}

void scaleOutput(vector<vector<T>> & output, T factor) {
  for (vector<T> & row : output) {
    for (T & value : row) { value = value * factor; }
  }
}

vector<T> components(const hemo::Array<T, 3> & v) {
  return vector<T>{v[0], v[1], v[2]};
}
} // namespace

void HemoCellParticleField::outputPositions(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                            std::string & name) {
  name = "Position";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(p.sv.position); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::dx); }
}

void HemoCellParticleField::outputVelocities(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                             std::string & name) {
  name = "Velocity";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(p.sv.v); });
  if (cellFields->hemocell.outputInSiUnits) {
    for (vector<T> & tf : output) {
      for (T & n : tf) { n = n * param::dx / param::dt; }
    }
  }
}

void HemoCellParticleField::outputForceBending(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                               std::string & name) {
  name = "Bending force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(*p.force_bending); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForceArea(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                            std::string & name) {
  name = "Area force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(*p.force_area); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForceLink(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                            std::string & name) {
  name = "Link force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(*p.force_link); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForceInnerLink(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                                 std::string & name) {
  name = "Inner link force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(*p.force_inner_link); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForceVolume(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                              std::string & name) {
  name = "Volume force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(*p.force_volume); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForceVisc(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                            std::string & name) {
  name = "Viscous force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(*p.force_visc); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForceRepulsion(Box3D domain, vector<vector<T>> & output, pluint ctype,
                                                 std::string & name) {
  name = "Repulsion force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(p.sv.force_repulsion); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputForces(Box3D domain, vector<vector<T>> & output, pluint ctype, std::string & name) {
  name = "Total force";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return components(p.force_total); });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::df); }
}

void HemoCellParticleField::outputTriangles(Box3D domain, vector<vector<plint>> & output, pluint ctype,
                                            std::string & name) {
  name = "Triangles";
  output.clear();
  vector<int> cellIds;
  getOutputCells(ctype, cellIds);
  const HemoCellField & cellField = *(*cellFields)[ctype];
  plint counter = 0;
  for (pluint c = 0; c < cellIds.size(); c++) {
    for (const hemo::Array<plint, 3> & triangle : cellField.triangle_list) {
      output.push_back({triangle[0] + counter, triangle[1] + counter, triangle[2] + counter});
    }
    counter += cellField.numVertex;
  }
}

void HemoCellParticleField::outputInnerLinks(Box3D domain, vector<vector<plint>> & output, pluint ctype,
                                             std::string & name) {
  name = "InnerLinks";
  output.clear();
  vector<int> cellIds;
  getOutputCells(ctype, cellIds);
  const HemoCellField & cellField = *(*cellFields)[ctype];
  const auto & innerEdges = cellField.mechanics->cellConstants.inner_edge_list;
  plint counter = 0;
  for (pluint c = 0; c < cellIds.size(); c++) {
    for (pluint i = 0; i < innerEdges.size(); i++) {
      output.push_back({innerEdges[i][0] + counter, innerEdges[i][1] + counter});
    }
    counter += cellField.numVertex;
  }
}

void HemoCellParticleField::outputVertexId(Box3D domain, vector<vector<T>> & output, pluint ctype, std::string & name) {
  name = "Vertex Id";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return vector<T>{T(p.sv.vertexId)}; });
}

void HemoCellParticleField::outputCellId(Box3D domain, vector<vector<T>> & output, pluint ctype, std::string & name) {
  name = "Cell Id";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return vector<T>{T(p.sv.cellId)}; });
}

void HemoCellParticleField::outputResTime(Box3D domain, vector<vector<T>> & output, pluint ctype, std::string & name) {
  name = "Res Time";
  outputPerVertex(*this, ctype, output, [](const HemoCellParticle & p) { return vector<T>{T(p.sv.restime)}; });
  if (cellFields->hemocell.outputInSiUnits) { scaleOutput(output, param::dt); }
}

}
