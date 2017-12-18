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

void HemoCellParticleField::AddOutputMap() {
  outputFunctionMap[OUTPUT_POSITION] = &HemoCellParticleField::outputPositions;
  outputFunctionMap[OUTPUT_FORCE] = &HemoCellParticleField::outputForces;
  outputFunctionMap[OUTPUT_FORCE_VOLUME] = &HemoCellParticleField::outputForceVolume;
  outputFunctionMap[OUTPUT_FORCE_AREA] = &HemoCellParticleField::outputForceArea;
  outputFunctionMap[OUTPUT_FORCE_LINK] = &HemoCellParticleField::outputForceLink;
  outputFunctionMap[OUTPUT_FORCE_INNER_LINK] = &HemoCellParticleField::outputForceInnerLink;
  outputFunctionMap[OUTPUT_FORCE_BENDING] = &HemoCellParticleField::outputForceBending;
  outputFunctionMap[OUTPUT_FORCE_VISC] = &HemoCellParticleField::outputForceVisc;
  outputFunctionMap[OUTPUT_VERTEX_ID] = &HemoCellParticleField::outputVertexId;
  outputFunctionMap[OUTPUT_CELL_ID] = &HemoCellParticleField::outputCellId;

}

void HemoCellParticleField::passthroughpass(int type, Box3D domain, vector<vector<double>>& output, pluint ctype, std::string & name) {
  void (HemoCellParticleField::*badideapointer)(Box3D,vector<vector<double>>&, pluint, std::string&) = outputFunctionMap[type];
  (this->*badideapointer)(domain,output,ctype,name);
}

void HemoCellParticleField::outputPositions(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  deleteIncompleteCells(ctype);
  name = "Position";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      if (particles_per_cell.at(cellid)[i] == -1) { continue; }
      sparticle = &particles[particles_per_cell.at(cellid)[i]];

      vector<double> pbv;
      pbv.push_back(sparticle->position[0]);
      pbv.push_back(sparticle->position[1]);
      pbv.push_back(sparticle->position[2]);
      output.push_back(pbv); //TODO, memory copy

    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::dx;
      }
    }
  }
}

void HemoCellParticleField::outputForceBending(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Bending force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];

      vector<double> tf;
      tf.push_back((*sparticle->force_bending)[0]);
      tf.push_back((*sparticle->force_bending)[1]);
      tf.push_back((*sparticle->force_bending)[2]);
      output.push_back(tf);
    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::df;
      }
    }
  }
}

void HemoCellParticleField::outputForceArea(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Area force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];
 
      vector<double> tf;
      tf.push_back((*sparticle->force_area)[0]);
      tf.push_back((*sparticle->force_area)[1]);
      tf.push_back((*sparticle->force_area)[2]);
      output.push_back(tf);
    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::df;
      }
    }
  }
}

void HemoCellParticleField::outputForceLink(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Link force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];
 
      vector<double> tf;
      tf.push_back((*sparticle->force_link)[0]);
      tf.push_back((*sparticle->force_link)[1]);
      tf.push_back((*sparticle->force_link)[2]);
      output.push_back(tf);
    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::df;
      }
    }
  }
}

void HemoCellParticleField::outputForceInnerLink(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Inner link force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];
 
      vector<double> tf;
      tf.push_back((*sparticle->force_inner_link)[0]);
      tf.push_back((*sparticle->force_inner_link)[1]);
      tf.push_back((*sparticle->force_inner_link)[2]);
      output.push_back(tf);
    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::df;
      }
    }
  }
}

void HemoCellParticleField::outputForceVolume(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Volume force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];

      vector<double> tf;
      tf.push_back((*sparticle->force_volume)[0]);
      tf.push_back((*sparticle->force_volume)[1]);
      tf.push_back((*sparticle->force_volume)[2]);
      output.push_back(tf);
    }
  }
}

void HemoCellParticleField::outputForceVisc(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Viscous force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];

      vector<double> tf;
      tf.push_back((*sparticle->force_visc)[0]);
      tf.push_back((*sparticle->force_visc)[1]);
      tf.push_back((*sparticle->force_visc)[2]);
      output.push_back(tf);
    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::df;
      }
    }
  }
}

void HemoCellParticleField::outputForces(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Total force";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];
 
      vector<double> tf;
      tf.push_back(sparticle->force_total[0]);
      tf.push_back(sparticle->force_total[1]);
      tf.push_back(sparticle->force_total[2]);
      output.push_back(tf);
    }
  }
  if(cellFields->hemocell.outputInSiUnits) {
    for (vector<double> & tf : output) {
      for (double & n : tf) {
        n = n * param::df;
      }
    }
  }
}

void HemoCellParticleField::outputTriangles(Box3D domain, vector<vector<plint>>& output, vector<vector<double>> & positions, pluint ctype, std::string & name) {
  name = "Triangles";
  output.clear();
  int counter = 0;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < (*cellFields)[ctype]->triangle_list.size(); i++) {
      vector<plint> triangle = {(*cellFields)[ctype]->triangle_list[i][0] + counter,
                          (*cellFields)[ctype]->triangle_list[i][1] + counter,
                          (*cellFields)[ctype]->triangle_list[i][2] + counter};
    }
    counter += (*cellFields)[ctype]->numVertex;
  }
   
}

void HemoCellParticleField::outputLines(Box3D domain, vector<vector<plint>>& output, vector<vector<double>> & positions, plint ctype, std::string & name) {
  name = "Lines";
  output.clear();
  unsigned int counter = 0;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (cellFields->celltype_per_cell[cellid] != ctype) { continue; }
    const vector<int> & cell = particles_per_cell.at(cellid);
    unsigned int nElems = 0;
    for (unsigned int i = 0; i < cell.size() - 1 ; i++) {
      if (cell[i] == -1) { continue; }
      nElems++;
      if (cell[i+1] == -1 ) { i++; continue;}
      vector<plint> line = {i+counter,(i+1)+counter};
      output.push_back(line);
    }
    counter += nElems;
  }
   
}

void HemoCellParticleField::outputVertexId(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Vertex Id";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];
      vector<double> tf;
      tf.push_back((sparticle->vertexId));
      output.push_back(tf);
    }
  }
}

void HemoCellParticleField::outputCellId(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Cell Id";
  output.clear();
  HemoCellParticle * sparticle;
  const map<int,bool> & lpc = get_lpc();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (particles_per_cell.at(cellid)[0] == -1) { continue; }
    if (ctype != particles[particles_per_cell.at(cellid)[0]].celltype) continue;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++) {
      sparticle = &particles[particles_per_cell.at(cellid)[i]];
      vector<double> tf;
      tf.push_back((sparticle->cellId));
      output.push_back(tf);
    }
  }
}