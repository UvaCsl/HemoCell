#ifndef HEMOCELL_PARTICLE_FIELD_OUTPUT_FUNCTIONS_CPP
#define HEMOCELL_PARTICLE_FIELD_OUTPUT_FUNCTIONS_CPP

#include "hemoCellParticleField.h"

void HemoCellParticleField::AddOutputMap() {
  outputFunctionMap[OUTPUT_POSITION] = &HemoCellParticleField::outputPositions;
  outputFunctionMap[OUTPUT_FORCE] = &HemoCellParticleField::outputForces;
  outputFunctionMap[OUTPUT_FORCE_VOLUME] = &HemoCellParticleField::outputForceVolume;
  outputFunctionMap[OUTPUT_FORCE_AREA] = &HemoCellParticleField::outputForceArea;
  outputFunctionMap[OUTPUT_FORCE_INPLANE] = &HemoCellParticleField::outputForceInPlane;
  outputFunctionMap[OUTPUT_FORCE_BENDING] = &HemoCellParticleField::outputForceBending;
  
}

void HemoCellParticleField::passthroughpass(int type, Box3D domain, vector<vector<double>>& output, pluint ctype, std::string & name) {
  //Too Much c++ will give you cancer like this function
  void (HemoCellParticleField::*cancerpointer)(Box3D,vector<vector<double>>&, pluint, std::string&) = outputFunctionMap[type];
  (this->*cancerpointer)(domain,output,ctype,name);
}

void HemoCellParticleField::outputPositions(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  deleteIncompleteCells(ctype);
  name = "Position";
  output.clear();
  HemoCellParticle * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> pbv;
      pbv.push_back(sparticle->position[0]);
      pbv.push_back(sparticle->position[1]);
      pbv.push_back(sparticle->position[2]);
      output.push_back(pbv); //TODO, memory copy

    }
  }
}

void HemoCellParticleField::outputForceBending(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Bending Force";
  output.clear();
  HemoCellParticle * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_bending)[0]);
      tf.push_back((*sparticle->force_bending)[1]);
      tf.push_back((*sparticle->force_bending)[2]);
      output.push_back(tf);
    }
  }
}

void HemoCellParticleField::outputForceArea(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Area Force";
  output.clear();
  HemoCellParticle * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_area)[0]);
      tf.push_back((*sparticle->force_area)[1]);
      tf.push_back((*sparticle->force_area)[2]);
      output.push_back(tf);
    }
  }
}

void HemoCellParticleField::outputForceInPlane(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "In Plane Force";
  output.clear();
  HemoCellParticle * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_inplane)[0]);
      tf.push_back((*sparticle->force_inplane)[1]);
      tf.push_back((*sparticle->force_inplane)[2]);
      output.push_back(tf);
    }
  }
}
void HemoCellParticleField::outputForceVolume(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Volume Force";
  output.clear();
  HemoCellParticle * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back((*sparticle->force_volume)[0]);
      tf.push_back((*sparticle->force_volume)[1]);
      tf.push_back((*sparticle->force_volume)[2]);
      output.push_back(tf);
    }
  }
}

void HemoCellParticleField::outputForces(Box3D domain,vector<vector<double>>& output, pluint ctype, std::string & name) {
  name = "Total Force";
  output.clear();
  HemoCellParticle * sparticle;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;
    for (pluint i = 0; i < particles_per_cell[cellid].size(); i++) {
      sparticle = particles_per_cell[cellid][i];

      vector<double> tf;
      tf.push_back(sparticle->force_total[0]);
      tf.push_back(sparticle->force_total[1]);
      tf.push_back(sparticle->force_total[2]);
      output.push_back(tf);
    }
  }
}

void HemoCellParticleField::outputTriangles(Box3D domain, vector<vector<plint>>& output, vector<vector<double>> & positions, pluint ctype, std::string & name) {
  name = "Triangles";
  output.clear();
  int counter = 0;
  for ( const auto &lpc_it : lpc ) {
    int cellid = lpc_it.first;
    if (!particles_per_cell[cellid][0]) { continue; }
    if (ctype != particles_per_cell[cellid][0]->celltype) continue;

    for (pluint i = 0; i < (*cellFields)[ctype]->triangle_list.size(); i++) {
      vector<plint> triangle = {(*cellFields)[ctype]->triangle_list[i][0] + counter,
                          (*cellFields)[ctype]->triangle_list[i][1] + counter,
                          (*cellFields)[ctype]->triangle_list[i][2] + counter};

      //Do not add triangles over periodic boundaries
      bool toolarge = false;
      /*for (pluint x = 0; x < 3; x++) {
        for (pluint y = x +1; y < 3; y++) {
          if ((abs(positions[triangle[x]][0] - positions[triangle[y]][0]) > 2.0) ||
              (abs(positions[triangle[x]][1] - positions[triangle[y]][1]) > 2.0) ||
              (abs(positions[triangle[x]][2] - positions[triangle[y]][2]) > 2.0))
          {
            toolarge = true;
            break;
          }
        }
      }*/
      if (!toolarge) {
        output.push_back(triangle);
      }
    }
    counter += (*cellFields)[ctype]->numVertex;
  }
   
}

#endif