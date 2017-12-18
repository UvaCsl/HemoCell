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
#ifndef Hemo_CommonCellConstants_H
#define Hemo_CommonCellConstants_H

class CommonCellConstants;

#include "hemocell_internal.h"
#include "geometryUtils.h"
#include "config.h"

//Forward declaration:
class HemoCellField;

class CommonCellConstants {
  private:
  CommonCellConstants(HemoCellField & cellField_,
                      vector<hemo::Array<plint,3>> triangle_list_,
                      vector<hemo::Array<plint,2>> edge_list_,
                      vector<double> edge_length_eq_list_,
                      vector<double> edge_angle_eq_list_,
                      vector<double> surface_patch_center_dist_eq_list,
                      vector<hemo::Array<plint,2>> edge_bending_triangles_list_,
                      vector<hemo::Array<plint,2>> edge_bending_triangles_outer_points_,
                      vector<double> triangle_area_eq_list_,
                      vector<hemo::Array<plint,6>> vertex_vertexes_,
                      vector<hemo::Array<plint,6>> vertex_edges_,
                      vector<hemo::Array<signed int,6>> vertex_edges_sign_,
                      vector<unsigned int> vertex_n_vertexes_,
                      vector<hemo::Array<hemo::Array<plint,2>,6>> vertex_outer_edges_per_vertex_,
                      vector<hemo::Array<hemo::Array<signed int,2>,6>> vertex_outer_edges_per_vertex_sign_,
                      double volume_eq_, double area_mean_eq_,
                      double edge_mean_eq_, double angle_mean_eq_,
                      vector<hemo::Array<plint,2>> inner_edge_list_,
                      vector<double> inner_edge_length_eq_list_);
  public: 
  static CommonCellConstants CommonCellConstantsConstructor(HemoCellField &, Config & modelCfg_);


  HemoCellField & cellField;
  const vector<hemo::Array<plint,3>> triangle_list;
  const vector<hemo::Array<plint,2>> edge_list;

  const vector<double> edge_length_eq_list;
  const vector<double> edge_angle_eq_list;
  const vector<double> surface_patch_center_dist_eq_list;
  const vector<hemo::Array<plint,2>> edge_bending_triangles_list;
  const vector<hemo::Array<plint,2>> edge_bending_triangles_outer_points;
  const vector<double> triangle_area_eq_list;
  const vector<hemo::Array<plint,6>> vertex_vertexes;
  const vector<hemo::Array<plint,6>> vertex_edges;
  const vector<hemo::Array<signed int,6>> vertex_edges_sign;
  const vector<unsigned int> vertex_n_vertexes;
  const vector<hemo::Array<hemo::Array<plint,2>,6>> vertex_outer_edges_per_vertex;
  const vector<hemo::Array<hemo::Array<signed int,2>,6>> vertex_outer_edges_per_vertex_sign;

  const double volume_eq;
  const double area_mean_eq;
  const double edge_mean_eq;
  const double angle_mean_eq;
  const vector<hemo::Array<plint,2>> inner_edge_list;
  const vector<double> inner_edge_length_eq_list;

};
#endif
