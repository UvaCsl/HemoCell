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
#ifndef HEMO_XDMF_OUTPUT_H
#define HEMO_XDMF_OUTPUT_H

#include "core/globalDefs.h"

#include <string>
#include <vector>

namespace hemo {

/// A dataset of a block file that is listed as an XDMF attribute.
struct XdmfAttribute {
  std::string name;                 ///< HDF5 dataset name, also the attribute name.
  int components;                   ///< Values per vertex (cells) or per lattice node (fluid).
  std::string numberType = "Float"; ///< XDMF NumberType of the dataset ("Float" or "Int").
  int precision = 4;                ///< Bytes per value of the dataset.
};

/**
 * @brief What one atomic block wrote to its HDF5 file in one output step.
 *
 * Filled by the HDF5 writers on the rank that wrote the file and passed to
 * writeXdmf(). Structured blocks (fluid, CEPAC) store one value per lattice
 * node of the block's bulk; the XDMF file describes them as the cells of a
 * regular grid, so neighbouring blocks tile the domain without overlap. Cell
 * blocks store a surface mesh: the vertices of the cells the block writes.
 */
struct XdmfBlockFile {
  std::string identifier;  ///< Output name: cell type, "Fluid" or "CEPAC"; "_PRE" appended in the pre-inlet.
  std::string file;        ///< HDF5 file, relative to the output directory.
  long blockId = 0;        ///< Atomic block ID; blocks are listed in ascending ID.
  bool structured = false; ///< Regular grid (fluid) or surface mesh (cells).

  long nodes[3] = {0, 0, 0};        ///< Structured: number of nodes in z, y, x (dataset dimensions).
  double origin[3] = {0., 0., 0.};  ///< Structured: lower corner of the first node's cell, z, y, x.
  double spacing[3] = {1., 1., 1.}; ///< Structured: node spacing, z, y, x.

  long vertices = 0;   ///< Cells: rows of the "Position" dataset (0: no cells in this block).
  long triangles = 0;  ///< Cells: rows of the "Triangles" dataset (0: not written).
  long innerLinks = 0; ///< Cells: rows of the "InnerLinks" dataset (0: not written).

  std::vector<XdmfAttribute> attributes; ///< Datasets listed as attributes.
};

/**
 * @brief Write the XDMF (version 3) description of one output step. Collective.
 *
 * Rank 0 gathers the block files of all ranks and writes, for every
 * identifier, into the output directory:
 *  - `<identifier>.<iteration>.xmf`: this step, a spatial collection of the
 *    blocks with its time;
 *  - `<identifier>.xmf`: all steps of the identifier found in the output
 *    directory (including those of earlier runs restarted from a checkpoint),
 *    as a temporal collection that includes the step files.
 * Cell blocks without vertices are omitted; cell types without a "Position"
 * dataset cannot be described and must not be passed. Both files are written
 * to a temporary name first and then renamed, so a reader never sees a
 * partial file. Paths in the files are relative to the output directory.
 *
 * @param[in] blocks The block files written by this rank in this step.
 * @param[in] iter Iteration of the step (names the step file).
 * @param[in] time Time of the step: seconds for SI output, iterations otherwise.
 */
void writeXdmf(const std::vector<XdmfBlockFile> & blocks, plb::plint iter, double time);

} // namespace hemo
#endif
