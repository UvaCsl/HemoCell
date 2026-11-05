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
#ifndef READ_POSISIONS_OF_BLOOD_CELLS_H
#define READ_POSISIONS_OF_BLOOD_CELLS_H

#include "hemoCellFields.h"
#include "config.h"
#include "constantConversion.h"

#include "core/geometry3D.h"

#include <set>
#include "offLattice/triangularSurfaceMesh.hh"

namespace hemo {
/// One cell read from a .pos file.
struct InitialCellPlacement {
  hemo::Array<T, 3> position; ///< Cell centre in micrometres (shifted to the pre-inlet origin if one is active).
  hemo::Array<T, 3> angles;   ///< Rotation about X, Y, Z in radians, negated for the left-handed frame.
  plint cellId;               ///< Global cell id, sequential over all types in registration order.
};

/**
 * @brief Place the cells listed in `<type name>.pos` of every cell type into the particle fields.
 *
 * Rank 0 reads the files and broadcasts their text; every rank parses them once and places the
 * cells whose centre lies within `particleEnvelope` of each of its atomic blocks. Sets
 * `cellFields.number_of_cells`. Collective: every rank must call it.
 *
 * @param[in,out] cellFields All cell types; their particle fields receive the cells.
 * @param[in] dx Lattice spacing in metres (`param::dx`).
 * @param[in] cfg Configuration, `domain/particleEnvelope` is read.
 */
void readPositionsBloodCellField3D(HemoCellFields & cellFields, T dx, Config & cfg);

/// Sum of the cell counts in the first line of every `<type name>.pos` file (read on each rank).
int getTotalNumberOfCells(HemoCellFields & cellFields);

/// Cells centred in the domain that cross a periodic boundary while image
/// placement is off: those clear of the walls on a block (`candidates`) and
/// those that touch a wall on some block (`blocked`, they would not be placed
/// anyway).
struct CellsCrossingPeriodicBoundary {
  std::set<plint> candidates;
  std::set<plint> blocked;
};

class ReadPositionsBloodCellField3D : public plb::BoxProcessingFunctional3D {
public:
  /// @param cellFields_ All cell types (their meshes and particle fields).
  /// @param dx_ Lattice spacing in metres.
  /// @param cells_ Parsed cells per cell type; must outlive the functional (owned by the caller).
  /// @param particleEnvelope_ Cells are placed in a block if their centre is within this many lattice units of it.
  /// @param crossing_ If not null, receives the IDs of cells centred in the
  ///        domain that cross a periodic boundary (image placement off);
  ///        must outlive the functional.
  ReadPositionsBloodCellField3D(HemoCellFields & cellFields_, T dx_,
                                const std::vector<std::vector<InitialCellPlacement>> & cells_, plint particleEnvelope_,
                                CellsCrossingPeriodicBoundary * crossing_ = 0)
      : dx(dx_), cellFields(cellFields_), cells(cells_), particleEnvelope(particleEnvelope_), crossing(crossing_) {}
  /// Arguments: [0] fluid lattice, [1..] particle field of each cell type.
  virtual void processGenericBlocks(plb::Box3D domain, std::vector<plb::AtomicBlock3D *> fields);
  virtual ReadPositionsBloodCellField3D * clone() const;
  virtual void getTypeOfModification(std::vector<plb::modif::ModifT> & modified) const;
  void getModificationPattern(std::vector<bool> & isWritten) const;
  virtual plb::BlockDomain::DomainT appliesTo() const;
  T dx;
  HemoCellFields & cellFields;
  const std::vector<std::vector<InitialCellPlacement>> & cells;
  plint particleEnvelope;
  CellsCrossingPeriodicBoundary * crossing;
};
}
#endif
