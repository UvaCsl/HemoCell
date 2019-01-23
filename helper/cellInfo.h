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
#ifndef CELLINFORMATION_H
#define CELLINFORMATION_H

#include "hemocell.h"
#include "hemoCellFunctional.h"
#include "array.h"

/* THIS CLASS IS NOT THREAD SAFE!*/
/* Calculate and store Cell-Specific Information
 * 
 * call the static members for the information you need 
 * (eg. CellInformationFunctionals::getCellVolume(hemocell))
 * then the information will be available in  the map:
 *  CellInformationFunctionals::info_per_cell
 * Clean the old cell results afterwards with 
 *  CellInformationFunctionals::clear_list()
 *
 */
namespace hemo {
  
struct CellInformation {
  hemo::Array<T,3> position = {};
  hemo::Array<T,3> velocity = {};
  T volume = 0;
  T area = 0;
  T stretch = 0;
  hemo::Array<T,6> bbox = {};
  pluint blockId = UINTMAX_MAX;
  pluint cellType = UINTMAX_MAX;
  bool centerLocal = false;
  int base_cell_id = 0;
};

class CellInformationFunctionals {
  // Legacy interface, remove at some point
  static HemoCell * hemocell;
  
  class CellVolume: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellVolume * clone() const;
  };
  class CellArea: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellArea * clone() const;
  };
  class CellPosition: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellPosition * clone() const;
  };
  class CellStretch: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellStretch * clone() const;
  };
  class CellBoundingBox: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellBoundingBox * clone() const;
  };
  class CellAtomicBlock: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellAtomicBlock * clone() const;
  };
  class CellType: public HemoCellFunctional {
  void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    CellType * clone() const;
  };
  //end legacy interface
  
  class allCellInformation: public HemoCellFunctional {
    HemoCell * hemocell;
    map<int,CellInformation> & info_per_cell;
  public:
    allCellInformation(HemoCell * hemocell_, map<int,CellInformation> & info_per_cell_) :
    hemocell(hemocell_), info_per_cell(info_per_cell_) {}
  private:
    void processGenericBlocks(plb::Box3D, std::vector<plb::AtomicBlock3D*>);
    allCellInformation * clone() const;
  };
  
  
public:
  // Legacy interface, should become private at some point.
  static map<int,CellInformation> info_per_cell;
  static void clear_list();
  static void calculate_vol_pos_area(HemoCell *); /*This excludes Stretch, since it is compute intensive!*/
  static void calculateCellVolume(HemoCell *);
  static void calculateCellArea(HemoCell *);
  static void calculateCellPosition(HemoCell *);
  static void calculateCellStretch(HemoCell *);
  static void calculateCellBoundingBox(HemoCell *);
  static void calculateCellAtomicBlock(HemoCell *);
  static void calculateCellType(HemoCell *);
  static void calculateCellVelocity(HemoCell *);
  
  //Interface that should be used, less error prone at the cost of some extra computation and memory usage 
  static pluint getTotalNumberOfCells(HemoCell *);
  static pluint getNumberOfCellsFromType(HemoCell *, string type);
  /// Calculate all possible macroscopic information and return it within the reference.
  static void calculateCellInformation(HemoCell *, map<int, CellInformation> &);
  
};
}
#endif
