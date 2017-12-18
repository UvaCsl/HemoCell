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

#include "hemocell_internal.h"
#include "hemoCellFunctional.h"

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

struct CellInformation {
  hemo::Array<double,3> position;
  double volume;
  double area;
  double stretch;
  hemo::Array<double,6> bbox;
  pluint blockId;
  pluint cellType;
  bool centerLocal;
};

class CellInformationFunctionals {
  static HemoCell * hemocell;
  

  class CellVolume: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellVolume * clone() const;
  };
  class CellArea: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellArea * clone() const;
  };
  class CellPosition: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellPosition * clone() const;
  };
  class CellStretch: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellStretch * clone() const;
  };
  class CellBoundingBox: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellBoundingBox * clone() const;
  };
  class CellAtomicBlock: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellAtomicBlock * clone() const;
  };
  class CellType: public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
    CellType * clone() const;
  };
  
public:
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
  static pluint getTotalNumberOfCells(HemoCell *);
  static pluint getNumberOfCellsFromType(HemoCell *, string type);
};

#endif
