#ifndef CELLINFORMATION_H
#define CELLINFORMATION_H

#include "hemocell_internal.h"
#include "hemoCellFunctional.h"
#include "hemocell.h"
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
  Array<double,3> position;
  double volume;
  double area;
  double stretch;
  Array<double,6> bbox;
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
  static void getCellVolume(HemoCell *);
  static void getCellArea(HemoCell *);
  static void getCellPosition(HemoCell *);
  static void getCellStretch(HemoCell *);
  static void getCellBoundingBox(HemoCell *);
  static void getCellAtomicBlock(HemoCell *);
  static void getCellType(HemoCell *);
  static pluint getTotalNumberOfCells(HemoCell *);
};

#endif
