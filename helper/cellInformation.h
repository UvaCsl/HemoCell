#ifndef CELLINFORMATION_H
#define CELLINFORMATION_H

#include "hemocell_internal.h"
#include "hemoCellFunctional.h"
#include "hemocell.h"

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
  
public:
  static map<int,CellInformation> info_per_cell;
  static void clear_list();
  static void calculate_all(HemoCell *);
  static void getCellVolume(HemoCell *);
  static void getCellArea(HemoCell *);
  static void getCellPosition(HemoCell *);
};

#endif
