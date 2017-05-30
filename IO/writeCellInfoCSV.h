/* 
 * File:   writeCellInfoCSV.h
 * Author: vikko
 *
 * Created on 30 May 2017, 14:27
 */

#ifndef WRITECELLINFOCSV_H
#define WRITECELLINFOCSV_H

#include "hemocell_internal.h"
#include "hemocell.h"

void writeCellInfo_CSV(HemoCell *);

class WriteCellInfoCSV : public HemoCellFunctional {
  void processGenericBlocks(Box3D, std::vector<AtomicBlock3D*>);
  WriteCellInfoCSV * clone() const;
};

#endif /* WRITECELLINFOCSV_H */

