/* 
 * File:   writeCellInfoCSV.cpp
 * Author: vikko
 * 
 * Created on 30 May 2017, 14:27
 */

#include "writeCellInfoCSV.h"
#include "cellInformation.h"

void writeCellInfo_CSV(HemoCell * hemocell) {
  
  CellInformationFunctionals::clear_list();
  CellInformationFunctionals::calculate_vol_pos_area(hemocell);
  CellInformationFunctionals::getCellAtomicBlock(hemocell);
  CellInformationFunctionals::getCellType(hemocell);

  
  vector<std::string> fileNames = vector<std::string>(hemocell->cellfields->size());
  vector<ofstream> csvFiles = vector<ofstream>(fileNames.size());
  for (unsigned int i = 0 ; i < fileNames.size(); i++ ) {
    fileNames[i] = global::directories().getOutputDir() + "/csv/" + zeroPadNumber(hemocell->iter) + '/'  + createFileName("CellInfo_" + (*hemocell->cellfields)[i]->name + ".",global::mpi().getRank(),3) + ".csv";
    csvFiles[i].open(fileNames[i], ofstream::trunc);
    csvFiles[i] << "X,Y,Z,area,volume,atomic_block,cellId" << endl;
  }
  
  for (const auto & pair : CellInformationFunctionals::info_per_cell) {
    const pluint cid = pair.first;
    const CellInformation & cinfo = pair.second;
    
    csvFiles[cinfo.cellType] << cinfo.position[0] << "," << cinfo.position[1] << "," << cinfo.position[2] << ",";
    csvFiles[cinfo.cellType] << cinfo.area << "," << cinfo.volume << "," << cinfo.blockId << "," << cid << endl;
  }
  
  for (ofstream & file : csvFiles) {
    file.close();
  }
    
  CellInformationFunctionals::clear_list();
}