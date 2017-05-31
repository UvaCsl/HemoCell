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
  
  std::string fileName = global::directories().getOutputDir() + "/csv/" + zeroPadNumber(hemocell->iter) + '/'  + createFileName("CellInfo.",global::mpi().getRank(),3) + ".csv";
  ofstream csvFile;
  csvFile.open(fileName.c_str(), ofstream::trunc);
  
  csvFile << "X,Y,Z,area,volume,atomic_block,cellId" << endl;
  
  for (const auto & pair : CellInformationFunctionals::info_per_cell) {
    const pluint cid = pair.first;
    const CellInformation & cinfo = pair.second;
    
    csvFile << cinfo.position[0] << "," << cinfo.position[1] << "," << cinfo.position[2] << ",";
    csvFile << cinfo.area << "," << cinfo.volume << "," << cinfo.blockId << "," << cid << endl;
  }
  csvFile.close();
  
  CellInformationFunctionals::clear_list();
}