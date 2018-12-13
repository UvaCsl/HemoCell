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
#include "writeCellInfoCSV.h"
#include "cellInfo.h"
#include "hemocell.h"

namespace hemo {

void writeCellInfo_CSV(HemoCell & hemocell) {
  global.statistics.getCurrent()["writeCellCSVInfo"].start();

  map<int,CellInformation> info_per_cell;
  CellInformationFunctionals::calculateCellInformation(&hemocell,info_per_cell);
  for (auto it = info_per_cell.cbegin(); it != info_per_cell.cend() ;) 
  {
    if (!it->second.centerLocal) {
      it = info_per_cell.erase(it);
    } else {
      ++it;
    }
  }
  
  HemoCellGatheringFunctional<CellInformation>::gather(info_per_cell);
 
  if (!global::mpi().getRank()) {
    vector<std::string> fileNames = vector<std::string>(hemocell.cellfields->size());
    vector<ofstream> csvFiles = vector<ofstream>(fileNames.size());
    for (unsigned int i = 0 ; i < fileNames.size(); i++ ) {
      fileNames[i] = global::directories().getOutputDir() + "/csv/" +  (*hemocell.cellfields)[i]->name + "." + zeroPadNumber(hemocell.iter) + ".csv";
      csvFiles[i].open(fileNames[i], ofstream::trunc);
      csvFiles[i] << "X,Y,Z,area,volume,atomic_block,cellId,baseCellId,velocity_x,velocity_y,velocity_z" << endl;
    }

    for (auto & pair : info_per_cell) {
      const plint cid = pair.first;
      CellInformation & cinfo = pair.second;

      if (hemocell.outputInSiUnits) {
        cinfo.position *= param::dx;
        cinfo.area *= param::dx*param::dx;
        cinfo.volume *= param::dx*param::dx*param::dx;
      }

      if (cinfo.centerLocal) {
        csvFiles[cinfo.cellType] << cinfo.position[0] << "," << cinfo.position[1] << "," << cinfo.position[2] << ",";
        csvFiles[cinfo.cellType] << cinfo.area << "," << cinfo.volume << "," << cinfo.blockId << "," << cid  << "," << cinfo.base_cell_id << ",";
        csvFiles[cinfo.cellType] << cinfo.velocity[0] << "," << cinfo.velocity[1] << "," << cinfo.velocity[2] << endl;
      }
    }

    for (ofstream & file : csvFiles) {
      file.close();
    }
  }
  global.statistics.getCurrent().stop();
}

}