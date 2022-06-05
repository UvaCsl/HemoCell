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

#include <vector>
#include <getopt.h>
#include <limits.h>

#include "packing.h"

void PrintHelp() {
  std::cerr <<
          "USAGE: packCells sX sY sZ [OPTIONAL ARGUMENTS ...]\n"
          "\n"
          "OPTIONAL ARGUMENTS:\n"
          "  --hematocrit <0-1.0>                 -h The hematocrit of the solution (human blood only!)\n"
          "  --plt_ratio <ratio>                     The ratio of PLT per RBC, default=0.07\n"
          "  --rbc <n>                               Number of Red Blood Cells\n"
          "  --plt <n>                               Number of Platelets\n"
          "  --rbc_m <n>                             Number of mouse Red Blood Cells\n"
          "  --plt_m <n>                             Number of mouse Wilde-Type platelets\n"
          "  --plt_mko <n>                           Number of mouse knockout Platelets\n"
          "  --wbc <n>                               Number of White Blood Cells\n"
          "  --vrbc <n>                              Number of Stage V gametocytes\n"
          "  --cell <name> <n> <e1, e2, diameter>    Custom Celltype described by ellipsoid\n"
          "  --noRotate                              Disallow rotation of ellipsoids\n"
          "  --scale <ratio>                      -s Scales the neighbourhood grid (only change this if you know what you are doing!)\n"
          "  --maxiter <n>                           Maximum number of iterations\n"
          "  --help                                  Print this"
          "\n"
          "OUTPUT:\n"
          "  <Cell>.pos for every celltype. The first line contains the number of cells.\n"
          "  The rest of the lines are the cells in \"Location<X Y Z> Rotation<X Y Z>\" format."
          "\n"
          "  Cells.pov is a tool for visualization in povray\n"
          "\n"
          "NOTE:\n"
          "  sX, sY and sZ are the domain size in [µm]\n"
          "  sX, sY, sZ and output are in [µm], not lattie units!\n"
          "  --hematocrit and --RBC are mutually exclusive\n"
          "  --hematocrit and --PLT are mutually exclusive\n"
          "  --PLT_ratio does not work without --hematocrit\n"
          ; 
          
}

const char* const short_opts = "h:s";

const option long_opts[] = {
            {"hematocrit", 1, nullptr, 0},
            {"plt_ratio",  1, nullptr, 10},
            {"rbc",        1, nullptr, 2},
            {"plt",        1, nullptr, 3},
            {"wbc",        1, nullptr, 4},
            {"vrbc",       1, nullptr, 11},
            {"rbc_m",      1, nullptr, 12},
            {"plt_m",      1, nullptr, 13},
            {"plt_mko",    1, nullptr, 14},
            {"cell",       1, nullptr, 5},
            {"noRotate",   0, nullptr, 6},
            {"scale",      1, nullptr, 7},
            {"maxiter",    1, nullptr, 8},
            {"help",       0, nullptr, 9},
            {NULL, 0, 0, 0}
};

// Human cell types
double rbcA = 8.4, rbcB = 4.4, rbcC = 8.4; // Inreased to fully cover biconcave shape of RBCs
double pltA = 2.4 , pltB = 1.05, pltC = 2.4;
double wbcA = 8.4, wbcB = 8.4, wbcC = 8.4;
double vrbcA= 3.5, vrbcB=6.0, vrbcC= 11.0; // Set according to the reconstructed stl dimensions

// Mouse cell types
double rbc_mA = 5.8, rbc_mB = 3.4, rbc_mC = 5.8;
double plt_mA = 1.84, plt_mB = 1.05, plt_mC = 1.84; 
double plt_mkoA = 1.71, plt_mkoB = 1.71, plt_mkoC = 1.71; 


int main(int argc, char *argv[]) 
{
  
  int val;  
  double scale = -1.0;
  double plt_ratio = 0.07;
  double sX,sY,sZ;
  int maxIter = INT_MAX-100; //-100 is compatibility with packing code
  float hematocrit = 0.0;
  bool hematocrit_set = false;
  bool RBC_PLT_set = false;
  bool doRotate = true;
  vector<CellType> cellTypes;
  
  
  while ((val = getopt_long_only(argc,argv,short_opts,long_opts,nullptr)) > -1) {
    if (val == 'h') { val = 0; }
    if (val == 's') { val = 7; }
    
    switch(val){
      case(0):
        if (RBC_PLT_set) {
          cerr << "--hematocrit is mutually exclusive with --plt and --rbc, exiting..." << endl;
          PrintHelp();
          return 1;
        }
        hematocrit_set = true;
        hematocrit = atof(optarg);
        break;
      case(10):
        plt_ratio = atof(optarg);
        break;
      case(2):
        if (hematocrit_set) {
          cerr << "--hematocrit is mutually exclusive with --plt and --rbc, exiting..." << endl;
          PrintHelp();
          return 1;
        }
        RBC_PLT_set = true;
        cellTypes.push_back({"RBC",rbcA,rbcB,rbcC,atoi(optarg)});
        break;
      case(3):
        if (hematocrit_set) {
          cerr << "--hematocrit is mutually exclusive with --plt and --rbc, exiting..." << endl;
          PrintHelp();
          return 1;
        }
        RBC_PLT_set = true;
        cellTypes.push_back({"PLT",pltA,pltB,pltC,atoi(optarg)});
        break;
      case(4):
        cellTypes.push_back({"WBC",wbcA,wbcB,wbcC,atoi(optarg)});
        break;
      case(11):
        cellTypes.push_back({"vRBC",vrbcA,vrbcB,vrbcC,atoi(optarg)});
        break;  
      case(12):
        cellTypes.push_back({"RBC_m",rbc_mA,rbc_mB,rbc_mC,atoi(optarg)});
        break;
      case(13):
        cellTypes.push_back({"PLT_m",plt_mA,plt_mB,plt_mC,atoi(optarg)});
        break;
      case(14):
        cellTypes.push_back({"PLT_mko",plt_mkoA,plt_mkoB,plt_mkoC,atoi(optarg)});
        break;      
      case(5):
        if (optind+3 > argc) {
          PrintHelp();
          return 1;
        }
        cellTypes.push_back({string(optarg),atof(argv[optind+1]),atof(argv[optind+2]),atof(argv[optind+3]),atoi(argv[optind])});
        optind = optind+4;
        break;
      case(6):
        doRotate = false;
        break;
      case(7):
        scale = atof(optarg);
        break;
      case(8):
        maxIter = atoi(optarg);
        break;
      case(9):
      case('?'):
      default:
        PrintHelp();
        return 1;
    }
  }
  if (optind + 3 > argc) {
    cout << "Insufficient arguments." << endl << endl;
    PrintHelp();
    return 1;
  }
  sX = atoi(argv[optind]);
  sY = atoi(argv[optind+1]);
  sZ = atoi(argv[optind+2]);
  
  if (!hematocrit_set && cellTypes.size() == 0) {
    cout << "You need to specify at least a celltype (--RBC, --PLT, --WBC, --CELL, etc.) or the hematocrit (--hematocrit), exiting..." << endl;
    return 1;
  }
  if (hematocrit_set) {
    double domainVol = sX*sY*sZ;
    double rbcVolNominal = 97.0; // This is set to match the model used in HemoCell! 
    //TODO subtract wbc number? why?
    int nRBC = (int)round(hematocrit * domainVol / rbcVolNominal);
    cellTypes.push_back({"RBC",rbcA,rbcB,rbcC,nRBC});
    //TODO relative to ALL cells? why?
    int nPLT = round(nRBC*plt_ratio);
    if (nPLT > 0 ) {
      cellTypes.push_back({"PLT",pltA,pltB,pltC,nPLT});
    }
  }
  
  string povFileName = "cells.pov";

  if(scale <= 0.0) { // Calculate optimal scaling automatically
    // Get the largest radius
    // double maxCellD = 0.0; 
    // for(auto cell : cellTypes) {
    //   double max_loc = max(max(cell.dx, cell.dy), cell.dz);
    //   if (max_loc > maxCellD)
    //     maxCellD = max_loc;
    // }

    // // Get the smallest domain side
    // double minDomainD = min(min(sX, sY),sZ);

    // scale = maxCellD/minDomainD; // -> scaled R > 1.0

    cout << "Setting the resolution according to the most numerous species. Overlap between smaller cell types might occur! (e.g., If you have many platelets, but they are not the most numerous, you will need an increased scaling factor to avoid overlap between them)." << endl;
    double maxCellD = 0.0;
    int maxCellCount = 0;
    for(auto cell : cellTypes) {
      if (maxCellCount < cell.number) {
        maxCellCount = cell.number;
        maxCellD = max(max(cell.dx, cell.dy), cell.dz);
      }
    }

    scale = (1./maxCellD)*1.2; // The 1.2 is an added "safety" range 
  }

  cout.setf (ios::fixed, ios::floatfield);
  cout << "Parameters:" << endl;
  cout << "  Domain Size [µm]: ( " << sX << " , " << sY << " , " << sZ << " )" << endl; 
  cout << "  Maximum Iterations : " << maxIter << endl;
  cout << "  Scale              : " << scale << endl;
  cout << "  Rotation           : " << doRotate << endl;
  if (hematocrit_set) {
    cout << "  Hematocrit    : " << hematocrit << endl;
    cout << "  PLT/RBC Ratio : " << plt_ratio << endl;
  }
  cout << "The following cell types are defined:" << endl;
  for (CellType & cell : cellTypes) {
    cout << "  " << cell.name << endl;
    cout << "    No   : " << cell.number << endl;
    cout << "    Sizes: (" << cell.dx << " , " << cell.dy << " , " << cell.dz << " )" << endl;
    cout << "    Volume: " << 4.0/3.0*PI*cell.dx*cell.dy*cell.dz/8.0 << " [µm3]." << endl;
  }
  
  Packing pack;
  
  pack.setRndRotation(doRotate); // This needs to be set first!
  
  pack.initBlood(sX, sY, sZ, maxIter, scale, cellTypes);

  pack.execute();

  pack.saveBloodCellPositions();
  
  //pack.savePov(povFileName.c_str(), sX, sY, sZ, wbcNumber);

  return 0;
}
