#include <vector>
#include <getopt.h>
#include <limits.h>

#include "packing.h"

void PrintHelp() {
  std::cerr <<
          "USAGE: packCells sX sY sZ [OPTIONAL ARGUMENTS ...]\n"
          "\n"
          "OPTIONAL ARGUMENTS:\n"
          "  --hematocrit <0-1.0>                 -h The hematocrit of the solution\n"
          "  --plt_ratio <ratio>                     The ratio of PLT per RBC, default=0.07\n"
          "  --rbc <n>                               Number of Red Blood Cells\n"
          "  --plt <n>                               Number of Platelets\n"
          "  --wbc <n>                               Number of White Blood Cells\n"
          "  --cell <name> <n> <e1, e2, diameter>    Custom Celltype described by ellipsoid\n"
          "  --allowRotate                        -r Allow for rotation of ellipsoids\n"
          "  --scale <ratio>                         Scales the simulated annealing and can affect performance, default=0.3\n"
          "  --maxiter <n>                           Maximum number of iterations\n"
          "  --help                                  Print this"
          "\n"
          "OUTPUT:\n"
          "  <Cell>.pos for every celltype. First line is the number of cells.\n"
          "  The rest of the lines is the cells in \"Location<X Y Z> Rotation<X Y Z>\" format."
          "\n"
          "  Cells.pov for visualization in, for example, povray\n"
          "\n"
          "NOTE:\n"
          "  sX, sY and sZ are the domain size\n"
          "  sX, sY, sZ and output are in micrometers[µm]\n"
          "  --hematocrit and --RBC are mutually exclusive\n"
          "  --hematocrit and --PLT are mutually exclusive\n"
          "  --PLT-ratio is an No-Op without --hematocrit\n"
          ; 
          
}

const char* const short_opts = "h:r";

const option long_opts[] = {
            {"hematocrit", 1, nullptr, 0},
            {"plt_ratio",  1, nullptr, 10},
            {"rbc",        1, nullptr, 2},
            {"plt",        1, nullptr, 3},
            {"wbc",        1, nullptr, 4},
            {"cell",       1, nullptr, 5},
            {"allowRotate",0, nullptr, 6},
            {"scale",      1, nullptr, 7},
            {"maxiter",    1, nullptr, 8},
            {"help",       0, nullptr, 9},
            {NULL, 0, 0, 0}
};

double rbcA = 8.4, rbcB = 4.4, rbcC = 8.4; // Inreased to cover biconcave shape of RBCs
double pltA = 2.4 , pltB = 1.05, pltC = 2.4;
double wbcA = 8.4, wbcB = 8.4, wbcC = 8.4;


int main(int argc, char *argv[]) 
{
  
  int val;  
  double scale = 0.3;
  double plt_ratio = 0.07;
  double sX,sY,sZ;
  int maxIter = INT_MAX-100; //-100 is compatibility with packing code
  float hematocrit = 0.0;
  bool hematocrit_set = false;
  bool RBC_PLT_set = false;
  bool doRotate = false;
  vector<CellType> cellTypes;
  
  
  while ((val = getopt_long_only(argc,argv,short_opts,long_opts,nullptr)) > -1) {
    if (val == 'h') { val = 0; }
    if (val == 'r') { val = 6; }
    
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
      case(5):
        if (optind+3 > argc) {
          PrintHelp();
          return 1;
        }
        cellTypes.push_back({string(optarg),atof(argv[optind+1]),atof(argv[optind+2]),atof(argv[optind+3]),atoi(argv[optind])});
        optind = optind+4;
        break;
      case(6):
        doRotate = true;
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
    cout << "You need to specify at least a celltype (--RBC, --PLT, --WBC, --CELL) or the hematocrit (--hematocrit), exiting..." << endl;
    return 1;
  }
  if (hematocrit_set) {
    double domainVol = sX*sY*sZ;
    double rbcVolNominal = 90.0; // This is set to match the model used in ficsion! 
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

  cout.setf (ios::fixed, ios::floatfield);
  cout << "Loaded parameters, we found:" << endl;
  cout << "  Domain Size (µm): ( " << sX << " , " << sY << " , " << sZ << " )" << endl; 
  cout << "  Maximum Iterations : " << maxIter << endl;
  cout << "  Scale              : " << scale << endl;
  cout << "  Rotation           : " << doRotate << endl;
  if (hematocrit_set) {
    cout << "  Hematocrit    : " << hematocrit << endl;
    cout << "  PLT/RBC Ratio : " << plt_ratio << endl;
  }
  cout << "We have found the following Cells:" << endl;
  for (CellType & cell : cellTypes) {
    cout << "  " << cell.name << endl;
    cout << "    No   : " << cell.number << endl;
    cout << "    Sizes: (" << cell.dx << " , " << cell.dy << " , " << cell.dz << " )" << endl;
  }
  
  Packing pack;
  
  pack.setRndRotation(doRotate); // This needs to be set first!
  
  pack.initBlood(sX, sY, sZ, maxIter, scale, cellTypes);

  pack.execute();

  pack.saveBloodCellPositions();
  
  //pack.savePov(povFileName.c_str(), sX, sY, sZ, wbcNumber);

  return 0;
}
