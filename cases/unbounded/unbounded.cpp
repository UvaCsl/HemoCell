#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include <fenv.h>

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  pcout << "(unbounded) (Parameters) calculating flow parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();
  
  int nx, ny, nz;
  nx = ny = nz = (*cfg)["domain"]["refDirN"].read<int>() ;
  //nx = 2 * ny ; 
  
  pcout << "(unbounded) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>()),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));

  //hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
  //nx+2, ny+2, nz+2, new BGKdynamics<double,DESCRIPTOR>(1.0/param::tau) );

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.lattice->periodicity().toggle(1,true);
  hemocell.lattice->periodicity().toggle(2,true);
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.0,0.0,0.0));

  //pcout << getMultiBlockInfo(hemocell.lattice) << endl;
  //pcout << hemocell.lattice->getLocalInfo() << endl;

  //Driving Force
  /*
  pcout << "(PipeFlow) (Fluid) Setting up driving Force" << endl; 
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
  DESCRIPTOR<T>::ExternalField::forceBeginsAt,                                                                                    
  plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
  */

  hemocell.lattice->initialize();

  //Adding all the cells
  hemocell.initializeCellfield();

  pcout << " done ? readPositionsBloodCels " << endl;

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC_HO", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  // Turn on periodicity in the all directions
  hemocell.setSystemPeriodicity(0, true);
  hemocell.setSystemPeriodicity(1, true);
  hemocell.setSystemPeriodicity(2, true);

  //todo add statistics here
  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    //pcout << " OK, here 1 " << endl;
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
    }
  //pcout << " OK, here 2 " << endl;
  //Restructure atomic blocks on processors when possible
  //hemocell.doRestructure();
  
  if (hemocell.iter == 0) {
    pcout << "(unbounded) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  //unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  //unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  pcout << "(unbounded) Starting simulation..." << endl;
  pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC_HO") << endl;
  pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
  int ncells = CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC_HO");
  pcout << " | RBC Volume ratio [x100%]: " << ncells * 77.0 * 100 / (nx * ny * nz) << endl;

  //CellInformationFunctionals::calculateCellBoundingBox(&hemocell);
  //plb::Box3D arr = hemocell.cellfields->immersedParticles->getBoundingBox();
  //plb::Box3D arr = hemocell.cellfields->getParticleField3D().getBoundingBox();
  hemo::Array<double,6> arr = hemocell.cellfields->cellFields[0]->getOriginalBoundingBox();
  plb::Box3D arr2 = hemocell.lattice->getBoundingBox();
  //CellInformationFunctionals::calculateCellBoundingBox(&hemocell);
  //hemo::Array<double,6> arr = CellInformationFunctionals::info_per_cell[0].bbox;
  //pcout <<"( "<< arr[0]  << " , " << arr[1] << " , "<< arr[2] << " )  and ( " << arr[3] << " , " << arr[4] << " , " 
  //<< " "<< arr[5] << " ) " << endl;

  pcout <<"lattice size : ( "<< arr2.x0  << " , " << arr2.y0 << " , "<< arr2.z0 << " )  and ( " << arr2.x1 << " , " << arr2.y1 << " , "  
  << arr2.z1 << " ) " << endl;

  pcout << " ========================================================== " << endl;
  plint domainVol = nx * ny * nz ;
  pcout << "(main) Number of fluid cells: " << domainVol << " (" << domainVol * (param::dx * param::dx * param::dx * 1e18) << " um^3)" << endl;

  //for (pluint iCell = 0; iCell < hemocell.cellfields->size(); ++iCell) {
  //plint nCells = cellFields[iCell]->getNumberOfCells_Global();

  // Get the dynamic volume from the radius increased by the particle-particle force cutoff
  //Array<T,2> xRange, yRange, zRange;
  //meshes[iCell]->computeBoundingBox (xRange, yRange, zRange);
  plb::Array<double, 3> cellBounds (arr[3]-arr[0], arr[4]-arr[1], arr[5]-arr[2]);
  //Array<T, 3> cellRatio;
  T dynVolume = 0;

  for(int iDir = 0; iDir < 3; iDir++)
    dynVolume *= 1.0 + (/*kernel * */ 0.5 * 1e-6 / param::dx) / cellBounds[iDir];

  pcout << "(main)   Volume ratio [x100%]: " << ( ncells * 77.2 * 100 ) / domainVol << endl;
  pcout << "(main)   nCells (global) = " << ncells << endl ;
  pcout << ", Cell volume = " <<  dynVolume * (param::dx * param::dx * param::dx * 1e18) << "  um^3 in LU " << dynVolume << endl;

  pcout<< " ========================================================== " << endl;


    //}

  while (hemocell.iter < tmax ) {
    hemocell.iterate();
    /*
    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
		      plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    */

    // Only enable if PARMETIS build is available
    /*
     if (hemocell.iter % tbalance == 0) {
       if(hemocell.calculateFractionalLoadImbalance() > 3) {
         hemocell.doLoadBalance();
         hemocell.doRestructure();
       }
     }
   */
    if (hemocell.iter % tmeas == 0) {
      pcout << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      pcout << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC_HO");
      pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
      pcout << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
      

      hemocell.writeOutput();
    }
    /*if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
      }*/
  }

  pcout << "(main) Simulation finished :) " << endl;

  return 0;
}
