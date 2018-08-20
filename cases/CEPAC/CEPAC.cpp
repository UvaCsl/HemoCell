#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "wbcHighOrderModel.h"
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

  int nx, ny, nz;
  nx = ny = nz =(*cfg)["domain"]["refDirN"].read<int>() ;
//  ny = nz = (*cfg)["domain"]["refDirN"].read<int>() ; 
  hlog << "(unbounded) (Parameters) calculating flow parameters" << endl;
  param::lbm_pipe_parameters((*cfg),nx);
  param::printParameters();
  
  //TODO calculate either here or from diffusion config in Parameters
//  double Dp = 0; //0.25e-4;
//  double D  = Dp * param::dt / (param::dx * param::dx);
//  param::tau_CEPAC = ( 3. * D ) + 0.5; //param::tau; 
//  param::tau_CEPAC = param::tau;
  
  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx,ny,nz, (*cfg)["domain"]["fluidEnvelope"].read<int>()),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.lattice->periodicity().toggle(1,true);
  hemocell.latticeEquilibrium(1.,plb::Array<T, 3>(0.,0.,0.));

  //Driving Force
  double shear_rate = (*cfg)["parameters"]["shearRate"].read<double>(); //input shear rate s-1
  hlog << "shear_rate = " << shear_rate << endl;

  double velocity_max = (shear_rate*(nz*param::dx));
  hlog << "velocity_max = " << velocity_max << endl;

//  double velocity_max_lbm = velocity_max * (param::dx / param::dt);
  double velocity_max_lbm = velocity_max * (param::dt / param::dx);


  pcout << "velocity_max_lbm = " << velocity_max_lbm << endl;

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT",ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.enableSolidifyMechanics("PLT");

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC_HO", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);
  outputs = {OUTPUT_DENSITY};
  hemocell.setCEPACOutputs(outputs);

  //Boundary Conditions
  plb::initializeAtEquilibrium(*hemocell.cellfields->CEPACfield, (*hemocell.cellfields->CEPACfield).getBoundingBox(), 0.0, {0.0,0.0,0.0});

  Box3D topChannel( 0, nx-1, 0, ny-1, nz-1, nz-1);
  Box3D bottomChannel( 0, nx-1, 0, ny-1, 0, 0);
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
  OnLatticeAdvectionDiffusionBoundaryCondition3D<T,CEPAC_DESCRIPTOR>* ADboundaryCondition = createLocalAdvectionDiffusionBoundaryCondition3D<T, CEPAC_DESCRIPTOR>();

  boundaryCondition->setVelocityConditionOnBlockBoundaries (*hemocell.lattice, topChannel );
  setBoundaryVelocity(*hemocell.lattice, topChannel, plb::Array<T,3>(velocity_max_lbm,0,0));

  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );

//  ADboundaryCondition->addTemperatureBoundary2P(bottomChannel, *hemocell.cellfields->CEPACfield); 
//  setBoundaryDensity(*hemocell.cellfields->CEPACfield,bottomChannel,0.0);

//  ADboundaryCondition->addTemperatureBoundary2P(topChannel, *hemocell.cellfields->CEPACfield); 
//  setBoundaryDensity(*hemocell.cellfields->CEPACfield,topChannel,0.0);

  //3 by 3 source on bottom of channel
  Box3D CEPACsource(1, 4, (ny/2-2), (ny/2)-2, 4, 8); //((nx/2)-1,(nx/2)+1,(ny/2)-2,(ny/2)+2, 2, 2);
  ADboundaryCondition->addTemperatureBoundary2N(CEPACsource,*hemocell.cellfields->CEPACfield); 
  setBoundaryDensity(*hemocell.cellfields->CEPACfield,CEPACsource,0.05);
  
  hemocell.lattice->initialize();   
  hemocell.cellfields->CEPACfield->initialize();   

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  hlog << "(PipeFlow) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    hemocell.iterate();
    
    if (hemocell.iter % tmeas == 0) {
        hlog << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
        hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
        hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC_HO");
        hlog << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
        FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); T toMpS = param::dx / param::dt;
        hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
        ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); T topN = param::df * 1.0e12;
        hlog << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

        hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }
  hemo::global.statistics.outputStatistics();
  hlog << "(main) Simulation finished :) " << endl;

  return 0;
}
