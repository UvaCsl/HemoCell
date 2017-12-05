#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "wbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "preInlet.h"
#include <fenv.h>

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  pcout << "(PipeFlow) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl; 
  MultiScalarField3D<int> *flagMatrix = 0; 
  VoxelizedDomain3D<double> * voxelizedDomain = 0; 
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),  
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),  
                       (*cfg)["domain"]["refDirN"].read<int>(),  
                       (*cfg)["domain"]["refDir"].read<int>(),  
                       voxelizedDomain, flagMatrix,  
                       (*cfg)["domain"]["blockSize"].read<int>()); 
     
  pcout << "(PipeFlow) (Parameters) calculating flow parameters" << endl;
  param::lbm_pipe_parameters((*cfg),(*cfg)["domain"]["refDirN"].read<int>());
  param::printParameters();
  
  pcout << "(PipeFlow) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));


  pcout << "(PipeFlow) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl; 
  defineDynamics(*hemocell.lattice, *flagMatrix, (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  //Driving Force
  pcout << "(PipeFlow) (Fluid) Setting up driving Force" << endl; 
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;


  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<double>(), (*cfg)["domain"]["RepCutoff"].read<double>());
  //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC_HO", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  
  Box3D outlet(hemocell.lattice->getNx()-1,hemocell.lattice->getNx()-1,0,hemocell.lattice->getNy(),0,hemocell.lattice->getNz());
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
  boundary->setPressureConditionOnBlockBoundaries(*hemocell.lattice,outlet);
  setBoundaryDensity(*hemocell.lattice,outlet, 1.0);


  
  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }
 
  //Restructure atomic blocks on processors when possible
  //hemocell.doRestructure(false); // cause errors
  
  if (hemocell.iter == 0) {
    pcout << "(PipeFlow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  //Setting Preinlet creation
  Box3D preinletBox(0,20,0,hemocell.lattice->getBoundingBox().y1,0,hemocell.lattice->getBoundingBox().z1);
  PreInlet preinlet(preinletBox,"preinlet",(*cfg)["ibm"]["stepMaterialEvery"].read<int>(),Direction::Xpos,hemocell, true);
  
  
  pcout << "(PipeFlow) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    preinlet.update();
    hemocell.iterate();
    

    // Only enable if PARMETIS build is available
    /*
     if (hemocell.iter % tbalance == 0) {
       if(hemocell.calculateFractionalLoadImbalance() > (*cfg)["parameters"]["maxFlin"].read<double>()) {
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
      ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); double topN = param::df * 1.0e12;
      pcout << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

      // Additional useful stats, if needed
      //finfo = FluidInfo::calculateForceStatistics(&hemocell);
      //Set force as required after this function;
      // setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
      //           DESCRIPTOR<T>::ExternalField::forceBeginsAt,
      //           hemo::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
      // pcout << "Fluid force, Minimum: " << finfo.min << " Maximum: " << finfo.max << " Average: " << finfo.avg << endl;
      // ParticleStatistics pinfo = ParticleInfo::calculateVelocityStatistics(&hemocell);
      // pcout << "Particle velocity, Minimum: " << pinfo.min << " Maximum: " << pinfo.max << " Average: " << pinfo.avg << endl;

      hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  pcout << "(main) Simulation finished :) " << endl;

  return 0;
}
