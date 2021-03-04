/*
 * Name: Vuong Ngo
 * Student ID: 10800301
 * Description:
 * This code is based on the pipeflow with preinlet case.
 * Some minor changes have been made in order to implement the curved vessel geometry
 * in the simulation.
 * The positions of the inlet and outlet slice have been altered to match the curved
 * vessel simulation domain.
 */

#include "hemocell.h"
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "preInlet.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;


/*Main program*/
int main (int argc, char * argv[]) {

  /*Check for config.xml*/
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }


  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  /*Define parameters for lattice Boltzmann simulation, read from config.xml */
  hlog << "(Stl preinlet) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl;
  MultiScalarField3D<int> *flagMatrix = 0;
  VoxelizedDomain3D<double> * voxelizedDomain = 0;
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),
            (*cfg)["domain"]["fluidEnvelope"].read<int>(),
            (*cfg)["domain"]["refDirN"].read<int>(),
            (*cfg)["domain"]["refDir"].read<int>(),
            voxelizedDomain, flagMatrix,
            (*cfg)["domain"]["blockSize"].read<int>());

  hlog << "(Stl preinlet) (Parameters) setting lbm parameters" << endl;

  param::lbm_base_parameters((*cfg));
  param::printParameters();

  hlog << "(PreInlets) creating preInlet" << endl;
  hemocell.preInlet = new hemo::PreInlet(&hemocell,flagMatrix);

  Box3D slice = flagMatrix->getBoundingBox();

/* Some prints to check the values of the inlet slice

  hlog << "SLICE x0 = " << slice.x0 << endl;
  hlog << "SLICE x1 = " << slice.x1 << endl;
  hlog << "SLICE y0 = " << slice.y0 << endl;
  hlog << "SLICE y1 = " << slice.y1 << endl;
  hlog << "SLICE z0 = " << slice.z0 << endl;
  hlog << "SLICE z1 = " << slice.z1 << endl;

*/

  // You can put the inlet inside the domain, or any other side if needed
  //slice.x0 = slice.x1 - 10;

  // Now we just put it at the beginning of the domain
  slice.x0 = slice.x1;

  /*Slightly increase the height of the inlet to hopefully let the simulation run.*/
  slice.y1 = 42;

  // Direction:: -> define the inflow direction (preInlet is on the X negative side)
  hemocell.preInlet->preInletFromSlice(Direction::Xpos,slice);

  hlog << "(Stl preinlet) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.initializeLattice(voxelizedDomain->getMultiBlockManagement());

    if (!hemocell.partOfpreInlet) {
      hemocell.lattice->periodicity().toggleAll(false);
    }

    // Setting Preinlet creation
    hemocell.preInlet->initializePreInlet();

    hlog << "(Stl preinlet) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl;
    boundaryFromFlagMatrix(hemocell.lattice,flagMatrix,hemocell.partOfpreInlet);

    hemocell.preInlet->createBoundary();

    /*Disable statistics so it runs faster*/
    hemocell.lattice->toggleInternalStatistics(false);

    /*Equilibrate everything*/
    hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

    // Driving Force
  hemocell.preInlet->calculateDrivingForce();

  hemocell.lattice->initialize();

  // Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

    hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

    //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<double>(), (*cfg)["domain"]["RepCutoff"].read<double>());
    //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

    vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
    hemocell.setOutputs("RBC_HO", outputs);
    hemocell.setOutputs("PLT", outputs);

    outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE,OUTPUT_BOUNDARY};
    hemocell.setFluidOutputs(outputs);

    /* Printing out the values of the outlet slice to see

    Box3D test = hemocell.lattice->getBoundingBox();

    hlog << "TEST x0 = " << test.x0 << endl;
    hlog << "TEST x1 = " << test.x1 << endl;
    hlog << "TEST y0 = " << test.y0 << endl;
    hlog << "TEST y1 = " << test.y1 << endl;
    hlog << "TEST z0 = " << test.z0 << endl;
    hlog << "TEST z1 = " << test.z1 << endl;

    */

    // For the main simulation domain we have to define outlets
    if (!hemocell.partOfpreInlet) {
      Box3D bb = hemocell.lattice->getBoundingBox();
      //Box3D outlet(bb.x0,bb.x0+2,bb.y0,bb.y1,bb.z0,bb.z1);

      /*Used the x, y, z values of the outlet of the curved vessel*/
      Box3D outlet(100, 102, 80, 122, 0, 42);

      //outlet.y0 = (plint) ((2 * outlet.y1) / 3);
      //outlet.y0 = 80;
      //outlet.y1 = 122;

      OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = new BoundaryConditionInstantiator3D
        < T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T,DESCRIPTOR> > ();
      boundary->addPressureBoundary0N(outlet,*hemocell.lattice,boundary::density);
      setBoundaryDensity(*hemocell.lattice,outlet, 1.0);
    }

    //loading the cellfield
    if (not cfg->checkpointed) {
      hemocell.loadParticles();
      hemocell.writeOutput();
    } else {
      hemocell.loadCheckPoint();
    }

    //Restructure atomic blocks on processors when possible
    //hemocell.doRestructure(false); // cause errors(?)

    if (hemocell.iter == 0) {
      pcout << "(Stl preinlet) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
      for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
          hemocell.lattice->collideAndStream();
      }
    }

    unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
    unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
    unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
    unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();



    pcout << "(Stl preinlet) Starting simulation..." << endl;

    while (hemocell.iter < tmax ) {
      //preinlet.update();
      hemocell.iterate();

      if (hemocell.partOfpreInlet) {
        //Set driving force as required after each iteration
        hemocell.preInlet->setDrivingForce();
      }

      hemocell.preInlet->applyPreInlet();

      // Load-balancing! Only enable if PARMETIS build is available
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
         // finfo = FluidInfo::calculateForceStatistics(&hemocell);
          //Set force as required after this function;
          //setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                    // DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                   //  hemo::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
         // pcout << "Fluid force, Minimum: " << finfo.min << " Maximum: " << finfo.max << " Average: " << finfo.avg << endl;
          //ParticleStatistics pinfo = ParticleInfo::calculateVelocityStatistics(&hemocell);
          //pcout << "Particle velocity, Minimum: " << pinfo.min << " Maximum: " << pinfo.max << " Average: " << pinfo.avg << endl;

          hemocell.writeOutput();
      }

      if (hemocell.iter % tcheckpoint == 0) {
          hemocell.saveCheckPoint();
      }
    }

  pcout << "(main) Simulation finished :) " << endl;

  return 0;

}
