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
#include "hemocell.h"
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

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  int nx, ny, nz;
  nx = ny = nz = (*cfg)["domain"]["refDirN"].read<int>() ;
  nx = nx*2;
  hlog << "(Parameters) setting lbm parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

  double shear_rate = 1800; //input shear rate s-1
  pcout << "shear_rate = " << shear_rate << endl;

  double velocity_max = (shear_rate*(nz/1e6))/4;
  pcout << "velocity_max = " << velocity_max << endl;

  double velocity_max_lbm = velocity_max * ( (*cfg)["domain"]["dt"].read<double>() /(*cfg)["domain"]["dx"].read<double>() );
  pcout << "velocity_max_lbm = " << velocity_max_lbm << endl;



  hlog << "(preinlet_shear) (Fluid) Initializing Palabos Fluid Field" << endl;
  MultiBlockManagement3D management = defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>());
  MultiBlockLattice3D<double,DESCRIPTOR> * temporary_lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
          management,
          defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
          new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));
  hemocell.lattice = temporary_lattice;

  //Set up boundaries (to help preinlet creation)
  //Outlet
  Box3D bb = hemocell.lattice->getBoundingBox();
  Box3D outlet(bb.x1-2,bb.x1,bb.y0,bb.y1,bb.z0,bb.z1);
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = new BoundaryConditionInstantiator3D
          < T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T,DESCRIPTOR> > ();
  boundary->addPressureBoundary0P(outlet,*hemocell.lattice,boundary::density);
  //Top
  Box3D topChannel(bb.x0, bb.x1, bb.y0, bb.y1, bb.z0, bb.z0 );
  Box3D bottomChannel(bb.x0, bb.x1, bb.y0, bb.y1, bb.z1, bb.z1 );

  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR> );

  //TODO
  //Bottom
  //TODO



  hlog << "(PreInlets) creating preInlet" << endl;
  hemocell.preInlet = new hemo::PreInlet(&hemocell,management);

   // Setting Preinlet slice
  Box3D slice = hemocell.lattice->getBoundingBox();
  slice.x1 = slice.x0 = slice.x0+2;
  hemocell.preInlet->preInletFromSlice(Direction::Xneg,slice);

  delete temporary_lattice;
  hemocell.lattice = 0;

  hlog << "(Stl preinlet) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.initializeLattice(management);



  // Boundaries for the domain
  if (!hemocell.partOfpreInlet) {
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = new BoundaryConditionInstantiator3D
          < T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T,DESCRIPTOR> > ();
  boundary->addPressureBoundary0P(outlet,*hemocell.lattice,boundary::density);
  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
  boundaryCondition->setVelocityConditionOnBlockBoundaries (*hemocell.lattice, topChannel );
  setBoundaryVelocity(*hemocell.lattice, topChannel, plb::Array<T,3>(0.75*velocity_max_lbm,0,0));
  }

  hemocell.preInlet->initializePreInlet();

  hemocell.lattice->periodicity().toggle(1,true);

  hlog << "(Stl preinlet) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl;

  hemocell.preInlet->createBoundary();
  if (hemocell.partOfpreInlet) {
    Box3D bb = hemocell.lattice->getBoundingBox();
    Box3D domain(bb.x0, bb.x1, bb.y0, bb.y1, bb.z0, bb.z0 );
    for (int x  = domain.x0 ; x <= domain.x1 ; x++) {
      for (int y  = domain.y0 ; y <= domain.y1 ; y++) {
        for (int z  = domain.z0 ; z <= domain.z1 ; z++) {
          defineDynamics(*hemocell.lattice,x,y,z,hemocell.lattice->getBackgroundDynamics().clone());
        }
      }
    }
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    boundaryCondition->setVelocityConditionOnBlockBoundaries (*hemocell.lattice, domain );
    setBoundaryVelocity(*hemocell.lattice, topChannel, plb::Array<T,3>(0.75*velocity_max_lbm,0,0));
  }


  hemocell.lattice->toggleInternalStatistics(false);

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  hemocell.lattice->initialize();

  // Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE,OUTPUT_BOUNDARY};
  hemocell.setFluidOutputs(outputs);

  // For the main simulation domain we have to define outlets
  if (!hemocell.partOfpreInlet) {
    Box3D bb = hemocell.lattice->getBoundingBox();
    Box3D outlet(bb.x1-2,bb.x1,bb.y0,bb.y1,bb.z0,bb.z1);
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundary = new BoundaryConditionInstantiator3D
          < T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T,DESCRIPTOR> > ();
    boundary->addPressureBoundary0P(outlet,*hemocell.lattice,boundary::density);
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
      //hemocell.preInlet->setDrivingForce();
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
      pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
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
