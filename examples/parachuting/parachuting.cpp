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
#include <fenv.h>

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  hlogfile << "(PipeFlow) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl; 

  std::auto_ptr<MultiScalarField3D<int>> flagMatrix;
  std::auto_ptr<VoxelizedDomain3D<T>> voxelizedDomain; 
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),  
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),  
                       (*cfg)["domain"]["refDirN"].read<int>(),  
                       (*cfg)["domain"]["refDir"].read<int>(),  
                       voxelizedDomain, flagMatrix,  
                       (*cfg)["domain"]["blockSize"].read<int>(),
                       (*cfg)["domain"]["particleEnvelope"].read<int>()); 

  param::lbm_pipe_parameters((*cfg),flagMatrix.get());
  param::printParameters();
  
  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
            voxelizedDomain.get()->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

  defineDynamics(*hemocell.lattice, *flagMatrix.get(), (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1.,plb::Array<T, 3>(0.,0.,0.));

  //Driving Force
  T poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / param::pipe_radius / param::pipe_radius;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  //hemocell.setInitialMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
  //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC_HO", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);

  // Enable boundary particles
  //hemocell.enableBoundaryParticles((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["BRepCutoff"].read<T>(),(*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
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
    hlog << "(PipeFlow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  hlog << "(PipeFlow) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    hemocell.iterate();
    
    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    
    if (hemocell.iter % tmeas == 0) {
        hlog << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
        hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
        hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC_HO") <<
endl;
        FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); T toMpS = param::dx / param::dt;
        hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
        ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); T topN = param::df * 1.0e12;
        hlog << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

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

  hlog << "(main) Simulation finished :) " << endl;

  return 0;
}
