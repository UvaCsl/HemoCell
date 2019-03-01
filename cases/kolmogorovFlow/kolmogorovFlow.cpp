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
#include "helper/hemocellInit.hh"
#include "helper/cellInfo.h"
#include "helper/fluidInfo.h"
#include "helper/particleInfo.h"


int main(int argc, char* argv[])
{   
	if(argc < 2)
	{
			cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
			return -1;
	}

	HemoCell hemocell(argv[1], argc, argv);
	Config * cfg = hemocell.cfg;

	

// ----------------- Read in config file & calc. LBM parameters ---------------------------
	pcout << "(KolmogorovFlow) (Parameters) calculating flow parameters" << endl;
	plint domainSize = (*cfg)["domain"]["refDirN"].read<int>();
  plint nx = domainSize;
  plint ny = domainSize;
  plint nz = domainSize;
  param::lbm_pipe_parameters((*cfg),ny/4);
  param::printParameters();

	// ------------------------ Init lattice --------------------------------

	pcout << "(KolmogorovFlow) Initializing lattice: " << nx <<"x" << ny <<"x" << nz << " [lu]" << std::endl;

	plint extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with with 2.

	hemocell.lattice = new MultiBlockLattice3D<double,DESCRIPTOR>(
			defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
			defaultMultiBlockPolicy3D().getBlockCommunicator(),
			defaultMultiBlockPolicy3D().getCombinedStatistics(),
			defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
			new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

	pcout << "(KolmogorovFlow) Re corresponds to u_max = " << (param::re * param::nu_p)/(hemocell.lattice->getBoundingBox().getNy()/2.0*param::dx) << " [m/s]" << endl;

	//Driving Force
  pcout << "(KolmogorovFlow) (Fluid) Setting up driving force using parallel planes approximation" << endl; 
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/4.0;
  double poiseuilleForce =  16 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
	// -------------------------- Define boundary conditions ---------------------

	OnLatticeBoundaryCondition3D<double,DESCRIPTOR>* boundaryCondition
			= createLocalBoundaryCondition3D<double,DESCRIPTOR>();

	hemocell.lattice->toggleInternalStatistics(false);

	Box3D top   = Box3D(0, nx-1, 0,  (ny-1)/2,    0, nz-1);
  Box3D bottom  = Box3D(0, nx-1, (ny-1)/2+1, ny-1, 0, nz-1);
  if (ny%2 != 0) {
    top.y1--;
  }

  pcout << "TOP domain height: " << top.y1 - top.y0 <<  " BOTTOM domain height: " << bottom.y1 - bottom.y0 << endl;
	hemocell.lattice->initialize();


	// ----------------------- Init cell models --------------------------
	
	hemocell.initializeCellfield();
	hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
	hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
	vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA, OUTPUT_FORCE_VISC, OUTPUT_FORCE_INNER_LINK, OUTPUT_VERTEX_ID, OUTPUT_CELL_ID}; 
	hemocell.setOutputs("RBC", outputs);

	hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
	hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
	hemocell.setOutputs("PLT", outputs);

	hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

	outputs = {OUTPUT_VELOCITY,OUTPUT_FORCE};
	hemocell.setFluidOutputs(outputs);


	// Turn on periodicity in all directions
  	hemocell.setSystemPeriodicity(0, true);
  	hemocell.setSystemPeriodicity(1, true);
  	hemocell.setSystemPeriodicity(2, true);

// ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

	//loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    pcout << "(KolmogorovFlow) CHECKPOINT found!" << endl;
    hemocell.loadCheckPoint();
  }


  if (hemocell.iter == 0) { 
    pcout << "(KolmogorovFlow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl; 
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {  
      hemocell.lattice->collideAndStream();  
    } 
  }


  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();


  while (hemocell.iter < tmax ) {
    setExternalVector( *hemocell.lattice, top,
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, plb::Array<double,DESCRIPTOR<T>::d>(poiseuilleForce,0.0,0.0));

    setExternalVector( *hemocell.lattice, bottom,
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, plb::Array<double,DESCRIPTOR<T>::d>(-poiseuilleForce,0.0,0.0));
    
    hemocell.iterate();

    if (hemocell.iter % tmeas == 0) {
      pcout << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      pcout << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
      pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
      pcout << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
      ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); double topN = param::df * 1.0e12;
      pcout << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;
    setExternalVector( *hemocell.lattice, top,
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, plb::Array<double,DESCRIPTOR<T>::d>(poiseuilleForce,0.0,0.0));

    setExternalVector( *hemocell.lattice, bottom,
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, plb::Array<double,DESCRIPTOR<T>::d>(-poiseuilleForce,0.0,0.0));
      hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  pcout << "(KolmogorovFlow) Simulation finished :)" << std::endl;
  return 0;
}
