#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "helper/hemoCellStretch.h"
#include "helper/cellInformation.h"


int main(int argc, char* argv[])
{   
	if(argc < 2)
	{
			cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
			return -1;
	}

	HemoCell hemocell(argv[1], argc, argv);
	Config * cfg = hemocell.cfg;

// ---------------------------- Calc. LBM parameters -------------------------------------------------
	pcout << "(CellStretch) (Parameters) calculating shear flow parameters" << endl;
	double nxyz = 30*(1e-6/(*cfg)["domain"]["dx"].read<double>());
  (*cfg)["domain"]["dt"].read(param::dt);
  (*cfg)["domain"]["dx"].read(param::dx);
  (*cfg)["domain"]["nuP"].read(param::nu_p);
  (*cfg)["domain"]["kBT"].read(param::kBT_p);
  (*cfg)["domain"]["rhoP"].read(param::rho_p);
  param::tau = 3.0 * (param::nu_p * param::dt / (param::dx * param::dx)) + 0.5;
  param::nu_lbm = 1.0/3.0 * (param::tau - 0.5);
  param::dm = param::rho_p * (param::dx * param::dx * param::dx);
  param::df = param::dm * param::dx / (param::dt * param::dt);
  param::ef_lbm = (*cfg)["parameters"]["stretchForce"].read<double>()*1e-12 / param::df;
  param::kBT_lbm = param::kBT_p/(param::df*param::dx);
  param::printParameters();
  
	// ------------------------ Init lattice --------------------------------

	pcout << "(CellStretch) Initializing lattice: " << nxyz << "^3 [lu] cube" << std::endl;

	plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.

			hemocell.lattice = new MultiBlockLattice3D<double,DESCRIPTOR>(
					defaultMultiBlockPolicy3D().getMultiBlockManagement(nxyz, nxyz, nxyz, extendedEnvelopeWidth),
					defaultMultiBlockPolicy3D().getBlockCommunicator(),
					defaultMultiBlockPolicy3D().getCombinedStatistics(),
					defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
	#if HEMOCELL_CFD_DYNAMICS == 1
					new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));
	#elif HEMOCELL_CFD_DYNAMICS == 2
					new GuoExternalForceMRTdynamics<T, DESCRIPTOR>(1.0/param::tau)); // Use with MRT dynamics!
	#endif


	hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);

  OnLatticeBoundaryCondition3D<double,DESCRIPTOR>* boundaryCondition
          = createLocalBoundaryCondition3D<double,DESCRIPTOR>();
  Box3D x_n = Box3D(0,0,0,nxyz,0,nxyz); 
  boundaryCondition->addPressureBoundary0N(x_n, *hemocell.lattice);
  Box3D x_p = Box3D(nxyz,nxyz,0,nxyz,0,nxyz); 
  boundaryCondition->addPressureBoundary0N(x_p, *hemocell.lattice);
  Box3D y_n = Box3D(0,nxyz,0,0,0,nxyz); 
  boundaryCondition->addPressureBoundary0N(y_n, *hemocell.lattice);
  Box3D y_p = Box3D(0,nxyz,nxyz,nxyz,0,nxyz); 
  boundaryCondition->addPressureBoundary0N(y_p, *hemocell.lattice);
  Box3D z_n = Box3D(0,nxyz,0,nxyz,0,0); 
  boundaryCondition->addPressureBoundary0N(z_n, *hemocell.lattice);
  Box3D z_p = Box3D(0,nxyz,0,nxyz,nxyz,nxyz); 
  boundaryCondition->addPressureBoundary0N(z_p, *hemocell.lattice);
  
  setBoundaryDensity(*hemocell.lattice,x_n,1.);
  setBoundaryDensity(*hemocell.lattice,x_p,1.);
  setBoundaryDensity(*hemocell.lattice,y_n,1.);
  setBoundaryDensity(*hemocell.lattice,y_p,1.);
  setBoundaryDensity(*hemocell.lattice,z_n,1.);
  setBoundaryDensity(*hemocell.lattice,z_p,1.);


  hemocell.latticeEquilibrium(1., Array<double, 3>(0.,0.,0.));

	hemocell.lattice->initialize();

	// ----------------------- Init cell models --------------------------
	
	hemocell.initializeCellfield();
	hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
	vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA}; 
	hemocell.setOutputs("RBC_HO", outputs);

	outputs = {OUTPUT_VELOCITY};
	hemocell.setFluidOutputs(outputs);

// ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

	//loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles((*cfg)["sim"]["particlePosFile"].read<string>());
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

// Setting up the stretching
  unsigned int n_forced_lsps = 7;
  HemoCellStretch cellStretch(*(*hemocell.cellfields)["RBC_HO"],n_forced_lsps, param::ef_lbm);
  
  pcout << "(CellStretch) External stretching force [flb]: " << param::ef_lbm << endl;
    
  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  while (hemocell.iter < tmax ) {  

    hemocell.cellfields->applyConstitutiveModel();    // Calculate Force on Vertices

    cellStretch.applyForce(); //IMPORTANT, not done normally in hemocell.iterate()
    
    hemocell.cellfields->spreadParticleForce();
    hemocell.lattice->timedCollideAndStream();
    hemocell.cellfields->interpolateFluidVelocity();
    hemocell.cellfields->syncEnvelopes();
    hemocell.cellfields->advanceParticles();
    hemocell.iter++;
    
    if (hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();
      CellInformationFunctionals::getCellVolume(&hemocell);
      pcout << "Our Cell has a volume of: " << (CellInformationFunctionals::info_per_cell[0].volume)/pow(1e-6/param::dx,3) << " µm^3" <<endl;
      CellInformationFunctionals::getCellArea(&hemocell);
      pcout << "Our Cell has an area of:  " << (CellInformationFunctionals::info_per_cell[0].area)/pow(1e-6/param::dx,2) << " µm^2" <<endl;
      CellInformationFunctionals::getCellPosition(&hemocell);
      Array<double,3> position = CellInformationFunctionals::info_per_cell[0].position/(1e-6/param::dx);
      pcout << "Our cells center is at:   {"<<position[0]<<","<<position[1]<<","<<position[2] << "} µm" <<endl;
      CellInformationFunctionals::getCellStretch(&hemocell);
      pcout << "Our Cell has a stretch of:"<<(CellInformationFunctionals::info_per_cell[0].stretch)/(1e-6/param::dx) << " µm" <<endl;
      CellInformationFunctionals::getCellBoundingBox(&hemocell);
      Array<double,6> bbox = CellInformationFunctionals::info_per_cell[0].bbox/(1e-6/param::dx);
      pcout << "Our Cell has a bbox of:   x: {"<< bbox[0] << "," << bbox[1] << "}" << endl;
      pcout << "                          y: {"<< bbox[2] << "," << bbox[3] << "}" << endl;
      pcout << "                          z: {"<< bbox[4] << "," << bbox[5] << "} µm" <<endl;
      
      CellInformationFunctionals::clear_list();


    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  return 0;
}
