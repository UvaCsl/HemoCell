#include "hemocell.h"
#include "rbcHighOrderModel.h"
// #include "wbcHighOrderModel.h" // in case of WBC uncomment this
#include "helper/hemoCellStretch.h"
#include "helper/cellInfo.h"


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

  param::lbm_base_parameters(*cfg);

  param::ef_lbm = (*cfg)["parameters"]["stretchForce"].read<double>()*1e-12 / param::df;

  param::printParameters();
  
	// ------------------------ Init lattice --------------------------------

  plint nx = 50; //plint nx = 70; // in case of WBC use larger dimensions
  plint ny = 25; //plint ny = 35;
  plint nz = 25; //plint nz = 35;

	pcout << "(CellStretch) Initializing lattice: " << nx <<"x" << ny <<"x" << nz << " [lu]" << std::endl;

	plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.

	hemocell.lattice = new MultiBlockLattice3D<double,DESCRIPTOR>(
			defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
			defaultMultiBlockPolicy3D().getBlockCommunicator(),
			defaultMultiBlockPolicy3D().getCombinedStatistics(),
			defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
			new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));


	hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);

  OnLatticeBoundaryCondition3D<double,DESCRIPTOR>* boundaryCondition
          = createLocalBoundaryCondition3D<double,DESCRIPTOR>();

  boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice);
  setBoundaryVelocity(*hemocell.lattice, hemocell.lattice->getBoundingBox(), plb::Array<double,3>(0.,0.,0.) );
  // Box3D x_n = Box3D(0,0,0,ny,0,nz); 
  // boundaryCondition->addPressureBoundary0N(x_n, *hemocell.lattice);
  // Box3D x_p = Box3D(nx,nx,0,ny,0,nz); 
  // boundaryCondition->addPressureBoundary0N(x_p, *hemocell.lattice);
  // Box3D y_n = Box3D(0,nx,0,0,0,nz); 
  // boundaryCondition->addPressureBoundary0N(y_n, *hemocell.lattice);
  // Box3D y_p = Box3D(0,nx,ny,ny,0,nz); 
  // boundaryCondition->addPressureBoundary0N(y_p, *hemocell.lattice);
  // Box3D z_n = Box3D(0,nx,0,ny,0,0); 
  // boundaryCondition->addPressureBoundary0N(z_n, *hemocell.lattice);
  // Box3D z_p = Box3D(0,nx,0,ny,nz,nz); 
  // boundaryCondition->addPressureBoundary0N(z_p, *hemocell.lattice);
  
  // setBoundaryDensity(*hemocell.lattice,x_n,1.);
  // setBoundaryDensity(*hemocell.lattice,x_p,1.);
  // setBoundaryDensity(*hemocell.lattice,y_n,1.);
  // setBoundaryDensity(*hemocell.lattice,y_p,1.);
  // setBoundaryDensity(*hemocell.lattice,z_n,1.);
  // setBoundaryDensity(*hemocell.lattice,z_p,1.);


  hemocell.latticeEquilibrium(1., hemo::Array<double, 3>({0.,0.,0.}));

	hemocell.lattice->initialize();

	// ----------------------- Init cell models --------------------------
	
	hemocell.initializeCellfield();
	hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
	// hemocell.addCellType<WbcHighOrderModel>("WBC_HO", WBC_SPHERE);
	vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC}; 
	hemocell.setOutputs("RBC_HO", outputs);
	// hemocell.setOutputs("WBC_HO", outputs);

	outputs = {OUTPUT_VELOCITY,OUTPUT_FORCE};
	hemocell.setFluidOutputs(outputs);
        
  hemocell.outputInSiUnits = true; //HDF5 output in SI units (except location (so fluid location, particle location is still in LU)

// ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

	//loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    pcout << "(CellStretch) CHECKPOINT found!" << endl;
    hemocell.loadCheckPoint();
  }

// Setting up the stretching
  unsigned int n_forced_lsps = 1 + 6;// + 12;
  HemoCellStretch cellStretch(*(*hemocell.cellfields)["RBC_HO"],n_forced_lsps, param::ef_lbm);
  // HemoCellStretch cellStretch(*(*hemocell.cellfields)["WBC_HO"],n_forced_lsps, param::ef_lbm);
  
  pcout << "(CellStretch) External stretching force [pN(flb)]: " <<(*cfg)["parameters"]["stretchForce"].read<double>() << " (" << param::ef_lbm  << ")" << endl;
    
  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  // Get undeformed values
  double volume_lbm = (*hemocell.cellfields)["RBC_HO"]->meshmetric->getVolume();
  // double volume_lbm = (*hemocell.cellfields)["WBC_HO"]->meshmetric->getVolume();
  double surface_lbm = (*hemocell.cellfields)["RBC_HO"]->meshmetric->getSurface();
  // double surface_lbm = (*hemocell.cellfields)["WBC_HO"]->meshmetric->getSurface();
  double volume_eq = volume_lbm/pow(1e-6/param::dx,3);
  double surface_eq = surface_lbm/pow(1e-6/param::dx,2);

  hemo::Array<double,6> bb =  (*hemocell.cellfields)["RBC_HO"]->getOriginalBoundingBox();
  // hemo::Array<double,6> bb =  (*hemocell.cellfields)["WBC_HO"]->getOriginalBoundingBox();
  pcout << "Original Bounding box:" << endl;
  pcout << "\tx: " << bb[0] << " : " << bb[1] << endl;
  pcout << "\ty: " << bb[2] << " : " << bb[3] << endl;
  pcout << "\tz: " << bb[4] << " : " << bb[5] << endl;
  
  // Creating output log file
  plb_ofstream fOut;
  if(cfg->checkpointed)
    fOut.open("stretch.log", std::ofstream::app);
  else
    fOut.open("stretch.log");


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

      // Fill up the static info structure with desired data
      CellInformationFunctionals::calculateCellVolume(&hemocell);
      CellInformationFunctionals::calculateCellArea(&hemocell);
      CellInformationFunctionals::calculateCellPosition(&hemocell);
      CellInformationFunctionals::calculateCellStretch(&hemocell);
      CellInformationFunctionals::calculateCellBoundingBox(&hemocell);

      double volume = (CellInformationFunctionals::info_per_cell[0].volume)/pow(1e-6/param::dx,3);
      double surface = (CellInformationFunctionals::info_per_cell[0].area)/pow(1e-6/param::dx,2);
      hemo::Array<double,3> position = CellInformationFunctionals::info_per_cell[0].position/(1e-6/param::dx);
      hemo::Array<double,6> bbox = CellInformationFunctionals::info_per_cell[0].bbox/(1e-6/param::dx);
      double largest_diam = (CellInformationFunctionals::info_per_cell[0].stretch)/(1e-6/param::dx);

      pcout << "\t Cell center at: {" <<position[0]<<","<<position[1]<<","<<position[2] << "} µm" <<endl;  
      pcout << "\t Diameters: {" << bbox[1]-bbox[0] <<", " << bbox[3]-bbox[2] <<", " << bbox[5]-bbox[4] <<"}  µm" << endl;
      pcout << "\t Surface: " << surface << " µm^2" << " (" << surface / surface_eq * 100.0 << "%)" << "  Volume: " << volume << " µm^3" << " (" << volume / volume_eq * 100.0 << "%)"<<endl;
      pcout << "\t Largest diameter: " << largest_diam << " µm." << endl;
      
      fOut << hemocell.iter << " " << bbox[1]-bbox[0] << " " << bbox[3]-bbox[2] << " " << bbox[5]-bbox[4] << " " << volume / volume_eq * 100.0 << " " << surface / surface_eq * 100.0 << endl;

      CellInformationFunctionals::clear_list();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  fOut.close();
  pcout << "(CellStretch) Simulation finished :)" << std::endl;

  return 0;
}
