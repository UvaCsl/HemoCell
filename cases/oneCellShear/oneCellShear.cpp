#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "helper/hemocellInit.hh"
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

	


// ----------------- Read in config file & calc. LBM parameters ---------------------------
	pcout << "(OneCellShear) (Parameters) calculating shear flow parameters" << endl;
	//double nxyz = 20.0*(1e-6/(*cfg)["domain"]["dx"].read<double>());
  plint nx = 40;
  plint ny = 40;
  plint nz = 20;
  param::lbm_shear_parameters((*cfg),ny);
  param::printParameters();

	// ------------------------ Init lattice --------------------------------

	pcout << "(CellStretch) Initializing lattice: " << nx <<"x" << ny <<"x" << nz << " [lu]" << std::endl;

	plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.

	hemocell.lattice = new MultiBlockLattice3D<double,DESCRIPTOR>(
			defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
			defaultMultiBlockPolicy3D().getBlockCommunicator(),
			defaultMultiBlockPolicy3D().getCombinedStatistics(),
			defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
			new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

	pcout << "(OneCellShear) Re corresponds to u_max = " << (param::re * param::nu_p)/(hemocell.lattice->getBoundingBox().getNy()*param::dx) << " [m/s]" << endl;
	// -------------------------- Define boundary conditions ---------------------

	OnLatticeBoundaryCondition3D<double,DESCRIPTOR>* boundaryCondition
			= createLocalBoundaryCondition3D<double,DESCRIPTOR>();

	hemocell.lattice->toggleInternalStatistics(false);

	iniLatticeSquareCouette(*hemocell.lattice, nx, ny, nz, *boundaryCondition, param::shearrate_lbm);

	hemocell.lattice->initialize();

	// ----------------------- Init cell models --------------------------
	
	hemocell.initializeCellfield();
	hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
	vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA, OUTPUT_FORCE_VISC}; 
	hemocell.setOutputs("RBC_HO", outputs);

	outputs = {OUTPUT_VELOCITY};
	hemocell.setFluidOutputs(outputs);

// ---------------------- Initialise particle positions if it is not a checkpointed run ---------------

	//loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    pcout << "(OneCellShear) CHECKPOINT found!" << endl;
    hemocell.loadCheckPoint();
  }


  if (hemocell.iter == 0) { 
    pcout << "(OneCellShear) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl; 
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {  
      hemocell.lattice->collideAndStream();  
    } 
  }

  pcout << "(OneCellShea) Shear rate: " << (*cfg)["domain"]["shearrate"].read<double>() << " s^-1." << endl;

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  // Get undeformed cell values
  CellInformationFunctionals::calculateCellVolume(&hemocell);
  CellInformationFunctionals::calculateCellArea(&hemocell);
  double volume_eq = (CellInformationFunctionals::info_per_cell[0].volume)/pow(1e-6/param::dx,3);
  double surface_eq = (CellInformationFunctionals::info_per_cell[0].area)/pow(1e-6/param::dx,2);

  double D0 = 2.0 * (*cfg)["ibm"]["radius"].read<double>() * 1e6;

  // Creating output log file
  plb_ofstream fOut;
  if(cfg->checkpointed)
    fOut.open("stretch.log", std::ofstream::app);
  else
    fOut.open("stretch.log");


  while (hemocell.iter < tmax ) {
    
    hemocell.iterate();

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
      double rel_D2 = (largest_diam/D0)*(largest_diam/D0);
      double def_idx = (rel_D2 - 1.0) / (rel_D2 + 1.0) * 100.0;

      pcout << "\t Cell center at: {" <<position[0]<<","<<position[1]<<","<<position[2] << "} µm" << endl;  
      pcout << "\t Diameters: {" << bbox[1]-bbox[0] <<", " << bbox[3]-bbox[2] <<", " << bbox[5]-bbox[4] <<"}  µm" << endl;
      pcout << "\t Surface: " << surface << " µm^2" << " (" << surface / surface_eq * 100.0 << "%)" << "  Volume: " << volume << " µm^3" << " (" << volume / volume_eq * 100.0 << "%)"<< endl;
      pcout << "\t Largest diameter: " << largest_diam << " µm." << endl;
      pcout << "\t Deformation index: " << def_idx << " [%]" << endl;

      fOut << hemocell.iter << " " << bbox[1]-bbox[0] << " " << bbox[3]-bbox[2] << " " << bbox[5]-bbox[4] << " " << volume / volume_eq * 100.0 << " " << surface / surface_eq * 100.0 << " " << largest_diam << " " << def_idx << endl;

      CellInformationFunctionals::clear_list();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  fOut.close();
  pcout << "(OneCellShear) Simulation finished :)" << std::endl;
  return 0;
}