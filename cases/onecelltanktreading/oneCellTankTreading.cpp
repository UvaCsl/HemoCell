#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "helper/hemocellInit.hh"


int main(int argc, char* argv[])
{   
	if(argc < 2)
	{
			cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
			return -1;
	}

	HemoCell hemocell(argv[1], argc, argv);
	Config * cfg = hemocell.cfg;

	
	// ------------------------- Read in config file ------------------------------------------------

/*	std::vector<T> ea;
	document["parameters"]["eulerAngles"].read(ea);
	if (ea.size() != 3) { ea.resize(3, 0.0); }
	eulerAngles = Array<T,3>(ea[0], ea[1], ea[2]);
	for(int i =0; i < 3; i++)
			eulerAngles[i] *= pi/180.;
*/


// ---------------------------- Calc. LBM parameters -------------------------------------------------
	pcout << "(OneCellShear) (Parameters) calculating shear flow parameters" << endl;
	double nxyz = 20*(1e-6/(*cfg)["domain"]["dx"].read<double>());
  param::lbm_shear_parameters((*cfg),nxyz);

	// ------------------------ Init lattice --------------------------------

	pcout << "(OneCellShear) Initializing lattice: " << nxyz << "x, y and z (cube)" << std::endl;

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

	pcout << "(OneCellShear) dx = " << param::dx << ", " <<
			"dt = " << param::dt << ", " <<
			"dm = " << param::dm << ", " <<
			"dN = " << param::df << ", " <<
			"shear rate = " << param::shearrate_lbm << ", " <<
			std::endl;

	pcout << "(OneCellShear) tau = " << param::tau << " Re = " << param::re << " u_lb_max(based on Re) = " << param::u_lbm_max << " nu_lb = " << param::nu_lbm << endl;
	pcout << "(main) Re corresponds to u_max = " << (param::re * param::nu_p)/(hemocell.lattice->getBoundingBox().getNy()*param::dx) << " [m/s]" << endl;
	// -------------------------- Define boundary conditions ---------------------

	OnLatticeBoundaryCondition3D<double,DESCRIPTOR>* boundaryCondition
			= createLocalBoundaryCondition3D<double,DESCRIPTOR>();

	hemocell.lattice->toggleInternalStatistics(false);

	iniLatticeSquareCouette(*hemocell.lattice, nxyz, nxyz, nxyz, *boundaryCondition, param::shearrate_lbm);

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


  if (hemocell.iter == 0) { 
    pcout << "(PipeFlow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl; 
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {  
      hemocell.lattice->collideAndStream();  
    } 
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  while (hemocell.iter < tmax ) {
    
    hemocell.iterate();
    
    //Set force as required after each iteration
    //setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
    //    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
    //    Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    if (hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  return 0;
}
/*
    pcout << std::endl << "(main) Starting simulation i=" << initIter  << ", tmeas = " << tmeas << std::endl;

    global::timer("HDFOutput").start();
    writeHDF5(lattice, dx, dt, 0);
    writeCellField3D_HDF5(RBCField, dx, dt, 0);
    global::timer("HDFOutput").stop();

    SimpleFicsionProfiler simpleProfiler(tmeas);
    simpleProfiler.writeInitial(nx, ny, nz, nCells, numVerticesPerCell);

    plb_ofstream fOut("stretch.log");
    //Array<T,3> stretch = RBCField[0]->get3D(CCR_POSITION_MAX) - RBCField[0]->get3D(CCR_POSITION_MIN); // ***** Added by Britt
    //pcout << "dx: " << stretch[0] << "dy: " << stretch[1]<<  "dz: " << stretch[2] << endl; // *****
    //fOut << 1 << " " << stretch[0] << " " << stretch[1] << " " << stretch[2] << endl;


    // ------------------------ Starting main loop --------------------------
    

        
        if ((iter%tmeas)==0) {
            //global::profiler().writeReport();
            pcout << "(main) Iteration:" << iter << "; time/it: "<< dtIteration*1.0/tmeas << " physical time: " << iter*dt << " s;";
            //Array<T,3> stretch = RBCField[0]->get3D(CCR_POSITION_MAX) - RBCField[0]->get3D(CCR_POSITION_MIN); 
            //pcout << " dx: " << stretch[0] << " dy: " << stretch[1]<<  " dz: " << stretch[2] << endl; 
            //fOut << iter << " " << stretch[0] << " " << stretch[1] << " " << stretch[2] << endl;
            //RBCField[0]->saveMesh("stretchedCell.stl");

    }
    pcout << "(main) Simulation finished." << std::endl;
}*/
