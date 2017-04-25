#include "hemocell.h"

int main(int argc, char *argv[]) {
	if(argc < 2)
	{
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

  pcout << "(PipeFlow) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
#if HEMOCELL_CFD_DYNAMICS == 1
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));
#elif HEMOCELL_CFD_DYNAMICS == 2
            new GuoExternalForceMRTdynamics<double, DESCRIPTOR>(1.0/param::tau)); // Use with MRT dynamics!
#endif

  pcout << "(PipeFlow) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl; 
  defineDynamics(*hemocell.lattice, *flagMatrix, (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
	hemocell.latticeEquilibrium(1.,Array<double, 3>(0.,0.,0.));

  //Driving Force
  pcout << "(PipeFlow) (Fluid) Setting up driving Force" << endl; 
  param::lbm_parameters((*cfg),flagMatrix->getBoundingBox());
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
                    Array<double, DESCRIPTOR<double>::d>(poiseuilleForce, 0.0, 0.0));

  hemocell.lattice->initialize();	


	

  return 0;
}
