#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "preInlet.h"
#include <fenv.h>

/// A functional, used to instantiate bounce-back nodes at the locations of the sphere
template<typename T>
class SphereShapeDomain3D : public plb::DomainFunctional3D {
public:
    SphereShapeDomain3D(plint cx_, plint cy_, plint cz_, plint radius)
        : cx(cx_),
          cy(cy_),
          cz(cz_),
          radiusSqr(radius*radius)
    { }
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return (iX-cx) * (iX-cx) + (iY-cy) * (iY-cy) + (iZ-cz) * (iZ-cz) <= radiusSqr;
    }
    virtual SphereShapeDomain3D<T>* clone() const {
        return new SphereShapeDomain3D<T>(*this);
    }
private:
    plint cx, cy, cz;
    plint radiusSqr;
};

// ----------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;


//  pcout << "(PipeFlow) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl; 
//  MultiScalarField3D<int> *flagMatrix = 0; 
//  VoxelizedDomain3D<double> * voxelizedDomain = 0; 
// getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),  
//                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),  
//                       (*cfg)["domain"]["refDirN"].read<int>(),  
//                       (*cfg)["domain"]["refDir"].read<int>(),  
//                       voxelizedDomain, flagMatrix,  
//                       (*cfg)["domain"]["blockSize"].read<int>()); 
     
  pcout << "(Flowaroundsphere) (Parameters) calculating flow parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();
 
// ---------------------------- Create sphere geometry ------------------------------------------------

  pcout << "(main) Flowaroundsphere geometry ..." << std::endl;

  plint lengthChannel = 2*(*cfg)["domain"]["refDirN"].read<int>();

  plint nx = lengthChannel; //length preinlet 50um
  plint ny = lengthChannel;
  plint nz = lengthChannel;

  plint sphere_diameter = 15*2; // 2, 5, 9 or 15 um.

  plint sphere_x = 2*lengthChannel/4;
  plint sphere_y = lengthChannel/2;
  plint sphere_z = sphere_diameter/2;

// ----------------------------------------------------------------------------------------------------------------------------
  double shear_rate = 1800; //input shear rate s-1
  pcout << "shear_rate = " << shear_rate << endl;

  double velocity_max = (shear_rate*(lengthChannel/1e6))/4;
  pcout << "velocity_max = " << velocity_max << endl;

  double velocity_max_lbm = velocity_max * ( (*cfg)["domain"]["dt"].read<double>() /(*cfg)["domain"]["dx"].read<double>() );
  pcout << "velocity_max_lbm = " << velocity_max_lbm << endl;

// ---------------------------------------------------------------------------------------------
  plint extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with with 2.

  pcout << "(Flowaroundsphere) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),     
       //voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));

//  pcout << "(PipeFlow) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl; 
//  defineDynamics(*hemocell.lattice, *flagMatrix, (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

//--------------------------------- boundary conditions ---------------------------------------------------
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

//  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
//              new SphereShapeDomain3D<T>(sphere_x,sphere_y,sphere_z, sphere_diameter/2.0),
//              new BounceBack<T, DESCRIPTOR> ); //b

  Box3D topChannel(0, 2*lengthChannel-1, 0, lengthChannel-1, lengthChannel-1, lengthChannel-1 );
  Box3D bottomChannel( 0, 2*lengthChannel-1, 0, lengthChannel-1, 0, 0);

  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );

  boundaryCondition->setVelocityConditionOnBlockBoundaries (*hemocell.lattice, topChannel );
  setBoundaryVelocity(*hemocell.lattice, topChannel, plb::Array<T,3>(0.75*velocity_max_lbm,0,0)); // #calculated form: vmax=shearrate*h/32 = 0.04375, v_topwall,p = 0.75*vmax, v_top_lbm=v_topwall*(dt/dx)

//  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR> );

// -----------------------------------

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.lattice->periodicity().toggle(1,true); //insert if periodic boundery condtions in y-direction
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  //Driving Force
//  pcout << "(PipeFlow) (Fluid) Setting up driving Force" << endl; 
//  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
//  double poiseuilleForce = 0.0; // 8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
//  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
//                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
//                    plb::Array<double, DESCRIPTOR<double>::d>(poiseuilleForce, 0.0, 0.0));

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
 
  hemocell.setMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
 
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());
  
  //repulsion force aangezet...
//  hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<double>(), (*cfg)["domain"]["RepCutoff"].read<double>());
//  hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES};
  hemocell.setOutputs("RBC_HO", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY};
  hemocell.setFluidOutputs(outputs);

  //todo add statistics here
  
  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
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
  unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  //Setting Preinlet creation
  Box3D preinletBox(0,20,0,hemocell.lattice->getBoundingBox().y1,0,hemocell.lattice->getBoundingBox().z1);
  createPreInlet preinlet(preinletBox,"preinlet",(*cfg)["ibm"]["stepMaterialEvery"].read<int>(),Direction::Xpos,hemocell,tmax,true);

  pcout << "(Flowaroundsphere) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    preinlet.saveCurrent();
    hemocell.iterate();
//  pcout << hemocell.iter << endl;
    
//    //Set driving force as required after each iteration
//    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
//               DESCRIPTOR<T>::ExternalField::forceBeginsAt,
//               plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    
    // Only enable if PARMETIS build is available
    // if (hemocell.iter % tbalance == 0) {
    //   if(hemocell.calculateFractionalLoadImbalance() > 3) {
    //    // hemocell.doLoadBalance();
    //   }
    // }
   
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
      //           plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
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
