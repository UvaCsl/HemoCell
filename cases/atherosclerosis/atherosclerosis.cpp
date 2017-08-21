#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include <fenv.h>

/// A functional, used to instantiate bounce-back nodes at the locations of the sphere
template<typename T>
class StenosisShapeDomain3D : public plb::DomainFunctional3D {
public:
    StenosisShapeDomain3D(plint cx_, plint cy_, plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(radius*radius)
    { }
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return (iX-cx) * (iX-cx) + (iY-cy) * (iY-cy) <= radiusSqr;
    }
    virtual StenosisShapeDomain3D<T>* clone() const {
        return new StenosisShapeDomain3D<T>(*this);
    }
private:
    plint cx, cy;
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
//  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),  
//                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),  
//                       (*cfg)["domain"]["refDirN"].read<int>(),  
//                       (*cfg)["domain"]["refDir"].read<int>(),  
//                       voxelizedDomain, flagMatrix,  
 //                      (*cfg)["domain"]["blockSize"].read<int>()); 
     
//  pcout << "(PipeFlow) (Parameters) calculating flow parameters" << endl;
//  param::lbm_pipe_parameters((*cfg),flagMatrix->getBoundingBox());
//  param::printParameters();

  pcout << "(PipeFlow) (Parameters) calculating flow parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

// ---------------------------- Create geometry ------------------------------------------------

  pcout << "(main) Flowaroundsphere geometry ..." << std::endl;

  plint lengthChannel = 2*(*cfg)["domain"]["refDirN"].read<int>();
  plint heightChannel = 2*17; //52um

  plint nx = 6*lengthChannel;
  plint ny = lengthChannel;
  plint nz = heightChannel; 

  plint sphere_diameter = 0.8*2*lengthChannel; // 20, 40, 60, 80 % stenosis.

  plint sphere_x = 3*lengthChannel;
  plint sphere_y = lengthChannel;

// ----------------------------------------------------------------------------------------------------------------------------
  double shear_rate = 1000; //s-1
  pcout << "shear_rate = " << shear_rate << endl;
  
  double flowQ = (0.5e-6)/3600; // m^3/s //shear_rate * ((100e-6 * 52e-6 * 52e-6)/6);
  pcout << "flow = " << flowQ << endl;
//  double u_max = flowQ * (100e-6) * (52e-6);
//  pcout << "u_max = " << u_max << endl;
//  double u_max_lbm = u_max * (*cfg)["domain"]["dx"].read<double>() / (*cfg)["domain"]["dt"].read<double>();
//  pcout << "u_max_lbm = " << u_max_lbm << endl;
// ---------------------------------------------------------------------------------------------
  plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.
  
  pcout << "(PipeFlow) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),  //voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));


  pcout << "(PipeFlow) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl; 
//  defineDynamics(*hemocell.lattice, *flagMatrix, (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
              new StenosisShapeDomain3D<T>(sphere_x,sphere_y, sphere_diameter/2.0),
              new BounceBack<T, DESCRIPTOR> ); //b

  Box3D topChannel(0, 6*lengthChannel-1, 0, lengthChannel-1, heightChannel-1, heightChannel-1 );
  Box3D bottomChannel( 0, 6*lengthChannel-1, 0, lengthChannel-1, 0, 0);
  Box3D backChannel( 0, 6*lengthChannel-1, lengthChannel-1, lengthChannel-1, 0, heightChannel-1);
  Box3D frontChannel( 0, 6*lengthChannel-1, 0, 0, 0, heightChannel-1);

  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, frontChannel, new BounceBack<T, DESCRIPTOR> );

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  double dpdz = (flowQ*12*3.0e-3)/(17e-6*17e-6*17e-6*300e-6); //(shear_rate * 2e-3) / (4*52e-6);
  double dpdz_lbm = dpdz * ((*cfg)["domain"]["dx"].read<double>() * (*cfg)["domain"]["dx"].read<double>() * (*cfg)["domain"]["dt"].read<double>()*(*cfg)["domain"]["dt"].read<double>() /param::dm);
  pcout << "dpdz_lbm = " << dpdz_lbm << endl;

  //Driving Force
  pcout << "(PipeFlow) (Fluid) Setting up driving Force" << endl; 
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce = dpdz_lbm; //u_max_lbm * 8 * param::nu_lbm / (lengthChannel*heightChannel );//8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
                    plb::Array<double, DESCRIPTOR<double>::d>(poiseuilleForce, 0.0, 0.0));

  pcout << "poiseuilleForce = " << poiseuilleForce << endl;

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setMinimumDistanceFromSolid("RBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<double>(), (*cfg)["domain"]["RepCutoff"].read<double>());
  //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC_HO", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
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

  pcout << "(PipeFlow) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    hemocell.iterate();
    

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
