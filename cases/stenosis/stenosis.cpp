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
    StenosisShapeDomain3D(plint xbottomL_, plint xbottomR_, plint xtopL_, plint xtopR_, plint xcirc_, plint ycirc_, plint ybottom_, plint ytop_, plint radiusCyl_)
        : xbottomL(xbottomL_),
          xbottomR(xbottomR_),
          xtopL(xtopL_),
          xtopR(xtopR_),
          xcirc(xcirc_),
          ycirc(ycirc_),
          ybottom(ybottom_),
          ytop(ytop_),
          radiusCyl(radiusCyl_),
          radiusSqr(radiusCyl*radiusCyl)
    {}
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return ((iX-xcirc)*(iX-xcirc) + (iY-ycirc)*(iY-ycirc) <= radiusSqr) ||
               //(iX <= -(95/151)*(iZ-(9964/19))  && iX >= xbottomL && iZ <= ytop) ;
               (iX <= xtopR && iX >= xbottomL && iY <= ycirc ) ||
               //(iX <= ( iZ*(xbottomR-xtopR)/(-ycirc) + xbottomR ) && iX >= xtopR && iZ <= ycirc );
               (iX <= ((iY - 514.16683048)/-1.60677134525) && iX >= 127.73502714 && iY <= 308.92584909) ;
    }
    virtual StenosisShapeDomain3D<T>* clone() const {
        return new StenosisShapeDomain3D<T>(*this);
    }
private:
    plint xbottomL;
    plint xbottomR;
    plint xtopL;
    plint xtopR;
    plint xcirc;
    plint ycirc;
    plint ybottom;
    plint ytop;
    plint radiusCyl;
    plint radiusSqr;

};

//-------------------------------------------------------------------------------------

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
//                       (*cfg)["domain"]["blockSize"].read<int>()); 
     
//  pcout << "(PipeFlow) (Parameters) calculating flow parameters" << endl;
//  param::lbm_pipe_parameters((*cfg),flagMatrix->getBoundingBox());
//  param::printParameters();

  pcout << "(PipeFlow) (Parameters) calculating flow parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

// ---------------------------- Create geometry ------------------------------------------------


  pcout << "(main) setting dimensions ..." << std::endl;

  plint extendedEnvelopeWidth = 2;  // Depends on the requirements of the ibmKernel. 4 or even 2 might be enough (also depends on dx)

  plint lengthChannel = (*cfg)["domain"]["refDirN"].read<int>()*2;
  plint heightChannel = 2*174;
  plint depthChannel = 2*80;

  plint nx = 3*lengthChannel;
  plint ny = heightChannel;
  plint nz = depthChannel;

  plint radiusCyl = 2*7.5; //
  plint widthSt = 2*110.0;

  plint ytop = 2*158.0;
  plint xbottomL = 2*50.0;
  plint xbottomR = xbottomL + widthSt;
  plint xtopL = xbottomL;
  plint xtopR = xtopL + 2*radiusCyl;
  plint xcirc = xtopL + radiusCyl;
  plint ycirc = ytop - radiusCyl;
  //plint ytop = top - radiusCyl;
  plint ybottom = 0;

  pcout << "parameters Stenosis: " <<
  "ytopR = " << xtopR << ", " <<
  "xbottomR = " << xbottomR << ", " <<
  "ytopL = " << xtopL << ", " <<
  "radiuscyl = " << radiusCyl << ", " <<
  "ytop = " << ytop << ", " << endl;

// ----------------------------------------------------------------------------------------------------------------------------
  double shear_rate = 1800; //input shear rate s-1
  pcout << "shear_rate = " << shear_rate << endl;

  double flowQ = (shear_rate*130e-6*80e-6*80e-6)/6;
  pcout << "flow = " << flowQ << endl;
//  double u_max = flowQ * (100e-6) * (174e-6);
  //pcout << "u_max = " << u_max << endl;
  //double u_max_lbm = u_max * (*cfg)["domain"]["dx"].read<double>() / (*cfg)["domain"]["dt"].read<double>();
  //pcout << "u_max_lbm = " << u_max_lbm << endl;

// -------------------------------------------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------
//  plint extendedEnvelopeWidth = 1;  // Because we might use ibmKernel with with 2.
  
  pcout << "(PipeFlow) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),  //voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));

 // pcout << "(PipeFlow) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl; 
 // defineDynamics(*hemocell.lattice, *flagMatrix, (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new StenosisShapeDomain3D<T>(xbottomL, xbottomR, xtopL, xtopR, xcirc, ycirc, ybottom, ytop, radiusCyl),
                new BounceBack<T, DESCRIPTOR> );

  Box3D topChannel(0, nx-1, 0, ny-1, nz-1, nz-1);
  Box3D bottomChannel( 0, nx-1, 0, ny-1, 0, 0);
  Box3D backChannel( 0, nx-1, ny-1, ny-1, 0, nz-1);
  Box3D frontChannel( 0, nx-1, 0, 0, 0, nz-1);

  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, frontChannel, new BounceBack<T, DESCRIPTOR> );

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
 

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  double dpdz = (flowQ*12*3.0e-3)/(80e-6*80e-6*80e-6*130e-6); //(shear_rate * 2e-3) / (4*174e-6);
  double dpdz_lbm = dpdz * ((*cfg)["domain"]["dx"].read<double>() * (*cfg)["domain"]["dx"].read<double>() * (*cfg)["domain"]["dt"].read<double>()*(*cfg)["domain"]["dt"].read<double>() /param::dm);
  pcout << "dpdz_lbm = " << dpdz_lbm << endl;

  //Driving Force
  pcout << "(PipeFlow) (Fluid) Setting up driving Force" << endl; 
  double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce = dpdz_lbm; //u_max_lbm * 8 * param::nu_lbm / (lengthChannel*heightChannel ); //8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
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
    hemocell.iterate();
    
    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
    
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
