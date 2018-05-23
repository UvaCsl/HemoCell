#include "hemocell.h"
#include "wbcHighOrderModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include <fenv.h>

/// A functional, used to instantiate bounce-back nodes at the locations of the sphere
template<typename T>
class TriangleShapeDomain3D : public plb::DomainFunctional3D {
public:
    TriangleShapeDomain3D(plint xbottomL_, plint xbottomR_, plint ytop_, plint ybottom_, plint ny_, double slope_)
        : xbottomL(xbottomL_),
          xbottomR(xbottomR_),
          ytop(ytop_),
          ybottom(ybottom_),
          ny(ny_),
          slope(slope_)

    {}
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return (iX > xbottomL && iX <= xbottomR && 
        	( iY <= ybottom-slope*(iX-xbottomL) || iY >= ytop+slope*(iX-xbottomL) ) );

    }
    virtual TriangleShapeDomain3D<T>* clone() const {
        return new TriangleShapeDomain3D<T>(*this);
    }
private:
    plint xbottomL;
    plint xbottomR;
    plint ytop;
    plint ybottom;
    plint ny;
    double slope;

};

//-------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  pcout << "(Capillary_lykov) (Parameters) calculating flow parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

// ---------------------------- Create geometry ------------------------------------------------


  pcout << "(main) setting dimensions ..." << std::endl;

  plint extendedEnvelopeWidth = 2;  // Depends on the requirements of the ibmKernel. 4 or even 2 might be enough (also depends on dx)

  plint lengthChannel = (*cfg)["domain"]["refDirN"].read<int>();

  plint nx = lengthChannel ;
  plint ny = 36;
  plint nz = 36;

  
  plint triangleLength = 50;
  plint contractionGap = 12;

  plint ybottom = (ny-contractionGap)/2;
  plint ytop = (ny-1) - ybottom;
  plint xbottomL = (lengthChannel - triangleLength)/2;
  plint xbottomR = xbottomL + triangleLength;
  double slope = 3.4 / 25.0; // half triangle height decrease: 3.4um    trangle length: 25um

  pcout << "parameters of the triangle: " <<
  "slope = " << slope << ", " <<
  "ytop = " << ytop << ", " <<
  "ybottom = " << ybottom << ", " <<
  "xbottomL = " << xbottomL << ", " <<
  "xbottomR = " << xbottomR << ", " << endl;
 
  pcout << "(Capillary_lykov) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),  //voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));

  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new TriangleShapeDomain3D<T>(xbottomL, xbottomR, ytop, ybottom, ny, slope),
                new BounceBack<T, DESCRIPTOR> );

  hemocell.lattice->toggleInternalStatistics(false);
  //hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.lattice->periodicity().toggle(1,true);
  hemocell.lattice->periodicity().toggle(2,true);
 

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  double dpdz = 6.7e5; // 
  double dpdz_lbm = dpdz * ((*cfg)["domain"]["dx"].read<double>() * (*cfg)["domain"]["dx"].read<double>() * (*cfg)["domain"]["dt"].read<double>()*(*cfg)["domain"]["dt"].read<double>() /param::dm);
  pcout << "dpdz_lbm = " << dpdz_lbm << endl;

  //Driving Force
  pcout << "(Capillary_lykov) (Fluid) Setting up driving Force" << endl; 
  //double rPipe = (*cfg)["domain"]["refDirN"].read<int>()/2.0;
  double poiseuilleForce = dpdz_lbm; //u_max_lbm * 8 * param::nu_lbm / (lengthChannel*heightChannel ); //8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
                    plb::Array<double, DESCRIPTOR<double>::d>(poiseuilleForce, 0.0, 0.0));

  pcout << "poiseuilleForce = " << poiseuilleForce << endl;

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<WbcHighOrderModel>("WBC_HO", WBC_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("WBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setMinimumDistanceFromSolid("WBC_HO", 1); //Micrometer! not LU

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<double>(), (*cfg)["domain"]["RepCutoff"].read<double>());
  //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK, OUTPUT_FORCE_INNER_LINK, OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC, OUTPUT_VERTEX_ID, OUTPUT_CELL_ID};
  hemocell.setOutputs("WBC_HO", outputs);
  
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
    pcout << "(Capillary_lykov) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  pcout << "(Capillary_lykov) Starting simulation..." << endl;

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
      pcout << " | # of WBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "WBC_HO");
    
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
