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
#include "wbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include <fenv.h>
#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

/// A functional, used to instantiate bounce-back nodes at the locations of the sphere
template<typename T>
class StenosisShapeDomain3D : public plb::DomainFunctional3D {
public:
    StenosisShapeDomain3D( plint xtopL_, plint xtopR_, plint xcircL_, plint xcircR_, plint ycirc_, plint ytop_, plint radiusCyl_,double a_, double bL_, double bR_, double y_)
        : //xbottomL(xbottomL_),
          //xbottomR(xbottomR_),
          xtopL(xtopL_),
          xtopR(xtopR_),
          xcircL(xcircL_),
          xcircR(xcircR_),
          ycirc(ycirc_),
          //ybottom(ybottom_),
          ytop(ytop_),
          radiusCyl(radiusCyl_),
          radiusSqr(radiusCyl*radiusCyl),
          a(a_),
          bL(bL_),
          bR(bR_),
          y(y_)

    {}
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return ((iX-xcircL)*(iX-xcircL) + (iY-ycirc)*(iY-ycirc) <= radiusSqr) ||
               ((iX-xcircR)*(iX-xcircR) + (iY-ycirc)*(iY-ycirc) <= radiusSqr) ||
               (iX <= xcircR && iX >= xcircL && iY <= ytop ) ||
               (iX >= (iY-bL)/a && iX <= xcircL && iY <= y) ||
               (iX <= (iY-bR)/(-a) && iX >= xcircR && iY <= y );

    }
    virtual StenosisShapeDomain3D<T>* clone() const {
        return new StenosisShapeDomain3D<T>(*this);
    }
private:
    //plint xbottomL;
    //plint xbottomR;
    plint xtopL;
    plint xtopR;
    plint xcircL;
    plint xcircR;
    plint ycirc;
    //plint ybottom;
    plint ytop;
    plint radiusCyl;
    plint radiusSqr;
    double a;
    double bL;
    double bR;
    double y;

};



//-------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

//  pcout << "(Flowchamber) (Parameters) calculating flow parameters" << endl;
//  param::lbm_pipe_parameters((*cfg),ny);
//  param::printParameters();


//-------------------------------------------------------------------------------------

  hlogfile << "(main) setting dimensions .." << endl;

  plint widthSt = 2*(*cfg)["parameters"]["widthStenosis"].read<int>();; //length of stenosis
  plint radiusCyl = 2*5; // rounding stenosis
  double pi = std::acos(-1);
  plint c_angle_degrees = (*cfg)["parameters"]["angleStenosis"].read<int>(); //60; //angle of stenosis
  double percentageSt = (*cfg)["parameters"]["percentageStenosis"].read<double>(); //0.8;
  double angle = (90-c_angle_degrees) * pi/180; //in degrees
  double c_angle = c_angle_degrees *pi/180;
  //calculate raakpunten and b from y =ax+b
  double h = std::sin(angle)*radiusCyl;
  double w = std::cos(angle)*radiusCyl;
  double a = std::tan(c_angle);
  plint widthChannel = 2*(*cfg)["parameters"]["widthChannel"].read<int>(); //4*(*cfg)["domain"]["refDirN"].read<int>()+60;
  plint heightChannel = 2*(*cfg)["parameters"]["heightChannel"].read<int>(); //4*(*cfg)["domain"]["refDirN"].read<int>();
//  cout << "HeightChannel = " << heightChannel << endl;
//  cout << "percentageSt = " << percentageSt << endl;
//  cout << "som = " << heightChannel*percentageSt << endl;
//  cout << " widthCons = " << (heightChannel*percentageSt)/a << endl;
  plint widthConst =  (heightChannel*percentageSt)/a;
  plint lengthChannel = 4*(*cfg)["domain"]["refDirN"].read<int>()+widthSt+(2*widthConst);

  plint nx = lengthChannel;
  plint ny = heightChannel;
  plint nz = widthChannel;

  plint hydraulic_radius = (2*ny*nz)/(2*ny+2*nz);

  pcout << "(Flowchamber) (Parameters) calculating flow parameters" << endl;
  param::lbm_pipe_parameters((*cfg),hydraulic_radius);
  param::printParameters();

  plint ytop = heightChannel*percentageSt; //percentage of stenosis
  plint xtopL = nx/2 - widthSt/2;
  plint xtopR = nx/2 + widthSt/2;
//  plint xbottomL = xtopL-2*46;
//  plint xbottomR = xtopR+2*46;
  plint xcircL = xtopL + radiusCyl;
  plint xcircR = xtopR - radiusCyl;
  plint ycirc = ytop - radiusCyl;
  //plint ytop = top - radiusCyl;
//  plint ybottom = 0;

  double xL = xcircL-w;
  double y = ycirc+h;
  double xR = xcircR+w;

  double bL = y - a*xL;
  double bR = y + a*xR;

  pcout << "parameters Stenosis: " <<
  "pi = " << pi << "," <<
  "angle = " << angle << "," <<
  "a = " << a << "," <<
  "bL = " << bL << "," <<
  "bR = " << bR << "," <<
  "y = " << y << "," <<
  "ytopR = " << xtopR << ", " <<
//  "xbottomR = " << xbottomR << ", " <<
  "ytopL = " << xtopL << ", " <<
  "radiuscyl = " << radiusCyl << ", " <<
  "w = " << w << "," <<
  "h = " << h << "," <<
  "widthConst = " << widthConst << "," <<
  "ytop = " << ytop << ", " << endl;


//-------------------------------------------------------------------------------------

  hlog << "(Flowchamber) (Fluid) Initializing Palabos Fluid Field" << endl;

  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>()),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

//  Box3D topChannel(0, nx-1, 0, ny-1, nz-1, nz-1);
//  Box3D bottomChannel( 0, nx-1, 0, ny-1, 0, 0);
  Box3D backChannel( 0, nx-1, ny-1, ny-1, 0, nz-1);
  Box3D frontChannel( 0, nx-1, 0, 0, 0, nz-1);

//  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR>(1.) );
//  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR>(1.) );
  defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR>(1.) );
  defineDynamics(*hemocell.lattice, frontChannel, new BounceBack<T, DESCRIPTOR>(1.) );

  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new StenosisShapeDomain3D<T>(xtopL, xtopR, xcircL, xcircR, ycirc, ytop, radiusCyl, a, bL, bR, y),
                new BounceBack<T, DESCRIPTOR>(1.) );

 //  defineDynamics(*hemocell.lattice, *flagMatrix, (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
  hemocell.lattice->periodicity().toggle(2,true);

  hemocell.latticeEquilibrium(1.,plb::Array<T, 3>(0.,0.,0.));

// ---------------------------------------------------------------------------------------------

  //Driving Force
  T rPipe = hydraulic_radius;
  T poiseuilleForce =  4.5e-6;//8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

  hemocell.lattice->initialize();

// ---------------------------------------------------------------------------------------------

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 1); //Micrometer! not LU

  hemocell.addCellType<WbcHighOrderModel>("WBC_HO", WBC_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("WBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("WBC_HO", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
  //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

// ---------------------------------------------------------------------------------------------

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("WBC_HO", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE,OUTPUT_SHEAR_STRESS,OUTPUT_STRAIN_RATE,OUTPUT_SHEAR_RATE,OUTPUT_BOUNDARY};
  hemocell.setFluidOutputs(outputs);

// ---------------------------------------------------------------------------------------------

  // Enable boundary particles
  //hemocell.enableBoundaryParticles((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["BRepCutoff"].read<T>(),(*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  //Restructure atomic blocks on processors when possible
  //hemocell.doRestructure(false); // cause errors

  if (hemocell.iter == 0) {
    hlog << "(Flowchamber) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  hlog << "(Flowchamber) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    hemocell.iterate();

    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    // Only enable if PARMETIS build is available
    /*
     if (hemocell.iter % tbalance == 0) {
       if(hemocell.calculateFractionalLoadImbalance() > (*cfg)["parameters"]["maxFlin"].read<T>()) {
         hemocell.doLoadBalance();
         hemocell.doRestructure();
       }
     }
   */
    if (hemocell.iter % tmeas == 0) {
        hlog << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
        hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
        hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
        hlog << ", WBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "WBC_HO");
        hlog << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
        FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); T toMpS = param::dx / param::dt;
        hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
        ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); T topN = param::df * 1.0e12;
        hlog << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

        // Additional useful stats, if needed
        //finfo = FluidInfo::calculateForceStatistics(&hemocell);
        //Set force as required after this function;
        // setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
        //           DESCRIPTOR<T>::ExternalField::forceBeginsAt,
        //           hemo::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
        // pcout << "Fluid force, Minimum: " << finfo.min << " Maximum: " << finfo.max << " Average: " << finfo.avg << endl;
        // ParticleStatistics pinfo = ParticleInfo::calculateVelocityStatistics(&hemocell);
        // pcout << "Particle velocity, Minimum: " << pinfo.min << " Maximum: " << pinfo.max << " Average: " << pinfo.avg << endl;
        hemocell.writeOutput();
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  hlog << "(main) Simulation finished :) " << endl;

  return 0;
}
