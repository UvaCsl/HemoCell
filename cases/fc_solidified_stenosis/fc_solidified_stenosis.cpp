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
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "bindingField.h"
#include "preInlet.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;

/// A functional, used to instantiate bounce-back nodes at the locations of the sphere
template <typename T>
class StenosisShapeDomain3D : public plb::DomainFunctional3D {
public:
  StenosisShapeDomain3D(plint xtopL_, plint xtopR_, plint xcircL_, plint xcircR_, plint ycirc_, plint ytop_,
                        plint radiusCyl_, double a_, double bL_, double bR_, double y_)
      : //xbottomL(xbottomL_),
        //xbottomR(xbottomR_),
        xtopL(xtopL_), xtopR(xtopR_), xcircL(xcircL_), xcircR(xcircR_), ycirc(ycirc_),
        //ybottom(ybottom_),
        ytop(ytop_), radiusCyl(radiusCyl_), radiusSqr(radiusCyl * radiusCyl), a(a_), bL(bL_), bR(bR_), y(y_)

  {}
  virtual bool operator()(plint iX, plint iY, plint iZ) const {
    return ((iX - xcircL) * (iX - xcircL) + (iY - ycirc) * (iY - ycirc) <= radiusSqr) ||
           ((iX - xcircR) * (iX - xcircR) + (iY - ycirc) * (iY - ycirc) <= radiusSqr) ||
           (iX <= xcircR && iX >= xcircL && iY <= ytop) || (iX >= (iY - bL) / a && iX <= xcircL && iY <= y) ||
           (iX <= (iY - bR) / (-a) && iX >= xcircR && iY <= y);
  }
  virtual StenosisShapeDomain3D<T> * clone() const { return new StenosisShapeDomain3D<T>(*this); }

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

int main(int argc, char * argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  //-------------------------------------------------------------------------------------

  hlogfile << "(main) setting dimensions .." << endl;

  plint widthSt = 2 * (*cfg)["parameters"]["widthStenosis"].read<int>();
  ;                        //length of stenosis
  plint radiusCyl = 2 * 5; // rounding stenosis
  double pi = std::acos(-1);
  plint c_angle_degrees = (*cfg)["parameters"]["angleStenosis"].read<int>();       //60; //angle of stenosis
  double percentageSt = (*cfg)["parameters"]["percentageStenosis"].read<double>(); //0.8;
  double angle = (90 - c_angle_degrees) * pi / 180;                                //in degrees
  double c_angle = c_angle_degrees * pi / 180;
  //calculate raakpunten and b from y =ax+b
  double h = std::sin(angle) * radiusCyl;
  double w = std::cos(angle) * radiusCyl;
  double a = std::tan(c_angle);
  plint widthChannel =
      2 * (*cfg)["parameters"]["widthChannel"].read<int>(); //4*(*cfg)["domain"]["refDirN"].read<int>()+60;
  plint heightChannel =
      2 * (*cfg)["parameters"]["heightChannel"].read<int>(); //4*(*cfg)["domain"]["refDirN"].read<int>();
  //  cout << "HeightChannel = " << heightChannel << endl;
  //  cout << "percentageSt = " << percentageSt << endl;
  //  cout << "som = " << heightChannel*percentageSt << endl;
  //  cout << " widthCons = " << (heightChannel*percentageSt)/a << endl;
  plint widthConst = (heightChannel * percentageSt) / a;
  plint lengthChannel = 4 * (*cfg)["domain"]["refDirN"].read<int>() + widthSt + (2 * widthConst);

  plint nx = lengthChannel;
  plint ny = heightChannel;
  plint nz = widthChannel;

  plint ytop = heightChannel * percentageSt; //percentage of stenosis
  plint xtopL = nx / 2 - widthSt / 2;
  plint xtopR = nx / 2 + widthSt / 2;
  //  plint xbottomL = xtopL-2*46;
  //  plint xbottomR = xtopR+2*46;
  plint xcircL = xtopL + radiusCyl;
  plint xcircR = xtopR - radiusCyl;
  plint ycirc = ytop - radiusCyl;
  //plint ytop = top - radiusCyl;
  //  plint ybottom = 0;

  double xL = xcircL - w;
  double y = ycirc + h;
  double xR = xcircR + w;

  double bL = y - a * xL;
  double bR = y + a * xR;

  pcout << "parameters Stenosis: " << "pi = " << pi << "," << "angle = " << angle << "," << "a = " << a << ","
        << "bL = " << bL << "," << "bR = " << bR << "," << "y = " << y << "," << "ytopR = " << xtopR << ", " <<
      //  "xbottomR = " << xbottomR << ", " <<
      "ytopL = " << xtopL << ", " << "radiuscyl = " << radiusCyl << ", " << "w = " << w << "," << "h = " << h << ","
        << "widthConst = " << widthConst << "," << "ytop = " << ytop << ", " << endl;

  //-------------------------------------------------------------------------------------

  hlog << "(Parameters) setting lbm parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();


  hlog << "(fc_solidified_stenosis) (Fluid) Initializing Palabos Fluid Field" << endl;
  MultiBlockManagement3D management =
      defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>());
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
      management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
      defaultMultiBlockPolicy3D().getCombinedStatistics(),
      defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
      new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0 / param::tau));

  //Set up boundaries (to help preinlet creation)
  Box3D topChannel(0, nx - 1, ny - 1, ny - 1, 0, nz - 1);
  Box3D bottomChannel(0, nx - 1, 0, 0, 0, nz - 1);

  defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR>);
  defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR>);
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                 new StenosisShapeDomain3D<T>(xtopL, xtopR, xcircL, xcircR, ycirc, ytop, radiusCyl, a, bL, bR, y),
                 new BounceBack<T, DESCRIPTOR>(1.));

  hlog << "(PreInlets) creating preInlet" << endl;
  hemocell.preInlet = new hemo::PreInlet(&hemocell, management);

  // Setting Preinlet slice
  Box3D slice = hemocell.lattice->getBoundingBox();
  slice.x1 = slice.x0 = slice.x0 + 1;
  hemocell.preInlet->preInletFromSlice(Direction::Xneg, slice);

  hlog << "(Stl preinlet) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.initializeLattice(management);

  hemocell.preInlet->initializePreInlet();

  hemocell.lattice->periodicity().toggle(2, true);

  hlog << "(Stl preinlet) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl;

  hemocell.preInlet->createBoundary();

  hemocell.lattice->toggleInternalStatistics(false);

  hemocell.latticeEquilibrium(1., plb::Array<double, 3>(0., 0., 0.));

  hemocell.preInlet->calculateDrivingForce();

  hemocell.lattice->initialize();

  // Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.enableSolidifyMechanics("PLT");

  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,      OUTPUT_TRIANGLES,  OUTPUT_FORCE,      OUTPUT_FORCE_VOLUME,
                         OUTPUT_FORCE_BENDING, OUTPUT_FORCE_LINK, OUTPUT_FORCE_AREA, OUTPUT_FORCE_VISC};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,    OUTPUT_DENSITY,      OUTPUT_FORCE,    OUTPUT_SHEAR_RATE,
             OUTPUT_STRAIN_RATE, OUTPUT_SHEAR_STRESS, OUTPUT_BOUNDARY, OUTPUT_BINDING_SITES};
  hemocell.setFluidOutputs(outputs);

  //Define binding sites
  if (!hemocell.partOfpreInlet) {
    Box3D bb = hemocell.lattice->getBoundingBox();
    Box3D outlet(bb.x1 - 2, bb.x1, bb.y0, bb.y1, bb.z0, bb.z1);
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> * boundary =
        new BoundaryConditionInstantiator3D<T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T, DESCRIPTOR>>();
    boundary->addPressureBoundary0P(outlet, *hemocell.lattice, boundary::density);
    Box3D bindingbox = hemocell.lattice->getBoundingBox();
    hemocell.cellfields->populateBindingSites(&bindingbox);
  }


  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  //Restructure atomic blocks on processors when possible
  //hemocell.doRestructure(false); // cause errors(?)

  if (hemocell.iter == 0) {
    pcout << "(Stl preinlet) fresh start: warming up cell-free fluid domain for "
          << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();


  pcout << "(Stl preinlet) Starting simulation..." << endl;

  while (hemocell.iter < tmax) {
    //preinlet.update();
    hemocell.iterate();

    if (hemocell.partOfpreInlet) {
      //Set driving force as required after each iteration
      hemocell.preInlet->setDrivingForce();
    }

    hemocell.preInlet->applyPreInlet();

// Load-balancing! Only enable if PARMETIS build is available
#ifdef HEMO_PARMETIS
    if (hemocell.iter % tbalance == 0) {
      if (hemocell.calculateFractionalLoadImbalance() > (*cfg)["parameters"]["maxFlin"].read<double>()) {
        hemocell.doLoadBalance();
        hemocell.doRestructure();
      }
    }
#endif

    if (hemocell.iter % tmeas == 0) {
      pcout << "(main) Stats. @ " << hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      pcout << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
      pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell);
      double toMpS = param::dx / param::dt;
      pcout << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS
            << " m/s, rel. app. viscosity: " << (param::u_lbm_max * 0.5) / finfo.avg << endl;
      ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell);
      double topN = param::df * 1.0e12;
      pcout << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max
            << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

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
    if (hemocell.iter % tcheckpoint == 0) { hemocell.saveCheckPoint(); }
  }

  pcout << "(main) Simulation finished :) " << endl;

  return 0;
}
