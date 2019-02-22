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

/// A functional, used to instantiate bounce-back nodes at the locations
template<typename T>
class CapillaryEllipseDomain3D : public plb::DomainFunctional3D {
public:
    CapillaryEllipseDomain3D(double centerX_, double centerY_, double outerA_, double outerB_, double innerA_, double innerB_, double xFrom_, double xTo_, double yTop_, double yBottom_)
        : centerX(centerX_),            // Outer and inner ellipse with the same center points
          centerY(centerY_),            // 
          outerASqr(outerA_ * outerA_), // Outer ellipse major radius squared
          outerBSqr(outerB_ * outerB_), // Outer ellipse minor radius squared
          innerASqr(innerA_ * innerA_), // Inner ellipse major radius squared
          innerBSqr(innerB_ * innerB_), // Inner ellipse minor radius squared
          xFrom(xFrom_),                // Draw nodes from this x points
          xTo(xTo_),                    // Draw nodes up to 
          yTop(yTop_),
          yBottom(yBottom_)           // The endpoint of the outer (clipped) ellipse have to be continuous with the contracted domain
    {}
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return ( ( (iX < xTo && iX > xFrom ) && ( ((iX-centerX)*(iX-centerX)*outerBSqr + (iY-centerY) * (iY-centerY) * outerASqr ) >= outerBSqr*outerASqr) && (iY >= yTop || iY <= yBottom ) ) || // large ellipse
                 ( (iX <= xTo && iX >= xFrom ) && ((iX-centerX)*(iX-centerX)*innerBSqr + (iY-centerY) * (iY-centerY) * innerASqr ) <= innerASqr*innerBSqr) ); // small ellipse

    }
    virtual CapillaryEllipseDomain3D<T>* clone() const {
        return new CapillaryEllipseDomain3D<T>(*this);
    }
private:
    double centerX;
    double centerY;
    double outerASqr;
    double outerBSqr;
    double innerASqr;
    double innerBSqr;
    double xFrom;
    double xTo;
    double yTop;
    double yBottom;

};


/// A functional, used to instantiate bounce-back nodes at rectangle locations
template<typename T>
class CapillaryRectangleDomain3D : public plb::DomainFunctional3D {
public:
    CapillaryRectangleDomain3D(double xMin_, double xMax_, double yMin_, double yMax_)
        : xMin(xMin_),
          xMax(xMax_),
          yMin(yMin_),
          yMax(yMax_)

    {}
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return ( (iX <= xMax && iX >= xMin ) && (iY <= yMax || iY >= yMin ) ) ;

    }
    virtual CapillaryRectangleDomain3D<T>* clone() const {
        return new CapillaryRectangleDomain3D<T>(*this);
    }
private:
    double xMin;
    double xMax;
    double yMin;
    double yMax;

};

/// A functional, used to instantiate bounce-back nodes at the middle domain
template<typename T>
class CapillaryMiddleRectangleDomain3D : public plb::DomainFunctional3D {
public:
    CapillaryMiddleRectangleDomain3D(double xMin_, double xMax_, double yWall_, double ny_, double cD_)
        : xMin(xMin_),
          xMax(xMax_),
          yWall(yWall_),
          ny(ny_),
          cD(cD_)

    {}
    virtual bool operator() (plint iX, plint iY, plint iZ) const {
        return ( (iX <= xMax && iX >= xMin ) && ((iY <= yWall || iY >= ny-yWall ) ||  (iY >= yWall+cD && iY <= ny-yWall-cD)  ) ) ;

    }
    virtual CapillaryMiddleRectangleDomain3D<T>* clone() const {
        return new CapillaryMiddleRectangleDomain3D<T>(*this);
    }
private:
    double xMin;
    double xMax;
    double yWall;
    double ny;
    double cD;

};

//-------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  pcout << "(Capillary) (Parameters) calculating flow parameters" << endl;
  param::lbm_base_parameters((*cfg));
  param::printParameters();

// ---------------------------- Create geometry ------------------------------------------------


  pcout << "(Capillary) setting dimensions ..." << std::endl;

  plint extendedEnvelopeWidth = 2;  // Depends on the requirements of the ibmKernel. 4 or even 2 might be enough (also depends on dx)

  plint heightChannel = (*cfg)["domain"]["refDirN"].read<int>();
  double capillaryD = (*cfg)["domain"]["capillaryD"].read<double>() / param::dx ;
  
  plint nx = 8*heightChannel;
  plint ny = heightChannel;
  plint nz = heightChannel;
  
  pcout << "(Capillary) nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;
  pcout << "(Capillary) capillary diameter lbm: " << capillaryD << std::endl;

  double ellipseCutX = (double) 0.1875*nx; // If we have heightChannel=60, then nx = 480 and 0.1875*480 = 90 --> the straight capillary is 480-2*90= 300 cell long (150 um) 
  double wallNy = 2 ;    // number of wall nodes at the middle domain

  double outerMajorR = (double) ny-2*wallNy;
  double outerMinorR = (double) 0.5*(ny-2*wallNy); 
  
  double ellipseLineFrontIntersect = ellipseCutX-outerMajorR+capillaryD;
  double ellipseLineEndIntersect = (double) nx-ellipseCutX+outerMajorR-capillaryD;

  double innerMinorR = outerMinorR - capillaryD; 
  double innerMajorR = outerMajorR * innerMinorR / outerMinorR;
  // pcout << "(Capillary) oMajR: " << outerMajorR << " oMinR: " << outerMinorR << " inMInR: " << innerMinorR << " inMajR: " << innerMajorR << std::endl;
  if (innerMinorR < 1 || innerMajorR < 1) {
    pcout << "(Capillary) The generation of the capillary geometry is not possible";
    return -1;
  }

  double centerX = ellipseCutX;
  double centerY = (double) ny*0.5;
  // Solving a quadratic equation for ellipse and line intersection (cutting an ellipse with a vertical line)
  // a*x**2 + b*x + c = 0
  double a = 1.0;
  double b = -2.0*centerY;
  double c = centerY*centerY - outerMinorR*outerMinorR *(1 - (ellipseLineFrontIntersect-centerX)*(ellipseLineFrontIntersect-centerX) / (outerMajorR*outerMajorR));
  double d = b*b-4.0*a*c;
  double yTop = 0.0;
  double yBottom = 0.0;

  if (d<0.0){
    pcout << "(Capillary) Negative  discriminant: Check the geometry setup!" << endl;
    return -1;
  }
  else if (d == 0){
    yTop = -b / (2.0*a);
    yBottom = yTop;
  }
  else {
    yTop = (-b + std::sqrt(d)) /(2.0*a);
    yBottom = (-b - std::sqrt(d)) /(2.0*a);
  }

// -------------------------------------------------------------------------------------------------------------------------------
  
  pcout << "(Capillary) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),  //voxelizedDomain->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));

  // Front contraction with a rectangle
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new CapillaryRectangleDomain3D<T>(0.0, ellipseLineFrontIntersect+1, yTop+1, yBottom-1 ), new BounceBack<T, DESCRIPTOR> );
  
  // Left half ellipse
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new CapillaryEllipseDomain3D<T>(centerX, centerY, outerMajorR, outerMinorR, innerMajorR, innerMinorR,
                ellipseLineFrontIntersect+1, ellipseCutX, yTop, yBottom), new BounceBack<T, DESCRIPTOR> );
  
  // Middle rectangle
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new CapillaryMiddleRectangleDomain3D<T>(ellipseCutX, (double) nx-ellipseCutX, wallNy, (double) ny, capillaryD), new BounceBack<T, DESCRIPTOR> );

  // Right half ellipse
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new CapillaryEllipseDomain3D<T>(nx-centerX, centerY, outerMajorR, outerMinorR, innerMajorR, innerMinorR,
                (double) nx-ellipseCutX, ellipseLineEndIntersect-1, yTop, yBottom), new BounceBack<T, DESCRIPTOR> );

  // End contraction
  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                new CapillaryRectangleDomain3D<T>(ellipseLineEndIntersect-1, (double) nx , yTop+1, yBottom-1), new BounceBack<T, DESCRIPTOR> );

  hemocell.lattice->toggleInternalStatistics(false);
  //hemocell.lattice->periodicity().toggleAll(false);
  hemocell.lattice->periodicity().toggle(0,true);
  //hemocell.lattice->periodicity().toggle(2,true);

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

  double dForce = (*cfg)["domain"]["drivingForce"].read<double>();

  //Driving Force
  pcout << "(Capillary) (Fluid) Setting up driving Force" << endl; 
  double poiseuilleForce = dForce * (param::dx * param::dx * param::dt*param::dt / param::dm);
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
                    plb::Array<double, DESCRIPTOR<double>::d>(poiseuilleForce, 0.0, 0.0));

  pcout << "(Capillary) poiseuilleForce = " << poiseuilleForce << endl;

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<WbcHighOrderModel>("WBC_HO", WBC_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("WBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("WBC_HO", 1); //Micrometer! not LU

  //hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  //hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK, OUTPUT_FORCE_INNER_LINK, OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC, OUTPUT_VERTEX_ID, OUTPUT_CELL_ID};
  hemocell.setOutputs("WBC_HO", outputs);
  //hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  if (hemocell.iter == 0) {
    pcout << "(Capillary) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  // unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  pcout << "(Capillary) Starting simulation..." << endl;

  while (hemocell.iter < tmax ) {
    hemocell.iterate();
    
    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
   
    if (hemocell.iter % tmeas == 0) {
      pcout << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
      pcout << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      pcout << " | # of WBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "WBC_HO");
      //pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
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
