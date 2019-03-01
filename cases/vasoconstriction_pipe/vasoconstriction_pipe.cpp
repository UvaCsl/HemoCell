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
#include <fenv.h>


/// A functional, used to instantiate bounce-back nodes at the locations of the sphere
template<typename T>
class CylinderShapeDomain3D : public plb::DomainFunctional3D {
public:
	CylinderShapeDomain3D(plint Lxcirc_, plint Lycirc_, plint LradiusCyl_, plint xbegin_, plint xend_, double SradiusCyl_, plint L_constr_)
		:Lxcirc(Lxcirc_),
		 Lycirc(Lycirc_),
		 LradiusCyl(LradiusCyl_),
		 LradiusSqr(LradiusCyl*LradiusCyl),
		 xbegin(xbegin_), 
		 xend(xend_),
		 SradiusCyl(SradiusCyl_),
		 //SradiusSqr(SradiusCyl*SradiusCyl),
		 L_constr(L_constr_)
	{}
	
	virtual bool operator() (plint iX, plint iY, plint iZ) const {
		return ((iZ-Lxcirc)*(iZ-Lxcirc) + (iY-Lycirc)*(iY-Lycirc) >= LradiusSqr && iX >= 0 && iX < xbegin) ||
			//((iZ-Lxcirc)*(iZ-Lxcirc) + (iY-Lycirc)*(iY-Lycirc) >= (((LradiusCyl-SradiusCyl)/2)*(1-std::cos((2*iX)/((L_constr/10)*std::acos(-1))))+(SradiusCyl))*(((LradiusCyl-SradiusCyl)/2)*(1-std::cos((2*iX)/((L_constr/10)*std::acos(-1))))+(SradiusCyl)) && iX >= xbegin && iX < xend) || 
			//((iZ-Sxcirc)*(iZ-Sxcirc) + (iY-Sycirc)*(iY-Sycirc) <= SradiusSqr && iX >=100 && iX <=150) ||
			( (iZ-Lxcirc)*(iZ-Lxcirc) + (iY-Lycirc)*(iY-Lycirc) >= ( ((LradiusCyl-SradiusCyl)/2)*std::cos( (2*std::acos(-1)/L_constr)*iX - (L_constr*xbegin) ) + ( ((LradiusCyl-SradiusCyl)/2) + SradiusCyl ) ) * ( ((LradiusCyl-SradiusCyl)/2)*std::cos( (2*std::acos(-1)/L_constr)*iX - (L_constr*xbegin) ) + ( ((LradiusCyl-SradiusCyl)/2) + SradiusCyl ) )  && iX >= xbegin && iX < xend ) ||
			((iZ-Lxcirc)*(iZ-Lxcirc) + (iY-Lycirc)*(iY-Lycirc) >= LradiusSqr && iX >= xend && iX <= xend+xbegin);
			
	}

	virtual CylinderShapeDomain3D<T>* clone() const{
		return new CylinderShapeDomain3D<T>(*this);
	}
	
private:
	plint Lxcirc;
	plint Lycirc;
	plint LradiusCyl;
	plint LradiusSqr;
	plint xend;
	plint xbegin;
	double SradiusCyl;
	plint L_constr;
};

//-------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  pcout << "(vasoconstriction) (Parameters) calculating flow parameters" << endl;
  param::lbm_pipe_parameters((*cfg), (*cfg)["domain"]["refDirN"].read<int>()); 
  param::printParameters();

// ---------------------------- Create geometry ------------------------------------------------


  pcout << "(vasoconstriction) setting dimensions ..." << std::endl;

  plint Cfactor = 2;
  plint nx = 6*(*cfg)["domain"]["refDirN"].read<int>(); //100lu=50um 
  plint ny = (*cfg)["domain"]["refDirN"].read<int>()+Cfactor;
  plint nz = ny+Cfactor;

  plint LradiusCyl = (ny-Cfactor)/2;
  plint Lxcirc = ny/2;
  plint Lycirc = ny/2;

//  plint SradiusCyl = 20;
//  plint Sxcirc = ny/2;
//  plint Sycirc = ny/2;

  plint L_constr = 400;
  double perc_constr = 0.56; //28/50;
  plint xbegin = (nx-L_constr)/2;
  plint xend = xbegin + L_constr;
  
  double SradiusCyl = LradiusCyl * (1.0-perc_constr);
  
  pcout << "(vasoconstriction) (Parameters)";
  pcout << "xbegin= " << xbegin; 
  pcout << "xend= " << xend;
  pcout << "radius lare = " << LradiusCyl;
  pcout << "radius small= " << SradiusCyl << endl;
  
// ---------------------------------------------------------------------------------------------
  plint extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with with 2.
  
  pcout << "(vasoconstriction) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));


  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
	      new CylinderShapeDomain3D<T>(Lxcirc, Lycirc, LradiusCyl, xbegin, xend, SradiusCyl,L_constr),
	      new BounceBack<T, DESCRIPTOR> );  


  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(true);
 

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));


  //Driving Force
  pcout << "(vasoconstriction) (Fluid) Setting up driving Force" << endl; 
  double poiseuilleForce = 8 * param::nu_lbm * (param::u_lbm_max * 0.5) / LradiusCyl / LradiusCyl;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<double>::ExternalField::forceBeginsAt,
                    plb::Array<double, DESCRIPTOR<double>::d>(poiseuilleForce, 0.0, 0.0));

  pcout << "poiseuilleForce = " << poiseuilleForce << endl;

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 0.1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());


  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC,OUTPUT_VERTEX_ID,OUTPUT_CELL_ID};
  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  // store the average velocity without cells (compare cellular flow vel. to this to get viscosity)
  double cellFreeVel = 0;

  // if it is a fresh start do some warmup (note, viscosity will not be accurate at all when starting from a checkpoint).
  if (hemocell.iter == 0) {

    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    pcout << "(vasoconstriction) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
    FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
    cellFreeVel = finfo.avg;
    pcout << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS << endl;
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();

  pcout << "(vasoconstriction) Starting simulation..." << endl;

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
      pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
      pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
      pcout << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << cellFreeVel / finfo.avg << endl;
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

  pcout << "(vasoconstriction) Simulation finished :) " << endl;

  return 0;
}
