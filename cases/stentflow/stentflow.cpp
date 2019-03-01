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
	CylinderShapeDomain3D(plint Vxcirc_, plint Vycirc_, plint VradiusCyl_, plint H2xcirc_, plint H2zcirc_, plint H2radiusCyl_)
		:Vxcirc(Vxcirc_),
		 Vycirc(Vycirc_),
		 VradiusCyl(VradiusCyl_),
		 VradiusSqr(VradiusCyl*VradiusCyl),
		 H2xcirc(H2xcirc_), 
		 H2zcirc(H2zcirc_),
		 H2radiusCyl(H2radiusCyl_),
		 H2radiusSqr(H2radiusCyl*H2radiusCyl)
	{}
	
	virtual bool operator() (plint iX, plint iY, plint iZ) const {
		return ((iX-Vxcirc)*(iX-Vxcirc) + (iY-Vycirc)*(iY-Vycirc) <= VradiusSqr) ||
			((iX-H2xcirc)*(iX-H2xcirc) + (iZ-H2zcirc)*(iZ-H2zcirc) <= H2radiusSqr);
	}

	virtual CylinderShapeDomain3D<T>* clone() const{
		return new CylinderShapeDomain3D<T>(*this);
	}
	
private:
	plint Vxcirc;
	plint Vycirc;
	plint VradiusCyl;
	plint VradiusSqr;
	plint H2xcirc;
	plint H2zcirc;
	plint H2radiusCyl;
	plint H2radiusSqr;
};

//-------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  pcout << "(stentflow) (Parameters) calculating flow parameters" << endl;
  param::lbm_pipe_parameters((*cfg), (*cfg)["domain"]["refDirN"].read<int>()); 
  param::printParameters();

// ---------------------------- Create geometry ------------------------------------------------


  pcout << "(stentflow) setting dimensions ..." << std::endl;

  plint ny = (*cfg)["domain"]["refDirN"].read<int>(); 
  plint nz = ny;
  plint nx = 2 * ny;

  plint VradiusCyl = 20;
  plint Vxcirc = ny/2;
  plint Vycirc = ny/2;

  plint HradiusCyl = 20; 
  plint Hxcirc = ny/2;
  plint Hzcirc = ny/2;


// ---------------------------------------------------------------------------------------------
  plint extendedEnvelopeWidth = 2;  // Because we might use ibmKernel with with 2.
  
  pcout << "(stentflow) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.lattice = new MultiBlockLattice3D<double, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, extendedEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<double, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(1.0/param::tau));


  defineDynamics(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
	      new CylinderShapeDomain3D<T>(Vxcirc, Vycirc, VradiusCyl, Hxcirc, Hzcirc, HradiusCyl),
	      new BounceBack<T, DESCRIPTOR> );  


  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(true);
 

  hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));


  //Driving Force
  pcout << "(stentflow) (Fluid) Setting up driving Force" << endl; 
  double poiseuilleForce = param::u_lbm_max * 8 * param::nu_lbm / (nx*ny); //8 * param::nu_lbm * (param::u_lbm_max * 0.5) / rPipe / rPipe;
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

    pcout << "(stentflow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
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

  pcout << "(stentflow) Starting simulation..." << endl;

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

  pcout << "(stentflow) Simulation finished :) " << endl;

  return 0;
}
