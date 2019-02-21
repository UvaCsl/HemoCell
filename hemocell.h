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
#ifndef HEMOCELL_H
#define HEMOCELL_H

//Load Constants
#include "constant_defaults.h"
#include "config.h"

/* CORE libs */
#include "hemoCellFields.h"

/* IO */
#include "loadBalancer.h"
#include "profiler.h"

/* MECHANICS */
#include "cellMechanics.h"
#include "constantConversion.h"

/* Helpers */
#include "preInlet.h"

/* Always used palabos functions in case files*/
#ifndef COMPILING_HEMOCELL_LIBRARY
#include "voxelizeDomain.h"
#include "palabos3D.h"
#include "palabos3D.hh"
using namespace hemo;
using namespace plb;
#endif



namespace hemo { 

/*!
 * The HemoCell class contains all the information, data and methods to set up a
 * basic HemoCell simulation.
 *  
 * Most importantly:
 *    - It contains one Palabos FluidField
 *    - It contains one CellFields
 *
 * After everything is set HemoCell::iterate() is used to iterate the
 * simulation.
 */
class HemoCell {
  public:
  /**
   * Creates an hemocell object
   * 
   * @param configFileName the location of the main config file
   * 
   * Unfortunately, due to palabos regulations, it is required to pass the
   * commandline arguments
   */
  HemoCell(char * configFileName, int argc, char* argv[]);

  /*
   * Clean up hemocell
   */
  ~HemoCell();
  
  /**
   *  Set all the fluid nodes to these values
   * 
   *  @param rho the desired density in lbm units
   *  @param vel the desired macroscopic velocity of each node
   */
  void latticeEquilibrium(T rho, hemo::Array<T, 3> vel);

  /**
   * Initialice the cellfields structure (and thus also the particlefield)
   */
  void initializeCellfield();

 /**
  * Add a celltype
  * valid options for constructType are:
  * RBC_FROM_SPHERE <- RBC
  * ELLIPSOID_FROM_SPHERE <- platelet
  * STRING_FROM_VERTEXES ->von willibrand factor
  * use as addCelltype<RbcHO>("RBC", RBC_FROM_SPHERE) for example
  * Since it is a template, it must be in the header class, maybe move to .hh
  * file for readability ...
  */
  template<class Mechanics>
  void addCellType(string name, int constructType) {
    HemoCellField * cellfield = cellfields->addCellType(name,constructType);
    Mechanics * mechanics = new Mechanics(*cellfield->materialCfg, *cellfield);
    cellfield->mechanics = mechanics;
    cellfield->statistics();
  }

  /**
   * Set the output of a celltype
   * the outputs string should contain constants like VELOCITY_OUTPUT defined
   * in the constants_defaults file
   * 
   * @param outputs a vector of constants that define the desired output
   * @param name the name of the CellType ("RBC", "PLT")
   */
  void setOutputs(string name, vector<int> outputs);
  
  //Sets the repulsion constant and cutoff distance, also enables repulsion
  bool repulsionEnabled = false;
  bool boundaryRepulsionEnabled = false;
  void setRepulsion(T repulsionConstant, T repulsionCutoff);

  //Set the timescale separation of the particles of a particle type
  void setMaterialTimeScaleSeparation(string name, unsigned int separation);
  
  //Enable solidify mechanics of a celltype
  void enableSolidifyMechanics(string name) {
    hlog << "(HemoCell) Enabling Solidify Mechanics for " << name << " mechanical model" << endl;
    (*cellfields)[name]->doSolidifyMechanics = true;
  }
  
  //Set the separation of when velocity is interpolated to the particle
  void setParticleVelocityUpdateTimeScaleSeparation(unsigned int separation);

  //Set the timescale separation of the repulsion force for all particles
  void setRepulsionTimeScaleSeperation(unsigned int separation);

  void setSolidifyTimeScaleSeperation(unsigned int separation);

  //Set the timescale separation of the interior viscosity, in between update and raytracing (expensive) update
  void setInteriorViscosityTimeScaleSeperation(unsigned int separation, unsigned int separation_entire_grid);
  
  //Enable Boundary particles and set the boundary particle constants
  void enableBoundaryParticles(T boundaryRepulsionConstant, T boundaryRepulsionCutoff, unsigned int timestep = 1);
  
  //Set the minimum distance of the particles of a type to the solid, must be called BEFORE loadparticles
  void setInitialMinimumDistanceFromSolid(string name, T distance);
  
  //Set the output of the fluid field
  void setFluidOutputs(vector<int> outputs);
  
  //Set the output of the CEPAC field
  void setCEPACOutputs(vector<int> outputs);
  
  //Explicitly set the periodicity of the domain along the different axes
  void setSystemPeriodicity(unsigned int axis, bool bePeriodic);

  //Set the number or times a particle should be able to wrap around in a certain direction (both negative and positive (default: 100)
  void setSystemPeriodicityLimit(unsigned int axis, int limit);
  
  //Load the particles
private:
  bool loadParticlesIsCalled = false;
public:
  ///Load the particles from their .pos files
  void loadParticles();

  ///Load a checkpoint
  void loadCheckPoint();
  
  ///Save a checkpoint
  void saveCheckPoint();

  ///Specify whether the output is in SI or LBM units
  bool outputInSiUnits = false;
  
  ///Write the specified output to hdf5 files
  void writeOutput();
  
  /// Do an iteration, If the system is driven by an external vector, you must set 
  /// it again after calling iterate
  void iterate();

  /// Check if any exis signal was caught
  void checkExitSignals();

  //Load balancing library functions
  /// Calculate and return the fractional load imbalance 
  T calculateFractionalLoadImbalance();
  
  ///Load balance the domain (only necessary with nAtomic blocks > nMpi processors, also checkpoints
  void doLoadBalance();
  
  ///Restructure the grid, has an optional argument to specify whether a checkpoint from this iteration is available, default is YES!
  void doRestructure(bool checkpoint_avail = true);
  
  ///Initialize the fluid field with the given management, should be done after specifing the pre inlets and before initializing the cellfields
  void initializeLattice(MultiBlockManagement3D const & management);
 
  PreInlet * preInlet = 0;
  
  bool partOfpreInlet = false;
  
  map<plint,plint> BlockToMpi;
  
  LoadBalancer * loadBalancer = 0;
  ///The fluid lattice
  MultiBlockLattice3D<T, DESCRIPTOR> * lattice = 0;
  
  MultiBlockManagement3D * preinlet_management = 0, * lattice_management = 0;
  
  Config * cfg = 0;
  ///The cellfields contains the particle field and all celltypes
  HemoCellFields * cellfields = 0;
  unsigned int iter = 0;
  
  XMLreader * documentXML = 0; //Needed for legacy checkpoint reading TODO fix
  private:
  /// Store the last time (iteration) output occured
  unsigned int lastOutputAt = 0;
  std::chrono::high_resolution_clock::duration lastOutput = std::chrono::high_resolution_clock::duration::zero();
  
  /// To be run right before the first iteration, all checking should move here
  void sanityCheck();
  /// Checked in iteration, do sanity check when not yet done
  bool sanityCheckDone = false;
};
}
#endif // HEMOCELL_H
